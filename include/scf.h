#ifndef _SCF_H_INCLUDED_
#define _SCF_H_INCLUDED_

/* An SCF solver for given h and V in ORTHOGONAL basis.
 * Spin-restricted case is assumed.*/

#include <Eigen/Dense>
#include <iostream>
#include <cstdio>
#include "hf.h"
#include "hred.h"
#include "diis.h"

using namespace Eigen;
using namespace std;

#define cpind(i,j) (i>j)?(ioff[i]+j):(ioff[j]+i)
#define RCA_THR 1E-2    // thr to switch to diis
#define SCF_THR 1E-10
#define MAX_SCF_ITER 1000

class SCF
{
    private:
        void _print_ ();
        void _rca_ (const MatrixXd&, double);
    public:
        int K, N;               // basis set size, # of e^-
        int diis, rca;
        MatrixXd h, P, C, F, occ;
        VectorXd e;             // energy levels
        double *V;
        double E_scf;           // electronic energy
        double scf_error;       // ||[F, P]||

        SCF (MatrixXd&, double *, int, int);
                                // SCF setup, diis and rca are ON
        SCF (MatrixXd&, double *, int, int, int ,int);
                                // SCF setup, give diis and rca
        void _init_ ();
        void _fock_ ();
        void _scf_ ();
        double _get_hf_E_ ();
};

// default setup: diis and rca are both toggled
SCF::SCF (MatrixXd& hinp, double *Vinp, int Nbs, int Ne)
{
    K = Nbs;    N = Ne;
    h.setZero (K, K);   h = hinp;
    int K4 = K * K; K4 = K4 * K4;
    V = _darray_gen_ (K4);
    for (int i = 0; i < K4; i++)    V[i] = Vinp[i];

    // diis and rca default: ON
    diis = 1;   rca = 1;
}

// specify whether diis and rca are used or not
SCF::SCF (MatrixXd& hinp, double *Vinp, int Nbs, int Ne, int idiis, int irca)
{
    K = Nbs;    N = Ne;
    h.setZero (K, K);   h = hinp;
    int K4 = K * K; K4 = K4 * K4;
    V = _darray_gen_ (K4);
    for (int i = 0; i < K4; i++)    V[i] = Vinp[i];

    // scf acceleration
    diis = idiis;   rca = irca;
}

// SCF initialization
void SCF::_init_ ()
{
    int i, j, k, l;

    // setup matrix size
    P.setZero (K, K);    C.setZero (K, K);   F.setZero (K, K);
    occ.setZero (K, K);    occ.topLeftCorner (N, N).diagonal ().setOnes ();

    // core guess for F
    _eigh_ (h, C, e);
    P = C * occ * C.transpose ();
    _fock_ ();
}

void SCF::_fock_ ()
{
    MatrixXd G; G.setZero (K, K);
    int mu, nu, la, si;

	for (nu = 0; nu < K; nu++)
		for (mu = 0; mu < K; mu++)
			for (la = 0; la < K; la++)
				for (si = 0; si < K; si++)
					G(mu, nu) += P(si, la) * (2. *
							V[index4(mu,la,nu,si,K)] -
							V[index4(mu,la,si,nu,K)]);

	F = G + h;
}

void SCF::_rca_ (const MatrixXd& Pold, double dE)
{
	double dE_deriv = ((P - Pold) * F).sum();
	double alpha = - dE_deriv / (2. * (dE - dE_deriv));
	if (!(alpha < 1. && alpha > 0.))
		if (dE < 0)	alpha = 0.999;
		else if (dE > 0) alpha = 0.001;
	MatrixXd Ptemp = alpha * P + (1. - alpha) * Pold;
	P = Ptemp;
}

void SCF::_scf_ ()
{
    _init_ ();
    int iter = 1;
    const double mixing_beta = 0.3;
    double Eold;
    MatrixXd Pold;  Pold.setZero (K, K);
    MatrixXd err;   err.setZero (K, K);
    DIIS Diis;  Diis._diis_init_ (K);
    cout << "#iter\tscf error\n";

    // check
    cout << "h:\n" << h << "\n\n";
    cout << "P:\n" << P << "\n\n";

    while (iter < MAX_SCF_ITER)
    {
        Eold = _get_hf_E_ ();
        Pold = P;

        // forming new P
        _eigh_ (F, C, e);
        P = C * occ * C.transpose ();

        // optimized damping
        if (rca)    _rca_ (Pold, _get_hf_E_ () - Eold);

        // forming new F
        _fock_ ();

        // calculating error and doing rca/diis
        err = F * P - P * F;
        scf_error = err.norm () / K;
        if (diis && rca && scf_error < RCA_THR)
        {
            rca = 0; cout << "\nEnd RCA, Start DIIS...\n\n";
        }
        Diis.Fs.append (F);
        Diis.errs.append (err);
        if (diis && !rca)   F = Diis._next_diis_rhf_ ();

        // check convergence
        if (scf_error < SCF_THR)    break;

        // print and increase iter
        printf ("%4d\t%10.7e\n", iter, scf_error);
        iter ++;
    }
    if (iter >= MAX_SCF_ITER)
    {
        cout << "\nSCF subroutine fails to convege!\n\n";
        exit (1);
    }
    else
    {
        cout << "\nSCF converges in " << iter << " iterations!\n\n";
        _print_ ();
    }
}

double SCF::_get_hf_E_ ()
{
    return ((h + F) * P).trace();
}

void SCF::_print_ ()
{
    E_scf = _get_hf_E_ ();
    printf ("SCF energy: %18.16f\n\n", E_scf);
    cout << "Energy levels:\n" << e << "\n\n";
	cout << "Density matrix:\n" << P << "\n\n";
}
#endif
