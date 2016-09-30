#ifndef _SCF_H_INCLUDED_
#define _SCF_H_INCLUDED_

/* An SCF solver for given h and V in ORTHOGONAL basis.
 * Spin-restricted case is assumed.*/

#include <Eigen/Dense>
#include <iostream>
#include <cstdio>
#include "hf.h"
#include "diis.h"

using namespace Eigen;
using namespace std;

#define cpind(i,j) (i>j)?(ioff[i]+j):(ioff[j]+i)
#define RCA_THR 1E-2    // thr to switch to diis
#define MAX_SCF_ITER 1000

class SCF
{
    private:
        void _rca_ (const MatrixXd&, double);
        double SCF_THR;
    public:
        string scf_type;
        int K, N;               // basis set size, # of e^-
        int diis, rca;
        MatrixXd h, P, C, F, occ;
        VectorXd e;             // energy levels
        double *V, U;           // U if hubbard
        double E_scf;           // electronic energy
        double scf_error;       // ||[F, P]||

        void _init_ (const MatrixXd&, double *, const int, const int,
            const double, const string, const int, const int);
        void _init_hub_ (const MatrixXd&, const double, const int, const int,
            const double, const string, const int, const int);
                                // SCF setup for HUBBARD, diis and rca are ON
        void _guess_ (string);
        void _fock_ ();
        void _fock_hub_ ();
        void _scf_ ();
        double _get_hf_E_ ();
        void _print_ ();
};

// default setup: diis and rca are both toggled
/*void SCF::_init_ (const MatrixXd& hinp, double *Vinp, int Nbs, int Ne)
{
    scf_type = "gen";
    K = Nbs;    N = Ne;
    h.setZero (K, K);   h = hinp;
    int lenh = K * (K + 1) / 2;
    int lenV = lenh * (lenh + 1) / 2;
    V = _darray_gen_ (lenV);
    for (int i = 0; i < lenV; i++)    V[i] = Vinp[i];

    // diis and rca default: ON
    diis = 1;   rca = 1;

    // setup matrix size
    P.setZero (K, K);    C.setZero (K, K);   F.setZero (K, K);
    occ.setZero (K, K);    occ.topLeftCorner (N, N).diagonal ().setOnes ();
}*/

// specify whether diis and rca are used or not
void SCF::_init_ (const MatrixXd& hinp, double *Vinp, const int Nbs, const int Ne,
    const double scf_thresh = 1E-8, const string guess = "core",
    const int idiis = 1, const int irca = 1)
{
    scf_type = "gen";
    K = Nbs;    N = Ne;
    h.setZero (K, K);   h = hinp;
    int lenh = K * (K + 1) / 2;
    int lenV = lenh * (lenh + 1) / 2;
    V = _darray_gen_ (lenV);
    for (int i = 0; i < lenV; i++)    V[i] = Vinp[i];

    // scf acceleration
    diis = idiis;   rca = irca;

    // setup matrix size
    P.setZero (K, K);    C.setZero (K, K);   F.setZero (K, K);
    occ.setZero (K, K);    occ.topLeftCorner (N, N).diagonal ().setOnes ();

    // initial guess
    SCF_THR = scf_thresh;
    _guess_ (guess);
}

void SCF::_init_hub_ (const MatrixXd& hubh, const double hubU, const int Nbs,
    const int Ne, const double scf_thresh = 1E-8, const string guess = "core",
    const int idiis = 0, const int irca = 1)
{
    scf_type = "hub";
    K = Nbs;    N = Ne;
    h.setZero (K, K);   h = hubh;
    U = hubU;

    // scf acceleration
    diis = idiis;   rca = irca;

    // setup matrix size
    P.setZero (K, K);    C.setZero (K, K);   F.setZero (K, K);
    occ.setZero (K, K);    occ.topLeftCorner (N, N).diagonal ().setOnes ();

    // initial guess
    SCF_THR = scf_thresh;
    _guess_ (guess);
}

// SCF initialization
void SCF::_guess_ (string type)
{
    // core guess for F
    if (type == "core")
    {
        _eigh_ (h, C, e);
        P = C * occ * C.transpose ();
        (scf_type == "gen") ? (_fock_ ()) : (_fock_hub_ ());
    }
}

void SCF::_fock_ ()
{
    MatrixXd G; G.setZero (K, K);
    int mu, nu, la, si, mn, ms, ls, ln;

	for (nu = 0; nu < K; nu++) for (mu = 0; mu < K; mu++)
    {
        mn = cpind(mu,nu);
        for (si = 0; si < K; si++)
        {
            ms = cpind(mu,si);
            for (la = 0; la < K; la++)
            {
                ls = cpind(la,si);
                ln = cpind(la,nu);
				G(mu, nu) += P(si, la) * (2. * V[cpind(mn,ls)] - V[cpind(ms,ln)]);
            }
        }
    }
	F = G + h;
}

void SCF::_fock_hub_ ()
{
    VectorXd n = P.diagonal ();
    F = h;
    F += (U * n).asDiagonal ();
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
    int iter = 1;
    double Eold;
    MatrixXd Pold;  Pold.setZero (K, K);
    MatrixXd err;   err.setZero (K, K);
    DIIS Diis;  Diis._diis_init_ (K);

    //cout << "#iter\tscf error\n";
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
        (scf_type == "gen") ? (_fock_ ()) : (_fock_hub_ ());

        // calculating error and doing rca/diis
        err = F * P - P * F;
        scf_error = err.norm () / K;
        if (diis && rca && scf_error < RCA_THR)
        {
            rca = 0;
            //cout << "\nEnd RCA, Start DIIS...\n\n";
        }
        Diis.Fs.append (F);
        Diis.errs.append (err);
        if (diis && !rca)   F = Diis._next_diis_rhf_ ();

        // check convergence
        if (scf_error < SCF_THR)    break;

        // print and increase iter
        //printf ("%4d\t%10.7e\n", iter, scf_error);
        iter ++;
    }
    if (iter >= MAX_SCF_ITER)
    {
        cout << "Unfortunately, the SCF subroutine REFUSES to convege!\t"
            << scf_error << "\n\n";
        exit (1);
    }
    else
    {
    //    cout << "\nSCF converges in " << iter << " iterations!\n\n";
    //    _print_ ();
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

#undef RCA_THR
#endif
