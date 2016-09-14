#ifndef _DMET_H_INCLUDED_
#define _DMET_H_INCLUDED_

#include <iostream>
#include <cstdio>
#include <Eigen/Dense>
#include "hf.h"
#include "schmidt.h"
#include "hred.h"
#include "dfci.h"
#include "scf.h"
#include "read.h"

using namespace std;
using namespace Eigen;

class DMET
{
    public:
        HUBBARD hub;
        SCHMIDT sm;
        HRED hr;
        void _dmet_init_ (char *);
        void _dmet_iter_ ();
        void _dmet_check_ ();
        double _dmet_energy_ (MatrixXd&, double *, MatrixXd&, int);
        double _dmet_energy_ (MatrixXd&, double *, MatrixXd&, double *, int);
};

void DMET::_dmet_init_ (char *fname)
{
    // read parameters from the input file
    _read_ (fname, hub, sm);

    // set up ioff -- the lookup table
	int K = sm.Nimp * 2;
	_gen_ioff_ (K * (K + 1) / 2);

    // Hubbard Hartree-Fock calculation
	hub._hubbard_rhf_ ();
    //hub._print_ ();

    // Schmidt
	sm._schmidt_ (hub);
    //sm._print_ (hub);

	// Construct Hred
	hr._xform_ (hub, sm);
}

void DMET::_dmet_check_ ()
{
    SCF scf;
    scf._init_ (hr.h, hr.V, hr.Ni, hr.Ni / 2);
    scf._guess_ ("core");
    scf._scf_ ();
    printf ("HF-in-HF embedding energy: %18.16f\n\n",
		_dmet_energy_ (scf.h, scf.V, scf.P, scf.N));

    DFCI dfci;
    dfci._init_ (hr);
    cout << "FCI initialization succeeds!\n" << dfci.tot <<
		" alpha strings are generated!\n\n";
	dfci._dfci_ ();
	dfci._1PDM_ ();
	cout << "scf 1PDM:\n" << scf.P << "\n\n";
	cout << "dfci 1PDM:\n" << dfci.P << "\n\n";
    cout << "P_tot, fragment block:\n";
    for (int i = 0; i < sm.Nimp; i++)
    {
        for (int j = 0; j < sm.Nimp; j++)
            printf ("%10.7f\t", hub.P(sm.frag[i], sm.frag[j]));
        cout << "\n";
    }   cout << "\n";
    cout << "T^{dagger} P_tot T:\n" << sm.T.transpose () * hub.P * sm.T << "\n\n";
	//cout << "check idempotency:\n" << hred_scf.P * hred_scf.P << "\n\n";
	dfci._2PDM_ ();
	printf ("FCI-in-HF embedding energy: %18.16f\n\n",
			_dmet_energy_ (hr.h, scf.V, dfci.P, dfci.G, dfci.N));
}

void DMET::_dmet_iter_ ()
{

}

// Mean-field case: 2PDM is not needed
double DMET::_dmet_energy_ (MatrixXd& h, double *V, MatrixXd& P, int N)
{
    int i, mu, nu, la, si, mn, ls, mnls, K = h.rows ();
    double E = 0;

    for (mu = 0; mu < N; mu++)
        for (nu = 0; nu < K; nu++)
            E += h(mu, nu) * P(nu, mu);
    E *= 2.;

    // Szabo89book page 141
    /*int ml, ns;
    for (mu = 0; mu < N; mu++)
        for (nu = 0; nu < K; nu++)
        {
            mn = cpind(mu,nu);
            for (la = 0; la < K; la++)
            {
                ml = cpind(mu,la);
                for (si = 0; si < K; si++)
                {
                    ls = cpind(si,la);  ns = cpind(si,nu);
                    E += P(nu, mu) * P(la, si) *
                        (2. * V[cpind(mn,ls)] -
                        V[cpind(ml,ns)]);
                }
            }
        }
    */

    // CHECK: mean-field 2PDM
    int lenh = K * (K + 1) / 2;
    int lenV = lenh * (lenh + 1) / 2;
    double *G = _darray_gen_ (lenV);
    for (mu = 0; mu < K; mu++)
        for (nu = 0; nu < K; nu++)
        {
            mn = cpind(mu,nu);
            for (si = 0; si < K; si++)
                for (la = 0; la < K; la++)
                {
                    ls = cpind(la,si);
                    //if (ls <= mn)
                    G[cpind(mn,ls)] += 2. * P(nu, mu) * P(la, si) - P(la, mu) * P(nu, si);
                }
        }

    // CHECK
    /*int ml, ns;
    for (mu = 0; mu < K; mu++)
        for (nu = 0; nu < K; nu++)
            for (la = 0; la < K; la++)
            {
                ml = cpind(mu,la);
                for (si = 0; si < K; si++)
                {
                    ns = cpind(nu,si);
                    printf ("%d;%d;%d;%d;%18.16f\n", mu, nu, la, si, G[cpind(ml,ns)]);
                }
            }
    */

    for  (mu = 0; mu < K; mu++)
        for (nu = 0; nu < N; nu++)
        {
            mn = cpind(mu,nu);
            for (la = 0; la < K; la++)
                for (si = 0; si < K; si++)
                {
                    ls = cpind(la,si);  mnls = cpind(mn,ls);
                    E += G[mnls] * V[mnls];
                }
        }

    return E;
}


// for correlated case where 2PDM Gamma is needed too
double DMET::_dmet_energy_ (MatrixXd& h, double *V, MatrixXd& P, double *G, int N)
{
    int mu, nu, la, si, K = h.rows ();
    double E1 = 0, E2 = 0;

    for (mu = 0; mu < N; mu++)
        for (nu = 0; nu < K; nu++)
            E1 += h(mu, nu) * P(nu, mu);
    E1 /= (sm.Nimp / 2.);
    printf ("One electron part: %18.16f\n\n", E1);

    // old storage style of V and G
    /*for (mu = 0; mu < N; mu++)
        for (nu = 0; nu < K; nu++)
            for (la = 0; la < K; la++)
                for (si = 0; si < K; si++)
                {
                    //cout << G[index4(mu,nu,la,si,K)] * V[index4(mu,nu,la,si,K)] << "\n";
                    //printf ("%d;%d;%d;%d;%d\n", mu, nu, la, si, index4(mu,nu,la,si,K));
                    E += G[index4(mu,nu,la,si,K)] * V[index4(mu,nu,la,si,K)];
                }
    */

    // new storage style
    int mn, ls;
    for (mu = 0; mu < N; mu++)
        for (nu = 0; nu < K; nu++)
        {
            mn = cpind(mu,nu);
            for (la = 0; la < K; la++)
                for (si = 0; si < K; si++)
                {
                    ls = cpind(la,si);
                    E2 += G[cpind(mn,ls)] * V[cpind(mn,ls)];
                }
        }
    E2 /= sm.Nimp;
    printf ("Two electron part: %18.16f\n\n", E2);

    return E1 + E2;
}

#endif
