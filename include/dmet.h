#ifndef _DMET_H_INCLUDED_
#define _DMET_H_INCLUDED_

#include <iostream>
#include <cstdio>
#include <Eigen/Dense>
#include "hf.h"
#include "dfci.h"
#include "frag.h"
#include "bootstrap.h"
#include "scf.h"
#include "read.h"

using namespace std;
using namespace Eigen;

class DMET
{
    public:
        HUBBARD hub;
        BOOTSTRAP bs;
        void _dmet_init_ (char *);
        void _dmet_check_ ();
        double _dmet_energy_ (MatrixXd&, double *, MatrixXd&, int);
        double _dmet_energy_ (MatrixXd&, double *, MatrixXd&, double *, int);
};

void DMET::_dmet_init_ (char *fname)
{
    // read parameters from the input file
    _read_ (fname, hub, bs.frag);

    // set up ioff -- the lookup table
	int K = bs.frag.Nimp * 2;
	_gen_ioff_ (hub.K * (hub.K + 1) / 2);

    // Hubbard Hartree-Fock calculation
	hub._hubbard_rhf_ ();
    //hub._print_ ();

    // Bootstrap
    bs._init_ (hub);
    bs._bs_opt_ ();

    // DMET energy;
    _dmet_energy_ (bs.frag.h, bs.frag.V, bs.frag.dfci.P,
        bs.frag.dfci.G, bs.frag.dfci.N);
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
                    G[cpind(mn,ls)] += 2. * P(nu, mu) * P(la, si) - P(la, mu) * P(nu, si);
                }
        }

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
    int mu, nu, mn, la, si, ls, K = h.rows ();
    const FRAG& frag = bs.frag;

    // center
    double *Ec1, *Ec2;
    Ec1 = _darray_gen_ (frag.Ncenter);
    Ec2 = _darray_gen_ (frag.Ncenter);
    for (int i = 0; i < frag.Ncenter; i++)
    {
        int ind = 0;
        for (auto I = frag.fragments.begin (); I != frag.fragments.end (); I++)
            if (I->second == frag.center[i])    break;
            else    ind++;

        for (nu = 0; nu < K; nu++)
            Ec1[i] += h(ind, nu) * P(nu, ind);
        Ec1[i] *= 2.;

        for (nu = 0; nu < K; nu++)
        {
            mn = cpind(ind,nu);
            for (la = 0; la < K; la++)
                for (si = 0; si < K; si++)
                {
                    ls = cpind(la,si);
                    Ec2[i] += G[cpind(mn,ls)] * V[cpind(mn,ls)];
                }
        }
    }

    // final energy
    double E1 = 0, E2 = 0;
    for (mu = 0; mu < N; mu++)
        for (nu = 0; nu < K; nu++)
            E1 += h(mu, nu) * P(nu, mu);
    E1 /= (N / 2.);

    // new storage style
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
    E2 /= N;

    cout << "========================\n";
    cout << "|      DMET ENERGY     |\n";
    cout << "========================\n";

    printf ("     ");
    for (int i = 0; i < frag.Ncenter; i++)
        printf ("center %d    ", frag.center[i]);
    printf ("final\n");
    printf ("1P   ");
    for (int i = 0; i < frag.Ncenter; i++)
        printf ("%10.7f  ", Ec1[i]);
    printf ("%10.7f\n", E1);
    printf ("2P   ");
    for (int i = 0; i < frag.Ncenter; i++)
        printf ("%10.7f  ", Ec2[i]);
    printf ("%10.7f\n", E2);
    printf ("tot  ");
    for (int i = 0; i < frag.Ncenter; i++)
        printf ("%10.7f  ", Ec1[i] + Ec2[i]);
    printf ("%10.7f\n", E1 + E2);

    return E1 + E2;
}

void DMET::_dmet_check_ ()
{
    cout << "=======================" << endl;
    cout << "|      read check     |" << endl;
    cout << "=======================" << endl;

    printf ("hub.K = %d\nhub.N = %d\nhub.U = %f\n", hub.K, hub.N, hub.U);
    cout << "hub.BC = " << hub.BC << "\n\n";
    cout << "Nimp: " << bs.frag.Nimp << "\n";
    for (auto i = bs.frag.fragments.begin (); i != bs.frag.fragments.end (); i++)
        cout << i->first << ": " << i->second << "\n";
    cout << "\n";
    cout << "Ncenter: " << bs.frag.Ncenter << "\n";
    for (int i = 0; i < bs.frag.Ncenter; i++)
        cout << bs.frag.center[i] << "\n";
    cout << "\n";
    cout << "Npop: " << bs.frag.Npop << "\n";
    for (int i = 0; i < bs.frag.Npop; i++)
    {
        for (int j = 0; j < bs.frag.popcon[i].size (); j++)
            cout << i << ", " << j << ": " << bs.frag.popcon[i][j] << "\n";
        cout << endl;
    }
    cout << "N1e: " << bs.frag.N1e << "\n";
    for (int i = 0; i < bs.frag.N1e; i++)
    {
        for (int j = 0; j < bs.frag.bad_1econ[i].size (); j++)
            cout << "bad " << i << ", " << j << ": ("
                << bs.frag.bad_1econ[i][j][0]
                << bs.frag.bad_1econ[i][j][1] << ")\n";
        cout << "good " << i << ": ("
            << bs.frag.good_1econ[i][0]
            << bs.frag.good_1econ[i][1] << ")\n";
        cout << endl;
    }
    cout << "N2e: " << bs.frag.N2e << "\n";
    for (int i = 0; i < bs.frag.N2e; i++)
    {
        for (int j = 0; j < bs.frag.bad_2econ[i].size (); j++)
            cout << "bad " << i << ", " << j << ": " << bs.frag.bad_2econ[i][j] << "\n";
        cout << "good " << i << ": " << bs.frag.good_2econ[i] << "\n";
        cout << endl;
    }

    cout << "=======================" << endl;
    cout << "|       HF-in-HF      |" << endl;
    cout << "=======================" << endl;
    SCF scf;
    scf._init_ (bs.frag.h, bs.frag.V, 2 * bs.frag.Nimp, bs.frag.Nimp);
    scf._guess_ ("core");
    scf._scf_ ();
    printf ("HF-in-HF embedding energy: %18.16f\n\n",
		_dmet_energy_ (bs.frag.h, bs.frag.V, scf.P, scf.N));

    cout << "=======================" << endl;
    cout << "|      FCI-in-HF      |" << endl;
    cout << "=======================" << endl;
    DFCI dfci;
    dfci._init_ (bs.frag.h, bs.frag.V, 2 * bs.frag.Nimp, bs.frag.Nimp);
    cout << "FCI initialization succeeds!\n" << dfci.tot <<
		" alpha strings are generated!\n\n";
	dfci._dfci_ ();
	dfci._1PDM_ ();
	cout << "scf 1PDM:\n" << scf.P << "\n\n";
	cout << "dfci 1PDM:\n" << dfci.P << "\n\n";
    cout << "P_tot, fragment block:\n";
    for (int i = 0; i < bs.frag.Nimp; i++)
    {
        for (int j = 0; j < bs.frag.Nimp; j++)
            printf ("%10.7f\t", hub.P(bs.frag.sm.frag[i], bs.frag.sm.frag[j]));
        cout << "\n";
    }   cout << "\n";
    cout << "T^{dagger} P_tot T:\n" << bs.frag.sm.T.transpose () * hub.P * bs.frag.sm.T << "\n\n";
	//cout << "check idempotency:\n" << hred_scf.P * hred_scf.P << "\n\n";
	dfci._2PDM_ ();
	printf ("FCI-in-HF embedding energy: %18.16f\n\n",
			_dmet_energy_ (bs.frag.h, bs.frag.V, dfci.P, dfci.G, dfci.N));
}

#endif
