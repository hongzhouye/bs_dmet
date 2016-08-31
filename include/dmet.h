#ifndef _DMET_H_INCLUDED_
#define _DMET_H_INCLUDED_

#include <iostream>
#include <cstdio>
#include <Eigen/Dense>
#include "hf.h"
#include "schmidt.h"

using namespace std;
using namespace Eigen;

class DMET
{
    public:
        double _dmet_energy_ (MatrixXd&, double *, MatrixXd&, int);
        double _dmet_energy_ (MatrixXd&, double *, MatrixXd&, double *, int);
};

double DMET::_dmet_energy_ (MatrixXd& h, double *V, MatrixXd& P, int N)
{
    int i, mu, nu, la, si, K = h.rows ();
    double E = 0;

    for (mu = 0; mu < N; mu++)
        for (nu = 0; nu < K; nu++)
            E += h(mu, nu) * P(nu, mu);
    E *= 2.;

    // Szabo89book page 141
    /*for (mu = 0; mu < N; mu++)
        for (nu = 0; nu < K; nu++)
            for (la = 0; la < K; la++)
                for (si = 0; si < K; si++)
                    E += P(nu, mu) * P(la, si) *
                        (2. * V[index4(mu,si,nu,la,K)] -
                        V[index4(mu,si,la,nu,K)]);
    */

    // CHECK: mean-field 2PDM
    double *G = _darray_gen_ (K * K * K * K);
    for (mu = 0; mu < K; mu++)
        for (si = 0; si < K; si++)
            for (nu = 0; nu < K; nu++)
                for (la = 0; la < K; la++)
                    G[index4(mu,si,nu,la,K)] = 2. * P(nu, mu) * P(la, si) - P(la, mu) * P(nu, si);

    for (mu = 0; mu < N; mu++)
        for (nu = 0; nu < K; nu++)
            for (la = 0; la < K; la++)
                for (si = 0; si < K; si++)
                    E += G[index4(mu,nu,la,si,K)] * V[index4(mu,nu,la,si,K)];
    /*double sum = 0.;
    for (mu = 0; mu < K; mu++)
    	for (nu = 0; nu < K; nu++)
    		sum += G[index4(mu,nu,mu,nu,K)];
    printf ("sum = %f\n", sum);*/
    return E;
}


// for correlated case where 2PDM Gamma is needed too
double DMET::_dmet_energy_ (MatrixXd& h, double *V, MatrixXd& P, double *G, int N)
{
    int mu, nu, la, si, K = h.rows ();
    double E = 0;

    for (mu = 0; mu < N; mu++)
        for (nu = 0; nu < K; nu++)
            E += h(mu, nu) * P(nu, mu);
    E *= 2.;

    for (mu = 0; mu < N; mu++)
        for (nu = 0; nu < K; nu++)
            for (la = 0; la < K; la++)
                for (si = 0; si < K; si++)
                {
                    //cout << G[index4(mu,nu,la,si,K)] * V[index4(mu,nu,la,si,K)] << "\n";
                    //printf ("%d;%d;%d;%d;%d\n", mu, nu, la, si, index4(mu,nu,la,si,K));
                    E += G[index4(mu,nu,la,si,K)] * V[index4(mu,nu,la,si,K)];
                }
    return E;
}

#endif
