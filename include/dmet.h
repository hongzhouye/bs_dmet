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
    double *G = _darray_gen_ (K * K * K * K);
    for (mu = 0; mu < K; mu++)
        for (nu = 0; nu <= mu; nu++)
        {
            mn = cpind(mu,nu);
            for (si = 0; si < K; si++)
                for (la = 0; la <= si; la++)
                {
                    ls = cpind(la,si);
                    if (ls <= mn)   G[cpind(mn,ls)] = 2. * P(nu, mu) * P(la, si) - P(la, mu) * P(nu, si);
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

    /*double sum = 0.;    int mm, nn;
    for (mu = 0; mu < K; mu++)
    {
        mm = cpind(mu,mu);
    	for (nu = 0; nu < K; nu++)
        {
            nn = cpind(nu,nu);
    		sum += G[cpind(mm,nn)];
        }
    }
    printf ("sum = %f\n", sum);
    */
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
