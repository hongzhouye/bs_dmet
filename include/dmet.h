#ifndef _DMET_H_INCLUDED_
#define _DMET_H_INCLUDED_

#include <iostream>
#include <cstdio>
#include <Eigen/Dense>
#include "hf.h"

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
                    E += P(nu, mu) * P(la, si) *
                        (2. * V[index4(mu,si,nu,la,K)] -
                        V[index4(mu,si,la,nu,K)]);
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
                    E += G[index4(mu,nu,si,la,K)] * V[index4(mu,nu,la,si,K)];
    return E;
}

#endif
