#ifndef _HRED_H_INCLUDED_
#define _HRED_H_INCLUDED_

#include <Eigen/Dense>
#include <cstdio>
#include "hf.h"
#include "hubbard.h"
#include "schmidt.h"

// Reduced Hamiltonian
class HRED
{
	public:
		int Ni;			// 2 * Nimp
		//MatrixXd h;
		//double *V;
		void _xform_ (const HUBBARD&, const SCHMIDT&, MatrixXd&, double *);
		void _xform_ (const MatrixXd&, double *, const SCHMIDT&, MatrixXd& ,
			double *V);
};

// hacking _xform_ for HUBBARD model
void HRED::_xform_ (const HUBBARD& hub, const SCHMIDT& sm,
	MatrixXd& h, double *V)
{
	int i, j, k, l, mu, ij, kl, ijkl;
	const MatrixXd& T = sm.T;

	// Environment's contribution to himp
	MatrixXd hc = ((sm.TE * sm.TE.transpose ()).diagonal () * hub.U).asDiagonal ();

	// himp
	h = T.transpose () * (hub.h + hc) * T;

	// Vimp
	Ni = 2 * sm.Nimp;

	for (i = 0; i < Ni; i++)
		for (j = 0; j <= i; j++)
		{
			ij = cpind(i,j);
			for (k = 0; k < Ni; k++)
				for (l = 0; l <= k; l++)
				{
					kl = cpind(k,l);
					if (kl <= ij)
					{
						ijkl = cpind(ij,kl);
						for (mu = 0; mu < hub.K; mu++)	V[ijkl] += T(mu, i) * T(mu, k) * T(mu, j) * T(mu, l);
						V[ijkl] *= hub.U;
					}
				}
		}
}

// general case
void HRED::_xform_ (const MatrixXd& hinp, double *Vinp, const SCHMIDT& sm,
	MatrixXd& h, double *V)
{
	int i, j, k, l, ij, kl, ijkl, ik, jl;
	const MatrixXd& T = sm.T;

	// Environment's contribution to himp
	MatrixXd PE = sm.TE * sm.TE.transpose ();	int K = PE.rows ();
	MatrixXd hc (K, K);	hc.setZero ();
	for (i = 0; i < K; i++)	for (j = 0; j < K; j++)
		for (k = 0; k < K; k++)	for (l = 0; l < K; l++)
		{
			ij = cpind(i,j);	kl = cpind(k,l);
			ik = cpind(i,k);	jl = cpind(j,l);
			hc(i, j) += PE(k, l) * (2. * Vinp[cpind(ij,kl)] - Vinp[cpind(ik,jl)]);
		}

	// himp
	Ni = 2 * sm.Nimp;

	h = T.transpose () * (hinp + hc) * T;

	// Vimp
	int mu, nu, mn, la, si, ls;
	for (i = 0; i < Ni; i++)	for (j = 0; j <= i; j++)
	{
		ij = cpind(i,j);
		for (k = 0; k < Ni; k++)	for (l = 0; l <= k; l++)
		{
			kl = cpind(k,l);
			if (kl <= ij)
			{
				ijkl = cpind(ij,kl);
				for (mu = 0; mu < K; mu++)	for (nu = 0; nu < K; nu++)
				{
					mn = cpind(mu,nu);
					for (la = 0; la < K; la++)	for (si = 0; si < K; si++)
					{
						ls = cpind(la,si);
						V[ijkl] += sm.T(mu, i) * sm.T(nu, j) *
							sm.T(la, k) * sm.T(si, l) * Vinp[cpind(mn,ls)];
					}
				}
			}
		}
	}
}

#endif
