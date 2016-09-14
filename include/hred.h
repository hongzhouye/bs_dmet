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
};

void HRED::_xform_ (const HUBBARD& hub, const SCHMIDT& sm,
	MatrixXd& h, double *V)
{
	int i, j, k, l, mu, ij, kl, ijkl;
	MatrixXd T = sm.T;

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
#endif
