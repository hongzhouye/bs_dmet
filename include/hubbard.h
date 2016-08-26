#ifndef _HUBBARD_H_INCLUDED_
#define _HUBBARD_H_INCLUDED_

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <Eigen/Dense>
#include "hf.h"
//#include "diis.h"
//#include "matdeque.h"

#define jup(i) (i==K-1)?(0):(i+1)
#define jdn(i) (i==0)?(K-1):(i-1)

using namespace std;
using namespace Eigen;

class HUBBARD
{
	public:
		int N, K;
		double U;
		char BC;
		VectorXd n, e;
		MatrixXd h, F, occ, C, P;
		void _init_ ();
		void _build_F_ ();
		double _error_ ();
		void _hubbard_rhf_ ();
		double _get_E_ ();
		void _print_ ();
};

// Initialization
void HUBBARD::_init_ ()
{
	int mu;

	// set the dimension for all matrices/vectors
	n.setZero (K);
	e.setZero (K);
	h.setZero (K, K);
	F.setZero (K, K);
	occ.setZero (K, K);
	C.setZero (K, K);
	P.setZero (K, K);

	// setup occupation matrices
	occ.topLeftCorner (N, N).diagonal ().setConstant (1.);

	// CORE guess (Szabo89book page 148)
	n.setZero ();

	// setup h matrix
	for (mu = 0; mu < K - 1; mu ++)
		h (mu, mu + 1) = h (mu + 1, mu) = -1.;
	h (0, K - 1) = h (K - 1, 0) = (BC == 'a') ? (1.) : (-1.);
}

// build the Fock matrices for both spins
void HUBBARD::_build_F_ ()
{
	F = h;
	F += U * (n.asDiagonal ());
}

// calculate diis error
double HUBBARD::_error_ ()
{
	MatrixXd err;
	err = F * P - P * F;

	return err.norm () / K;
}

// solve the Hubbard model for the translational symmetric case
// restricted spin symmetry is assumed
void HUBBARD::_hubbard_rhf_ ()
{
	// initialization
	_init_ ();

	VectorXd nnew;
	double er = 100., erc, mix_beta = 0.7;
	int iter = 0;

	printf ("#iter\terror\t\t|[F, p]|\n");
	while (iter < SCF_ITER && er > SCF_CONV)
	{
		// build new F
		_build_F_ ();

		// erc = ||[F, P]||
		erc = _error_ ();

		// diagonalizing F
		_eigh_ (F, C, e);
		P = C * occ * C.transpose ();	nnew = P.diagonal ();

		er = (n - nnew).norm ();
		iter ++;
		printf ("%4d\t%.3e\t%.3e\n", iter, er, erc);

		if (er < SCF_CONV)
		{
			cout << "\nSCF converges in " << iter << " cycles.\n\n";
			break;
		}
		else
		{
			// mix densities
			n = mix_beta * nnew + (1. - mix_beta) * n;
		}
	}
	cout << "Hubbard C:\n" << C << "\n\n";
	cout << endl;
}

// get energy
double HUBBARD::_get_E_ ()
{
	int mu, nu;
	double Etot = 0.;
	/*for (mu = 0; mu < K - 1; mu++)
		Etot += -2. * (P (mu, mu + 1) + P (mu + 1, mu)) + U * n (mu) * n (mu);
	Etot += -2. * (P (0, K - 1) + P (K - 1, 0)) + U * n (K - 1) * n (K - 1);
	*/

	/*for (mu = 0; mu < N; mu++)	Etot += e(mu);
	Etot *= 2.;
	Etot -= U * (n.cwiseProduct (n)).sum ();*/

	Etot = ((h + F) * P).trace ();

	return Etot;
}

void HUBBARD::_print_ ()
{
	cout << "energy levels:\n" << e << "\n\n";
	cout << "Occupation:\n" << n << "\n\n";
	cout << "Density matrix:\n" << P << "\n\n";
	printf ("Total Energy: %18.16f\n\n", _get_E_ ());
}

#endif
