#ifndef _HUBBARD_H_INCLUDED_
#define _HUBBARD_H_INCLUDED_

/* A class for the Hubbard model setup.
 * h, P, C, U etc. are for the mean-field bath.
 * In frag.h, we will define Hamiltonian's for
 * each fragment.
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <Eigen/Dense>
#include "scf.h"
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
		string BC;
		VectorXd n, e;
		MatrixXd h, F, occ, C, P;
		void _init_ ();
		void _hubbard_rhf_ ();
		double _get_E_ ();
		void _print_ ();
		MatrixXd _PFrag_ (int *, int);
		MatrixXd _PFrag_ (int *, int, const MatrixXd&);
};

// Initialization
void HUBBARD::_init_ ()
{
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
	for (int mu = 0; mu < K - 1; mu ++)
		h (mu, mu + 1) = h (mu + 1, mu) = -1.;
	h (0, K - 1) = h (K - 1, 0) = (BC == "a") ? (1.) : (-1.);
}

// solve the Hubbard model for the translational symmetric case
// restricted spin symmetry is assumed
void HUBBARD::_hubbard_rhf_ ()
{
	_init_ ();
	SCF scf;
	scf._init_hub_ (h, U, K, N);
	scf._scf_ ();
	P = scf.P;	C = scf.C;	F = scf.F;	e = scf.e;	n = P.diagonal ();
}

// get energy
double HUBBARD::_get_E_ ()
{
	int mu, nu;
	double Etot = 0.;

	Etot = ((h + F) * P).trace ();

	return Etot;
}

void HUBBARD::_print_ ()
{
	cout << "energy levels:\n" << e << "\n\n";
	cout << "Occupation:\n" << n << "\n\n";
	cout << "Density matrix:\n" << P << "\n\n";
	printf ("Initial SCF Energy: %18.16f\n\n", _get_E_ ());
}

MatrixXd HUBBARD::_PFrag_ (int *frag, int Nimp)
{
	MatrixXd PF;	PF.setZero (Nimp, Nimp);
	for (int mu = 0; mu < Nimp; mu++)
        for (int nu = 0; nu < Nimp; nu++)
            PF(mu, nu) = P(frag[mu], frag[nu]);
	return PF;
}

MatrixXd HUBBARD::_PFrag_ (int *frag, int Nimp, const MatrixXd& Ptot)
{
	MatrixXd PF;	PF.setZero (Nimp, Nimp);
	for (int mu = 0; mu < Nimp; mu++)
        for (int nu = 0; nu < Nimp; nu++)
            PF(mu, nu) = Ptot(frag[mu], frag[nu]);
	return PF;
}

#endif
