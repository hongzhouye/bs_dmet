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

// eigen solver for Hermitian matrix
void _eigh_ (const MatrixXd& A, MatrixXd& U, VectorXd& D)
{
	esXd es;
	es.compute (A);
	D = es.eigenvalues ();
	U = es.eigenvectors ();
}

class HUBBARD
{
	public:
		int Nup, Ndn, K;
		double U;
		char BC;
		VectorXd nup, ndn, eup, edn;
		MatrixXd Fup, Fdn, occ, occb, Cup, Cdn, Pup, Pdn;
		//Matdeq Fsup, Fsdn, errsup, errsdn;
		HUBBARD (char *);
		void _init_ ();
		void _build_F_ ();
		double _error_ ();
		void _hubbard_general_ ();
		double _get_E_ ();
		void _print_ ();
};

// read parameters
HUBBARD::HUBBARD (char * name)
{
	FILE * pw = fopen (name, "r");
	fscanf (pw, "%d", &K);	fscanf (pw, "\n");
	fscanf (pw, "%d", &Nup);	fscanf (pw, "\n");
	fscanf (pw, "%d", &Ndn);	fscanf (pw, "\n");
	fscanf (pw, "%c", &BC);	fscanf (pw, "\n");
	fscanf (pw, "%lf", &U);
	fclose (pw);
}

// Initialization
void HUBBARD::_init_ ()
{
	// set the dimension for all matrices
	nup.setZero (K); ndn.setZero (K); eup.setZero (K); edn.setZero (K);
	Fup.setZero (K, K); Fdn.setZero (K, K); 
	occ.setZero (K, K); occb.setZero (K, K); 
	Cup.setZero (K, K); Cdn.setZero (K, K); 
	Pup.setZero (K, K); Pdn.setZero (K, K);

	// setup occupation matrices
	occ.topLeftCorner (Nup, Nup).diagonal ().setConstant (1.);
	occb.topLeftCorner (Ndn, Ndn).diagonal ().setConstant (1.);

	// setup nup/ndn vectors
	nup.setZero ();
	ndn.setZero ();
}

// build the Fock matrices for both spins
void HUBBARD::_build_F_ ()
{
	Fup.setZero ();
	Fdn.setZero ();

	int i;
	for (i = 0; i < K; i++)
	{
		Fup (i, i) = U * ndn (i);
		Fdn (i, i) = U * nup (i);
		Fup (i, jup(i)) = Fup (i, jdn(i)) = -1.;
		Fdn (i, jup(i)) = Fdn (i, jdn(i)) = -1.;
	}
	if (BC == 'a')
		Fup (0, K - 1) = Fup (K - 1, 0) = 
			Fdn (0, K - 1) = Fdn (K - 1, 0) = 1;
}

// calculate diis error
double HUBBARD::_error_ ()
{
	MatrixXd errup, errdn;
	errup = Fup * Pup - Pup * Fup;		
	errdn = Fdn * Pdn - Pdn * Fdn;		

	return (errup.norm () + errdn.norm ()) / (2. * K);
}

// solve the Hubbard model for the translational symmetric case
void HUBBARD::_hubbard_general_ ()
{
	// initialization 
	_init_ ();

	VectorXd nupnew, ndnnew;
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
		_eigh_ (Fup, Cup, eup);	
		_eigh_ (Fdn, Cdn, edn);	
		Pup = Cup * occ * Cup.transpose ();	nupnew = Pup.diagonal ();
		Pdn = Cdn * occb * Cdn.transpose ();	ndnnew = Pdn.diagonal ();

		er = ((nup - nupnew).norm () + (ndn - ndnnew).norm ()) / 2.;
		iter ++;
		printf ("%4d\t%.3e\t%.3e\n", iter, er, erc);	

		if (er < SCF_CONV)
		{
			cout << "SCF converges in " << iter << " cycles.\n\n";
			break;
		}
		else
		{
			// mix densities
			nup = mix_beta * nupnew + (1. - mix_beta) * nup;
			ndn = mix_beta * ndnnew + (1. - mix_beta) * ndn;
		}
	}
	cout << endl;
}

// get energy
double HUBBARD::_get_E_ ()
{
	int mu, nu;
	double Etot = 0.;
	for (mu = 0; mu < K; mu++)
		Etot += - (Pup (jup(mu), mu) + Pup (jdn(mu), mu) +
				Pdn (jup(mu), mu) + Pdn (jdn(mu), mu)) + 
			U * nup (mu) * ndn (mu) ;

	return Etot;
}

void HUBBARD::_print_ ()
{
	MatrixXd n_collect (K, 2);
	n_collect.block (0, 0, K, 1) = nup;
	n_collect.block (0, 1, K, 1) = ndn;
	MatrixXd e_collect (K, 2);
	e_collect.block (0, 0, K, 1) = eup;
	e_collect.block (0, 1, K, 1) = edn;
	cout << "Occupation:\n" << n_collect << "\n\n";
	cout << "energy levels:\n" << e_collect << "\n\n";
	cout << "Density matrix for spin alpha:\n" << Pup << "\n\n";
	cout << "Density matrix for spin beta:\n" << Pdn << "\n\n";
	cout << "Total Energy: " << _get_E_ () << "\n\n";
}

#endif
