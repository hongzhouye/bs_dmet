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
		int Nup, Ndn, K;
		double U;
		char BC;
		VectorXd nup, ndn, eup, edn;
		MatrixXd hup, hdn, Fup, Fdn, occ, occb, Cup, Cdn, Pup, Pdn;
		void _init_ ();
		void _build_F_ ();
		double _error_ ();
		void _hubbard_general_ ();
		double _get_E_ ();
		void _print_ ();
};

// Initialization
void HUBBARD::_init_ ()
{
	int mu;

	// set the dimension for all matrices
	nup.setZero (K); ndn.setZero (K); 
	eup.setZero (K); edn.setZero (K);
	hup.setZero (K, K); hdn.setZero (K, K);
	Fup.setZero (K, K); Fdn.setZero (K, K); 
	occ.setZero (K, K); occb.setZero (K, K); 
	Cup.setZero (K, K); Cdn.setZero (K, K); 
	Pup.setZero (K, K); Pdn.setZero (K, K);

	// setup occupation matrices
	occ.topLeftCorner (Nup, Nup).diagonal ().setConstant (1.);
	occb.topLeftCorner (Ndn, Ndn).diagonal ().setConstant (1.);

	// CORE guess (Szabo89book page 148)
	nup.setZero ();
	ndn.setZero ();

	// setup h matrix
	for (mu = 0; mu < K - 1; mu ++)
		hup (mu, mu + 1) = hup (mu + 1, mu) = 
			hdn (mu, mu + 1) = hdn (mu + 1, mu) = -1.;
	hup (0, K - 1) = hup (K - 1, 0) =
		hdn (0, K - 1) = hdn (K - 1, 0) = (BC == 'a') ? (1.) : (-1.);
}

// build the Fock matrices for both spins
void HUBBARD::_build_F_ ()
{
	Fup = hup;
	Fdn = hdn;
	Fup += U * (ndn.asDiagonal ());
	Fdn += U * (nup.asDiagonal ());
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
