#ifndef _HUBBARD_H_INCLUDED_
#define _HUBBARD_H_INCLUDED_

#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include "hf.h"
#include "diis.h"
#include "matdeque.h"

#define jup(i) (i==K-1)?(0):(i+1)
#define jdn(i) (i==0)?(K-1):(i-1)

using namespace std;
using namespace Eigen;

// eigen solver for Hermitian matrix
void _eigh_ (const MatrixKd& A, MatrixKd& U, VectorKd& D)
{
	esKd es;
	es.compute (A);
	D = es.eigenvalues ();
	U = es.eigenvectors ();
}

class HUBBARD
{
	public:
		double U;
		VectorKd nup, ndn, eup, edn;
		MatrixKd Fup, Fdn, occ, occb, Cup, Cdn, Pup, Pdn;
		Matdeq Fsup, Fsdn, errsup, errsdn;
		HUBBARD (double, string);
		void _build_F_ ();
		void _diis_ ();
		double _error_ ();
		void _hubbard_general_ ();
		double _get_E_ ();
		void _print_ ();
};

// read U and mode
HUBBARD::HUBBARD (double Ui, string mode)
{
	U = Ui;

	// setup occupation matrices
	occ.setZero ();	occ.topLeftCorner (Nup, Nup).diagonal ().setConstant (1.);
	occb.setZero (); occb.topLeftCorner (Ndn, Ndn).diagonal ().setConstant (1.);

	// setup nup/ndn vectors
	if (mode == "core")
	{
		nup.setZero ();
		ndn.setZero ();
		_build_F_ ();
		_eigh_ (Fup, Cup, eup);	
		_eigh_ (Fdn, Cdn, edn);	
		Pup = Cup * occ * Cup.transpose ();	nup = Pup.diagonal ();
		Pdn = Cdn * occb * Cdn.transpose ();	ndn = Pdn.diagonal ();
		_build_F_ ();	
	}
	/*else if (mode == "random")
	{
		nup.setRandom ().cwiseAbs ();
		ndn.setRandom ().cwiseAbs ();
	}*/
	else
	{
		cout << "Please provide a correct value for *mode*!\n";
		exit (0);
	}

	// setup for diis
	Fsup.max_size = MAX_DIIS;
	Fsdn.max_size = MAX_DIIS;
	errsup.max_size = MAX_DIIS;
	errsdn.max_size = MAX_DIIS;
	Fsup.allocate ();	Fsdn.allocate ();
	errsup.allocate ();	errsdn.allocate ();
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
}

// calculate diis error
double HUBBARD::_error_ ()
{
	MatrixKd errup, errdn;
	errup = Fup * Pup - Pup * Fup;		
	errdn = Fdn * Pdn - Pdn * Fdn;		
	errsup.append (errup);	errsdn.append (errdn);

	//cout << "errup:\n" << errup << "\n\n";
	//cout << "errdn:\n" << errdn << "\n\n";

	return (errup.norm () + errdn.norm ()) / (2. * K);
}

// solve the Hubbard model for the translational symmetric case
void HUBBARD::_hubbard_general_ ()
{
	VectorKd nupnew = nup, ndnnew = ndn;
	double er = 100., mix_beta = 0.7;
	int iter = 0;

	printf ("#iter\terror\n");
	while (iter < SCF_ITER && er > SCF_CONV)
	{
		if (iter < 100)
		{
			nup = mix_beta * nupnew + (1. - mix_beta) * nup;
			ndn = mix_beta * ndnnew + (1. - mix_beta) * ndn;
		}
		_eigh_ (Fup, Cup, eup);	
		_eigh_ (Fdn, Cdn, edn);	
		
		Pup = Cup * occ * Cup.transpose ();	nup = Pup.diagonal ();
		Pdn = Cdn * occb * Cdn.transpose ();	ndn = Pdn.diagonal ();

		_build_F_ ();
		Fsup.append (Fup);	Fsdn.append (Fdn);
		//cout << "Fup before diis:\n" << Fup << "\n\n";

		er = _error_ ();
		if (iter < 100)
		{
			_eigh_ (Fup, Cup, eup);	
			_eigh_ (Fdn, Cdn, edn);	
			Pup = Cup * occ * Cup.transpose ();	nupnew = Pup.diagonal ();
			Pdn = Cdn * occb * Cdn.transpose ();	ndnnew = Pdn.diagonal ();
		}
		else
		{
			_next_diis_ (errsup, Fsup, errsdn, Fsdn);
			Fup = Fsup.element[Fsup.now];
			Fdn = Fsdn.element[Fsdn.now];
		}
		//cout << "Fup after diis:\n" << Fup << "\n\n";

		iter ++;
		printf ("%4d\t%.3e\n", iter, er);	
	}
	cout << endl;
}

// get energy
double HUBBARD::_get_E_ ()
{
	int mu, nu;
	double Etot = 0.;
	for (mu = 0; mu < K; mu++)
		Etot += - (Pup (jup(mu), mu) + Pdn (jdn(mu), mu)) + U * nup (mu) * ndn (mu) ;

	return Etot;
}

void HUBBARD::_print_ ()
{
	Matrix<double, K, 2> n_collect;
	n_collect.block (0, 0, K, 1) = nup;
	n_collect.block (0, 1, K, 1) = ndn;
	cout << "Occupation:\n" << n_collect << "\n\n";
	cout << "Density matrix for spin alpha:\n" << Pup << "\n\n";
	cout << "Density matrix for spin beta:\n" << Pdn << "\n\n";
	cout << "Total Energy: " << _get_E_ () << "\n\n";
}

#endif
