#ifndef _SCHMIDT_H_INCLUDED_
#define _SCHMIDT_H_INCLUDED_

#include <Eigen/Dense>
#include <iostream>
#include "hubbard.h"
#include "hf.h"

using namespace std;
using namespace Eigen;

class SCHMIDT
{
	private:
		void _form_xform_mat_ ();
	public:
		int K, N, Nimp;
		int *frag;
		MatrixXd W;				// M = W * d * W^T
		VectorXd d;				// Singular values
		MatrixXd C;				// W-transformed C's
		MatrixXd T, TE;
		//SCHMIDT (char*, HUBBARD&);
		void _schmidt_ (HUBBARD&);
		MatrixXd _make_P_frag_ (HUBBARD&);
		void _print_ ();
};

void SCHMIDT::_schmidt_ (HUBBARD& hub)
{
	int i;

	// prefix 'r' for raw (i.e. untransformed)
	MatrixXd rCF;	rCF.setZero (Nimp, N);
	for (i = 0; i < Nimp; i++)
		rCF.row (i) = hub.C.block (frag[i], 0, 1, N);

	MatrixXd M;
	M = rCF.transpose () * rCF;

	_eigh_ (M, W, d);	_revert_ (W, d);

	C = hub.C.leftCols (N) * W;
	//cout << "Schmidt decomposed C:\n" << C << "\n\n";

	_form_xform_mat_ ();
}

void SCHMIDT::_form_xform_mat_ ()
{
	int i = 0;
	T.setZero (K, 2 * Nimp);
	T.rightCols (Nimp) = C.leftCols (Nimp);

	// normalization
	for (i = 0; i < Nimp; i++)	T.col (i + Nimp) /= sqrt (1. - d(i));
	for (i = 0; i < Nimp; i++)
	{
		T.block (frag[i], 0, 1, 2 * Nimp).setZero ();
		T(frag[i], i) = 1.;		// set C^{F} to eye(Nimp)
								// see Knizia13JCTC, Bulik14PRB
	}
	TE = C.rightCols (N - Nimp);
}

MatrixXd SCHMIDT::_make_P_frag_ (HUBBARD& hub)
{
	return T.transpose () * hub.P * T;
}

void SCHMIDT::_print_ ()
{
	cout << "Eigenvalues:\t" << d.transpose () << "\n\n";
	cout << "T\n" << T << "\n\n";
	cout << "TE\n" << TE << "\n\n";
}

#endif
