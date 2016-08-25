#ifndef _SCHMIDT_H_INCLUDED_
#define _SCHMIDT_H_INCLUDED_

#include <Eigen/Dense>
#include "hubbard.h"
#include "hf.h"

using namespace Eigen;

class SCHMIDT
{
	public:
		int K, N, Nimp;
		int start_row;
		MatrixXd W;				// M = W * d * W^T
		VectorXd d;				// Singular values
		MatrixXd C;				// W-transformed C's
		MatrixXd T, TE;
		//SCHMIDT (char*, HUBBARD&);
		void _schmidt_ (HUBBARD&);
		void _form_xform_mat_ ();
};

void SCHMIDT::_schmidt_ (HUBBARD& hub)
{
	int i;
	start_row = 0;

	// prefix 'r' for raw (i.e. untransformed)
	MatrixXd rCF;
	rCF = hub.C.block (start_row, 0, Nimp, N);

	MatrixXd M;
	M = rCF.transpose () * rCF;

	_eigh_ (M, W, d);	_revert_ (W, d);
	cout << "M:\n" << M << "\n\n";
	cout << "W:\n" << W << "\n\n";
	cout << "d:\n" << d << "\n\n";

	C = hub.C.leftCols (N) * W;
	cout << "Schmidt decomposed C:\n" << C << "\n\n";
	for (i = 0; i < Nimp; i++)
	{
		C.block (start_row, i, Nimp, 1) /= sqrt(d (i));
		C.block (0, i, start_row, 1) /= sqrt (1. - d (i));
		C.block (start_row + Nimp, i, K - Nimp - start_row, 1) /= sqrt (1. - d (i));
	}

	_form_xform_mat_ ();
}

void SCHMIDT::_form_xform_mat_ ()
{
	T.setZero (K, 2 * Nimp);
	T.block (start_row, 0, Nimp, Nimp) = C.block (start_row, 0, Nimp, Nimp);
	T.block (0, Nimp, start_row, Nimp) = C.block (0, 0, start_row, Nimp);
	T.block (start_row + Nimp, Nimp, K - Nimp - start_row, Nimp) = C.block (start_row + Nimp, 0, K - Nimp - start_row, Nimp);
	TE = C.topRightCorner (K, N - Nimp);

	cout << "C:\n" << C << "\n\n";
	cout << "T:\n" << T << "\n\n";
	cout << "TE:\n" << TE << "\n\n";
}

#endif
