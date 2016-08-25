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

	// prefix 'r' for raw (i.e. untransformed)
	MatrixXd rCF;		
	rCF = hub.C.block (0, 0, Nimp, N);

	MatrixXd M;
	M = rCF.transpose () * rCF;

	_eigh_ (M, W, d);	_revert_ (W, d);
	cout << "M:\n" << M << "\n\n";
	cout << "W:\n" << W << "\n\n";
	cout << "d:\n" << d << "\n\n";

	C = hub.C.block (0, 0, K, N) * W;

	for (i = 0; i < Nimp; i++)
	{
		C.block (0, i, Nimp, 1) /= sqrt(d (i));
		C.block (Nimp, i, K - Nimp, 1) /= sqrt (1. - d (i));
	}

	_form_xform_mat_ ();
}

void SCHMIDT::_form_xform_mat_ ()
{
	T.setZero (K, 2 * Nimp);
	T.topLeftCorner (Nimp, Nimp) = C.topLeftCorner (Nimp, Nimp);
	T.bottomRightCorner (K - Nimp, Nimp) = C.bottomLeftCorner (K - Nimp, Nimp);
	TE = C.topRightCorner (K, N - Nimp);

	cout << "C:\n" << C << "\n\n";
	cout << "T:\n" << T << "\n\n";
	cout << "TE:\n" << TE << "\n\n";
}

#endif
