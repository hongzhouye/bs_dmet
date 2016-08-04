#ifndef _SCHMIDT_H_INCLUDED_
#define _SCHMIDT_H_INCLUDED_

#include <Eigen/Dense>
#include "hubbard.h"
#include "hf.h"

using namespace Eigen;

class SCHMIDT
{
	public:
		int K, Nup, Ndn, Nimp;
		MatrixXd Wup, Wdn;				// M = W * d * W^T
		VectorXd dup, ddn;				// Singular values
		MatrixXd Cup, Cdn;				// W-transformed C's
		MatrixXd Tup, Tdn, TEup, TEdn;
		//SCHMIDT (char*, HUBBARD&);
		void _schmidt_ (HUBBARD&);
		void _form_xform_mat_ ();
};

void SCHMIDT::_schmidt_ (HUBBARD& hub)
{
	// prefix 'r' for raw (i.e. untransformed)
	MatrixXd rCFup, rCFdn;		
	rCFup = hub.Cup.block (0, 0, Nimp, Nup);
	rCFdn = hub.Cdn.block (0, 0, Nimp, Ndn);

	MatrixXd Mup, Mdn;
	Mup = rCFup.transpose () * rCFup;
	Mdn = rCFdn.transpose () * rCFdn;

	_eigh_ (Mup, Wup, dup);	_revert_ (Wup, dup);
	_eigh_ (Mdn, Wdn, ddn); _revert_ (Wdn, ddn);

	Cup = hub.Cup.block (0, 0, K, Nup) * Wup;
	Cdn = hub.Cdn.block (0, 0, K, Ndn) * Wdn;

	_form_xform_mat_ ();
}

void SCHMIDT::_form_xform_mat_ ()
{
	Tup.setZero (K, 2 * Nimp);
	Tdn.setZero (K, 2 * Nimp);
	Tup.topLeftCorner (Nimp, Nimp) = Cup.topLeftCorner (Nimp, Nimp);
	Tdn.topLeftCorner (Nimp, Nimp) = Cdn.topLeftCorner (Nimp, Nimp);
	Tup.bottomRightCorner (K - Nimp, Nimp) = Cup.bottomLeftCorner (K - Nimp, Nimp);
	Tdn.bottomRightCorner (K - Nimp, Nimp) = Cdn.bottomLeftCorner (K - Nimp, Nimp);
	TEup = Cup.topRightCorner (K, Nup - Nimp);
	TEdn = Cdn.topRightCorner (K, Ndn - Nimp);

	cout << "Cdn:\n" << Cdn << "\n\n";
	cout << "Tdn:\n" << Tdn << "\n\n";
	cout << "TEdn:\n" << TEdn << "\n\n";
}

#endif
