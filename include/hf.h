#ifndef _HF_INCLUDED_
#define _HF_INCLUDED_

#include <Eigen/Dense>

using namespace Eigen;

#define SCF_ITER 1E3
#define SCF_CONV 1E-6
#define MAX_DIIS 10
#define PINV_TOL 1E-10
#define LINSOLVER_TOL 1E-4
#define mixing_beta_HL 0.2
#define ZERO 1E-14

//typedef Matrix<double, K, K> MatrixKd;
//typedef Matrix<double, K, 1> VectorKd;
#define esXd SelfAdjointEigenSolver<MatrixXd> 

// eigen solver for Hermitian matrix
void _eigh_ (const MatrixXd& A, MatrixXd& U, VectorXd& D)
{
	esXd es;
	es.compute (A);
	D = es.eigenvalues ();
	U = es.eigenvectors ();
}

// revert the order of eigen-vectors/values
// since Eigen by default gives them in an ascending order
// this will return them with a DESCENDING style
void _revert_ (MatrixXd& U, VectorXd& d)
{
	int i, n = d.size ();
	MatrixXd Up;	Up.setZero (n, n);
	VectorXd dp;	dp.setZero (n);
	for (i = 0; i < n; i++)
	{
		dp (i) = (d (n - 1 - i) < ZERO) ? (0.) : (d (n - 1 - i));
		Up.col (i) = U.col (n - 1 - i);
	}

	U = Up;	d = dp;
}

#endif
