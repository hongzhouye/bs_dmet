#ifndef _HF_INCLUDED_
#define _HF_INCLUDED_

#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

#define SCF_ITER 1E3
#define SCF_CONV 1E-6
#define MAX_DIIS 10
#define PINV_TOL 1E-10
#define LINSOLVER_TOL 1E-4
#define mixing_beta_HL 0.2
#define ZERO 1E-14

#define esXd SelfAdjointEigenSolver<MatrixXd> 
#define index4(i,j,k,l,K) i*K*K*K+j*K*K+k*K+l

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

// generate double type array
double * _darray_gen_ (int size)
{
	int i;
	double * p = new double[size];
	if (p == NULL)
	{
		cout << "Allocate memory error!\n";
		exit (1);
	}
	for (i = 0; i < size; i++)	p[i] = 0.;
	return p;
}

// n choose k function
long int _nchoosek_ (int n, int k)
{
	int i;
	double prod = 1.;
	if (k == 0)	return 1;
	for (i = 1; i <= k; i++)
		prod *= (double) (n - k + i) / i;
	return (long int) prod;
}

#endif
