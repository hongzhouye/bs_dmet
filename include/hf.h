#ifndef _HF_INCLUDED_
#define _HF_INCLUDED_

#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <vector>
#include <map>

using namespace Eigen;
using namespace std;

#define SCF_ITER 1E3
#define SCF_CONV 1E-6
#define mixing_beta_HL 0.2
#define ZERO 1E-14

#define esXd SelfAdjointEigenSolver<MatrixXd>
#define index4(i,j,k,l,K) i*K*K*K+j*K*K+k*K+l
#define cpind(i,j) (i>j)?(ioff[i]+j):(ioff[j]+i)

typedef map<string, int> fmap;
typedef vector<string> vs;
typedef vector<int> vi;
typedef vector<vector<int> > vvi;

int *ioff;		// lookup table for compound indices

void _gen_ioff_ (int len)
// set up the lookup table 'ioff'
{
	int i;
	ioff = new int [len + 1];
	if (ioff == NULL)
	{
		cout << "failed to malloc memory for array ioff!\n";
		exit (1);
	}
	ioff[0] = 0;
	for (i = 1; i < len + 1; i++)	ioff[i] = ioff[i - 1] + i;
}

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

// generate int type array
int *_iarray_gen_ (int size)
{
	int i;
	int * p = new int[size];
	if (p == NULL)
	{
		cout << "Allocate memory error!\n";
		exit (1);
	}
	for (i = 0; i < size; i++)	p[i] = 0;
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
