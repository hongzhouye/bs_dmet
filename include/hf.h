#ifndef _HF_INCLUDED_
#define _HF_INCLUDED_

#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <vector>
#include <valarray>
#include <map>
#include <chrono>
//#include "easylogging++.h"

using namespace Eigen;
using namespace std;

#define SCF_ITER 1E3
#define SCF_CONV 1E-6
#define mixing_beta_HL 0.2
#define ZERO 1E-14

#define esXd SelfAdjointEigenSolver<MatrixXd>
#define index4(i,j,k,l,K) i*K*K*K+j*K*K+k*K+l
#define cpind(i,j) (i>j)?(ioff[i]+j):(ioff[j]+i)
#define SWAP(a,b,c) (a)=(b);(b)=(c);(c)=(a);

IOFormat HeavyFmt(12, 0, "\t", ";\n");
IOFormat Short(8, 0, ", ", "\n", "[", "]");

typedef map<string, int> msi;
typedef vector<string> vs;
typedef vector<int> vi;
typedef vector<vector<int> > vvi;
typedef vector<vector<int *> > vvis;
typedef vector<int *> vis;
typedef vector<long int*> vlis;
typedef vector<MatrixXd> vMatrixXd;

//typedef valarray<int> vi;
typedef valarray<double> vd;

typedef vector<int> iv1;
typedef vector<iv1> iv2;
typedef vector<iv2> iv3;
typedef vector<double> dv1;
typedef vector<dv1> dv2;
typedef vector<dv2> dv3;

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

// eigen solver for 2-by-2 Hermitian matrix
void _eigh2_ (const Matrix2d& A, Matrix2d& U, Vector2d& D)
{
	SelfAdjointEigenSolver<Matrix2d> es;
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

// generate int type array, with input value
int *_iarray_gen_ (int size, int init)
{
	int i;
	int * p = new int[size];
	if (p == NULL)
	{
		cout << "Allocate memory error!\n";
		exit (1);
	}
	for (i = 0; i < size; i++)	p[i] = init;
	return p;
}

// generate a 2D int type array
int **_iarray2_gen_ (int N1, int N2)
{
	int i, j;
	int **p = new int* [N1];
	for (i = 0; i < N1; i++)
	{
		p[i] = new int [N2];
		for (j = 0; j < N2; j++)	p[i][j] = 0;
	}
	return p;
}

// generate a 3D int type array
int ***_iarray3_gen_ (int N1, int N2, int N3)
{
	int i, j, k;
	int ***p = new int** [N1];
	for (i = 0; i < N1; i++)
	{
		p[i] = new int* [N2];
		for (j = 0; j < N2; j++)
		{
			p[i][j] = new int [N3];
			for (k = 0; k < N3; k++)	p[i][j][k] = 0;
		}
	}
	return p;
}

// copy array 1 to array 2 (size N given)
template <typename T>
void _copy_array_ (T *a1, T *a2, int N)
{
	for (int i = 0; i < N; i++)	a2[i] = a1[i];
}

// copy array 1 to array 2 (3D)
template <typename T>
void _copy_array3_ (T ***a1, T ***a2, int N1, int N2, int N3)
{
	int i,j,k;
	for (i = 0; i < N1; i++)
        for (j = 0; j < N2; j++)
            for (k = 0; k < N3; k++)
    			a2[i][j][k] = a1[i][j][k];
}

// take absolute values of a 3D array
template <typename T>
void _abs_array3_ (T ***a, int N1, int N2, int N3)
{
    int i,j,k;
    for (i = 0; i < N1; i++)
        for (j = 0; j < N2; j++)
            for (k = 0; k < N3; k++)
                if (a[i][j][k] < T(0))  a[i][j][k] = -a[i][j][k];
}

// n choose k function
long int _nchoosek_ (int n, int k)
{
	int i;
	double prod = 1.;
	if (k == 0)	return 1;
	else if (n < 0 || k < 0)	return 0;
	for (i = 1; i <= k; i++)
		prod *= (double) (n - k + i) / i;
	return (long int) prod;
}

#endif
