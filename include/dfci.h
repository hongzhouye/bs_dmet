#ifndef _DFCI_H_INCLUDED_
#define _DFCI_H_INCLUDED_

#include <Eigen/Dense>
#include <iostream>
#include <cstdio>
#include "hf.h"
#include "hred.h"
#include "ab_string.h"

using namespace Eigen;
using namespace std;

#define cpind(i,j) (i>j)?(ioff[i]+j):(ioff[j]+i)

class DFCI
{
	private:
		int *ioff;		// lookup table for compound indices
		int _str_match_ (short *, short *, short *);
		void _str_gen_ (AB_STRING *, int, int, long int *);
	public:
		int K;			// basis set size
		int N;			// # of e^-
		long int tot;	// total # of configurations
		double *h;		// reduced one-electron integrals
		double *V;		// two electron integrals
		MatrixXd G;		// Davidson matrix
		AB_STRING *astr;// only one is enough for Ms = 0
						// and spin-restricted case
		void _init_ (HRED&);
		void _dfci_ ();
};

int DFCI::_str_match_ (short *a, short *b, short *sgn)
// match two strings, if # of diff == 1, return position
{
	*sgn = 1;
	int i, j, count = 0, flag;
	short temp;
	for (i = 0; i < N; i++)
	{
		flag = 1;
		if (a[i] == b[i])	continue;
		else
			for (j = 0; j < N; j++)
				if (a[i] == b[j])
				{
					// swap b_i and b_j
					temp = b[i];	b[i] = b[j];	b[j] = temp;
					// and change sign of the b string
					*sgn = - *sgn;
					flag = 0;
					break;
				}
		if (flag)	count ++;
		if (count >= 2)	break;
	}
	if (count == 1)
	{
		for (i = 0; i < N; i++)
			if (a[i] != b[i])	return i;
	}
	else	return -1;
}

void DFCI::_init_ (HRED& hr)
{
	int i, j, k, l, ij, kl, ijkl;

	K = hr.Ni;	N = K / 2;
	int lenh = K * (K + 1) / 2;
	int lenV = lenh * (lenh + 1) / 2;
	h = _darray_gen_ (lenh);
	V = _darray_gen_ (lenV);

	// set up the lookup table 'ioff'
	ioff = new int [lenh + 1];
	if (ioff == NULL)
	{
		cout << "failed to malloc memory for array ioff!\n";
		exit (1);
	}
	ioff[0] = 0;
	for (i = 1; i < lenh + 1; i++)	ioff[i] = ioff[i - 1] + i;

	// set up h and V
	for (i = 0; i < K; i++)
		for (j = 0; j <= i; j++)
		{
			ij = cpind(i,j);
			h[ij] = hr.h (i, j);
			for (k = 0; k < K; k++) h[ij] -= V[index4(i,k,k,j,K)] / 2.;
			for (k = 0; k < K; k++)
				for (l = 0; l <= k; l++)
				{
					kl = cpind(k,l);
					ijkl = cpind(ij,kl);
					// physicists' notation to chemists' notation!
					V[ijkl] = hr.V[index4(i,k,j,l,K)];
				}
		}

	// malloc memory for 'astr'
	tot = _nchoosek_ (K, N);
	astr = new AB_STRING [tot];
	if (astr == NULL)
	{
		cout << "failed to malloc memory for ab strings!\n";
		exit (1);
	}

	// initialize 'astr'
	long int I, J, count = 0;
	for (I = 0; I < tot; I++)	astr[I]._init_ (K, N);
	_str_gen_ (astr, 0, N, &count);
	
	for (I = 0; I < tot; I++)
	{
		cout << "#" << I << ":\t";
		for (i = 0; i < N; i++)
			cout << astr[I].str[i] << "  ";
		cout << "\n";
	}	

	// set up info of interacting strings for each 'astr'
	short b[N], sgn;
	int pos;
	for (I = 0; I < tot; I++)
	{
		count = 0;
		for (J = 0; J < tot; J++)
			if (J == I)	continue;
			else
			{
				for (k = 0; k < N; k++)	b[k] = astr[J].str[k];
				pos = _str_match_ (astr[I].str, b, &sgn);

				if (pos == -1)	continue;
				else
				{
					i = b[pos], j = astr[I].str[pos];
					astr[I].cmpind[count] = cpind(i,j);
					astr[I].istr[count] = J;
					astr[I].sgn[count] = sgn;
					count ++;
				}
			}
	}

	// print out check
	/*for (I = 0; I < tot; I++)
	{
		for (i = 0; i < N; i++)	cout << astr[I].str[i] << "  ";
		cout << "\n";
		cout << "( istr, sgn, cmpind )\n";
		for (i = 0; i < astr[I].itot; i++)
			printf ("( %ld, %d, %d )\n", astr[I].istr[i], 
					astr[I].sgn[i], astr[I].cmpind[i]);
		cout << "\n\n";
	}*/
}

void DFCI::_str_gen_ (AB_STRING *as, int init, int n, long int *cnt)
// recursively generate N-electron strings from K orbitals
{
	int i, j;
	if (n == 0)
	{
		*cnt = *cnt + 1;
		return;
	}
	else
		for (i = init; i <= K - n; i++)
		{
			as[*cnt].str[N - n] = i;
			_str_gen_ (as, i + 1, n - 1, cnt);
			if (i != K - n && *cnt < tot)
				for (j = 0; j < N - n; j++)	as[*cnt].str[j] = as[*cnt - 1].str[j];
		}
}

void DFCI::_dfci_ ()
{
	//double diagH = _darray_gen_ (tot);
}

#endif
