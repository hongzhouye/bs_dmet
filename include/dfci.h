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
#define MAX_DVDS_ITER 100
#define DVDS_CONV 1e-8

// FCI ref: J. Olsen and B. O. Roos, J. Chem. Phys. 89 (4), 1988
// Davidson ref: E. Davidson, J. Comput. Phys. 17, 87 - 94 (1975)

class DFCI
{
	private:
		int *ioff;		// lookup table for compound indices
		int _str_match_ (short *, short *, short *);
		void _str_gen_ (AB_STRING *, int, int, long int *);
		void _dfci_initguess_ (MatrixXd&, string);
		VectorXd _diagH_ ();
		VectorXd _Hx_ (MatrixXd&, int);
	public:
		int K;			// basis set size
		int N;			// # of e^-
		long int tot;	// total # of configurations
		double *h;		// reduced one-electron integrals
		double *V;		// two electron integrals
		MatrixXd G;		// Davidson matrix
		AB_STRING *astr;// only one is enough for Ms = 0
						// and spin-restricted case
						// see the FCI ref above
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
	int max = 0;
	for (i = 0; i < K; i++)
		for (j = 0; j <= i; j++)
		{
			ij = cpind(i,j);
			h[ij] = hr.h (i, j);
			for (k = 0; k < K; k++)
			{
				h[ij] -= hr.V[index4(i,k,k,j,K)] / 2.;
				for (l = 0; l <= k; l++)
				{
					kl = cpind(k,l);
					if (kl > ij)	continue;
					ijkl = cpind(ij,kl);
					// physicists' notation to chemists' notation!
					V[ijkl] = hr.V[index4(i,k,j,l,K)];
				}
			printf("%d;%d;%10.7f\n", i, j, h[ij]);
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

	// CHECK: strings generated
	/*for (I = 0; I < tot; I++)
	{
		cout << "#" << I << ":\t";
		for (i = 0; i < N; i++)
			cout << astr[I].str[i] << "  ";
		cout << "\n";
	}*/

	// set up info of interacting strings for each 'astr'
	short b[N], sgn;
	int pos;
	for (I = 0; I < tot; I++)
	{
		count = 0;
		for (J = 0; J < tot; J++)
			if (J == I)
				for (j = 0; j < N; j++)
				{
					k = astr[I].itot - j - 1;
					i = astr[I].str[j];
					astr[I].cmpind[k] = cpind(i,i);
					astr[I].istr[k] = I;
					astr[I].sgn[k] = 1;
				}
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

	// CHECK: interacting strings
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

void DFCI::_dfci_initguess_ (MatrixXd& b, string mode)
{
	if (mode == "random")		b.setRandom ();
	else if (mode == "ones")	b.setOnes ();
	else if (mode == "major")	{b.setZero ();	b(0, 0) = 1.;}
	b.col(0).normalize ();
}

VectorXd DFCI::_diagH_ ()
// return H(I, J)
{
	int Ia, Ib, i, k, pos, ii, kk, ij;
	MatrixXd dgH (tot, tot);	dgH.setZero ();
	VectorXd F (tot);
	AB_STRING cstr, cstr2;

	// one-body term
	F.setZero ();
	for (Ib = 0; Ib < tot; Ib++)
	{
		cstr = astr[Ib];
		// loop over last N interacting strings
		for (i = 0; i < N; i++)
		{
			pos = cstr.itot - i - 1;
			ii = cstr.cmpind[pos];
			F(Ib) += h[ii];
			for (k = 0; k < N; k++)
			{
				pos = cstr.itot - k - 1;
				kk = cstr.cmpind[pos];
				F(Ib) += V[cpind(ii,kk)] / 2.;
			}
		}
		for (i = 0; i < cstr.itot - N; i++)
		{
			ij = cstr.cmpind[i];
			F(Ib) += V[cpind(ij,ij)] / 2.;
		}
	}
	for (Ia = 0; Ia < tot; Ia++)
	{
		dgH.col(Ia) += F;
		dgH.row(Ia) += F.transpose ();
	}

	// two-body cross term
	for (Ia = 0; Ia < tot; Ia++)
	{
		cstr = astr[Ia];
		for (Ib = 0; Ib < tot; Ib++)
		{
			cstr2 = astr[Ib];
			// loop over last N interacting strings for astr[Ia]
			for (i = 0; i < N; i++)
			{
				pos = cstr.itot - i - 1;
				ii = cstr.cmpind[pos];
				// loop over last N interacting strings for astr[Ib]
				for (k = 0; k < N; k++)
				{
					pos = cstr2.itot - k - 1;
					kk = cstr2.cmpind[pos];
					dgH(Ia, Ib) += V[cpind(ii,kk)];
				}
			}
		}
	}

	// convert matrix dgH to vector diagH
	VectorXd diagH (tot * tot);
	for (Ia = 0; Ia < tot; Ia++)
		for (Ib = 0; Ib < tot; Ib++)
			diagH(Ia * tot + Ib) = dgH(Ia, Ib);

	return diagH;
}

VectorXd DFCI::_Hx_ (MatrixXd& b, int col)
// see FCI ref above
{
	int Ia, Ib, Ja, Jb, Kb;
	int ij, kl;

	// convert vector b into matrix C
	MatrixXd C (tot, tot);
	for (Ia = 0; Ia < tot; Ia++)
		for (Ib = 0; Ib < tot; Ib++)
			C(Ia, Ib) = b(Ia * tot + Ib, col);

	// forming sigma_1
	VectorXd F (tot);
	MatrixXd sigma_1 (tot, tot), sigma_3 (tot, tot);
	sigma_1.setZero ();
	AB_STRING cstr, cstr2;		// current string in the 1st and 2nd loops
	for (Ib = 0; Ib < tot; Ib++)
	{
		F.setZero ();	cstr = astr[Ib];
		// loop over interacting strings of cstr
		for (Kb = 0; Kb < cstr.itot; Kb++)
		{
			kl = cstr.cmpind[Kb];
			F(cstr.istr[Kb]) += cstr.sgn[Kb] * h[kl];
			cstr2 = astr[cstr.istr[Kb]];
			// loop over interacting strings of cstr2
			for (Jb = 0; Jb < cstr2.itot; Jb++)
			{
				ij = cstr2.cmpind[Jb];
				F(cstr2.istr[Jb]) += cstr.sgn[Kb] * cstr2.sgn[Jb] * V[cpind(ij,kl)] / 2.;
			}
		}
		for (Ia = 0; Ia < tot; Ia++)	
		{
			sigma_1(Ia, Ib) += C.row(Ia) * F;
			sigma_1(Ib, Ia) += F.transpose () * C.col(Ia);
		}
	}

	// forming sigma_3
	MatrixXd G (tot, tot);
	for (Ia = 0; Ia < tot; Ia++)
	{
		cstr = astr[Ia];
		for (Ib = 0; Ib < tot; Ib++)
		{
			G.setZero ();
			cstr2 = astr[Ib];
			for (Ja = 0; Ja < cstr.itot; Ja++)
			{
				ij = cstr.cmpind[Ja];
				for (Jb = 0; Jb < cstr2.itot; Jb++)
				{
					kl = cstr2.cmpind[Jb];
					G(cstr.istr[Ja], cstr2.istr[Jb]) += cstr.sgn[Ja] * cstr2.sgn[Jb] * V[cpind(ij,kl)];
				}
			}
			sigma_3(Ia, Ib) = (G.cwiseProduct (C)).sum ();
		}
	}

	// convert matrix sigma to vector s
	MatrixXd sigma = sigma_1 + sigma_3;
	VectorXd s (tot * tot);
	for (Ia = 0; Ia < tot; Ia++)
		for (Ib = 0; Ib < tot; Ib++)
			s(Ia * tot + Ib) = sigma(Ia, Ib);

	return s;
}

void DFCI::_dfci_ ()
{
	// initial setup
	MatrixXd b (tot * tot, 1);
	_dfci_initguess_ (b, "major");
	MatrixXd s (tot * tot, 1);
	s.col(0) = _Hx_ (b, 0);
	MatrixXd At = b.transpose () * s;
	double lambda = At(0, 0);
	VectorXd alpha (1);	alpha (0) = 1.;
	VectorXd diagH = _diagH_ ();
	cout << "diagH:\n" << diagH << endl << endl;

	// Davidson iteration (see Davidson ref above)
	int iter = 1, i;
	VectorXd q (tot * tot);
	double temp;
	cout << "iter\terror" << endl;
	while (iter < MAX_DVDS_ITER)
	{
		q = (s - lambda * b) * alpha;
		if (q.norm () < DVDS_CONV)	break;
		else	cout << iter << "\t" << q.norm () << "\n";

		q = q.cwiseQuotient (VectorXd::Constant (tot * tot, lambda) - diagH);

		// Schmidt orthonormalization
		for (i = 0; i < iter; i++)
		{
			temp = b.col(i).transpose () * q;
			q = q - temp * b.col(i);
		}
		q.normalize ();
		b.conservativeResize (tot * tot, iter + 1);	b.col(iter) = q;
		s.conservativeResize (tot * tot, iter + 1);	s.col(iter) = _Hx_ (b, iter);

		//check
		//cout << "( " << b.rows() << ", " << b.cols() << " )" << "\n\n";

		// form tilde{A}
		At.conservativeResize (iter + 1, iter + 1);
		for (i = 0; i < iter; i++)	At(i, iter) = At(iter, i) = b.col(i).transpose () * s.col(iter);
		At(iter, iter) = q.transpose () * s.col(iter);

		// diagonalize tilde{A}
		MatrixXd U (iter + 1, iter + 1);
		VectorXd D (iter + 1);
		_eigh_ (At, U, D);
		alpha = U.col(0);
		lambda = D(0);

		// increase iter by 1
		iter++;
	}
	if (iter >= MAX_DVDS_ITER)	// if not converged
	{
		cout << "\nDavidson diagonalization failed to converge!\n";
		return;
	}

	cout << "\nDesired accuracy is reached after " << iter << "iterations!\n\n";
	MatrixXd v = b * alpha;
	cout << "Estimated error ||H v - lambda * v|| is " << (_Hx_ (v, 0) - lambda * v).norm () << "\n\n";
	printf ("FCI energy is %18.16f\n\n", lambda);

	iter ++;
}

#endif
