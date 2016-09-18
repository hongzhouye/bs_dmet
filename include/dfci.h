#ifndef _DFCI_H_INCLUDED_
#define _DFCI_H_INCLUDED_

#include <Eigen/Dense>
#include <iostream>
#include <cstdio>
#include "hf.h"
#include "ab_string.h"

using namespace Eigen;
using namespace std;

#define MAX_DVDS_ITER 100
#define DVDS_CONV 1e-8
#define TROY_CONV 1e-13

// FCI ref: J. Olsen and B. O. Roos, J. Chem. Phys. 89 (4), 1988
// Davidson ref: E. Davidson, J. Comput. Phys. 17, 87 - 94 (1975)

class DFCI
{
	private:
		int lenh, lenV;
		bool *ifzero;
		// for dfci
		int _str_match_ (short *, short *, short *);
		void _str_gen_ (AB_STRING *, int, int, long int *);
		void _dfci_initguess_ (MatrixXd&);
		void _set_h_V_ (const MatrixXd&, double *);
		VectorXd _diagH_ ();
		template<typename Derived>
		VectorXd _Hx_ (const MatrixBase<Derived>&, int);
		void _GS_ (VectorXd&, const MatrixXd&);
		// for troy fci
		void _get_istr_ (const VectorXd&, VectorXi&, VectorXi&, int);
		void _convert_istr_ (const VectorXi&, vlis&, int);
		MatrixXd _get_H0_ (const VectorXd&, const vlis&, int);
	public:
		string mode;
		int K;			// basis set size
		int N;			// # of e^-
		long int tot;	// total # of configurations
		double E_fci;
		MatrixXd C_fci;	// expansion coefficients
		double *h;		// reduced one-electron integrals
		double *V;		// two electron integrals
		MatrixXd P;		// 1PDM
		double *G;		// 2PDM
		AB_STRING *astr;// only one is enough for Ms = 0
						// and spin-restricted case
						// see the FCI ref above
		void _init_ (int, int);
		void _dfci_ (const MatrixXd&, double *);
		void _troyfci_ (const MatrixXd&, double *);	// modified from Troy's code
		void _1PDM_ ();
		void _2PDM_ ();
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

void DFCI::_init_ (int Nbs, int Ne)
{
	int i, j, k, l, ij, kl, ijkl;

	K = Nbs;	N = Ne;
	lenh = K * (K + 1) / 2;
	lenV = lenh * (lenh + 1) / 2;
	h = _darray_gen_ (lenh);
	V = _darray_gen_ (lenV);
	//ifzero = new bool [lenV];

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

	// CHECK: interacting strings -- simplified version
	/*for (I = 0; I < tot; I++)
	{
		cout << I << ": [";
		for (i = 0; i < astr[I].itot - 1; i++)
			cout << astr[I].istr[i] << ", ";
		cout << astr[I].istr[i] << "]" << endl << endl;
	}*/
}

void DFCI::_set_h_V_ (const MatrixXd& hinp, double *Vinp)
// convert input h and V
{
	int i, j, k, ij, ik, kj;
	for (i = 0; i < lenV; i++)	V[i] = Vinp[i];
	for (i = 0; i < K; i++)
		for (j = 0; j <= i; j++)
		{
			ij = cpind(i,j);
			h[ij] = hinp(i, j);
			for (k = 0; k < K; k++)
			{
				ik = cpind(i,k);
				kj = cpind(k,j);
				h[ij] -= V[cpind(ik,kj)] / 2.;
			}
		}
}

// davidson
void DFCI::_dfci_initguess_ (MatrixXd& b)
{
	if (mode == "random")		b.setRandom ();
	else if (mode == "ones")	b.setOnes ();
	else if (mode == "major")	{b.setZero ();	b(0, 0) = 1.;}
	else if (mode == "read")	b = C_fci;
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

template<typename Derived>
VectorXd DFCI::_Hx_ (const MatrixBase<Derived>& b, int col)
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

void DFCI::_GS_ (VectorXd& b, const MatrixXd& A)
// Schmidt orthonormalize b to columns of A
{
	int n = A.cols ();
	double temp;
	for (int i = 0; i < n; i++)
	{
		const VectorXd& q = A.col (i);
		temp = q.transpose () * b;
		b = b - temp * q;
	}
	b.normalize ();
}

void DFCI::_dfci_ (const MatrixXd& hinp, double *Vinp)
{
	// setup h and V
	_set_h_V_ (hinp, Vinp);
	// initial setup
	MatrixXd b (tot * tot, 1);
	_dfci_initguess_ (b);
	MatrixXd s (tot * tot, 1);
	s.col(0) = _Hx_ (b, 0);
	MatrixXd At = b.transpose () * s;
	double lambda = At(0, 0);
	VectorXd alpha (1);	alpha (0) = 1.;
	VectorXd diagH = _diagH_ ();

	// Davidson iteration (see Davidson ref above)
	int iter = 1, i;
	VectorXd q (tot * tot), precond (tot * tot);
	double error, temp, offset, beta = 1.;
	//cout << "iter\terror" << endl;
	while (iter < MAX_DVDS_ITER)
	{
		C_fci = b * alpha;
		q = s * alpha - lambda * C_fci;
		error = q.norm ();
		if (error < DVDS_CONV)	break;
		else	//cout << iter << "\t" << error << "\n";

		// offset preconditioner to avoid singularity
		if (iter < 3)	offset = 1e-14;	else	offset = 0.;
		// Davidson preconditioner
		//precond = VectorXd::Constant (tot * tot, lambda + offset) - diagH;
		// Bendazzoli preconditioner
		precond = VectorXd::Constant (tot * tot, lambda) - diagH
			+ C_fci.cwiseProduct (2. * q - beta * C_fci);
		q = q.cwiseQuotient (precond);

		// Schmidt orthonormalization
		_GS_ (q, b);
		b.conservativeResize (tot * tot, iter + 1);	b.col(iter) = q;
		s.conservativeResize (tot * tot, iter + 1);	s.col(iter) = _Hx_ (b, iter);

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
		cout << "\nDavidson diagonalization failed to converge!\t"
			<< error << "\n\n";
		return;
	}
	//CHECK
	//cout << "guess mode = " << mode << "\titer = " << iter << endl << endl;
	//cout << "\nDesired accuracy is reached after " << iter << " iterations!\n\n";
	MatrixXd v = b * alpha;	C_fci = v;
	//cout << "Estimated error ||H v - lambda * v|| is " << (_Hx_ (v, 0) - lambda * v).norm () << "\n\n";
	//printf ("FCI energy is %18.16f\n\n", lambda);	// no need to print out this value

	iter ++;
}

void DFCI::_get_istr_ (const VectorXd& Hd, VectorXi& istr0,
		VectorXi& istr1, int N0)
/* get str for Olson's precond */
{
	long int TOT = tot * tot;
	istr0.setZero (N0);	istr1.setZero (TOT - N0);
    double maxval;
	int maxind = 0, i, temp;
	bool flag;
	// first assume the first N0 elements are lowest
    for (int i = 0; i < N0; i++)	istr0(i) = i;
    for (int i = N0; i < TOT; i++)	istr1(i - N0) = i;

	do
	{
        flag = false;
        maxval = -500.;
		// find maxval in the first N0 elements
		for (int i = 0; i < N0; i++)
			if (Hd(istr0(i)) > maxval)
			{
				maxval = Hd(istr0(i));	maxind = i;
			}
		// if some element is less than maxval, break
		for (i = N0; i < TOT; i++)
		{
			if (Hd(istr1(i - N0)) < maxval)
			{
                SWAP(temp,istr0(maxind),istr1(i - N0));
				flag = true;
			}
			if (flag)	break;
		}
	}
	while (flag);
}

void DFCI::_convert_istr_ (const VectorXi& istr0tmp, vlis& istr0, int N0)
{
	int rem;
	long int Ia, Ib;
	for (int i = 0; i < N0; i++)
	{
		Ia = istr0tmp(i) / tot;
		Ib = istr0tmp(i) % tot;
		long int *tmp = new long int [3];
		tmp[0] = Ia;	tmp[1] = Ib;	tmp[2] = istr0tmp(i);
		istr0.push_back (tmp);
	}
}

MatrixXd DFCI::_get_H0_ (const VectorXd& Hd, const vlis& istr0, int N0)
{
	MatrixXd H0;	H0.setZero (N0, N0);
	int i, j, k, ij, kl, ijkl, sgnij, sgnkl, cnnct, cnnct2;
	long int Ia, Ib, Ja, Jb, Ka, Kb;
	AB_STRING cstra, cstrb, cstr2a, cstr2b;

/*** off-diagonal elements ***/
// 1P contribution
	for (j = 0; j < N0; j++)
	{
		Ja = istr0[j][0];	Jb = istr0[j][1];
		cstra = astr[Ja];	cstrb = astr[Jb];
		for (i = 0; i < N0; i++)
		{
			if (i != j)
			{
				Ib = istr0[i][1];	Ia = istr0[i][0];
/*** alpha contribution ***/
// if Ib != Jb, alpha contribution is zero
				if (Ib == Jb)
				{
// check whether Ia is connected to Ja;
					cnnct = cstra._check_connect_ (Ia);
					if (cnnct > 0)
					{
						ij = cstra.cmpind[cnnct];
						sgnij = cstra.sgn[cnnct];
						H0(i, j) += h[ij] * sgnij;
					}
				}
/*** beta contribution ***/
// if Ia != Ja, beta contribution is zero
				if (Ia == Ja)
				{
// check whether Ib is connected to Jb;
					cnnct = cstrb._check_connect_ (Ib);
					if (cnnct > 0)
					{
						ij = cstrb.cmpind[cnnct];
						sgnij = cstrb.sgn[cnnct];
						H0(i, j) += h[ij] * sgnij;
					}
				}
			}
		}
	}

// 2P same spin
	for (j = 0; j < N0; j++)
	{
		Ja = istr0[j][0];	Jb = istr0[j][1];
		cstra = astr[Ja];	cstrb = astr[Jb];
		for (i = 0; i < N0; i++)
		{
			if (i != j)
			{
				Ia = istr0[i][0];	Ib = istr0[i][1];
				cstr2a = astr[Ia];	cstr2b = astr[Ib];
/*** alpha contribution ***/
				if (Ib == Jb)
				{
					for (Ka = 0; Ka < tot; Ka++)
						if (Ka == Ja)
						{
							cnnct2 = cstr2a._check_connect_ (Ka);
							if (cnnct2 < 0)	continue;
							for (int p = 0; p < N; p++)
							{
								int P = cstra.itot - p - 1;
								ij = cstra.cmpind[P];
								kl = cstr2a.cmpind[cnnct2];
								sgnkl = cstr2a.sgn[cnnct2];
								H0(i, j) += V[cpind(ij,kl)] * sgnkl / 2.;
							}
						}
						else if (Ka == Ia)
						{
							cnnct = cstra._check_connect_ (Ka);
							if (cnnct < 0)	continue;
							for (int p = 0; p < N; p++)
							{
								ij = cstra.cmpind[cnnct];
								sgnij = cstra.sgn[cnnct];
								int P = cstr2a.itot - p - 1;
								kl = cstr2a.cmpind[P];
								H0(i, j) += V[cpind(ij,kl)] * sgnij / 2.;
							}
						}
						else
						{
							cnnct = cstra._check_connect_ (Ka);
							if (cnnct < 0)	continue;
							cnnct2 = cstr2a._check_connect_ (Ka);
							if (cnnct2 < 0)	continue;
							ij = cstra.cmpind[cnnct];
							sgnij = cstra.sgn[cnnct];
							kl = cstr2a.cmpind[cnnct2];
							sgnkl = cstr2a.sgn[cnnct2];
							H0(i, j) += V[cpind(ij,kl)] * sgnij * sgnkl / 2.;
						}
				}
/*** beta contribution ***/
				if (Ia == Ja)
				{
					for (Kb = 0; Kb < tot; Kb++)
					{
						cnnct = cstrb._check_connect_ (Kb);
						if (cnnct < 0)	continue;
						cnnct2 = cstr2b._check_connect_ (Kb);
						if (cnnct2 < 0)	continue;
						ij = cstrb.cmpind[cnnct];
						sgnij = cstrb.sgn[cnnct];
						kl = cstr2b.cmpind[cnnct2];
						sgnkl = cstr2b.sgn[cnnct2];
						H0(i, j) += V[cpind(ij,kl)] * sgnij * sgnkl / 2.;
					}
				}
			}
		}
	}

// 2P opposite spin
	for (j = 0; j < N0; j++)
	{
		Ja = istr0[j][0];	Jb = istr0[j][1];
		cstra = astr[Ja];	cstrb = astr[Jb];
		for (i = 0; i < N0; i++)
			if (i != j)
			{
				Ia = istr0[i][0];	Ib = istr0[i][1];
				cnnct = cstra._check_connect_ (Ia);
				if (cnnct < 0)	continue;
				cnnct2 = cstrb._check_connect_ (Ib);
				if (cnnct2 < 0)	continue;
				kl = cstra.cmpind[cnnct];
				sgnkl = cstra.sgn[cnnct];
				ij = cstrb.cmpind[cnnct2];
				sgnij = cstrb.sgn[cnnct2];
				ijkl = cpind(ij,kl);
				H0(i, j) += V[ijkl] * sgnij * sgnkl;
			}
	}

/*** diagonal elements ***/
	for (i = 0; i < N0; i++)
		H0(i, i) = Hd(istr0[i][2]);

/*** SANITY CHECK: diagonal elements ***/
/*	VectorXd H0d;	H0d.setZero (N0);
	for (i = 0; i < N0; i++)
	{
// 1P
		Ia = istr0[i][0];	Ib = istr0[i][1];
		cstra = astr[Ia];	cstrb = astr[Ib];
		for (int pos = 0; pos < N; pos++)
		{
			int POS = cstra.itot - pos - 1;
			ij = cstra.cmpind[POS];	//alpha
			H0d(i) += h[ij];
			ij = cstrb.cmpind[POS];	//beta
			H0d(i) += h[ij];
		}
// 2P same spin
		for (int p = 0; p < cstra.itot - N; p++)
		{
			ij = cstra.cmpind[p];
			if (ifzero[cpind(ij,ij)])	H0d(i) += V[cpind(ij,ij)] / 2.;
			ij = cstrb.cmpind[p];
			if (ifzero[cpind(ij,ij)])	H0d(i) += V[cpind(ij,ij)] / 2.;
		}
		for (int p1 = 0; p1 < N; p1++)
		{
			int P1 = cstra.itot - p1 - 1;
			ij = cstra.cmpind[P1];
			for (int p2 = 0; p2 < N; p2++)
			{
				int P2 = cstra.itot - p2 - 1;
				kl = cstra.cmpind[P2];
				if (ifzero[cpind(ij,kl)])	H0d(i) += V[cpind(ij,kl)] / 2.;
			}
			ij = cstrb.cmpind[P1];
			for (int p2 = 0; p2 < N; p2++)
			{
				int P2 = cstrb.itot - p2 - 1;
				kl = cstrb.cmpind[P2];
				if (ifzero[cpind(ij,kl)])	H0d(i) += V[cpind(ij,kl)] / 2.;
			}
		}
// 2P opposite spin
		for (int p1 = 0; p1 < N; p1++)
		{
			int P1 = cstra.itot - p1 - 1;
			ij = cstra.cmpind[P1];
			for (int p2 = 0; p2 < N; p2++)
			{
				int P2 = cstrb.itot - p2 - 1;
				kl = cstrb.cmpind[P2];
				if (ifzero[cpind(ij,kl)])	H0d(i) += V[cpind(ij,kl)];
			}
		}
	}
	cout << "CHECK:: H0d:\n" << H0d << "\n\n";
*/

	return H0;
}

void DFCI::_troyfci_ (const MatrixXd& hinp, double *Vinp)
{
// set up h and V
	_set_h_V_ (hinp, Vinp);

	const int N0 = 10;	// dim of precond
	VectorXd Hd = _diagH_ ();

// strings in/out precond
	//VectorXi istr0, istr1;
	//_get_istr_ (Hd, istr0, istr1, N0);

// convert precond strings in (Ia, Ib, Ia * tot + Ib) form
	//vlis istr0tmp;	_convert_istr_ (istr0, istr0tmp, N0);

// get H0
	//MatrixXd H0 = _get_H0_ (Hd, istr0tmp, N0);

// diagonalize H0
	//MatrixXd U0;	VectorXd D0;
	//_eigh_ (H0, U0, D0);

// initial guess
	VectorXd b (tot * tot);	b.setZero ();
	//if (mode == "major")	b(0) = 1.;
	//else if (mode == "read")	b = C_fci;
	b(0) = 1.;

	VectorXd s = _Hx_ (b, 0);
	double lambda = b.transpose () * s;
	VectorXd q (tot * tot);
	VectorXd prcd = q, bi = q, qi = q;

// fci iter
	int iter = 0;
	double error = 10., fac, beta = 1.;
	Matrix2d Hm, Um;	Vector2d Dm;
	do
	{
		iter ++;
// Davidson update
		q = s - lambda * b;
// Davidson's preconditioner
		fac = (iter < 3) ? (1.E-10) : (0);
		prcd = VectorXd::Constant (tot * tot, fac + lambda) - Hd;
// Bendazzoli's preconditioner
		//prcd = VectorXd::Constant (tot * tot, lambda) - Hd
		//	+ b.cwiseProduct (2. * q - beta * b);
		qi = q.cwiseQuotient (prcd);
// correct if qi is parallel to q
		q.normalize ();
		fac = 0.1;
		if (fabs(qi.transpose () * q) < fac)
		{
			qi += 2. * fac * q;
			fac = qi.transpose () * b;
	        qi -= fac * b;   qi.normalize ();
			fac = qi.transpose () * b;
	        qi -= fac * b;   qi.normalize ();
		}
		q = _Hx_ (qi, 0);
// construct Hm
		Hm(0, 0) = b.transpose () * s;
		Hm(1, 1) = qi.transpose () * q;
		Hm(0, 1) = b.transpose () * q;
		Hm(1, 0) = Hm(0, 1);
// diagonalize Hm
		SelfAdjointEigenSolver<Matrix2d> es2;
		es2.compute (Hm);
		Um = es2.eigenvectors ();
		Dm = es2.eigenvalues ();
// Olsen update
		if (iter < 50 || iter % 10)
		{
			error = Um(1, 0);
			b = Um(0, 0) * b + Um(1, 0) * qi;
			s = Um(0, 0) * s + Um(1, 0) * q;
		}
		else
		{
			error = Um(0, 1) / 2.;
			fac = 1. / sqrt (1. + error * error);
			b = fac * (b + error * qi);
			s = fac * (s + error * q);
		}
		fac = b.norm ();	b /= fac;	s /= fac;
// avoid accumulated error
		if (iter % 4 == 0)	s = _Hx_ (b, 0);
// update energy
		lambda = b.transpose () * s;
		lambda /= b.squaredNorm ();
		printf ("%4d\t%10.7e\n", iter, error);
	}
	while (fabs(error) > TROY_CONV && iter < MAX_DVDS_ITER);
	cout << "estimated error: " << (s - lambda * b).norm () << endl;
	if (fabs(error) <= TROY_CONV)	C_fci = b;
	else	cout << "TroyFCI fails to converge.\n\n";
}

// post-process
void DFCI::_1PDM_ ()
{
	P.setZero (K, K);

	AB_STRING cstr;

	int mu, nu, pos, mn, mnp;
	long int Ia, Ja, Ib, iIa, iJa, iIb;
	double Ca, Cb;
	for (Ja = 0; Ja < tot; Ja++)
	{
		iJa = Ja * tot;
		cstr = astr[Ja];
		for (pos = 0; pos < cstr.itot; pos++)
		{
			Ia = cstr.istr[pos];
			iIa = Ia * tot;
			Ca = Cb = 0;
			for (Ib = 0; Ib < tot; Ib++)
			{
				Ca += C_fci(iIa + Ib) * C_fci(iJa + Ib);
				//iIb = Ib * tot;
				//Cb += C_fci(iIb + Ia) * C_fci(iIb + Ja);
			}
			mn = cstr.cmpind[pos];
			for (mu = 0; mu < K; mu++)
				for (nu = 0; nu < K; nu++)
				{
					mnp = cpind(mu, nu);
					if (mn == mnp)
					{
						P(mu, nu) += (Ca
							//+ Cb
						) * cstr.sgn[pos];
					}
				}
		}
	}
	P.diagonal () *= 2.;	P /= 2.;
}

void DFCI::_2PDM_ ()
{
	int lenh = K * (K + 1) / 2;
	int lenV = lenh * (lenh + 1) / 2;
	G = _darray_gen_ (lenV);

	AB_STRING cstr, cstr2;
	int mu, nu, la, si, pos, pos2;
	int ns, ml, ms, mn, ls;		// compound indices
	long int Ia, Ib, Ja, Jb, Ka;

	// Gamma^{alpha,alpha}_1, Gamma^{beta,beta}_1
	double Ca = 0, Cb = 0;
	long int iIa, iJa, iIb;
	for (Ja = 0; Ja < tot; Ja++)
	{
		iJa = Ja * tot;
		cstr = astr[Ja];
		for (pos = 0; pos < cstr.itot; pos++)
		{
			Ka = cstr.istr[pos];
			ns = cstr.cmpind[pos];
			cstr2 = astr[Ka];
			for (pos2 = 0; pos2 < cstr2.itot; pos2++)
			{
				Ia = cstr2.istr[pos2];
				iIa = Ia * tot;
				Ca = Cb = 0;
				for (Ib = 0; Ib < tot; Ib++)
				{
					Ca += C_fci(iIa + Ib) * C_fci(iJa + Ib);
					//iIb = Ib * tot;
					//Cb += C_fci(iIb + Ia) * C_fci(iIb + Ja);
				}
				ml = cstr2.cmpind[pos2];
				G[cpind(ml,ns)] += cstr.sgn[pos] * cstr2.sgn[pos2]
					* (Ca
				//	+ Cb
					)
				//	/ 2.
				;
			}
		}
	}

	// Gamma^{alpha,beta}_1
	for (Ja = 0; Ja < tot; Ja++)
	{
		cstr = astr[Ja];
		for (pos = 0; pos < cstr.itot; pos++)
		{
			Ia = cstr.istr[pos];
			ml = cstr.cmpind[pos];
			for (Jb = 0; Jb < tot; Jb++)
			{
				cstr2 = astr[Jb];
				for (pos2 = 0; pos2 < cstr2.itot; pos2++)
				{
					Ib = cstr2.istr[pos2];
					ns = cstr2.cmpind[pos2];
					G[cpind(ml,ns)] += cstr.sgn[pos] * cstr2.sgn[pos2]
						* C_fci(Ia * tot + Ib) * C_fci(Ja * tot + Jb);
				}
			}
		}
	}

	// Gamma^{alpha/beta}_2
	for (Ja = 0; Ja < tot; Ja++)
	{
		cstr = astr[Ja];
		iJa = Ja * tot;
		for (pos = 0; pos < cstr.itot; pos++)
		{
			Ia = cstr.istr[pos];
			iIa = Ia * tot;
			ms = cstr.cmpind[pos];
			Ca = Cb = 0;
			for (Ib = 0; Ib < tot; Ib++)
			{
				Ca += C_fci(iIa + Ib) * C_fci(iJa + Ib);
				//Cb += C_fci(Ib * tot + Ia) * C_fci(Ib * tot + Ja);
			}
			for (mu = 0; mu < K; mu++)
				for (si = 0; si < K; si++)
					if (cpind(mu,si) == ms)
						for (nu = 0; nu < K; nu++)
						{
							mn = cpind(mu,nu);	ns = cpind(nu,si);
							G[cpind(mn,ns)] -= cstr.sgn[pos] * Ca;
							//G[cpind(mn,ns)] -= cstr.sgn[pos] * (Ca + Cb) / 2.;
						}
		}
	}

	// CHECK: sum rule
	/*double sum = 0;
	for (mu = 0; mu < K; mu++)
	{
		mn = cpind(mu,mu);
		for (la = 0; la < K; la++)
		{
			ls = cpind(la,la);
			sum += G[cpind(mn,ls)];
		}
	}
	cout << "sum = " << sum << "\n\n";
	*/
}
#endif
