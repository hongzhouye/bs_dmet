#ifndef _FCI_H_INCLUDED_
#define _FCI_H_INCLUDED_

#include <Eigen/Dense>
#include <iostream>
#include "hred.h"
#include "hubbard.h"
#include "ab_string.h"

using namespace Eigen;
using namespace std;

class FCI
{
	public:
		int K;				// basis set size
		int N;				// # of e^-
		long int TOT;
		double *V;
		MatrixXd h;
		AB_STRING *astr, *bstr;
		FCI (HRED&);
		void _fci_ ();
	private:
		int _match_str_ (short *, short *, short *, int, int *);
};

int FCI::_match_str_ (short * a, short * b, short * sign_b, int Ne, int * pos)
// put string b in such a form that maximally overlap with string a
{
	int i, j, count = 0, flag;
	short temp;
	for (i = 0; i < Ne; i++)
	{
		flag = 1;
		if (a[i] == b[i])	continue;
		else
			for (j = 0; j < Ne; j++)
				if (a[i] == b[j])
				{
					// swap b_i and b_j
					temp = b[i];	b[i] = b[j];	b[j] = temp;
					// and change sign of the b string
					*sign_b = - *sign_b;
					flag = 0;
					break;
				}
		if (flag)	count ++;
	}
	j = 0;
	if (count <= 2)
		for (i = 0; i < Ne; i++)
			if (a[i] != b[i])	{pos[j] = i;	j++;}
	return count;	
}

FCI::FCI (HRED& hr)
{
	long int count, tot = 0; 
	int i;
	K = hr.Ni;	N = K / 2;
	h = hr.h;	V = hr.V;
	astr = new AB_STRING[N + 1];
	bstr = new AB_STRING[N + 1];
	MatrixXi config_summary (N + 1, 2);
	for (i = 0; i <= N; i++)
	{
		count = 0;
		astr[i]._init_ (K, N);
		astr[i]._str_new_ (i);
		astr[i]._str_gen_ (0, N, i, &count);
		tot += count;
		config_summary (i, 0) = i;
		config_summary (i, 1) = (int) count;

		count = 0;
		bstr[i]._init_ (K, N);
		bstr[i]._str_new_ (i);
		bstr[i]._str_gen_ (0, N, i, &count);
	}
	TOT = tot * tot;
	cout << "Nexct Nconf\n" << config_summary << "\n\n";
}

void FCI:: _fci_ ()
{
	IOFormat Full(FullPrecision, 0, "", "\n", "", "");

	int Na = N, Nb = N;
	int Ia, Ja, ia, ja, Ib, Jb, ib, jb;
	int i, mu, nu, la, si;
	long int row, col, tot_Ia, tot_Ja, tot_Ib, tot_Jb;
	int diff_a, diff_b, pos_a[2], pos_b[2], id1, id2;
	short *arow, *acol, *brow, *bcol, sign_a, sign_b;
	double tempEa, tempEa0, tempEb, tempEab;
	
	cout << "Total number of unique configurations: " << TOT << endl << endl;

	MatrixXd FCI(TOT, TOT);	FCI.setZero ();
	
	row = col = 0;
	for (Ia = 0; Ia <= Na; Ia++)
	{
		tot_Ia = astr[Ia].tot;
		for (ia = 0; ia < tot_Ia; ia++)
		{
			arow = astr[Ia].str[ia];

			tempEa0 = 0;
			for (mu = 0; mu < Na; mu++)
			{
				tempEa0 += h(arow[mu], arow[mu]);
				for (nu = mu + 1; nu < Na; nu++)
					tempEa0 += V[index4(arow[mu],arow[nu],arow[mu],arow[nu],K)] - 
						V[index4(arow[mu],arow[nu],arow[nu],arow[mu],K)];
			}

			for (Ib = 0; Ib <= Nb; Ib++)
			{
				tot_Ib = (Nb)?(bstr[Ib].tot):(1);
				for (ib = 0; ib < tot_Ib; ib++)
				{
					brow = (Nb)?(bstr[Ib].str[ib]):(NULL);

					tempEb = 0;
					for (mu = 0; mu < Nb; mu++)
					{
						tempEb += h(brow[mu], brow[mu]);
						for (nu = mu + 1; nu < Nb; nu++)
							tempEb += V[index4(brow[mu],brow[nu],brow[mu],brow[nu],K)] - 
								V[index4(brow[mu],brow[nu],brow[nu],brow[mu],K)];
					}

					tempEab = 0;
					if (Nb)
						for (mu = 0; mu < Na; mu++)
							for (nu = 0; nu < Nb; nu++)
								tempEab += V[index4(arow[mu],brow[nu],arow[mu],brow[nu],K)];

					FCI (row, row) = tempEa0 + tempEb + tempEab;		// diagonal
					col = 0;	// don't forget to reset col!!!
					
					for (Ja = 0; Ja <= Na; Ja++)
					{
						tot_Ja = astr[Ja].tot;
						for (ja = 0; ja < tot_Ja; ja++)
						{
							acol = astr[Ja].str[ja];	sign_a = astr[Ja].s[ja];

							// string matching number for alpha
							if (abs (Ja - Ia) > 2)
								diff_a = 3;
							else
							{
								diff_a = _match_str_ (arow, acol, &sign_a, Na, pos_a);
								astr[Ja].s[ja] = sign_a;
								sign_a *=  astr[Ia].s[ia];
							}

							for (Jb = 0; Jb <= Nb; Jb++)
							{
								tot_Jb = (Nb)?(bstr[Jb].tot):(1);
								for (jb = 0; jb < tot_Jb; jb++)
								{
									if (col <= row)	{col++;	continue;}

									if (Nb)	{bcol = bstr[Jb].str[jb];	sign_b = bstr[Jb].s[jb];}
									else	{bcol = NULL;	sign_b = sign_a;}
									// string matching number for beta
									if (Nb)
										if (abs (Jb - Ib) > 2)
											diff_b = 3;
										else
										{
											diff_b = _match_str_ (brow, bcol, &sign_b, Nb, pos_b);
											bstr[Jb].s[jb] = sign_b;
											sign_b *= bstr[Ib].s[ib];
										}
									else	diff_b = 0;

									// differ too much ==> zero
									if (diff_a + diff_b > 2)	{col++;	continue;}
									
									// differ Not that much
									if (diff_a + diff_b == 2)
										if (Nb)
											if (diff_a == 2)
											{
												id1 = pos_a[0];	id2 = pos_a[1];
												tempEa = 
													(V[index4(arow[id1],arow[id2],acol[id1],acol[id2],K)]
													- V[index4(arow[id1],arow[id2],acol[id2],acol[id1],K)]);
												tempEb = 0;
												tempEab = 0;
											}
											else if (diff_a == 1)
											{
												id1 = pos_a[0];	id2 = pos_b[0];
												tempEa = tempEb = 0;
												tempEab =
													V[index4(arow[id1],brow[id2],acol[id1],bcol[id2],K)];
											}
											else
											{
												id1 = pos_b[0];	id2 = pos_b[1];
												tempEb = 
													(V[index4(brow[id1],brow[id2],bcol[id1],bcol[id2],K)]
													- V[index4(brow[id1],brow[id2],bcol[id2],bcol[id1],K)]);
												tempEa = 0;
												tempEab = 0;
											}
										else
										{
											id1 = pos_a[0];	id2 = pos_a[1];
											tempEa =
												V[index4(arow[id1],arow[id2],acol[id1],acol[id2],K)]
												- V[index4(arow[id1],arow[id2],acol[id2],acol[id1],K)];
											tempEb = tempEab = 0;
										}
									else if (diff_a + diff_b == 1)
										if (Nb)
											if (diff_a == 1)
											{
												id1 = pos_a[0];
												tempEa = h(arow[id1], acol[id1]);
												for (mu = 0; mu < Na; mu++)
													if (mu != id1)	tempEa += 
														V[index4(arow[mu],arow[id1],arow[mu],acol[id1],K)] -
														V[index4(arow[mu],arow[id1],acol[id1],arow[mu],K)];
												tempEb = tempEab = 0;
												for (mu = 0; mu < Nb; mu++)
													tempEab +=
														V[index4(brow[mu],arow[id1],brow[mu],acol[id1],K)];
											}
											else
											{
												id1 = pos_b[0];
												tempEb = h(brow[id1], bcol[id1]);
												for (mu = 0; mu < Na; mu++)
													if (mu != id1)	tempEb += 
														V[index4(brow[mu],brow[id1],brow[mu],bcol[id1],K)] -
														V[index4(brow[mu],brow[id1],bcol[id1],brow[mu],K)];
												tempEa = tempEab = 0;
												for (mu = 0; mu < Nb; mu++)
													tempEab +=
														V[index4(arow[mu],brow[id1],arow[mu],bcol[id1],K)];
											}
										else
										{
											id1 = pos_a[0];
											tempEa = h (arow[id1], acol[id1]);
											for (mu = 0; mu < Na; mu++)
												if (mu != id1)	tempEa +=
													V[index4(arow[mu],arow[id1],arow[mu],acol[id1],K)]
													- V[index4(arow[mu],arow[id1],acol[id1],arow[mu],K)];
											tempEb = tempEab = 0;
										}
									else
										cout << "zhen de you ei!\n";
									if (Nb)
										FCI (row, col) = FCI (col, row) = sign_a * sign_b * (tempEa + tempEb + tempEab);
									else
										FCI (row, col) = FCI (col, row) = sign_a * tempEa;

									col++;
								}
							}
						}
					}
					row++;
				}
			}
		}
	}

	cout << "FCI matrix:" << endl;
	if (TOT <= 10)	
		cout << FCI << endl << endl;
	else
		cout << FCI.topLeftCorner (10, 10) << endl << endl;

	// diagonalize FCI matrix
	SelfAdjointEigenSolver<MatrixXd> ess;
	ess.compute (FCI);
	VectorXd Es = ess.eigenvalues ();	// Eigenvalues
	Es = Es.array();
	//MatrixXd Us = ess.eigenvectors ();		// Eigenvectors
	//_init_Sspin_ ();
	//Es.col(1) = (Us.transpose() * Sspins * Us).diagonal();
	//for (i = 0; i < K2; i++)	Es (i, 1) = round (Es (i, 1));
	cout << "FCI energies for (Na, Nb) = (" << Na << ", " << Nb << "):" << endl << 
		Es.format(Full) << endl << endl;

}

#endif
