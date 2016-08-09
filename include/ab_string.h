#ifndef _AB_STRING_H_INCLUDED_
#define _AB_STRING_H_INCLUDED_

#include <cmath>
#include <iostream>

using namespace std;

long int _nchoosek_ (int n, int k)
{
	int i;
	double prod = 1.;
	if (k == 0)	return 1;
	for (i = 1; i <= k; i++)
		prod *= (double) (n - k + i) / i;
	return (long int) prod;
}

class AB_STRING
{
  public:
	  short **str;		// alpha/beta strings, stored in 2D arrays
	  short *s;
	  int K;			// basis set size
	  int Ne;			// # of alpha/beta occupied states
	  long int tot;		// total # of alpha/beta configurations
	  void _init_ (int, int);
	  void _str_new_ (int);
	  void _str_delete_ ();
	  void _str_gen_ (int, int, int, long int *);
};

void AB_STRING::_init_ (int Nbs, int Nel)
{
	K = Nbs; Ne = Nel;
}

void AB_STRING::_str_new_ (int Nexct)
{
	int i, j;

	// total # of configurations for a given Ne
	// and degree of excitation Nexct
	tot = _nchoosek_ (Ne, Nexct) * _nchoosek_ (K - Ne, Nexct);

	// alpha/beta strings are stored as short type 2D arrays
	// i.e. alpha[tot][Ne], beta[tot_b][Nb]
	str = new short * [tot];
	for (i = 0; i < tot; i++)	str[i] = new short[Ne];
	s = new short [tot];

	// initialize alpha/beta strings, don't forget sign
	for (i = 0; i < tot; i++)
	{
		s[i] = 1;
		for (j = 0; j < Ne; j++)
			str[i][j] = j;
	}
}

void AB_STRING::_str_delete_ ()
{
	int i;

	// delete alpha/beta strings
	for(i = 0; i < tot; i++)	delete [] str[i];
	delete [] str;
}

void AB_STRING::_str_gen_ (int occ_init, int virt_init, int Nexct, long int * count)
{
	int i, j;

	if (Nexct == 0)
	{
		/*cout << "count = " << *count;
		for (i = 0; i < Ne; i++)
			cout << "\t" << str[*count][i];
		cout << endl;
		*/
		*count += 1;
		return;
	}
	else
	{
		for (i = occ_init; i <= Ne - Nexct; i++)
			for (j = virt_init; j <= K - Nexct; j++)
			{
				str[*count][i] = j;
				_str_gen_ (i + 1, j + 1, Nexct - 1, count);
				if ((i != Ne - Nexct || j != K - Nexct) && *count < tot)
					for (int k = 0; k < i; k++)
						str[*count][k] = str[*count - 1][k];
			}
	}
}

#endif
