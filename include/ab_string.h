#ifndef _AB_STRING_H_INCLUDED_
#define _AB_STRING_H_INCLUDED_

#include <cmath>
#include <iostream>

using namespace std;

// ref: J. Olsen and B. O. Roos, J. Chem. Phys. 89 (4), 1988

class AB_STRING
{
  public:
	  short *str;		// alpha/beta strings, stored in 1D array
	  long int *istr;	// indices of interacting strings
	  short *sgn;		// sign of interacting strings
	  int *cmpind;		// compound index kl for \hat{E}_{kl}
	  int K;			// basis set size
	  int Ne;			// # of alpha/beta occupied states
	  int itot;			// total # of interacting strings
	  void _init_ (int, int);
      int _check_connect_ (long int);
};

void AB_STRING::_init_ (int Nbs, int Nel)
// Nbs: # of basis set functions
// Nel: # of electrons
{
	K = Nbs; Ne = Nel;
	itot = Ne * (K - Ne);	// only single excitation
    itot += Ne;             // connections to itself

	str = new short [Ne];
	sgn = new short [itot];
	istr = new long int [itot];
	cmpind = new int [itot];

	// initialize alpha/beta strings, don't forget sign
	int i;
	for (i = 0; i < Ne; i++)	str[i] = i;
}

int AB_STRING::_check_connect_ (const long int I)
{
    int low = 0, high = itot - 1, ind = -1;

    // cheat: check two ends
    if (istr[low] == I)  return low;
    else if (istr[high] == I)    return high;

    bool flag;
    do
    {   // largest integer less than (low + high) / 2.
        ind = (low + high) / 2;
        if (istr[ind] == I) break;
        else if (I > istr[ind]) low = ind;
        else    high = ind;

        if (abs (low - high) <= 1)
        {
            ind = -1;   break;
        }
    }
    while (true);

    return ind;
}

#endif
