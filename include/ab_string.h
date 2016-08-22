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

#endif
