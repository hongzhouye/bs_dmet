#include <iostream>
#include <cstdlib>
#include "include/hf.h"
#include "include/hubbard.h"
#include "include/schmidt.h"
#include "include/hred.h"
#include "include/dfci.h"
#include "include/read.h"

using namespace std;

int main (int argc, char * argv[])
{
	HUBBARD hub;
	SCHMIDT sm;

	// read parameters from the input file
	if (argc == 1)	_read_ ("input", hub, sm);
	else if (argc == 2)	_read_ (argv[1], hub, sm);
	else	{cout << "Too many input files!\n";	exit (1);}

	// Hubbard Hartree-Fock calculation
	hub._hubbard_general_ ();
	hub._print_ ();

	// Schmidt
	sm._schmidt_ (hub);

	// Construct Hred
	HRED hr;
	hr._xform_ (hub, sm);

	// FCI on fragment
	DFCI dfci;
	dfci._init_ (hr);
	cout << "FCI initialization succeeds!\n";
	dfci._dfci_ ();
	return 0;
}
