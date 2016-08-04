#include <iostream>
#include <cstdlib>
#include "include/hf.h"
#include "include/hubbard.h"
#include "include/schmidt.h"
#include "include/read.h"

using namespace std;

int main (int argc, char * argv[])
{
	if (argc < 2)
	{
		cout << "Please provide an input file!\n";
		exit (1);
	}

	HUBBARD hub;
	SCHMIDT sm;

	// read parameters from the input file
	_read_ (argv[1], hub, sm);

	// Hubbard Hartree-Fock calculation
	hub._hubbard_general_ ();
	hub._print_ ();

	// Schmidt
	sm._schmidt_ (hub);
	return 0;
}
