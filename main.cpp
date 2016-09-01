#include <iostream>
#include <cstdlib>
#include <Eigen/Dense>
#include "include/hf.h"
#include "include/dmet.h"
#include "include/hubbard.h"
#include "include/schmidt.h"
#include "include/hred.h"
#include "include/scf.h"
#include "include/dfci.h"
#include "include/read.h"

using namespace std;
using namespace Eigen;

int main (int argc, char * argv[])
{
	HUBBARD hub;
	SCHMIDT sm;

	// read parameters from the input file
	if (argc == 1)	_read_ ("input", hub, sm);
	else if (argc == 2)	_read_ (argv[1], hub, sm);
	else	{cout << "Too many input files!\n";	exit (1);}

	// set up ioff -- the lookup table
	int K = sm.Nimp * 2;
	_gen_ioff_ (K * (K + 1) / 2);

	// Hubbard Hartree-Fock calculation
	hub._hubbard_rhf_ ();
	hub._print_ ();

	// Schmidt
	sm._schmidt_ (hub);

	// Construct Hred
	HRED hr;
	hr._xform_ (hub, sm);
	//hr._write_ ();

	// SCF on Hred, for CHECK
	SCF hred_scf (hr.h, hr.V, hr.Ni, hr.Ni / 2);
	hred_scf._scf_ ();

	DMET dmet;
	printf ("HF-in-HF embedding energy: %18.16f\n\n",
		dmet._dmet_energy_ (hred_scf.h, hred_scf.V, hred_scf.P, hred_scf.N));

	// FCI on fragment
	DFCI dfci;
	dfci._init_ (hr);
	cout << "FCI initialization succeeds!\n" << dfci.tot <<
		" alpha strings are generated!\n\n";
	dfci._dfci_ ();
	dfci._1PDM_ ();
	cout << "scf 1PDM:\n" << hred_scf.P << "\n\n";
	cout << "dfci 1PDM:\n" << dfci.P << "\n\n";
	//cout << "check idempotency:\n" << hred_scf.P * hred_scf.P << "\n\n";
	dfci._2PDM_ ();
	printf ("FCI-in-HF embedding energy: %18.16f\n\n",
			dmet._dmet_energy_ (hr.h, hred_scf.V, dfci.P, dfci.G, dfci.N));
	
	return 0;
}
