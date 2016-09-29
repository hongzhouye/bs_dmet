#include <iostream>
#include <cstdlib>
#include <Eigen/Dense>
#include "include/hf.h"
#include "include/dmet.h"

using namespace std;
using namespace Eigen;
using namespace std::chrono;

int main (int argc, char * argv[])
{
	DMET dmet;

	// initialization
	if (argc == 1)	dmet._dmet_init_ ("input");
	else if (argc == 2)	dmet._dmet_init_ (argv[1]);
	else	{cout << "Too many input files!\n";	exit (1);}

	// bootstrap
	high_resolution_clock::time_point t1 = high_resolution_clock::now ();
	dmet._bs_dmet_ ();
	high_resolution_clock::time_point t2 = high_resolution_clock::now ();
	duration<double> dt_bs = duration_cast<duration<double> >(t2 - t1);
	cout << "=========================\n";
	cout << "Total wall time: " << dt_bs.count () << " sec" << endl << endl;

	// check
	//dmet._dmet_check_ ();

	// FCI check
	//dmet._fci_check_ ();

	return 0;
}
