#include <iostream>
#include <cstdlib>
#include <Eigen/Dense>
#include "include/hf.h"
#include "include/dmet.h"

using namespace std;
using namespace Eigen;

int main (int argc, char * argv[])
{
	DMET dmet;

	// initialization
	if (argc == 1)	dmet._dmet_init_ ("input");
	else if (argc == 2)	dmet._dmet_init_ (argv[1]);
	else	{cout << "Too many input files!\n";	exit (1);}

	// check
	dmet._dmet_check_ ();

	return 0;
}
