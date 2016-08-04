#include <iostream>
#include <cstdlib>
#include "include/hf.h"
#include "include/hubbard.h"

using namespace std;

int main (int argc, char * argv[])
{
	if (argc < 2)
	{
		cout << "Please provide an input file!\n";
		exit (1);
	}

	HUBBARD hub (argv[1]);
	hub._hubbard_general_ ();
	hub._print_ ();
	return 0;
}
