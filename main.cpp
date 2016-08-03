#include <iostream>
#include "include/hf.h"
#include "include/hubbard.h"

using namespace std;

int main ()
{
	HUBBARD hub (4, "core");
	hub._hubbard_general_ ();
	hub._print_ ();
	return 0;
}
