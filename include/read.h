#ifndef _READ_H_INCLUDED_
#define _READ_H_INCLUDED_

#include <iostream>
#include <fstream>
#include <algorithm>
#include <cctype>
#include <functional>
#include <string>
#include "hubbard.h"
#include "schmidt.h"

using namespace std;

string _uppercase_ (const string& s)
{
	string result (s.length (), ' ');
	transform (
			s.begin (),
			s.end (),
			result.begin (),
			ptr_fun <int, int> (toupper)
			);
	return result;
}

void _read_ (char * fname, HUBBARD& hub, SCHMIDT& sm)
{
	string line;
	ifstream input (fname);
	if (input)	// same as: if (input.good ())
	{
		while (getline (input, line))	// same as ".good ()"
		{
			if (_uppercase_ (line) == "HUBBARD")
			{
				getline (input, line);	hub.K = stoi (line);
				getline (input, line);	hub.N = stoi (line);
				input.get (hub.BC);		getline (input, line);	// eat '\n'
				getline (input, line);	hub.U = stod (line);
			}
			else if (_uppercase_ (line) == "SCHMIDT")
			{
				getline (input, line);	sm.Nimp = stoi (line);
				sm.K = hub.K; sm.N = hub.N;
			}
		}
		input.close ();
	}
	else	cout << "Hey man please give me a valid INPUT file, okay?\n";
}

#endif
