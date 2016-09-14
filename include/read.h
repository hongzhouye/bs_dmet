#ifndef _READ_H_INCLUDED_
#define _READ_H_INCLUDED_

#include <iostream>
#include <string>

/* FOR split */
#include <sstream>
#include <vector>

/* FOR _uppercase_ */
#include <fstream>
#include <algorithm>
#include <cctype>
#include <functional>

/* FOR dmet */
#include "hubbard.h"
#include "frag.h"
#include "hf.h"

using namespace std;

/* These two functions split a string with given delimeter */
void _split_ (const string &s, char delim, vector<string> &elems) {
    stringstream ss;
    ss.str(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

vector<string> _split_ (const string &s, char delim) {
    vector<string> elems;
    _split_ (s, delim, elems);
    return elems;
}

/* This function turns a string into all uppercase */
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

/* split a string like
 *      "N = 5" or "N=5" or "N= 5" or "N =5"
 * into a vector<string> with elements "N" and "5"
 */
vs _split_eq_ (string line)
{
    vs temp = _split_ (line, '=');

    vs namestr = _split_ (temp[0], ' ');
	string name = *(namestr.begin ());

	vs valuestr = _split_ (temp[1], ' ');
	string value = valuestr[valuestr.size () - 1];

    vs nv;
    nv.push_back (name);    nv.push_back (value);

    return nv;
}

void _read_ (char *fname, HUBBARD& hub, FRAG& frag)
{
	int i;
	string line;
	ifstream input (fname);
	if (input)	// same as: if (input.good ())
	{
		while (getline (input, line))	// same as ".good ()"
		{
            if (line.empty ())  continue;
			else if (_uppercase_ (line) == "&HUBBARD")
            /*** HUBBARD card ***/
            {
                while (getline (input, line))
    			{
                    if (line.empty ())  continue;
                    else if (_uppercase_ (_split_ (line, ' ')[0]) == "K")
    				    hub.K = stoi (_split_eq_ (line)[1]);
                    else if (_uppercase_ (_split_ (line, ' ')[0]) == "N")
                        hub.N = stoi (_split_eq_ (line)[1]);
                    else if (_uppercase_ (_split_ (line, ' ')[0]) == "BC")
                        hub.BC = _split_eq_ (line)[1];
                    else if (_uppercase_ (_split_ (line, ' ')[0]) == "U")
                        hub.U = stod (_split_eq_ (line)[1]);
                    else if (_uppercase_ (line) == "&END HUBBARD")  break;
                    else
                    {
                        cout << "Input file error: check the HUBBARD card.\n\n";
                        exit (1);
                    }
    			}
            }
			else if (_uppercase_ (line) == "&FRAGMENT")
            /*** FRAGMENT card ***/
			{
                while (getline (input, line))
                {
                    if (line.empty ())  continue;
                    // Nimp & frag sites
                    else if (_uppercase_ (_split_ (line, ' ')[0]) == "NIMP")
                    {
                        // Nimp
        				frag.Nimp = stoi (_split_eq_ (line)[1]);
                        // frag sites
                        for (int i = 0; i < frag.Nimp; i++)
                        {
                            getline (input, line);
                            vs fragstr = _split_ (line, ' ');
                            frag.fragments[fragstr[0]] = stoi (fragstr[1]);
                        }
                    }
                    // center site(s)
                    else if (_uppercase_ (_split_ (line, ' ')[0]) == "CENTER")
                    {
                        vs csitestr = _split_ (_split_eq_ (line)[1], ';');
                        frag.Ncenter = csitestr.size ();     // # of center sites
                        frag.center = new int [frag.Ncenter];
                        for (int i = 0; i < frag.Ncenter; i++)
                        // convert names to indices
                            frag.center[i] = frag.fragments[csitestr[i]];
                    }
                    // population constraints
                    else if (_uppercase_ (_split_ (line, ' ')[0]) == "NPOP")
                    {
                        frag.Npop = stoi (_split_eq_ (line)[1]);
                        for (int i = 0; i < frag.Npop; i++)
                        {
                            getline (input, line);
                            vs templine = _split_ (line, ';');
                            vi tempvi;
                            for (int j = 0; j < templine.size (); j++)
                                tempvi.push_back (frag.fragments[templine[j]]);
                            frag.popcon.push_back (tempvi);
                        }
                    }
                    else if (_uppercase_ (_split_ (line, ' ')[0]) == "N2E")
                    {
                        frag.N2e = stoi (_split_eq_ (line)[1]);
                        for (int i = 0; i < frag.N2e; i++)
                        {
                            getline (input, line);
                            vs templine = _split_ (line, ' ');
                            // bad sites
                            vs badstr = _split_ (templine[0], ';');
                            vi tempbad;
                            for (int j = 0; j < badstr.size (); j++)
                                tempbad.push_back (frag.fragments[badstr[j]]);
                            frag.bad_2econ.push_back (tempbad);
                            // good sites
                            vs goodstr = _split_ (templine[1], ';');
                            vi tempgood;
                            for (int j = 0; j < goodstr.size (); j++)
                                tempgood.push_back (frag.fragments[goodstr[j]]);
                            frag.good_2econ.push_back (tempgood);
                        }
                    }
                    else if (_uppercase_ (line) == "&END FRAGMENT") break;
                    else
                    {
                        cout << "Input file error: check the FRAGMENT card.\n\n";
                        exit (1);
                    }
                }
			}
		}
		input.close ();
	}
	else	cout << "Hey man please give me a valid INPUT file, okay?\n";
}

#endif
