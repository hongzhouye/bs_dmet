#ifndef _FRAG_H_INCLUDED_
#define _FRAG_H_INCLUDED_

#include <iostream>
#include <string>
#include <map>
#include "hf.h"

class FRAG
{
    public:
        int Nimp;           // # of sites
        msi fragments;      // map<string, int>

        int Ncenter;        // # of center sites
        int *center;        // index of center site

        int Npop;           // # of pop constraints
        vvi popcon;         // constr[i][j]: the i-th constr
                            // is applied to pop of site j

        int N2e;
        vvi bad_2econ;
        vvi good_2econ;

        MatrixXd h;         // 1P part of Hamiltonian
        double *V;          // 2P part of Hamiltonian
};

#endif
