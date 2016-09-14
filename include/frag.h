#ifndef _FRAG_H_INCLUDED_
#define _FRAG_H_INCLUDED_

#include "hubbard.h"
#include "schmidt.h"
#include "hred.h"
#include "hf.h"

class FRAG
{
    public:
        int Nimp;           // # of sites
        msi fragments;      // map<string, int>

        int Ncenter;        // # of center sites
        int *center;        // index of center site

        int Npop;           // # of pop constraints
        vvi popcon;         // popcon[i][j]: the i-th constr
                            // is applied to pop of site j

        int N2e;            // # of 2e on-top constraints
        vvi bad_2econ;      // to be matched
        vvi good_2econ;     // matching target

        MatrixXd h;         // 1P part of Hamiltonian
        double *V;          // 2P part of Hamiltonian

        SCHMIDT sm;
        HRED hr;

        void _init_ (const HUBBARD&);
};

void FRAG::_init_ (const HUBBARD& hub)
{
    // set up h and V
    int Ni = 2 * Nimp;
    h.setZero (Ni, Ni);
    int lenh = Ni * (Ni + 1) / 2, lenV = lenh * (lenh + 1) / 2;
	V = _darray_gen_ (lenV);

    // set up sm
    sm._init_ (hub, Nimp);
    int count = 0;
	for (auto i = fragments.begin (); i != fragments.end (); i++)
		sm.frag[count++] = i->second;

    // Schmidt decomposition
    sm._schmidt_ (hub);

    // make hred
    hr._xform_ (hub, sm, h, V);
}
#endif
