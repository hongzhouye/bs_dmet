#ifndef _FRAG_H_INCLUDED_
#define _FRAG_H_INCLUDED_

#include "hubbard.h"
#include "schmidt.h"
#include "hred.h"
#include "fciwrap.h"
#include "hf.h"

class FRAG
{
    public:
        int Nopt;           // tot # of constraints
        double target_filling;

        int Nimp;           // # of sites
        msi fragments;      // map<string, int>

        int Ncenter;        // # of center sites
        int *center;        // index of center site

        int Npop;           // # of pop constraints
        vvi popcon;         // popcon[i][j]: the i-th constr
                            // is applied to pop of site j

        int N1e;            // # of 1e coherence constraints
        vvis bad_1econ;     // to be matched
        vis good_1econ;     // matching target

        int N2e;            // # of 2e on-top constraints
        vvi bad_2econ;      // to be matched
        vi good_2econ;      // matching target

        MatrixXd h, hpot;   // 1P part of Hamiltonian
        double *V, *Vpot, *Vtot;
                            // 2P part of Hamiltonian
        int lenV;

        SCHMIDT sm;
        HRED hr;
        //DFCI dfci;
        FCIWRAP dfci;

        void _init_ (const HUBBARD&);
        void _solver_ (bool);
};

void FRAG::_init_ (const HUBBARD& hub)
{
    // set up h and V
    int Ni = 2 * Nimp;
    h.setZero (Ni, Ni);
    hpot.setZero (Ni, Ni);
    int lenh = Ni * (Ni + 1) / 2;
    lenV = lenh * (lenh + 1) / 2;
	V = _darray_gen_ (lenV);
    Vpot = _darray_gen_ (lenV);
    Vtot = _darray_gen_ (lenV);

    // set up sm
    sm._init_ (hub, Nimp);
    int count = 0;
	for (auto i = fragments.begin (); i != fragments.end (); i++)
		sm.frag[count++] = i->second;

    // Schmidt decomposition
    sm._schmidt_ (hub);

    // make hred
    hr._xform_ (hub, sm, h, V);

    // solver
    //dfci._init_ (2 * Nimp, Nimp);
    //dfci.mode = "major";
    dfci._init_ (2 * Nimp, Nimp, 12);
    dfci.guess_read = false;

    // tot # of constraints
    Nopt = Npop + N1e + N2e;

    // target_filling
    target_filling = (double) hub.N / hub.K;
}

void FRAG::_solver_ (bool _2PDM_flag_)
{
    for (int i = 0; i < lenV; i++)  Vtot[i] = V[i] + Vpot[i];
    for (int i = 0; i < 2 * Nimp; i++)
        for (int j = 0; j < 2 * Nimp; j++)
            for (int k = 0; k < 2 * Nimp; k++)
                for (int l = 0; l < 2 * Nimp; l++)
                {
                    int ik = cpind(i,k), jl = cpind(j,l);
                    if (fabs (V[cpind(ik,jl)]) > 1.E-10)
                        printf ("%3d %3d %3d %3d %18.16f\n", i,j,k,l,V[cpind(ik,jl)]);
                }
    dfci._solve_ (h + hpot, Vtot);
    //dfci._1PDM_ ();
    //if (_2PDM_flag_)    dfci._2PDM_ ();
}
#endif
