#ifndef _FCIWRAP_H_INCLUDED_
#define _FCIWRAP_H_INCLUDED_

#include "troyfci.h"
#include "hf.h"
#include <cstdio>

class FCIWRAP
/*********************************************************************
*   N       - Number of Basis Functions
*   No      - Number of Alpha electrons (=Number of Beta Electrons)
*   N2      - N * (N - 1) / 2
*   Nstr    - Number of Alpha CI Strings (=Number of Beta Strings)
*   N0      - Number of strings in preconditioner space
*   NS      - Number of eigenstates desired
*
*   h       - One electron Hamilontian
*   V       - Two electron Hamiltonian
*
*   Xi      - Input: NS Approximate guess eigenvectors
*           - Output: Lowest NS Eigenvectors
*   P       - One particle alpha transition density matrix
*   G       - 2 particle density matrix. Currently only for ground state
*   Ei      - Lowest NS eigenvalues
*   H1s     - Full hamiltonian for 1-site embedding
**********************************************************************/
{
    private:
        void _get_h_ (const MatrixXd&);
    public:
        int N, No, Nstr, N0MAX;
        double *h, *G, Ei, Sym;
        MatrixXd P;
        vMatrixXd Xi;
        bool guess_read;
        void _init_ (int, int, int);
        void _solve_ (const MatrixXd&, double *);
};

void FCIWRAP::_get_h_ (const MatrixXd& hinp)
{
    int i, j, ij;
    for (i = 0; i < N; i++) for (j = 0; j <= i; j++)
        h[cpind(i,j)] = hinp(i, j);
}

void FCIWRAP::_init_ (int Nbs, int Ne, int N0max)
{
//  gets values from input
    N = Nbs;    No = Ne;    N0MAX = N0max;
    Nstr = _nchoosek_ (N, No);

//  allocate memory for all matrices
    int lenh = N * (N + 1) / 2, lenV = lenh * (lenh + 1) / 2;
    h = _darray_gen_ (lenh);
    P.setZero (N, N);           G = _darray_gen_ (lenV);
}

void FCIWRAP::_solve_ (const MatrixXd& hinp, double *Vinp)
{
//  get h from input
    _get_h_ (hinp);

    //void _write_ (const MatrixXd&, double *);
    //if (Xi.size () == 0)    _write_ (hinp, Vinp);

//  allocate memory for NS-dependent quantities
    if (!guess_read)
    {
        if (Xi.size ()) Xi.clear ();
        Xi.push_back (MatrixXd::Zero (Nstr, Nstr));
    }

//  No Time To Explain! Get On The Car!
    _FCIman_ (N, No, N0MAX, 1, h, Vinp, Xi, &Ei, &Sym, P, G);
}

void _write_ (const MatrixXd& h, double *V)
{
    int N = h.rows ();
    int i,j,k,l;
    FILE *ph = fopen ("h", "w+");
    FILE *pV = fopen ("V", "w+");
    for (i = 0; i < N; i++) for (j = 0; j < N; j++)
    {
        fprintf (ph, "%d;%d;%18.16f\n", i,j,h(i, j));
        for (k = 0; k < N; k++) for (l = 0; l < N; l++)
        {
            int ik = cpind(i,k);    int jl = cpind(j,l);
            fprintf (pV, "%d;%d;%d;%d;%18.16f\n", i,j,k,l,V[cpind(ik,jl)]);
        }
    }
    fclose (ph);
    fclose (pV);
}

#endif
