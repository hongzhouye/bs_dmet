#ifndef _BOOTSTRAP_H_INCLUDED_
#define _BOOTSTRAP_H_INCLUDED_

#include "hf.h"
#include "hubbard.h"
#include "frag.h"
#include "dfci.h"
#include "nropt.h"

class BOOTSTRAP
{
    public:
        int Ni;
        FRAG frag;

        void _init_ (const HUBBARD&);
        void _bs_opt_ ();
};

void BOOTSTRAP::_init_ (const HUBBARD& hub)
{
    frag._init_ (hub);
    Ni = frag.Nimp * 2;

    frag.hpot.setZero (Ni, Ni);
    int lenh = Ni * (Ni + 1) / 2, lenV = lenh * (lenh + 1) / 2;
    frag.Vpot = _darray_gen_ (lenV);
}

void BOOTSTRAP::_bs_opt_ ()
{
    VectorXd _dmet_iter_ (VectorXd&, FRAG&);

    VectorXd u; u.setZero (frag.Nopt);
    u = _bfgs_opt_ (_dmet_iter_, u, frag);
}

VectorXd _dmet_iter_ (VectorXd& u, FRAG& frag)
{
    void _make_pot_ (const VectorXd&, FRAG&);
    VectorXd _calc_obj_ (const FRAG&);

    _make_pot_ (u, frag);
    /*cout << "u = " << u.transpose ().format (Short) << "\n";
    cout << "hpot:\n" << frag.hpot.format (Short) << "\n\n";*/
    frag._solver_ (true);
    return _calc_obj_ (frag);
}

VectorXd _calc_obj_ (const FRAG& frag)
{
    VectorXd obj;   obj.setZero (frag.Nopt);
    int shift = 0;

    // filling
    for (int i = 0; i < frag.Npop; i++)
    {
        // when degenerate, only consider one of them
        int j = frag.popcon[i][0];
        obj[shift] = frag.dfci.P(j, j) - frag.target_filling;
        shift++;
    }

    // 1e coherence
    for (int i = 0; i < frag.N1e; i++)
    {
        int goodi = frag.good_1econ[i][0];
        int goodj = frag.good_1econ[i][1];
        // when degenerate, only consider one of them
        int badi = frag.bad_1econ[i][0][0];
        int badj = frag.bad_1econ[i][0][1];
        obj[shift] = frag.dfci.P(badi, badj) - frag.dfci.P(goodi, goodj);
        shift++;
    }

    // 2e on-top
    for (int i = 0; i < frag.N2e; i++)
    {
        int goodj = frag.good_2econ[i];
        // when degenerate, only consider one of them
        int badj = frag.bad_2econ[i][0];
        int gjj, gjjjj, bjj, bjjjj;
        gjj =  cpind(goodj,goodj);  gjjjj = cpind(gjj,gjj);
        bjj =  cpind(badj,badj);    bjjjj = cpind(bjj,bjj);
        obj[shift] = frag.dfci.G[bjjjj] - frag.dfci.G[gjjjj];
        shift++;
    }

    // check
    if (shift != frag.Nopt)
    {
        cout << "ERROR in BS::_calc_obj_: # of opt param's does not check!\n\n";
        exit (0);
    }

    return obj;
}

void _make_pot_ (const VectorXd& u, FRAG& frag)
{
    int shift = 0;

    // filling
    for (int i = 0; i < frag.Npop; i++)
    {
        vi pc = frag.popcon[i];
        for (int j = 0; j < pc.size (); j++)
            frag.hpot(pc[j], pc[j]) = u[shift];
        shift++;
    }

    // 1e coherence
    for (int i = 0; i < frag.N1e; i++)
    {
        vis bad_1e = frag.bad_1econ[i];
        for (int j = 0; j < bad_1e.size (); j++)
        {
            int indi = bad_1e[j][0];
            int indj = bad_1e[j][1];
            frag.hpot(indi, indj) = u[shift];
            frag.hpot(indj, indi) = u[shift];
        }
        shift++;
    }

    // 2e on-top
    for (int i = 0; i < frag.N2e; i++)
    {
        // when making pots, it has nothing to do with good
        vi bad_2e = frag.bad_2econ[i];
        for (int j = 0; j < bad_2e.size (); j++)
        {
            int jj = cpind(bad_2e[j],bad_2e[j]);
            int jjjj = cpind(jj,jj);
            frag.Vpot[jjjj] = u[shift];
        }
        shift++;
    }

    // check
    if (shift != frag.Nopt)
    {
        cout << "ERROR in BS::_make_pot_: # of opt param's does not check!\n\n";
        exit (0);
    }
}

#endif
