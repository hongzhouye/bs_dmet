#ifndef _TROYFCI_H_INCLUDED_
#define _TROYFCI_H_INCLUDED_

#include <iostream>
#include <cstdio>
#include <cmath>
#include <Eigen/Dense>
#include "hf.h"

using namespace std;
using namespace Eigen;

#define THRESH 1.E-9
#define MAXITER 5000

// A global index matrix
int **Zindex;

void _Isort_ (int N, iv1& X, int *Isign)
//  Sorts an Integer arrary X by straight insertion
//  Isign is the sign of the Permutation needed to bring it in order
{
    int i, j, tmp;
    bool flag;

    *Isign = 1;
    for (i = 1; i < N; i++)
    {
        flag = true;
        tmp = X[i];
        for (j = i - 1; j >= 0; j--)
        {
            if (X[j] <= tmp)    {flag = false;   break;}
            else    {X[j + 1] = X[j];   *Isign = -*Isign;}
        }
        if (flag)  j = -1;
        X[j + 1] = tmp;
    }
}

void _RecurEx1_ (int N, int No, int i, int a, int Max1,
    iv1& Iocc, iv3& Ex1, int IRecur, int *m, int *Ind)
{
    int ii,jj,IsiO,IsaO,aind,iind,Isign,iimin;
    iv1 Isubst(No);

    if (!IRecur)    iimin = 0;
    else    iimin = Iocc[IRecur - 1] + 1;

    for (ii = iimin; ii < N; ii++)
    {
        Iocc[IRecur] = ii;

        if (IRecur == No - 1)
//      End Recursion Condition
        {
            IsiO = 0; IsaO = 0;
            for (jj = 0; jj < No; jj++)
                if (Iocc[jj] == i)  IsiO = 1;
                else if (Iocc[jj] == a) IsaO = 1;
            if (IsiO == 1 && IsaO == 0)
            {
                //_copy_array_ (Iocc, Isubst, No);
                Isubst = Iocc;
                for (jj = 0; jj < No; jj++)
                    if (i == Isubst[jj])    Isubst[jj] = a;
                _Isort_ (No, Isubst, &Isign);

                Ex1[*m][a][i] = Isign * (*Ind);
                *m += 1;
            }
            *Ind += 1;
        }
        else
            _RecurEx1_ (N, No, i, a, Max1, Iocc, Ex1,
                IRecur + 1, m, Ind);
    }
    //delete Isubst;
    return;
}

void _RecurEx2_ (int N, int No, int N2, int i, int j, int a, int b,
    int ij, int ab, int Max2, iv1& Iocc, iv3& Ex2, int IRecur,
    int *m, int *Ind)
{
    int ii,jj;
    bool IsiO,IsjO,IsaO,IsbO;
    int iind,jind,aind,bind,Isign,iimin;
    iv1 Isubst(No);

    if (!IRecur)    iimin = 0;
    else    iimin = Iocc[IRecur - 1] + 1;

    for (ii = iimin; ii < N; ii++)
    {
        Iocc[IRecur] = ii;

        if (IRecur == No - 1)
//      End Recursion Condition
        {
            IsiO = IsjO = IsaO = IsbO = false;

//          Is abj*i* |Iocc> Non-Zero?
            for (jj = 0; jj < No; jj++)
                if (Iocc[jj] == i)  IsiO = true;
                else if (Iocc[jj] == j)  IsjO = true;
                else if (Iocc[jj] == a)  IsaO = true;
                else if (Iocc[jj] == b)  IsbO = true;
            if (IsiO && IsjO && !IsaO && !IsbO)
            {
                //_copy_array_ (Iocc, Isubst, No);
                Isubst = Iocc;
                for (jj = 0; jj < No; jj++)
                    if (Isubst[jj] == i)    Isubst[jj] = a;
                    else if (Isubst[jj] == j)    Isubst[jj] = b;

                _Isort_ (No, Isubst, &Isign);

                Ex2[*m][ab][ij] = Isign * (*Ind);
                *m += 1;
            }
            *Ind += 1;
        }
        else
            _RecurEx2_ (N, No, N2, i, j, a, b, ij, ab,
                Max2, Iocc, Ex2, IRecur + 1, m, Ind);
    }
    //delete Isubst;
    return;
}

void _IString_ (int N, int No, int N2,
    int Max1, int Max2, iv3& Ex1, iv3& Ex2)
{
    int i, j, k, a, b, ij, ab, IRecur, m, Ind;
    iv1 Iocc(No);

//  Find Strings that differ by a Single Excitation
    for (i = 0; i < N; i++)
        for (a = 0; a < N; a++)
        {
            m = 0; Ind = 1; /* !!!!! INDEX STARTS FROM 1 !!!!! */
            _RecurEx1_ (N, No, i, a, Max1, Iocc, Ex1, 0, &m, &Ind);
        }

//  Find Strings that differ by a Double Excitation
    ij = -1;
    for (i = 0; i < N; i++)
        for (j = i + 1; j < N; j++)
        {
            ij++; ab = -1;
            for (a = 0; a < N; a++)
                for (b = a + 1; b < N; b++)
                {
                    ab++;  m = 0;
                    Ind = 1;    /* !!!!! INDEX STARTS FROM 1 !!!!! */
                    _RecurEx2_ (N, No, N2, i, j, a, b, ij, ab,
                        Max2, Iocc, Ex2, 0, &m, &Ind);
                }
        }
    //delete Iocc;
    return;
}

int _Index_ (int No, const iv1& Iocc)
{
    int i, Isign;
    int index = 0;  /* !!!!! INDEX STARTS FROM 1 !!!!! */

    for (i = 0; i < No; i++)    index += Zindex[i][Iocc[i]];
    return index;
}

void _GetHd_ (int N, int No, int Nstr, double *h, double *V,
    iv1& Iocca, iv1& Ioccb, MatrixXd& Hd, int IRecur)
{
    int i,j,imin,jmin,k,ka,kb,l,la,lb,Isigna,Isignb;
    int kka, kkb, lla, llb, kla, klb;
    double tmp;

    if (IRecur == 1)    imin = jmin = 0;
    else
    {
        imin = Iocca[IRecur - 2] + 1;
        jmin = Ioccb[IRecur - 2] + 1;
    }

    for (i = imin; i < N; i++)
    {
        Iocca[IRecur - 1] = i;
        for (j = jmin; j < N; j++)
        {
            Ioccb[IRecur - 1] = j;
            if (IRecur == No)
            {
                tmp = 0.;
                for (k = 0; k < No; k++)
                {
                    ka = Iocca[k];  kb = Ioccb[k];
                    kka = cpind(ka,ka); kkb = cpind(kb,kb);
//  Spin Contaminated Elements
                    tmp += h[kka] + h[kkb];
                    for (l = 0; l < No; l++)
                    {
                        la = Iocca[l];  lb = Ioccb[l];
                        lla = cpind(la,la); llb = cpind(lb,lb);
                        kla = cpind(ka,la); klb = cpind(kb,lb);
                        tmp += 0.5 * (
                            V[cpind(kka,lla)] - V[cpind(kla,kla)] +
                            V[cpind(kka,llb)] + V[cpind(kkb,lla)] +
                            V[cpind(kkb,llb)] - V[cpind(klb,klb)]);
                    }
                }
                Hd (_Index_(No, Iocca), _Index_ (No, Ioccb)) = tmp;
            }
            else
                _GetHd_ (N, No, Nstr, h, V, Iocca, Ioccb, Hd, IRecur + 1);
        }
    }
    return;
}

void _GetIstr_ (int N, int No, int Nstr, int N0,
    MatrixXd& Hd, iv1& Istr)
{
    int itmp, a, b, i;
    double min;

    for (a = 0; a < N0; a++)
    {
        min = 500.; itmp = 0;
        for (i = 0; i < Nstr; i++)
            if (Hd(i, i) < min) {min = Hd(i, i);    itmp = i;}

        Istr[a] = itmp;
        Hd(Istr[a], Istr[a]) += 1000.;
    }
    for (a = 0; a < N0; a++)
        Hd(Istr[a], Istr[a]) -= 1000.;
}

void _GetH0_ (int N, int No, int N0, int N2, int Max1, int Max2,
    int Nstr, const iv3& Ex1, const iv3& Ex2, const iv1& Istr,
    double *h, double *V, MatrixXd& H0)
{
    iv1 RIstr(Nstr), IY;   IY.assign (N0, -1);
    iv2 ifzero(N, iv1(N));
    iv3 AEx1(Max1, iv2(N, iv1(N))), AEx2(Max2, iv2(N2, iv1(N2)));
    int i,j,k,l,ij,kl,ii,jj,I1,I2,iiMax,J1,J2,jjMax,ik,jl,m,If0,mmax;
    double Vtmp,VS,VSS,htmp,hS;
    dv1 stmp(N0);
    MatrixXd H0tmp;
    bool flag;
    const double zero = 0., one = 1., two = 2., three = 3., four = 4.;

    for (i = 0; i < Max1; i++)  for (j = 0; j < N; j++) for (k = 0; k < N; k++)
        AEx1[i][j][k] = abs (Ex1[i][j][k]);
    for (i = 0; i < Max2; i++)  for (j = 0; j < N2; j++) for (k = 0; k < N2; k++)
        AEx2[i][j][k] = abs (Ex2[i][j][k]);
//  Reverse String Ordering
    for (m = 0; m < N0; m++)    RIstr[Istr[m]] = m;

//  Remove terms that are not in H0
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            for (ii = 0; ii < Max1; ii++)
            {
                If0 = false;
                for (m = 0; m < N0; m++)
                    if (Istr[m] + 1 == AEx1[ii][i][j])   If0 = true;
                if (!If0)   AEx1[ii][i][j] = 0;
            }

    for (i = 0; i < N2; i++)
        for (j = 0; j < N2; j++)
            for (ii = 0; ii < Max2; ii++)
            {
                If0 = false;
                for (m = 0; m < N0; m++)
                    if (Istr[m] + 1 == AEx2[ii][i][j])   If0 = true;
                if (!If0)   AEx2[ii][i][j] = 0;
            }

//  Check for Zero blocks in V
    for (i = 0; i < N; i++)
        for (k = 0; k < N; k++)
        {
            ik = cpind(i,k);    flag = true;
            for (j = 0; j < N && flag; j++)
                for (l = 0; l < N && flag; l++)
                {
                    jl = cpind(j,l);
                    if (fabs(V[cpind(ik,jl)]) > 1.E-10) flag = false;
                }
            if (flag)   ifzero[i][k] = 1;
        }

//  One Electron Part
    for (i = 0; i < N; i++) for (j = 0; j < N; j++)
    {
        htmp = h[cpind(i,j)];
        if (fabs (htmp) < 1.E-10)   continue;
        if (i == j) iiMax = _nchoosek_ (N - 1, No - 1);
        else    iiMax = _nchoosek_ (N - 2, No - 1);
        for (ii = 0; ii < iiMax; ii++)
        {
            I1 = AEx1[ii][i][j];    I2 = AEx1[ii][j][i];
            if (I1 == 0 || I2 == 0) continue;
            hS = (I1 == Ex1[ii][i][j]) ? (htmp) : (-htmp);
            for (m = 0; m < N0; m++)
                H0(m + RIstr[I2 - 1] * N0, m + RIstr[I1 - 1] * N0) += hS;
        }
    }
/*  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*   In Fortran, M(N*N, N*N) can also be viewd as Ms(N,N,N,N).
*   When viewd in this way, the conversion follows rule
*       Ms(a,b,c,d) --> M(a+b*N, c+d*N)
*   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */

//  Same Spin Two-Electron Part
    int IK, IL, JK, JL;
    ij = -1;
    for (i = 0; i < N; i++) for (j = i + 1; j < N; j++)
    {
        ij++;   kl = -1;
        for (k = 0; k < N; k++) for (l = k + 1; l < N; l++)
        {
            kl++;
            IK = cpind(i,k);    IL = cpind(i,l);
            JK = cpind(j,k);    JL = cpind(j,l);
            Vtmp = (V[cpind(IK,JL)] - V[cpind(IL,JK)] -
                V[cpind(JK,IL)] + V[cpind(JL,IK)]) / two;
            if (fabs (Vtmp) > 1.E-10)
            {
                if (i == k && j == l)   iiMax = _nchoosek_ (N - 2, No - 2);
                else if (i == k)    iiMax = _nchoosek_ (N - 3, No - 2);
                else if (i == l)    iiMax = _nchoosek_ (N - 3, No - 2);
                else if (j == k)    iiMax = _nchoosek_ (N - 3, No - 2);
                else if (j == l)    iiMax = _nchoosek_ (N - 3, No - 2);
                else    iiMax = _nchoosek_ (N - 4, No - 2);
                for (ii = 0; ii < iiMax; ii++)
                {
                    I1 = AEx2[ii][ij][kl];  I2 = AEx2[ii][kl][ij];
                    if (I1 == 0 || I2 == 0) continue;
                    VS = (I1 == Ex2[ii][ij][kl]) ? (Vtmp) : (-Vtmp);
                    for (m = 0; m < N0; m++)
                        H0(m + RIstr[I2 - 1] * N0, m + RIstr[I1 - 1] * N0) += VS;
                }
            }
        }
    }

//  Opposite Spin Two-Electron Part
    ik = -1;
    for (i = 0; i < N; i++) for (k = 0; k < N; k++)
    {
        ik++;   jl = -1;    IK = cpind(i,k);
        if (ifzero[i][k])   continue;
        if (i == k) iiMax = _nchoosek_ (N - 1, No - 1);
        else    iiMax = _nchoosek_ (N - 2, No - 1);

//  Gather together phases in Stmp
        stmp.assign (N0, 0);    IY.assign (N0, -1);
        for (ii = 0; ii < iiMax; ii++)
        {
            I1 = AEx1[ii][i][k];   I2 = AEx1[ii][k][i];
            if (I1 == 0 || I2 == 0) continue;
            VS = one;
            if (I1 != Ex1[ii][i][k])    VS = -VS;
            IY[RIstr[I1 - 1]] = RIstr[I2 - 1];
            stmp[RIstr[I1 - 1]] = VS;
        }
        mmax = 0;
        for (m = 0; m < N0; m++)
            if (IY[m] < 0) IY[m] = 0;
            else    mmax = m;

//  Collect Elements of H0
        for (j = 0; j < N; j++) for (l = 0; l < N; l++)
        {
            jl++;   JL = cpind(j,l);
            if (ik < jl)    continue;
            if (ik == jl)   Vtmp = V[cpind(IK,JL)] / two;
            if (ik > jl)    Vtmp = V[cpind(IK,JL)];
            if (fabs (Vtmp) < 1.E-10)   continue;

            if (j == l) jjMax = _nchoosek_ (N - 1, No - 1);
            else    jjMax = _nchoosek_ (N - 2, No - 1);

            for (jj = 0; jj < jjMax; jj++)
            {
                J1 = AEx1[jj][j][l];    J2 = AEx1[jj][l][j];
                if (J1 == 0 || J2 == 0) continue;
                VS = Vtmp;
                if (J1 != Ex1[jj][j][l])    VS = -VS;
                for (m = 0; m <= mmax; m++)
                {
                    H0(IY[m] + RIstr[J2 - 1] * N0, m + RIstr[J1 - 1] * N0) +=
                        VS * stmp[m];
                }
            }
        }
    }

    H0tmp.setZero (N0 * N0, N0 * N0);
    for (i = 0; i < N0; i++) for (j = 0; j < N0; j++)
        for (k = 0; k < N0; k++) for (l = 0; l < N0; l++)
        {
            H0tmp(i + j * N0, k + l * N0)
                = H0(i + j * N0, k + l * N0) +
                H0(j + i * N0, l + k * N0);
        }
    H0 = H0tmp;

    /*delete RIstr;   delete IY;
    for (i = 0; i < N; i++) delete [] ifzero[i];
    delete [] ifzero;
    for (i = 0; i < Max1; i++)
    {
        for (j = 0; j < N; j++) delete [] AEx1[i][j];
        delete [] AEx1[i];
    }
    delete [] AEx1;
    for (i = 0; i < Max2; i++)
    {
        for (j = 0; j < N2; j++) delete [] AEx2[i][j];
        delete [] AEx2[i];
    }
    delete [] AEx2;*/
}

double _MDOT_ (const MatrixXd& A, const MatrixXd& B)
{
    return (A.cwiseProduct (B)).sum ();
}

void _GS_ (MatrixXd& X, vMatrixXd Xi, int iS)
{
    double dum; int i;

//  Gramm-Schmidt Once
    for (i = 0; i < iS; i++)
    {
        dum = _MDOT_ (X, Xi[i]);
        X -= dum * Xi[i];
    }
    X.normalize (); // normalize w.r.t. its squared sum

//  Gramm-Schmidt Twice (for Stability)
    for (i = 0; i < iS; i++)
    {
        dum = _MDOT_ (X, Xi[i]);
        X -= dum * Xi[i];
    }
    X.normalize (); // normalize w.r.t. its squared sum
}

void _HX_ (int N, int No, int N2, int Max1, int Max2, int Nstr,
    const iv3& Ex1, const iv3& Ex2, const MatrixXd& X, double *h,
    double *V, MatrixXd& Y)
{
    int i,j,k,l,ij,kl,ii,jj,I1,I2,iiMax,J1,J2,jjMax,ik,jl;
    iv2 ifzero(N, iv1(N));  iv3 AEx1(Max1, iv2(N, iv1(N)));
    int IfSym;
    double Vtmp,VS,VSS,htmp,hS,Tmp,Spin;
    MatrixXd Xtmp, Ytmp;
    bool flag;
    const double zero = 0., one = 1., two = 2., three = 3., four = 4.;

    Y.setZero (Nstr, Nstr);
    for (i = 0; i < Max1; i++)  for (j = 0; j < N; j++) for (k = 0; k < N; k++)
        AEx1[i][j][k] = abs (Ex1[i][j][k]);

//  Check Spin Symmetry of X
    Spin = one; Y = X;  Y -= X.transpose ();    Tmp = Y.squaredNorm ();
    if (Tmp > 1.E-1)    Spin = -one;
    Y.setZero ();

//  Check for Zero Blocks in V
    for (i = 0; i < N; i++) for (k = 0; k < N; k++)
    {
        ik = cpind(i,k);    flag = true;
        for (j = 0; j < N && flag; j++) for (l = 0; l < N && flag; l++)
        {
            jl = cpind(j,l);
            if (fabs (V[cpind(ik,jl)]) > 1.E-10)    flag = false;
        }
        if (flag)   ifzero[i][k] = 1;
    }

//  Check Symmetry of V
// Don't need this in current version, as we always assume
//      <ij|kl> == <ji|lk>, i.e. e1 and e2 are exchangeable
/*    IfSym = 1;  flag = true;
    for (i = 0; i < N && flag; i++) for (k = 0; k < N && flag; k++)
    {
        ik = cpind(i,k);
        for (j = 0; j < N && flag; j++) for (l = 0; l < N && flag; l++)
        {
            jl = cpind(j,l);
            if (fabs (V[cpind(ik,jl)] - V[cpind(jl,ik)]) > 1.E-10)
            {   flag = false;   IfSym = 0;  }
        }
    }
*/

// One Electron Part
    for (i = 0; i < N; i++) for (j = 0; j < N; j++)
    {
        htmp = h[cpind(i,j)];
        if (fabs (htmp) < 1.E-10)   continue;
        if (i == j) iiMax = _nchoosek_ (N - 1, No - 1);
        else    iiMax = _nchoosek_ (N - 2, No - 1);
        for (ii = 0; ii < iiMax; ii++)
        {
            I1 = AEx1[ii][i][j];    I2 = AEx1[ii][j][i];
            hS = (I1 == Ex1[ii][i][j]) ? (htmp) : (-htmp);
            Y.col (I2 - 1) += X.col (I1 - 1) * hS;
        }
    }

//  Same Spin Two-Electron Part
    int IK, IL, JK, JL;
    ij = -1;
    for (i = 0; i < N; i++) for (j = i + 1; j < N; j++)
    {
        ij++;   kl = -1;
        for (k = 0; k < N; k++) for (l = k + 1; l < N; l++)
        {
            kl++;
            IK = cpind(i,k);    IL = cpind(i,l);
            JK = cpind(j,k);    JL = cpind(j,l);
            Vtmp = (V[cpind(IK,JL)] - V[cpind(IL,JK)] -
                V[cpind(JK,IL)] + V[cpind(JL,IK)]) / two;
            if (fabs (Vtmp) > 1.E-10)
            {
                if (i == k && j == l)   iiMax = _nchoosek_ (N - 2, No - 2);
                else if (i == k)    iiMax = _nchoosek_ (N - 3, No - 2);
                else if (i == l)    iiMax = _nchoosek_ (N - 3, No - 2);
                else if (j == k)    iiMax = _nchoosek_ (N - 3, No - 2);
                else if (j == l)    iiMax = _nchoosek_ (N - 3, No - 2);
                else    iiMax = _nchoosek_ (N - 4, No - 2);
                for (ii = 0; ii < iiMax; ii++)
                {
                    I1 = abs (Ex2[ii][ij][kl]);  I2 = abs (Ex2[ii][kl][ij]);
                    VS = (I1 == Ex2[ii][ij][kl]) ? (Vtmp) : (-Vtmp);
                    Y.col (I2 - 1) += VS * X.col (I1 - 1);
                }
            }
        }
    }

//  Opposite Spin Two-Electron Part
    ik = -1;
    Xtmp.setZero (Max1, Nstr);
    for (i = 0; i < N; i++) for (k = 0; k < N; k++)
    {
        ik++;   jl = -1;    IK = cpind(i,k);
        if (ifzero[i][k])   continue;
        if (i == k) iiMax = _nchoosek_ (N - 1, No - 1);
        else    iiMax = _nchoosek_ (N - 2, No - 1);

//  Gather together phases in Xtmp
        for (ii = 0; ii < iiMax; ii++)
        {
            I1 = AEx1[ii][i][k];
            VS = (I1 == Ex1[ii][i][k]) ? (one) : (-one);
            for (jj = 0; jj < Nstr; jj++)   Xtmp(ii, jj) = X(I1 - 1, jj) * VS;
        }

//  Collect Elements of Ytmp
        Ytmp.setZero (Max1, Nstr);
        for (j = 0; j < N; j++) for (l = 0; l < N; l++)
        {
            jl++;   JL = cpind(j,l);
            //if (IfSym)    // As mentioned above, sym is not needed
            if (ik < jl)    continue;
            if (ik == jl)   Vtmp = V[cpind(IK,JL)] / two;
            if (ik > jl)    Vtmp = V[cpind(IK,JL)];
            if (fabs (Vtmp) < 1.E-10)   continue;

            if (j == l) jjMax = _nchoosek_ (N - 1, No - 1);
            else    jjMax = _nchoosek_ (N - 2, No - 1);

            for (jj = 0; jj < jjMax; jj++)
            {
                J1 = AEx1[jj][j][l];    J2 = AEx1[jj][l][j];
                VS = (J1 == Ex1[jj][j][l]) ? (Vtmp) : (-Vtmp);
                Ytmp.col (J2 - 1) += Xtmp.col (J1 - 1) * VS;
            }
        }

//  Scatter Elements of Y
        for (ii = 0; ii < iiMax; ii++)
        {
            I1 = AEx1[ii][k][i];
            for (jj = 0; jj < Nstr; jj++)
                Y(I1 - 1, jj) += Ytmp(ii, jj);
        }
    }
//  Enforce MS = 0
    MatrixXd Ytemp = Spin * Y.transpose ();
    Y += Ytemp;

    /*for (i = 0; i < N; i++) delete [] ifzero[i];
    delete [] ifzero;
    for (i = 0; i < Max1; i++)
    {
        for (j = 0; j < N; j++) delete [] AEx1[i][j];
        delete [] AEx1[i];
    }
    delete [] AEx1;*/
}

void _RPDM_ (int N, int No, int N2, int Max1, int Max2, int Nstr,
    const iv3& Ex1, const iv3& Ex2, vMatrixXd& Xi, MatrixXd& P, double *P2)
{
    int i, j, k, l, ij, kl, ijkl, lenh = N*(N+1)/2, lenV = lenh*(lenh+1)/2;
    MatrixXd XH (Nstr, Nstr), X (Nstr, Nstr);
    double *T1 = _darray_gen_ (lenh), *T2 = _darray_gen_ (lenV);

//  One Particle Density Matrix (1PDM) for spin alpha
    P.setZero (N, N);   X = Xi[0];
    for (i = 0; i < N; i++) for (j = i; j < N; j++)
    {
        ij = cpind(i,j); T1[ij] += (i == j) ? (0.5) : (0.25);
        _HX_ (N, No, N2, Max1, Max2, Nstr, Ex1, Ex2, X, T1, T2, XH);
        P(i, j) = _MDOT_ (X, XH);   P(j, i) = P(i, j);
        T1[ij] = 0.;
    }

//  One Particle Density Matrix (2PDM)
    for (i = 0, ij = 0; i < N; i++) for (j = 0; j <= i; j++, ij++)
    for (k = 0, kl = 0; k < N; k++) for (l = 0; l <= k; l++, kl++)
    {
        if (ij < kl)    continue;
        ijkl = cpind(ij,kl);
        T2[ijkl] = 1.;
        _HX_ (N, No, N2, Max1, Max2, Nstr, Ex1, Ex2, X, T1, T2, XH);
        P2[ijkl] = _MDOT_ (X, XH);
        T2[ijkl] = 0.;
    }
    delete T1;  delete T2;
}

void _FCIman_ (int N, int No, int N0MAX, int NS, double *h, double *V,
    vMatrixXd& Xi, double *Ei, double *Sym, MatrixXd& P, double *P2)
/*********************************************************************
*   N       - Number of Basis Functions
*   No      - Number of Alpha electrons (=Number of Beta Electrons)
*   N2      - N * (N - 1) / 2
*   Nstr    - Number of Alpha CI Strings (=Number of Beta Strings)
*   Max1    - Max. # of states a 1 particle operator can connect to
*   Max2    - Max. # of states a 2 particle operator can connect to
*   N0      - Number of strings in preconditioner space
*   NS      - Number of eigenstates desired
*
*   h       - One electron Hamilontian
*   V       - Two electron Hamiltonian
*
*   Xi      - Input: NS Approximate guess eigenvectors
*           - Output: Lowest NS Eigenvectors
*   P       - One particle alpha transition density matrix
*   P2      - 2 particle density matrix. Currently only for ground state
*   Ei      - Lowest NS eigenvalues
*   H1s     - Full hamiltonian for 1-site embedding
**********************************************************************/
{
    const int N2 = _nchoosek_ (N, 2);
    const int Nstr = _nchoosek_ (N, No);
    const int N0 = (N0MAX < Nstr) ? (N0MAX) : (Nstr);
    const int Max1 = _nchoosek_ (N - 1, No - 1);
    const int Max2 = _nchoosek_ (N - 2, No - 2);
    const int lenh = N * (N + 1) / 2;
    int i,j,k,l,m,a,b,ij,kl,ab,ii,jj,kk,ll,ierr,iter,Info;
    int iX,iS,iSpin,iSym;
    iv1 Iocc(No), Isubst(No), Istr(N0);
    iv3 Ex1(Max1, iv2(N, iv1(N))), Ex2(Max2, iv2(N2, iv1(N2)));
    double Energy, DE, fac, norm, eps, *Uncertainty;
    MatrixXd Hd, H0, U0, X, X1, XH, X1H, Xtmp;
    VectorXd E0;    Matrix2d Hm, U; Vector2d Eig;
    const double zero = 0., one = 1., two = 2., three = 3., four = 4.;

    Energy = zero; fac = one;
//  Build indexing array for future use
    Zindex = _iarray2_gen_ (N, N);
    for (k = 0; k < No; k++)
        for (l = k; l <= N - No + k; l++)
            if (k == No - 1)    Zindex[k][l] = l - k;
            else
                for (m = N - l; m < N - k; m++)
                    Zindex[k][l] += _nchoosek_ (m , No - k - 1) -
                        _nchoosek_ (m - 1, No - k - 2);

//  Determine which strings are connected by various operators
    //Ex1 = _iarray3_gen_ (Max1, N, N);
    //Ex2 = _iarray3_gen_ (Max2, N2, N2);
    _IString_ (N, No, N2, Max1, Max2, Ex1, Ex2);

//  Build Diagonal part of H
    Hd.setZero (Nstr, Nstr);
    //Iocc = _iarray_gen_ (No);   Isubst = _iarray_gen_ (No);
    _GetHd_ (N, No, Nstr, h, V, Iocc, Isubst, Hd, 1);

//  Get N0 lowest strings
    //Istr = _iarray_gen_ (N0);
    _GetIstr_ (N, No, Nstr, N0, Hd, Istr);
    _Isort_ (N0, Istr, &Info);

//  Build + Diagonalize H0
    H0.setZero (N0 * N0, N0 * N0);
    _GetH0_ (N, No, N0, N2, Max1, Max2, Nstr, Ex1, Ex2,
        Istr, h, V, H0);
    _eigh_ (H0, U0, E0);

//  Big loop over states
    iX = -1;    Uncertainty = _darray_gen_ (NS);
    for (iS = 0; iS < NS; iS++)
//  Initial vector for this state (ensure it is a singlet for now)
//  This restarts from the input Xi vector if Xi is nonzero
    {
        X = Xi[iS]; iSpin = 1;  norm = X.squaredNorm ();
        if (norm < 1.E-2)
        {
            iSpin = -1;
            while (iSpin == -1)
            {
                iX++;
                X.setZero (Nstr, Nstr); ij = -1;
                for (i = 0; i < N0; i++)    for (j = 0; j < N0; j++)
                {
//  Initialize X with eigenvectors of H0
                    ij++;   X(Istr[j], Istr[i]) = U0(ij, iX);
                }
                X1 = X; X1 -= X.transpose ();
                norm = X1.squaredNorm () / X.squaredNorm ();
                if (norm < 1.E-2)   iSpin = 1;
            }
        }

//  Check if the state has even (+1) or odd (-1) S
//  Even is singlet, quintet...  | Odd is triplet, sestet...
        X1 = X; X1 -= X.transpose ();
        norm = X1.squaredNorm ();
        if (norm <= 1.E-1)  iSym = 1;   else    iSym = -1;
        Sym[iS] = iSym;

//  Fix up broken spin symmetry (if it exists)
        _GS_ (X, Xi, iS);

//  Compute Initial guess energy
        XH.setZero (Nstr, Nstr);
        _HX_ (N, No, N2, Max1, Max2, Nstr, Ex1, Ex2, X, h, V, XH);
        Energy = _MDOT_ (X, XH);
        //printf ("Guess energy: %18.16f\n", Energy);

//  Iteratively determine eigenvalue #iS of the Hamiltonian
        fac = 1.;  iter = 0;
        while (fabs (fac) > THRESH && iter < MAXITER)
        {
            iter++; DE = Energy;
//  Make the (orthogonal component of the) Davidson Update
            X1 = -(XH - Energy * X).cwiseQuotient (Hd -
                MatrixXd::Constant (Nstr, Nstr, Energy));

//  Build (H0-Energy)^-1 (excluding eigenvalues that might blow up)
            H0.setZero ();
            for (k = 0; k < N0 * N0; k++)
            {
                fac = E0(k) - Energy;
                if (fabs (fac) > 1.E-2)  fac = 1. / fac;
                else    fac = 0.;
                for (i = 0; i < N0 * N0; i++)   for (j = 0; j < N0 * N0; j++)
                    H0(i, j) += U0(i, k) * U0(j, k) * fac;
            }

//  Build the Davidson Update using H0
            ij = -1;
            for (i = 0; i < N0; i++)    for (j = 0; j < N0; j++)
            {
                ij++;   kl = -1;
                X1(Istr[j], Istr[i]) = zero;
                for (k = 0; k < N0; k++)    for (l = 0; l < N0; l++)
                {
                    kl++;
                    X1(Istr[j], Istr[i]) -= H0(ij, kl) *
                        (XH(Istr[l], Istr[k]) - Energy * X(Istr[l], Istr[k]));
                }
            }

//  Apply Spin Symmetry and Gramm-Schmidt to Update
            Xtmp = iSym * X1.transpose ();  X1 = 0.5 * (X1 + Xtmp);
            X1.normalize ();
            norm = _MDOT_ (X1, X);  X1 -= norm * X; X1.normalize ();
            norm = _MDOT_ (X1, X);  X1 -= norm * X; X1.normalize ();

//  Correct if Davidson has given us a bad vector.
//  X1 should not be orthogonal to (H-E)X
//  If it is (nearly) orthogonal, add a bit of (H-E)X to X1
//  If we don't do this, it occasionally gives false convergence
//  This should happen rarely and eps controls how often
//  this correction is invoked.
//  A more elegant fix might be to just have a better
//  preconditioner.
            X1H = XH - Energy * X;  _GS_ (X1H, Xi, iS); X1H.normalize ();
            eps = 1.E-1;    fac = fabs (_MDOT_ (X1, X1H));
            if (fac < eps)
            {
                X1 += 2. * eps * X1H;
                norm = _MDOT_ (X1, X);  X1 -= norm * X; X1.normalize ();
                norm = _MDOT_ (X1, X);  X1 -= norm * X; X1.normalize ();
            }

//  Make X1 orthogonal to lower eigenvectors
            _GS_ (X1, Xi, iS);

//  Act H on Davidson Vector
            _HX_ (N, No, N2, Max1, Max2, Nstr, Ex1, Ex2, X1, h, V, X1H);

//  Build Hm
            Hm(0, 0) = _MDOT_ (X, XH);
            Hm(1, 1) = _MDOT_ (X1, X1H);
            Hm(0, 1) = Hm(1, 0) = _MDOT_ (X, X1H);

//  Diagonalize Hm
            _eigh2_ (Hm, U, Eig);

//  Keep Lowest Eigenvector
            if (iter < 50 || iter % 10)
            {
                fac = U(1, 0);
                X = U(0, 0) * X + U(1, 0) * X1;
                XH = U(0, 0) * XH + U(1, 0) * X1H;
                //if (iter == 51) printf ("%18.16e  %18.16e\n", U(0,0),U(1,0));
            }
            else
//  If convergence is slow, sometimes it is because we are
//  Swiching back and forth between two near-solutions. To
//  fix that, every once in a while, we take half the step.
//  A more elegant fix might be to use more than one Davidson
//  vectors in the space.
            {
                fac = fabs(0.5 * U(1, 0));    norm = one / sqrt (one + fac * fac);
                X = norm * (X + fac * X1);
                XH = norm * (XH + fac * X1H);
            }
//  Normalize X and XH
            norm = X.norm ();   X /= norm;  XH /= norm;

//  Avoid accumulating roundoff error
            if (iter % 4 == 0)
                _HX_ (N, No, N2, Max1, Max2, Nstr, Ex1, Ex2, X, h, V, XH);

            Energy = _MDOT_ (X, XH) / X.squaredNorm ();

//  End of the Loop for FCI iteration
        }
//        printf ("Done after %4d iterations.\n", iter);
//  Energy, Uncertainty for THIS State
        _HX_ (N, No, N2, Max1, Max2, Nstr, Ex1, Ex2, X, h, V, XH);
        Energy = _MDOT_(X,XH)/X.squaredNorm ();
        Uncertainty[iS] = _MDOT_(XH,XH)/_MDOT_(XH,X)/_MDOT_ (XH,X) - 1.;
        Xi[iS] = X; Ei[iS] = Energy;
//  End of the BIG Loop over States
    }
/*    cout << "=============================\n";
    cout << "|        FCI Summary        |\n";
    cout << "=============================\n";
    printf ("State\tEnergy\t\t\tSpin Sym\tUncertainty\n");
    for (i = 0; i < NS; i++)
        printf ("%2d\t%18.16f\t%7.5f\t\t%18.16e\n", i + 1, Ei[i], Sym[i], Uncertainty[i]);
*/
//  Get the (Ground Staate) 1PDM and 2PDM
    _RPDM_ (N, No, N2, Max1, Max2, Nstr, Ex1, Ex2, Xi, P, P2);

    for (i = 0; i < N; i++) delete [] Zindex[i];
    delete [] Zindex;
}

#undef THRESH
#undef MAXITER

void _frank_hack_ ()
{

}

#endif
