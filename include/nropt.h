#ifndef _BFGS_H_INCLUDED_
#define _BFGS_H_INCLUDED_

#include "hf.h"
#include "frag.h"
#include "brent.h"

#define LINSOL_TOL1 1.E-8
#define LINSOL_TOL2 1.E-4
#define BRENT_TOL 1.E-4
#define DX_TOL 1.E-9
#define FX_TOL 1.E-5
#define MAX_STEP 10.
#define JAC_EPS 1.E-4

MatrixXd _fd_Jac_ (VectorXd (*func) (VectorXd&, FRAG&), VectorXd& u, FRAG& frag)
{
    int n = frag.Nopt;
    MatrixXd J; J.setZero (n, n);
    VectorXd d, ff, fff, fb, fbb;
    for (int i = 0; i < n; i++)
    {
        VectorXd tempu;
        d.setZero (n);
        d(i) = 1.;
        tempu = u + JAC_EPS * d;
        ff  = (*func) (tempu, frag);
        tempu = u + 2. * JAC_EPS * d;
        fff  = (*func) (tempu, frag);
        tempu = u - JAC_EPS * d;
        fb  = (*func) (tempu, frag);
        tempu = u - 2. * JAC_EPS * d;
        fbb  = (*func) (tempu, frag);
        J.col(i) = (-fff + 8. * ff - 8. * fb + fbb) / (12. * JAC_EPS);
        // DEBUG
        /*printf ("u[%d] = %10.7f\n", i, u(i));
        cout << "ff = " << ff.transpose ().format (Short) << endl;
        cout << "fff = " << fff.transpose ().format (Short) << endl;
        cout << "fb = " << fb.transpose ().format (Short) << endl;
        cout << "fbb = " << fbb.transpose ().format (Short) << endl;
        cout << "J(:, i) = " << J.col(i).transpose ().format (Short) << endl;*/
    }
    return J;
}

VectorXd _bfgs_opt_ (VectorXd (*func) (VectorXd&, FRAG&), VectorXd& u,
    FRAG& frag)
{
    bool converge = false;
    VectorXd fx;    fx.setZero (frag.Nopt);
    VectorXd dx;    fx.setZero (frag.Nopt);
    int iter = 0, siter;
    double scale, obj_norm;
    MatrixXd J;

    clock_t t1, t2, dt_Jac = 0, dt_brent = 0;

    while (converge != true)
    {
        iter++;
        fx = (*func) (u, frag);
        // switch FCI guess mode to read
        //frag.dfci.mode = "read";
        frag.dfci.guess_read = true;
        t1 = clock ();
        J = _fd_Jac_ (func, u, frag);
        t2 = clock ();  dt_Jac += t2 - t1;

        dx = J.colPivHouseholderQr (). solve (-fx);
        double linsol_err = (J * dx + fx).norm ();
        if (linsol_err > LINSOL_TOL2)
        {
            cout << "Error in _bfgs_opt_: linear solver fails.\n\n";
            exit (1);
        }
        else if (linsol_err > LINSOL_TOL1)
            cout << "Warning: accuracy of the linsolver is low.\n\n";

        t1 = clock ();
        obj_norm = _brent_ (-1., 1., 2., func, &siter, BRENT_TOL,
            &scale, u, dx, frag);
        t2 = clock ();  dt_brent += t2 - t1;

        if (scale < 0)  cout << "Warning: stepped backwards.\n\n";
        //frag.dfci.guess_read = false;

        if (dx.norm () < DX_TOL || fx.norm () < FX_TOL) converge = true;

        if (scale * dx.norm () > MAX_STEP)
        {
            cout << "Warning: limiting step size to max_step\n\n";
            scale = MAX_STEP / dx.norm ();
        }

        u += scale * dx;

        printf ("Iteration: %d\n", iter);
        // DEBUG
        //cout << "Jacobian:\n" << J << "\n\n";
        cout << "x:  " << u.transpose ().format (Short) << "\t" << u.norm () << endl;
        cout << "dx: " << dx.transpose ().format (Short) << "\t" << dx.norm () << endl;
        cout << "fx: " << fx.transpose ().format (Short) << "\t" << fx.norm () << endl;
        printf ("scale: %10.7f\n", scale);
        printf ("scale iterations: %d\n\n", siter); fflush (0);
    }

    cout << "=========================\n";
    cout << "|     TIME SUMMARY      |\n";
    cout << "=========================\n";
    cout << "JAC Time: " << dt_Jac /
        (double)(CLOCKS_PER_SEC / 1000) << " ms" << endl;
    cout << "BRT Time: " << dt_brent /
        (double)(CLOCKS_PER_SEC / 1000) << " ms" << endl;

    return u;
}

#undef LINSOL_TOL1
#undef LINSOL_TOL2
#undef BRENT_TOL
#undef DX_TOL
#undef FX_TOL
#undef MAX_STEP
#undef JAC_EPS

#endif
