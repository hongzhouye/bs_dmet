#ifndef _BRENT_H_INCLUDED_
#define _BRENT_H_INCLUDED_

#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include "hf.h"
#include "frag.h"

using namespace std;
using namespace Eigen;

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

double _obj_func_ (VectorXd (*func) (VectorXd&, FRAG&),
	VectorXd& u, FRAG& frag)
{
	VectorXd obj = (*func) (u, frag);
	return obj.norm ();
}

double _brent_ (double ax, double bx, double cx,
	VectorXd (*func) (VectorXd&, FRAG&), int *itr, double tol,
	double *xmin, VectorXd& u0, VectorXd& dx, FRAG& frag)
{
	int iter;
	double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
	double e = 0.0;

	a = (ax < cx ? ax : cx);
	b = (ax > cx ? ax : cx);
	x = w = v = bx;
	VectorXd tempu = u0 + x * dx;
	fw = fv = fx = _obj_func_ (func, tempu, frag);
	for (iter = 1; iter <= ITMAX; iter++) {
		*itr = iter;
		xm = 0.5 * (a + b);
		tol2 = 2.0 * (tol1 = tol * fabs (x) + ZEPS);
		if (fabs (x - xm) <= (tol2 - 0.5 * (b - a))) {
			*xmin = x;
			return fx;
		}
		if (fabs(e) > tol1) {
			r = (x - w) * (fx - fv);
			q = (x - v) * (fx - fw);
			p = (x - v) * q - (x - w) * r;
			q = 2.0 * (q - r);
			if (q > 0.0) p = -p;
			q = fabs (q);
			etemp = e;
			e = d;
			if (fabs (p) >= fabs (0.5 * q * etemp)
				|| p <= q * (a - x) || p >= q * (b - x))
				d = CGOLD * (e = (x >= xm ? a - x : b - x));
			else {
				d = p / q;
				u = x + d;
				if (u - a < tol2 || b - u < tol2)
					d = SIGN(tol1,xm - x);
			}
		} else {
			d = CGOLD * (e = (x >= xm ? a - x : b - x));
		}
		u = (fabs (d) >= tol1 ? x+d : x + SIGN(tol1,d));
		tempu = u0 + u * dx;
		fu = _obj_func_ (func, tempu, frag);
		if (fu <= fx) {
			if (u >= x) a = x; else b = x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a = u; else b = u;
			if (fu <= fw || w == x) {
				v = w;
				w = u;
				fv = fw;
				fw = fu;
			} else if (fu <= fv || v == x || v == w) {
				v = u;
				fv = fu;
			}
		}
	}
	cout << "Too many iterations in brent";
	*xmin = x;
	return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef SIGN
/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */

#endif
