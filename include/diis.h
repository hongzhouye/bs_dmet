#ifndef _DIIS_H_INCLUDED_
#define _DIIS_H_INCLUDED_

#include <Eigen/Dense>
#include <iostream>
#include "hf.h"
#include "matdeque.h"

using namespace std;
using namespace Eigen;

// MP pseudo-inverse solver for a real symmetric matrix B
// by first eigenvalue-decomposing it into U * S * U'
// then deleting small eigen-components i if
//
//			S(i, i) < PINV_TOL * S_max
// then,
//
//			 pinv(B) = u * (s)^(-1) * u'
//
// where u and s are truncated U and S.
MatrixXd _pinv_ (MatrixXd& B)
{
	int i, j;
	SelfAdjointEigenSolver<MatrixXd> es;
	es.compute(B);
	MatrixXd U0 = es.eigenvectors();
	VectorXd S0 = es.eigenvalues();
	//cout << "eigenvalues:" << endl << S0 << endl << endl;
	
	double max, min, temp;
	int flag = 1, nrow = S0.size(), last;
	VectorXd u_tmp (nrow);
	while (flag)
	{
		VectorXd::Index maxS, minS;
		max = S0.array().abs().matrix().maxCoeff (&maxS);
		min = S0.array().abs().matrix().minCoeff (&minS);
		last = S0.size() - 1;
		//cout << "(maxS, minS) = (" << maxS << ", " << minS << ")\n\n";
		if (min / max > PINV_TOL)	flag = 0;
		else if (minS != last)
		{
			temp = S0(minS);
			S0(minS) = S0(last);
			S0.conservativeResize(last);
			u_tmp = U0.col(minS);
			U0.col(minS) = U0.col(last);
			U0.conservativeResize(nrow, last);
		}
		else
		{
			S0.conservativeResize(last);
			U0.conservativeResize(nrow, last);
		}
	}
	//cout << "U0:" << endl << U0 << endl << endl;
	//cout << "S0:" << endl << S0 << endl << endl;

	MatrixXd Bi = U0 * (1. / S0.array()).matrix().asDiagonal() * U0.transpose();
	return Bi;
}

// A combination of four common linear solvers, solving
//
//			B * x = b
// for x and returning the relative error |B * x - b| / |b|.
double _lin_solver_ (MatrixXd& B, const MatrixXd& b, 
		VectorXd& x, string method)
{
	if (method == "cpqr")
		x = B.colPivHouseholderQr().solve(b);
	else if (method == "fpqr")
		x = B.fullPivHouseholderQr().solve(b);
	else if (method == "fplu")
		x = B.fullPivLu().solve(b);
	else if (method == "pinv")
	{
		MatrixXd Bi = _pinv_ (B);
		x = Bi * b;
		//cout << "Bi" << endl << Bi << endl << endl;
		//cout << "x" << endl << x << endl << endl;
		//cout << "error = " << (B * x - b).norm() / b.norm() << "\n\n";
	}
	else
	{
		cout << "Error when calling Linear Solver!\n";
		exit (1);
	}
	return (B * x - b).norm() / b.norm();
}

// diis
void _next_diis_ (Matdeq& errs, Matdeq& Fs, 
		Matdeq& errbs, Matdeq& Fbs)
{
	int n = (Fs.full_flag) ? (Fs.max_size) : (Fs.now);	
	double error;

	// Construct B matrix and b rhs vector (Helgaker 10.6.28)
	MatrixXd B = MatrixXd::Constant (n+1, n+1, -1.);
	int i, j;
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
			B(i, j) = (errs.element[i].cwiseProduct(errs.element[j])).sum()
					+ (errbs.element[i].cwiseProduct(errbs.element[j])).sum();
	B(n, n) = 0;
	VectorXd b = VectorXd::Constant (n + 1, 0);
	b(n) = -1.;

	// Solve B * x = b for x
	VectorXd x;
	error = _lin_solver_ (B, b, x, "cpqr");
	if (error > LINSOLVER_TOL)
	{	
		cout << "CPQR error = " << error << endl;
		cout << "CPQR fails, switch to FPQR ..." << endl;
		error = _lin_solver_ (B, b, x, "fpqr");
		if (error > LINSOLVER_TOL)
		{	
			cout << "FPQR error = " << error << endl;
			cout << "FPQR fails, switch to FPLU ..." << endl;
			error = _lin_solver_ (B, b, x, "fplu");
			if (error > LINSOLVER_TOL)
			{	
				cout << "FPLU error = " << error << endl;
				cout << "FPLU fails, switch to PINV ..." << endl;
				error = _lin_solver_ (B, b, x, "pinv");
				if (error > LINSOLVER_TOL)
				{
					cout << "PINV error = " << error << endl;
					cout << "B matrix is drastically singular!\n" <<
						"Linear Solver does not work!\n";
					exit (1);
				}
			}
		}
	}

	// Construct new F
	x.conservativeResize(n);
	Fs.element[Fs.now].setZero ();
	Fbs.element[Fbs.now].setZero ();
	for (i = 0; i < n; i++)
	{
		Fs.element[Fs.now] += x(i) * Fs.element[i];
		Fbs.element[Fbs.now] += x(i) * Fbs.element[i];
	}
}

#endif
