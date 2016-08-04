#ifndef _HF_INCLUDED_
#define _HF_INCLUDED_

#include <Eigen/Dense>

using namespace Eigen;

#define SCF_ITER 1E3
#define SCF_CONV 1E-6
#define MAX_DIIS 10
#define PINV_TOL 1E-10
#define LINSOLVER_TOL 1E-4
#define mixing_beta_HL 0.2

//typedef Matrix<double, K, K> MatrixKd;
//typedef Matrix<double, K, 1> VectorKd;
#define esXd SelfAdjointEigenSolver<MatrixXd> 

#endif
