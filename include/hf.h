#ifndef _HF_INCLUDED_
#define _HF_INCLUDED_

#include <Eigen/Dense>

using namespace Eigen;

#define Nup 5
#define Ndn 5
#define K 10 

#define SCF_ITER 1E3
#define SCF_CONV 1E-5
#define MAX_DIIS 10
#define PINV_TOL 1E-10
#define LINSOLVER_TOL 1E-4

typedef Matrix<double, K, K> MatrixKd;
typedef Matrix<double, K, 1> VectorKd;
#define mapKd(p) Map<MatrixKd>(p)
#define esKd SelfAdjointEigenSolver<MatrixKd> 

#endif
