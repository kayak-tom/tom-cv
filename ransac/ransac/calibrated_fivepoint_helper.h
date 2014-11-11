/* This code is for academic use only. More details at http://www.hilandtom.com/tombotterill/code */

#pragma once

//#include <lapackpp.h> Lapack++ implementation is much slower (3x)
#include <Eigen/Core>

//void fivepointSetupMatricesHelper(const LaGenMatDouble & Emat, LaGenMatDouble & Amat); Lapack++ implementation is much slower (3x)
void fivepointSetupMatricesHelper(const Eigen::Matrix<double, 9, 4, 0, 9, 4> & Emat, Eigen::Matrix<double, 10, 20, 0, 10, 20> & Amat);
