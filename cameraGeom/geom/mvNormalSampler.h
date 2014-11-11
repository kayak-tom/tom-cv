/*
 * mvNormalSampler.h
 *
 *  Created on: 24 Jan 2011
 *      Author: tom
 */

#ifndef MVNORMALSAMPLER_H_
#define MVNORMALSAMPLER_H_

#include <util/exception.h>
#include <util/optimisation_attributes.h>
#include <Eigen/Core>

template<class MatrixType>
void checkIsPosDef(const MatrixType & information) HOT;

template<class MatrixType>
MatrixType matrixSqrt(const MatrixType & cov) HOT;

template<class MatrixType>
class CMVNormalSampler
{
public:
    typedef Eigen::Matrix<double, MatrixType::RowsAtCompileTime, 1> VecType;
private:
    const VecType mean;
    VecType mvNormal;
    MatrixType covSqrt;
public:
    CMVNormalSampler(const VecType & mean, const MatrixType & cov) HOT;

    VecType sample() HOT;
};

#endif /* MVNORMALSAMPLER_H_ */
