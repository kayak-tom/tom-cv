/*
 * alignPoints.cpp
 *
 *  Created on: 10/12/2009
 *      Author: tom
 */

#include "alignPoints.h"
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/LU>

#include "util/opencv.h"

#include <boost/smart_ptr.hpp>

using namespace std;
using namespace Eigen;

//Also duplicated in scaleOptimiser.cpp
void svd(Eigen::Matrix<double, 4, Dynamic, RowMajor + AutoAlign> & A, Eigen::Matrix<double, 4, 4, RowMajor + AutoAlign> & U, Eigen::Matrix<double, 4, 1, ColMajor + AutoAlign> & D, Eigen::Matrix<double, Dynamic, Dynamic, RowMajor + AutoAlign> & V)
{
    CvMat cvA = cvMat((int)A.rows(), (int)A.cols(), CV_64F, A.data());
    CvMat cvU = cvMat((int)U.rows(), (int)U.cols(), CV_64F, U.data());
    CvMat cvV = cvMat((int)V.rows(), (int)V.cols(), CV_64F, V.data());
    CvMat cvD = cvMat((int)D.rows(), (int)D.cols(), CV_64F, D.data());
    cvSVD(&cvA, &cvD, &cvU, &cvV, CV_SVD_MODIFY_A);
}

void svd(Eigen::Matrix<double, 3, Dynamic, RowMajor + AutoAlign> & A, Eigen::Matrix<double, 3, 3, RowMajor + AutoAlign> & U, Eigen::Matrix<double, 3, 1, ColMajor + AutoAlign> & D, Eigen::Matrix<double, Dynamic, Dynamic, RowMajor + AutoAlign> & V)
{
    CvMat cvA = cvMat((int)A.rows(), (int)A.cols(), CV_64F, A.data());
    CvMat cvU = cvMat((int)U.rows(), (int)U.cols(), CV_64F, U.data());
    CvMat cvV = cvMat((int)V.rows(), (int)V.cols(), CV_64F, V.data());
    CvMat cvD = cvMat((int)D.rows(), (int)D.cols(), CV_64F, D.data());
    cvSVD(&cvA, &cvD, &cvU, &cvV /*, CV_SVD_MODIFY_A*/);
}

CAlignPoints::CAlignPoints()
{
    // TODO Auto-generated constructor stub

}

CAlignPoints::~CAlignPoints()
{
    // TODO Auto-generated destructor stub
}

bool CAlignPointsUmeyama::alignPoints(const T3dPointMatchVector & vPointMatches, C3dRotation & R, C3dPoint & t, double & s)
{
    // OK but needs NEWMAT or re-writing
    // DEBUGONLY(testForPlanarPoints(vPointMatches));

    C3dPoint mean1, mean2;

    const double sizeInv = 1.0 / vPointMatches.size();

    for(T3dPointMatchVector::const_iterator pPointPair = vPointMatches.begin(); pPointPair != vPointMatches.end(); pPointPair++)
    {
        const C3dPoint & p1 = pPointPair->p1();
        const C3dPoint & p2 = pPointPair->p2();
        mean1 += p1;
        mean2 += p2;
    }
    mean1 *= sizeInv;
    mean2 *= sizeInv;

    double var1 = 0, var2 = 0;
    Matrix3d Sigma;
    Sigma.setZero();

    for(T3dPointMatchVector::const_iterator pPointPair = vPointMatches.begin(); pPointPair != vPointMatches.end(); pPointPair++)
    {
        const C3dPoint & centredPoint1 = (pPointPair->p1() - mean1);
        const C3dPoint & centredPoint2 = (pPointPair->p2() - mean2);
        var1 += centredPoint1.sum_square();
        var2 += centredPoint2.sum_square();

        Vector3d v1,v2;
        centredPoint1.asVector(v1);
        centredPoint2.asVector(v2);

        Sigma += v2 * v1.transpose();
    }

    var1 *= sizeInv;
    var2 *= sizeInv;
    Sigma *= sizeInv;

    cout << Sigma<< "=Sigma" << endl;

    Matrix3d S; S.setIdentity();
    if(Sigma.determinant() <= 0){
        //return false;
        cout << "Warning: getTandR: large errors (neg determinant)...";
        S(2, 2) = -1;
    }
#if EIGEN_VERSION_AT_LEAST(2,90,0)
    JacobiSVD< Eigen::Matrix3d > svd(Sigma, ComputeFullU | ComputeFullV);
#else
    SVD< Eigen::Matrix3d > svd(Sigma);
#endif

    //cout << svd.matrixU() << "=U\n";
    //cout << svd.matrixV() << "=V\n";

    const Matrix3d & R_temp = svd.matrixU() * S * svd.matrixV().transpose();

    /*{
        Matrix3d diag; diag.setZero(); diag.diagonal() = svd.singularValues();
        cout << svd.matrixU() * diag * svd.matrixV().transpose() << "= sigma\n";
    }*/

    if(!isRotMat(R_temp)) return false;

    R=R_temp;
    //Was c = (D * S).trace() / var1; //Done: swap to get rid of 1 div
    //c = 1. / c;
    Vector3d SVs = S*svd.singularValues();
    //cout << SVs << "=sv's" << endl;

    double dTrace = SVs.sum();
    if(zero(dTrace)) return false;
    s = var1 / dTrace;
    t = mean2 * s - R * mean1;
    //PRINTMAT(mean2)
    //PRINTMAT(R * mean1)
    return true;
}

bool CAlignPointsUmeyama_Weighted::alignPoints(const T3dPointMatchVector & vPointMatches, C3dRotation & R, C3dPoint & t, double & s)
{
    C3dPoint mean1, mean2;

    int nSize = vPointMatches.size();

    Matrix<double, 3, Dynamic, RowMajor + AutoAlign> X(3, nSize);
    Matrix<double, 3, Dynamic> Y(3, nSize);

    int i=0;
    double dTotalLambda = 0;
    ARRAY(double, adLambda, nSize);
    for(T3dPointMatchVector::const_iterator pPointPair = vPointMatches.begin(); pPointPair != vPointMatches.end(); pPointPair++, i++)
    {
        const C3dPoint & centredPoint1 = (pPointPair->p1());
        const C3dPoint & centredPoint2 = (pPointPair->p2());

        CCamera P;
        //adLambda[i] = 1;
        adLambda[i] = 1.0/max<double>(centredPoint1.depth(P), centredPoint2.depth(P)); //Max better than min
        adLambda[i] = sqr(adLambda[i]);
        //adLambda[i] = 2/(centredPoint1.depth(P) + centredPoint2.depth(P)); //Max better than min
        //cout << centredPoint1.depth(P) << ' ' << centredPoint2.depth(P) << endl;
        //const double dLambda_inv = 1.0/adLambda[i];
        dTotalLambda += adLambda[i];
        mean1 += adLambda[i]*centredPoint1;
        mean2 += adLambda[i]*centredPoint2;
    }

    mean1 /= dTotalLambda;
    mean2 /= dTotalLambda;

    i=0;
    double var1 = 0, var2 = 0;
    Matrix3d Sigma;
    Sigma.setZero();

    for(T3dPointMatchVector::const_iterator pPointPair = vPointMatches.begin(); pPointPair != vPointMatches.end(); pPointPair++, i++)
    {
        const C3dPoint & centredPoint1 = (pPointPair->p1() - mean1);
        const C3dPoint & centredPoint2 = (pPointPair->p2() - mean2);

        Vector3d v1,v2;
        centredPoint1.asVector(v1);
        centredPoint2.asVector(v2);
        v1 *= adLambda[i];
        v2 *= adLambda[i];

        Sigma += v2 * v1.transpose();

        var1 += centredPoint1.sum_square();
        var2 += centredPoint2.sum_square();
    }
    var1 /= dTotalLambda;
    var2 /= dTotalLambda;

    Sigma /= dTotalLambda;

    cout << Sigma<< "=Sigma" << endl;

    Matrix3d S; S.setIdentity();
    if(Sigma.determinant() <= 0){
        //return false;
        cout << "Warning: getTandR: large errors (neg determinant)...";
        S(2, 2) = -1;
    }

#if EIGEN_VERSION_AT_LEAST(2,90,0)
    JacobiSVD< Eigen::Matrix3d > svd(Sigma, ComputeFullU | ComputeFullV);
#else
    SVD< Eigen::Matrix3d > svd(Sigma);
#endif

    const Matrix3d & R_temp = svd.matrixU() * S * svd.matrixV().transpose();

    if(!isRotMat(R_temp)) return false;

    R=R_temp;
    //Was c = (D * S).trace() / var1; //Done: swap to get rid of 1 div
    //c = 1. / c;
    Vector3d SVs = S*svd.singularValues();
    //cout << SVs << "=sv's" << endl;

    double dTrace = SVs.sum();
    if(zero(dTrace)) return false;
    s = var1 / dTrace;

    s=1.0;

    t = mean2 * s - R * mean1;
    //PRINTMAT(mean2)
    //PRINTMAT(R * mean1)
    return true;
}

#ifndef NEWER_NVG
bool CAlignPointsNew::alignPoints(const T3dPointMatchVector & vPointMatches, C3dRotation & R, C3dPoint & t, double & s)
{
    C3dPoint mean1, mean2;

    int nSize = vPointMatches.size();

    //const double sizeInv = 1.0 / nSize;

    Matrix<double, 3, Dynamic, RowMajor + AutoAlign> X(3, nSize);
    Matrix<double, 3, Dynamic> Y(3, nSize);

    /*for(T3dPointMatchVector::const_iterator pPointPair = vPointMatches.begin(); pPointPair != vPointMatches.end(); pPointPair++)
    {
        const C3dPoint & p1 = pPointPair->p1();
        const C3dPoint & p2 = pPointPair->p2();
        mean1 += p1;
        mean2 += p2;
    }
    mean1 *= sizeInv;
    mean2 *= sizeInv;*/

    int i=0;
    double dTotalLambda = 0;
    ARRAY(double, adLambda, nSize);
    for(T3dPointMatchVector::const_iterator pPointPair = vPointMatches.begin(); pPointPair != vPointMatches.end(); pPointPair++, i++)
    {
        const C3dPoint & centredPoint1 = (pPointPair->p1());
        const C3dPoint & centredPoint2 = (pPointPair->p2());

        CCamera P;
        //adLambda[i] = 1;
        adLambda[i] = 1.0/max<double>(centredPoint1.depth(P), centredPoint2.depth(P)); //Max better than min
        //adLambda[i] = 2/(centredPoint1.depth(P) + centredPoint2.depth(P)); //Max better than min
        //cout << centredPoint1.depth(P) << ' ' << centredPoint2.depth(P) << endl;
        //const double dLambda_inv = 1.0/adLambda[i];
        dTotalLambda += adLambda[i];
        mean1 += adLambda[i]*centredPoint1;
        mean2 += adLambda[i]*centredPoint2;
    }

    mean1 /= dTotalLambda;
    mean2 /= dTotalLambda;

    i=0;
    for(T3dPointMatchVector::const_iterator pPointPair = vPointMatches.begin(); pPointPair != vPointMatches.end(); pPointPair++, i++)
    {
        const C3dPoint & centredPoint1 = (pPointPair->p1() - mean1);
        const C3dPoint & centredPoint2 = (pPointPair->p2() - mean2);

        Vector3d v1,v2;
        centredPoint1.asVector(v1);
        centredPoint2.asVector(v2);
        v1 *= adLambda[i];
        v2 *= adLambda[i];

        X.block<3,1>(0, i) = v1;
        //X(3,i) =  adLambda[i];
        Y.col(i) = v2;
    }

    Eigen::Matrix<double, 3, 3, RowMajor + AutoAlign> U;
    Eigen::Matrix<double, 3, 1, ColMajor + AutoAlign> D;
    Eigen::Matrix<double, Dynamic, Dynamic, RowMajor + AutoAlign> V(nSize, nSize);

    //cout << "X=" << endl << X << endl ;

    svd(X, U, D, V);

    Eigen::Matrix<double, 3, Dynamic, RowMajor + AutoAlign> Diag(3, nSize);
    Diag.setZero();
    Diag.diagonal() = D;

/*    cout << "U=" << endl << U<< endl ;
    cout << "D=" << endl << D<< endl ;
    cout << "V=" << endl << V<< endl ;

    cout << "U*D*Vt=" << endl << U * Diag * V.transpose()<< endl ;*/

    for(int i=0;i<3;i++)
        Diag(i,i) = 1.0/Diag(i,i);

    Eigen::Matrix<double, Dynamic, 3 > pseudoInv = V * Diag.transpose() * U.transpose();

    //cout << "I matrix:\n" << X * pseudoInv << endl;
    Eigen::Matrix<double, 3, 3> A = Y * pseudoInv;
    //cout << "A matrix:\n" << A << endl;
    Matrix3d rot = A.block<3, 3>(0,0);

#if EIGEN_VERSION_AT_LEAST(2,90,0)
    JacobiSVD< Eigen::Matrix3d > svd2(rot, ComputeFullU | ComputeFullV);
#else
    SVD< Eigen::Matrix3d > svd2(rot);
#endif

    //cout << "SVs 2: " << svd2.singularValues().transpose();
    Matrix3d S; S.setIdentity();
    if(rot.determinant() <= 0){
        cout << "Warning: CAlignPointsNew::alignPoints: large errors (neg determinant)...";
        S(2, 2) = -1;
    }

    rot = svd2.matrixU() * S * svd2.matrixV().transpose();
    if(IS_DEBUG) CHECK(!isRotMat(rot), "Failed to extract a valid rot mat");

    R = C3dRotation(rot);
    //t = C3dPoint(A.col(3));
    t = mean2 - R * mean1;

    //if(!isRotMat(R_temp)) return false;

    //R=R_temp;
    //Was c = (D * S).trace() / var1; //Done: swap to get rid of 1 div
    //c = 1. / c;
    //VectorXd SVs = svd.singularValues();
    //cout << SVs << "=sv's" << endl;

    double dTrace = D.sum();
    if(zero(dTrace)) return false;
    //s = var1 / dTrace;
    //t = mean2 * s - R * mean1;
    //PRINTMAT(mean2)
    //PRINTMAT(R * mean1)
    return true;
}
#elif defined(OLD)
bool CAlignPointsNew::alignPoints(const T3dPointMatchVector & vPointMatches, C3dRotation & R, C3dPoint & t, double & s)
{
    C3dPoint mean1, mean2;

    int nSize = vPointMatches.size();

    //const double sizeInv = 1.0 / nSize;

    Matrix<double, 4, Dynamic, RowMajor + AutoAlign> X(4, nSize);
    Matrix<double, 3, Dynamic> Y(3, nSize);

    /*for(T3dPointMatchVector::const_iterator pPointPair = vPointMatches.begin(); pPointPair != vPointMatches.end(); pPointPair++)
    {
        const C3dPoint & p1 = pPointPair->p1();
        const C3dPoint & p2 = pPointPair->p2();
        mean1 += p1;
        mean2 += p2;
    }
    mean1 *= sizeInv;
    mean2 *= sizeInv;*/

    int i=0;
    ARRAY(double, adLambda, nSize);
    for(T3dPointMatchVector::const_iterator pPointPair = vPointMatches.begin(); pPointPair != vPointMatches.end(); pPointPair++, i++)
    {
        const C3dPoint & centredPoint1 = (pPointPair->p1());
        const C3dPoint & centredPoint2 = (pPointPair->p2());

        CCamera P;
        adLambda[i] = max<double>(centredPoint1.depth(P), centredPoint2.depth(P));
        cout << centredPoint1.depth(P) << ' ' << centredPoint2.depth(P) << endl;
        const double dLambda_inv = 1.0/adLambda[i];

        Vector3d v1,v2;
        centredPoint1.asVector(v1);
        centredPoint2.asVector(v2);
        v1 *= dLambda_inv;
        v2 *= dLambda_inv;

        X.block<3,1>(0, i) = v1;
        X(3,i) = dLambda_inv;
        Y.col(i) = v2;
    }

    Eigen::Matrix<double, 4, 4, RowMajor + AutoAlign> U;
    Eigen::Matrix<double, 4, 1, RowMajor + AutoAlign> D; 
    Eigen::Matrix<double, Dynamic, Dynamic, RowMajor + AutoAlign> V(nSize, nSize);

    //cout << "X=" << endl << X << endl ;

    svd(X, U, D, V);

    Eigen::Matrix<double, 4, Dynamic, RowMajor + AutoAlign> Diag(4, nSize);
    Diag.setZero();
    Diag.diagonal() = D;

/*    cout << "U=" << endl << U<< endl ;
    cout << "D=" << endl << D<< endl ;
    cout << "V=" << endl << V<< endl ;

    cout << "U*D*Vt=" << endl << U * Diag * V.transpose()<< endl ;*/

    for(int i=0;i<4;i++)
        Diag(i,i) = 1.0/Diag(i,i);

    Eigen::Matrix<double, Dynamic, 4 > pseudoInv = V * Diag.transpose() * U.transpose();

    //cout << "I matrix:\n" << X * pseudoInv << endl;
    Eigen::Matrix<double, 3, 4> A = Y * pseudoInv;
    //cout << "A matrix:\n" << A << endl;
    Matrix3d rot = A.block<3, 3>(0,0);

    Eigen::JacobiSVD<Matrix3d> svd2(rot, ComputeFullU | ComputeFullV);

    //cout << "SVs 2: " << svd2.singularValues().transpose();
    Matrix3d S; S.setIdentity();
    if(rot.determinant() <= 0){
        cout << "Warning: CAlignPointsNew::alignPoints: large errors (neg determinant)...";
        S(2, 2) = -1;
    }

    rot = svd2.matrixU() * S * svd2.matrixV().transpose();
    if(IS_DEBUG) CHECK(!isRotMat(rot), "Failed to extract a valid rot mat");

    R = C3dRotation(rot);
    t = C3dPoint(A.col(3));

    //if(!isRotMat(R_temp)) return false;

    //R=R_temp;
    //Was c = (D * S).trace() / var1; //Done: swap to get rid of 1 div
    //c = 1. / c;
    //VectorXd SVs = svd.singularValues();
    //cout << SVs << "=sv's" << endl;

    double dTrace = D.sum();
    if(zero(dTrace)) return false;
    //s = var1 / dTrace;
    //t = mean2 * s - R * mean1;
    //PRINTMAT(mean2)
    //PRINTMAT(R * mean1)
    return true;
}
#endif
