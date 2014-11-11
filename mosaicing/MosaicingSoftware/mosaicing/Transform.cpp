#include "Transform.h"
#include "cvGeom.h"
#include "GRCException.h"
#include "Warp.h"
#ifndef USE_OLD_OPENCV
#include "opencv2/opencv.hpp"
#endif

namespace grc
{

Enums::eTransformType Transform::eTransType_s = (Enums::eTransformType)-1;
int Transform::warpMethod_s = -1;

//!Transform factory--should this be in Transform class? 
Transform * Transform::newTransform()
{
    switch( eTransType_s )
    {
    case Enums::ePerspective:
        return new PerspectiveTransform;
    case Enums::eSimilarity:
        return new SimilarityTransform;
    case Enums::eAffine:
        return new AffineTransform;
    }
    throw new GRCException("Transform::newTransform: Type not recognised");
}

Transform * Transform::copyTransform(const Transform * T)
{
    switch( eTransType_s )
    {
    case Enums::ePerspective:
        return new PerspectiveTransform(*dynamic_cast<const PerspectiveTransform * >(T));
    case Enums::eSimilarity:
        return new SimilarityTransform(*dynamic_cast<const SimilarityTransform * >(T));
    case Enums::eAffine:
        return new AffineTransform(*dynamic_cast<const AffineTransform * >(T));
    }
    throw new GRCException("Transform::copyTransform: Type not recognised");
}

double Transform::operator()(int x, int y) const 
{
    return cvmGet(this, x, y);
}

void Transform::shift(double x, double y)
{
    (*this)(0,2) += x;
    (*this)(1,2) += y;    
}

//// PERSPECTIVE
PerspectiveTransform::PerspectiveTransform()
{
    **this = cvMat(ROWS, COLS, CV_64FC1, data_);
    cvSetIdentity((CvMat *)this);
}

PerspectiveTransform::PerspectiveTransform(const PerspectiveTransform & T) 
{
    memcpy(data_, T.data_, sizeof(data_));
    **this = cvMat(ROWS, COLS, CV_64FC1, data_);
}

CvPoint2D64f PerspectiveTransform::translation() const
{
    double dHomoScale = 1. / cvmGet(*this, 2,2);
    return cvPoint2D64f(cvmGet(*this, 0, 2)*dHomoScale, cvmGet(*this, 1, 2)*dHomoScale);
};

//! Add all the tunable parameters for this transform to the parameter vector
void PerspectiveTransform::addParams(TParamLocations * paramVector)
{
    for(int idx=0; idx<8; idx++)
    {
        double weight = 1.0;
        if(idx==2 || idx == 5) weight = 1000.0;
        if(idx == 6 || idx == 7) weight = 0.01;

        paramVector->push_back(param(data_ + idx, weight));
    }
};

CvPoint2D64f PerspectiveTransform::applyToPoint(CvPoint2D64f p) const
{
/*    double pointData[COLS], newPointData[ROWS];
    CvMat pointAsMat = cvMat(COLS, 1, CV_64FC1, pointData);
    pointData[0] = p.x; pointData[1] = p.y; pointData[2] = 1.0;
    CvMat transPointAsMat = cvMat(ROWS, 1, CV_64FC1, newPointData);
    cvMatMul((CvMat *)this, &pointAsMat, &transPointAsMat);
    double dHomogScale = 1.0/newPointData[2];
    return cvPoint2D64f(dHomogScale*newPointData[0], dHomogScale*newPointData[1]);
*/
    double x = data_[0]*p.x + data_[1]*p.y + data_[2];
    double y = data_[3]*p.x + data_[4]*p.y + data_[5];
    double t = data_[6]*p.x + data_[7]*p.y + data_[8];
    double t_inv = 1.0/t;

    return cvPoint2D64f(t_inv*x, t_inv*y);
}

double & PerspectiveTransform::operator()(int x, int y)
{
    return data.db[x*COLS+y];
}

void matMul(const double * s1, const double * s2, double * d, int I, int J, int K)
{
    //s1 is i by j, s2 is j by k, dest is i by k
    for(int i=0; i<I; i++)
    {
        const double * s1_i = s1 + i*J;
        for(int k=0; k<K; k++)
        {
            double val = 0;
            for(int j=0; j<J; j++)
                val += s1_i[j]*s2[j * K + k];

            d[i*K + k] = val;
        }
    }
}

PerspectiveTransform * PerspectiveTransform::accumulate(const Transform * T2) const
{
    PerspectiveTransform * m = new PerspectiveTransform;
    //cvMatMul(*this, *T2, *m); //don't need casts as opencv will check dims agree.
    matMul(data.db, T2->data.db, m->data.db, 3, 3, 3);
    return m;
}

PerspectiveTransform * PerspectiveTransform::accumulateInverse(const Transform * T2) const
{
    PerspectiveTransform * m = new PerspectiveTransform;
    PerspectiveTransform T2_inv(*dynamic_cast<const PerspectiveTransform *>(T2)); //make a copy

    cvInvert(T2_inv, T2_inv);

    cvMatMul(T2_inv, *this, *m); //don't need casts as opencv will check dims agree.
    matMul(T2_inv.data.db, data.db, m->data.db, 3, 3, 3);
    return m;
}

//!Apply to a matrix of points
void PerspectiveTransform::applyToPoints(const CvMat * positions, CvMat * newPositions) const
{
	cvMatMul(*this, positions, newPositions);

    for(int i=0; i<newPositions->cols; i++)
    {
        double x=cvmGet(newPositions, 0, i);
        double y=cvmGet(newPositions, 1, i);
        double t_inv= 1. / cvmGet(newPositions, 2, i);

        cvmSet(newPositions, 0, i, x*t_inv);
        cvmSet(newPositions, 1, i, y*t_inv);
        cvmSet(newPositions, 2, i, 1.0);
    }
}

//!Estimate transform from a set of points
void PerspectiveTransform::estimateFromPoints(const CvMat * points1, const CvMat * points2)
{
    cvFindHomography(points1, points2, (CvMat*)*this);
}

void PerspectiveTransform::applyToImage(const IplImage * sourceIm, IplImage * destIm) const
{
    WarpPerspective(sourceIm, destIm, *this, warpMethod_s);
}

//// Affine
AffineTransform::AffineTransform()
{
    **this = cvMat(ROWS, COLS, CV_64FC1, data_);
    cvSetIdentity(*this);
}
AffineTransform::AffineTransform(const AffineTransform & T)
{
    memcpy(data_, T.data_, sizeof(data_));
    **this = cvMat(ROWS, COLS, CV_64FC1, data_);
}

CvPoint2D64f AffineTransform::translation() const
{
    return cvPoint2D64f(cvmGet(*this, 0, 2), cvmGet(*this, 1, 2));
};

//! Add all the tunable parameters for this transform to the parameter vector
void AffineTransform::addParams(TParamLocations * paramVector)
{
    for(int idx=0; idx<6; idx++)
    {
        double weight = 1.0;
        if(idx==2 || idx == 5) weight = 1000.0;

        paramVector->push_back(param(data_ + idx, weight));
    }
};

CvPoint2D64f AffineTransform::applyToPoint(CvPoint2D64f p) const
{
/*    double pointData[COLS], newPointData[ROWS];
    CvMat pointAsMat = cvMat(COLS, 1, CV_64FC1, pointData);
    pointData[0] = p.x; pointData[1] = p.y; pointData[2] = 1.0;
    CvMat transPointAsMat = cvMat(ROWS, 1, CV_64FC1, newPointData);
    cvMatMul(*this, &pointAsMat, &transPointAsMat);
    return cvPoint2D64f(newPointData[0], newPointData[1]);*/

    double x = data_[0]*p.x + data_[1]*p.y + data_[2];
    double y = data_[3]*p.x + data_[4]*p.y + data_[5];

    return cvPoint2D64f(x, y);
}

double & AffineTransform::operator()(int x, int y)
{
    return data.db[x*COLS+y];
}

AffineTransform * AffineTransform::accumulate(const Transform * T2) const
{
    return accumulateInt(T2, false);
}

AffineTransform * AffineTransform::accumulateInt(const Transform * T2, bool bInvert) const
{
    AffineTransform * m = new AffineTransform;

    // Just expand mat. multiplication to avoid converting to 3x3
    const double * s1 = 0, * s2 = 0;
    double invMatData[3*3]; //here so stays in scope when pointer copied
    if(bInvert)
    {
        //Invert 3x3 mat
        double T2MatData[3*3];
        CvMat invMat = cvMat(3, 3, CV_64FC1, invMatData);
        CvMat T2Mat = cvMat(3, 3, CV_64FC1, T2MatData);

        const double * affineT2Data = dynamic_cast<const AffineTransform *>(T2)->data_;
        for(int i=0; i<6; i++)
            T2MatData[i] = affineT2Data[i];
        
        //Last row of affine mat
        T2MatData[6] = 0;
        T2MatData[7] = 0;
        T2MatData[8] = 1;

        cvInvert(&T2Mat, &invMat);

        s1 = invMatData;
        s2 = this->data_;
    }
    else
    {
        s1 = this->data_;
        s2 = dynamic_cast<const AffineTransform *>(T2)->data_;
    }

    for(int row=0; row<2; row++)
    {
        (*m)(row,0) = s1[0]*s2[0] + s1[1]*s2[0+3];
        (*m)(row,1) = s1[0]*s2[1] + s1[1]*s2[1+3];
        (*m)(row,2) = s1[0]*s2[2] + s1[1]*s2[2+3] + s1[2];
        s1 += COLS; //skip to second 
    }
    return m;
}
AffineTransform * AffineTransform::accumulateInverse(const Transform * T2) const
{
    return accumulateInt(T2, true);
}

//!Apply to a matrix of points
void AffineTransform::applyToPoints(const CvMat * positions, CvMat * newPositions) const
{
    CvMat newPositions2d;
    cvGetSubRect(newPositions, &newPositions2d, cvRect(0,0,newPositions->cols, 2));
	cvMatMul(*this, positions, &newPositions2d);

    for(int i=0; i<newPositions->cols; i++)
    {
        cvmSet(newPositions, 0, i, cvmGet(&newPositions2d, 0, i));
        cvmSet(newPositions, 1, i, cvmGet(&newPositions2d, 1, i));
        cvmSet(newPositions, 2, i, 1.0);
    }
}

//!Estimate transform from a set of points
void AffineTransform::estimateFromPoints(const CvMat * points1, const CvMat * points2)
{
    double perspectiveMatData[3*3];
    CvMat matPerspective = cvMat(3, 3, CV_64FC1, perspectiveMatData);
    cvFindHomography(points1, points2, &matPerspective);
//Todo: adapt cvGetAffineTransform--it's in Warp.cpp
    cvConvertScale(&matPerspective, &matPerspective, 1. / cvmGet(&matPerspective, 2, 2));
    for(int x=0; x < ROWS; x++)
        for(int y=0; y < COLS; y++)
        {
            cvmSet(*this, x, y, cvmGet(&matPerspective, x, y));
        }
}

void AffineTransform::applyToImage(const IplImage * sourceIm, IplImage * destIm) const
{
    WarpAffine(sourceIm, destIm, *this, warpMethod_s);
}

//// Similarity

//! Add all the tunable parameters for this transform to the parameter vector
void SimilarityTransform::addParams(TParamLocations * paramVector)
{
    paramVector->push_back(param(&scale_, 1));
    paramVector->push_back(param(&theta_, 1));
    paramVector->push_back(param(&translation_.x, 100));
    paramVector->push_back(param(&translation_.y, 100));
};

CvPoint2D64f SimilarityTransform::translation() const
{
    return translation_;
};

void SimilarityTransform::applyToImage(const IplImage * sourceIm, IplImage * destIm) const
{
    AffineTransform * pAT = getAffineTransform();
    pAT->applyToImage(sourceIm, destIm);
    delete pAT;
}

void SimilarityTransform::applyToPoints(const CvMat * positions, CvMat * newPositions) const
{
    AffineTransform * pAT = getAffineTransform();
    pAT->applyToPoints(positions, newPositions);
    delete pAT;
}

inline double sqr(double x) { return x*x; };

//!Estimate transform from a set of points
void SimilarityTransform::estimateFromPoints(const CvMat * points1, const CvMat * points2)
{
    //const CvMat * temp;
    //CV_SWAP(points1, points2, temp);

/*    AffineTransform * pAT = getAffineTransform();
    pAT->estimateFromPoints(points1, points2);
    delete pAT;*/
    //Umeyama's algorithm:
    //Find mean and s.d.
    int numPoints = points1->cols;
    double meanP1Data[2];
    CvMat meanP1 = cvMat(2, 1, CV_64FC1, meanP1Data);
    double meanP2Data[2];
    CvMat meanP2 = cvMat(2, 1, CV_64FC1, meanP2Data);
    cvSetZero(&meanP1);
    cvSetZero(&meanP2);
    
    for(int i = 0; i<numPoints; i++)
    {
        meanP1Data[0] += cvmGet(points1, 0, i);
        meanP1Data[1] += cvmGet(points1, 1, i);
        meanP2Data[0] += cvmGet(points2, 0, i);
        meanP2Data[1] += cvmGet(points2, 1, i);
    }

    double numPoints_inv = 1.0/numPoints;
    meanP1Data[0] *= numPoints_inv;
    meanP1Data[1] *= numPoints_inv;
    meanP2Data[0] *= numPoints_inv;
    meanP2Data[1] *= numPoints_inv;

    //Now calculate variance
    double varP1, varP2;
    varP1 = 0;
    varP2 = 0;

    double SIGMAData[4];
    CvMat SIGMA = cvMat(2, 2, CV_64FC1, SIGMAData);
    cvSetZero(&SIGMA);

    for(int i = 0; i<numPoints; i++)
    {
        double x1 = cvmGet(points1, 0, i) - meanP1Data[0];
        double y1 = cvmGet(points1, 1, i) - meanP1Data[1];
        double x2 = cvmGet(points2, 0, i) - meanP2Data[0];
        double y2 = cvmGet(points2, 1, i) - meanP2Data[1];

        varP1 += sqr(x1) + sqr(y1);
        varP2 += sqr(x2) + sqr(y2);

        SIGMAData[0] += x1*x2;
        SIGMAData[1] += x1*y2;
        SIGMAData[2] += y1*x2;
        SIGMAData[3] += y1*y2;
    }
    varP1 *= numPoints_inv;
    varP2 *= numPoints_inv;

    SIGMAData[0] *= numPoints_inv;
    SIGMAData[1] *= numPoints_inv;
    SIGMAData[2] *= numPoints_inv;
    SIGMAData[3] *= numPoints_inv;

    double DData[4];
    CvMat D = cvMat(2, 2, CV_64FC1, DData);
    cvSetZero(&D);
    double UData[4];
    CvMat U = cvMat(2, 2, CV_64FC1, UData);
    cvSetZero(&U);
    double VData[4];
    CvMat V = cvMat(2, 2, CV_64FC1, VData);
    cvSetZero(&V);
    double RotationData[4];
    CvMat Rotation = cvMat(2, 2, CV_64FC1, RotationData);
    cvSetZero(&Rotation);

    cvSVD(&SIGMA, &D, &U, &V);
    cvGEMM(&U, &V, 1, 0, 0, &Rotation, CV_GEMM_B_T);

//    theta_ = acos(RotationData[0]);
    theta_ = asin(RotationData[1]);

    scale_ = (1.0/varP1) * cvTrace(&D).val[0];

    double transData[2];
    CvMat trans = cvMat(2, 1, CV_64FC1, transData);
    cvSetZero(&trans);

    cvGEMM( &Rotation, &meanP1, -scale_, &meanP2, 1, &trans);

    translation_ = cvPoint2D64f(transData[0], transData[1]);

    //Applying this to p1 gives us p2
}


SimilarityTransform::SimilarityTransform(double theta, double scale, const CvPoint2D64f & translation)
    : theta_(theta), scale_(scale), translation_(translation)
{
}

SimilarityTransform::SimilarityTransform()
{
    //**this = cvMat(ROWS, COLS, CV_64FC1, data_);
    //cvSetIdentity(this);
    scale_ = 1.0;
    theta_ = 0.0;
    translation_ = cvPoint2D64f(0,0);
}

SimilarityTransform::SimilarityTransform(const SimilarityTransform & T)
{
    scale_ = T.scale_;
    theta_ = T.theta_;
    translation_ = T.translation_;
}

CvPoint2D64f SimilarityTransform::applyToPoint(CvPoint2D64f p) const
{
    AffineTransform * pAT = getAffineTransform();
    CvPoint2D64f transPoint = pAT->applyToPoint(p);
    delete pAT;
    return transPoint;
}

double & SimilarityTransform::operator()(int x, int y)
{
    throw new GRCException("SimilarityTransform: Not Implemented as matrix");
}

SimilarityTransform * SimilarityTransform::accumulate(const Transform * T2) const
{
    AffineTransform * T1aff = getAffineTransform();
    AffineTransform * T2aff = dynamic_cast<const SimilarityTransform *>(T2)->getAffineTransform();
    AffineTransform * aff = dynamic_cast<AffineTransform *>(T1aff->accumulate(T2aff));

    SimilarityTransform * pSim = aff->getSimilarityTransform();
    
    delete T1aff;
    delete T2aff;
    delete aff;

    return pSim;
}

SimilarityTransform * SimilarityTransform::accumulateInverse(const Transform * T2) const
{
    AffineTransform * T1aff = getAffineTransform();
    AffineTransform * T2aff = dynamic_cast<const SimilarityTransform *>(T2)->getAffineTransform();
    AffineTransform * aff = T1aff->accumulateInverse(T2aff);
    SimilarityTransform * pSim = aff->getSimilarityTransform();
//    m->scale_ = scale_ * T2sim->scale_;
//    m->translation_ = translation_ + T2sim->translation_; //Do we need to scale here?? No
    
    delete T1aff;
    delete T2aff;
    delete aff;

    return pSim;
}

AffineTransform * SimilarityTransform::getAffineTransform() const
{
    AffineTransform * pm = new AffineTransform;
    AffineTransform &m = *pm;
    //for(int x=0; x<2; x++)
    //    for(int y=0; y<2; y++)
    //        m(x,y) = scale_ * cvmGet(*this, x, y);
    double scaledCosTheta = scale_ * cos(theta_);
    double scaledSinTheta = scale_ * sin(theta_);

    m(0,0) = scaledCosTheta;
    m(1,1) = scaledCosTheta;
    m(0,1) = scaledSinTheta;
    m(1,0) = -scaledSinTheta;

    m(0,2) = translation_.x;
    m(1,2) = translation_.y;

    return pm;
}

SimilarityTransform * AffineTransform::getSimilarityTransform() const
{
    double subMatData[4];
    CvMat subMat = cvMat(2, 2, CV_64FC1, subMatData);
    subMatData[0] = data_[0];
    subMatData[1] = data_[1];
    subMatData[2] = data_[3];
    subMatData[3] = data_[4];
    
    double scale = sqrt(cvDet(&subMat));
    //double theta = acos((*(const Transform *)this)(0,0)/scale); //Todo: do a LS fit?
    double theta = asin((*(const Transform *)this)(0,1)/scale); //Todo: do a LS fit?
    CvPoint2D64f translation = cvPoint2D64f((*(const Transform *)this)(0,2), (*(const Transform *)this)(1,2));
    SimilarityTransform * pm = new SimilarityTransform(theta, scale, translation);

    //for(int x=0; x<2; x++)
    //    for(int y=0; y<2; y++)
    //        m(x,y) = scale_ * cvmGet(*this, x, y);

    return pm;
}
void SimilarityTransform::shift(double x, double y)
{
    translation_.x += x;
    translation_.y += y;    
}


}
