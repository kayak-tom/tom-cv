/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once
#ifndef GEOM_H
#define GEOM_H

#include "util/location.h"
#include "util/dynArray.h"
#include <iostream>
#include <Eigen/Core>
#include "util/calibration.h"

template<int I, int J, int K> void matMul(const double * s1, const double * s2, double * d) {
    //s1 is i by j, s2 is j by k, dest is i by k
    double * dest = d;
    double temp_dest[I * K];
    if (s1 == d || s2 == d) {
        dest = temp_dest;
    }
    for (int i = 0; i < I; i++) {
        const double * s1_i = s1 + i*J;
        for (int k = 0; k < K; k++) {
            double val = 0;
            for (int j = 0; j < J; j++)
                val += s1_i[j] * s2[j * K + k];

            dest[i * K + k] = val;
        }
    }
    if (s1 == d || s2 == d) {
        for (int i = I * K; i > 0; i--) {
            *d = *dest;
            d++;
            dest++;
        }
    }
}

class C3dPoint;
class CCamera;
//class Matrix;
//class ColumnVector;

class C3dRotationMat // this is a Rotation mat
{
    double R[9];
public:

    C3dRotationMat() {
        R[0] = R[4] = R[8] = 1;
        R[1] = R[2] = R[3] = R[5] = R[6] = R[7] = 0;
    };

    double operator[](int i) const {
        return R[i];
    };

    C3dRotationMat operator*(const C3dRotationMat & R2) const {
        C3dRotationMat multRot;
        matMul < 3, 3, 3 > (R, R2.R, multRot.R);
        return multRot;
    };

    C3dPoint operator*(const C3dPoint & p) const;

    C3dRotationMat t() const //'transpose' -- inverse rotation
    {
        C3dRotationMat invertedRot;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                invertedRot.R[i + 3 * j] = R[3 * i + j];

        return invertedRot;
    };

    CCamera operator|(const C3dPoint & p) const;

    double angle() const;

    bool operator==(const C3dRotationMat & p) const;

    C3dPoint headingVector() const;
};

class C3dRotationQuat {
    friend class C3dPoint;
    friend std::ostream& operator<<(std::ostream& s, const C3dRotationQuat& X);

    void lengthOk() const; //just for consistency checking
    void fromMat(const Eigen::Matrix3d & rot);
protected:
    double Q[4];

public:

    double operator[](int i) const {
        return Q[i];
    };
    double & operator[](int i) {
        return Q[i];
    };
    void normalise();

    C3dRotationQuat(const CCamera & cam);

    C3dRotationQuat() {
        Q[0] = Q[1] = Q[2] = 0;
        Q[3] = 1;
    }

    C3dRotationQuat(double x, double y, double z, double w) {
        Q[0] = x;
        Q[1] = y;
        Q[2] = z;
        Q[3] = w;
        normalise();
    }
    C3dRotationQuat(C3dPoint axis, double angle);
    C3dRotationQuat(const Eigen::Matrix3d & cam);

    C3dRotationQuat(const Eigen::Vector4d & v) {
        Q[0] = v(0);
        Q[1] = v(1);
        Q[2] = v(2);
        Q[3] = v(3);
    }

    C3dRotationQuat(const Eigen::Vector4f & v) {
        Q[0] = (double) v(0);
        Q[1] = (double) v(1);
        Q[2] = (double) v(2);
        Q[3] = (double) v(3);
    }

    C3dRotationQuat operator*(const C3dRotationQuat & R2) const;

    C3dPoint operator*(const C3dPoint & p) const;

    C3dRotationQuat t() const //'transpose' -- inverse rotation
    {
        C3dRotationQuat invertedRot;

        for (int i = 0; i < 3; i++)
            invertedRot.Q[i] = -Q[i];

        invertedRot.Q[3] = Q[3];

        return invertedRot;
    };

    void asMat(Eigen::Matrix3d & R) const;

    CCamera operator|(const C3dPoint & p) const;

    double angle() const {
        double cosAng = fabs(Q[3]);
        if (cosAng > 1.0) cosAng = 1.0;
        double ang = 2 * acos(cosAng);
        if(IS_DEBUG) CHECK(ang < 0, "C3dRotationQuat::angle(): acos returning val less than 0");
        if(IS_DEBUG) CHECK(std::isnan(ang), "C3dRotationQuat::angle(): acos returning nan");
        if (ang > M_PI) ang -= 2 * M_PI;
        return ang;
    };

    bool operator==(const C3dRotationQuat & p) const;

    bool operator!=(const C3dRotationQuat & p) const {
        return !(*this == p);
    };

    C3dPoint headingVector() const;

    void toPhiThetaPsi(double & phi, double & theta, double & psi) const;
    void toPhiThetaPsi2(double & phi, double & theta, double & psi) const;

    static void quatMult(const C3dRotationQuat & Q1, const C3dRotationQuat & Q2, C3dRotationQuat & Qres) HOT;
    static void quatMultByVec(const C3dRotationQuat & Q1, const C3dPoint & vec, C3dRotationQuat & Qres) HOT;
    static void conjQuatMultVec(const C3dRotationQuat & Q1, const C3dPoint & vec, C3dRotationQuat & Qres) HOT;
    static void quatMultConj(const C3dRotationQuat & Q1, const C3dRotationQuat & Q2, C3dRotationQuat & Qres) HOT;
    static void quatMultConjIntoVec(const C3dRotationQuat & Q1, const C3dRotationQuat & Q2, C3dPoint & Qres) HOT;
    static void quatMultIntoVec(const C3dRotationQuat & Q1, const C3dRotationQuat & Q2, C3dPoint & Qres) HOT;

    C3dRotationQuat conj() const {
        return C3dRotationQuat(-Q[0], -Q[1], -Q[2], Q[3]);
    }
    
    void setRandom(double dAngle);
    void setRandom();
};

inline double diff(const C3dRotationQuat & q1, const C3dRotationQuat & q2) {
    return (q1 * q2.t()).angle();
}

typedef C3dRotationQuat C3dRotation;

class C2dPoint;
class CCamCalibMatrix;

class C3dPoint {
protected:
    double x, y, z;
    friend class C2dPoint;
    friend class C3dRotationQuat;
public:

    C3dPoint() : x(0), y(0), z(0) {
        if(IS_DEBUG) CHECK(std::isnan(x), "C3dPoint: NaN");
    };

    C3dPoint(double x, double y, double z) : x(x), y(y), z(z) {
        if(IS_DEBUG) CHECK(std::isnan(x + y + z), "C3dPoint: NaN");
    };

    C3dPoint(const C3dPoint & p) : x(p.x), y(p.y), z(p.z) {
        if(IS_DEBUG) CHECK(std::isnan(x + y + z), "C3dPoint: NaN");
    };
    C3dPoint(const C2dPoint & p, const double depth);
    C3dPoint(const Eigen::Vector3d & vec);
    C3dPoint(const Eigen::Vector3f & vec);
    Eigen::Vector3d asVector() const { return Eigen::Vector3d(x, y, z); }

    void operator*=(double s) {
        x *= s;
        y *= s;
        z *= s;
    };

    void operator/=(double s) {
        if(IS_DEBUG) CHECK(s == 0, "C3dPoint: Divide by 0");
        double s_inv = 1.0 / s;
        x *= s_inv;
        y *= s_inv;
        z *= s_inv;
    };

    void operator+=(const C3dPoint &s) {
        x += s.x;
        y += s.y;
        z += s.z;
    };

    void operator-=(const C3dPoint &s) {
        x -= s.x;
        y -= s.y;
        z -= s.z;
    };

    C3dPoint operator*(double p) const {
        C3dPoint scaled(x*p, y*p, z * p);
        return scaled;
    };

    C3dPoint operator+(const C3dPoint & p) const {
        C3dPoint added(x + p.x, y + p.y, z + p.z);
        return added;
    };

    C3dPoint operator-(const C3dPoint & p) const {
        C3dPoint subtracted(x - p.x, y - p.y, z - p.z);
        return subtracted;
    };

    C3dPoint operator-() const {
        C3dPoint reversed(-x, -y, -z);
        return reversed;
    };

    //Element-wise (for variances)

    C3dPoint operator*(const C3dPoint & p) const {
        C3dPoint elwiseProd(x * p.x, y * p.y, z * p.z);
        return elwiseProd;
    }

    void rotate(const C3dRotationMat & R) //Todo: should't rotate in place--need to rename
    {
        double xNew = R[0] * x + R[1] * y + R[2] * z;
        double yNew = R[3] * x + R[4] * y + R[5] * z;
        double zNew = R[6] * x + R[7] * y + R[8] * z;
        x = xNew;
        y = yNew;
        z = zNew;
    };

    void rotateInv(const C3dRotationMat & R) //Todo: should't rotate in place--need to rename (?)
    {
        double xNew = R[0] * x + R[3] * y + R[6] * z;
        double yNew = R[1] * x + R[4] * y + R[7] * z;
        double zNew = R[2] * x + R[5] * y + R[8] * z;
        x = xNew;
        y = yNew;
        z = zNew;
    };
    void rotate(const C3dRotationQuat & R); //Todo: should't rotate in place--need to rename
    void rotateInv(const C3dRotationQuat & R); //Todo: should't rotate in place--need to rename

    void operator=(const C3dPoint & p) {
        x = p.x;
        y = p.y;
        z = p.z;
    };

    void operator=(double d) {
        x = d;
        y = d;
        z = d;
    };

    void operator=(int d) {
        x = d;
        y = d;
        z = d;
    };

    /*void operator=(const ColumnVector p);

    ColumnVector asCV() const;*/

    C2dPoint photo(const CCamera & P) const;
    bool testInFront(const CCamera & P) const;

    bool testInFront(const CCamera & P, const CCamera & Pp) const {
        return testInFront(P) && testInFront(Pp);
    };
    double depth(const CCamera & P) const;

    bool operator==(const C3dPoint & p) const;

    double sum_square() const {
        return x * x + y * y + z*z;
    };

    double length() const {
        return sqrt(sum_square());
    };

    void normalise() //normalise to length 1
    {
        double dLength = length();
        if (fabs(dLength - 1) > 0.0000001)
            * this /= dLength;
    }

    double getX() const {
        return x;
    };

    double getY() const {
        return y;
    };

    double getZ() const {
        return z;
    };

    void asVector(Eigen::Vector3d & vec) const;

    void addNoise(double dSD);

    void setRandomNormal();
    void setRandom(double dDepth);
    void setRandomPlanar(double dDepth);
};

C3dPoint operator*(double p, const C3dPoint & T);
C3dPoint crossproduct(const C3dPoint & p1, const C3dPoint & p2);
double dotproduct(const C2dPoint & p1, const C2dPoint & p2);
double dotproduct(const C3dPoint & p1, const C3dPoint & p2);
//Angle in radians between 2 unit vectors

inline double angle(const C3dPoint & p1, const C3dPoint & p2) {
    return fabs(acos(std::min<double>(dotproduct(p1, p2), 1.0)));
}

class C2dPoint : public CSimple2dPoint {
    //double x, y;
public:

    C2dPoint(const C3dPoint & p3) : CSimple2dPoint(0, 0) {
        double dHomoInv = 1.0 / p3.getZ();
        x = p3.getX() * dHomoInv;
        y = p3.getY() * dHomoInv;
    };

    C2dPoint() : CSimple2dPoint(0, 0) {
    };

    C2dPoint(double x, double y) : CSimple2dPoint(x, y) {
    };

    C2dPoint(CLocation loc) : CSimple2dPoint(loc.dx(), loc.dy()) {
    };
    C2dPoint(const Eigen::Vector2d & v);
    C2dPoint(const Eigen::Vector3d & v);

    C2dPoint(const CSimple2dPoint & p) : CSimple2dPoint(p) {
    }

    //double getX() const { return CSimple2dPoint::getX(); };
    //double getY() const { return CSimple2dPoint::getY(); };

    void operator*=(double s) {
        x *= s;
        y *= s;
    };

    void operator=(const C2dPoint & p) {
        x = p.x;
        y = p.y;
    };
    //void operator=(C2dPoint & p) { x=p.x;y=p.y; }

    void operator+=(const C2dPoint & p) {
        x += p.x;
        y += p.y;
    };

    void operator-=(const C2dPoint & p) {
        x -= p.x;
        y -= p.y;
    };

    bool operator==(const C2dPoint & p) const {
        return x == p.x && y == p.y;
    };

    C2dPoint operator+(const C2dPoint & p) const {
        return C2dPoint(x + p.x, y + p.y);
    };
    C2dPoint operator-(const C2dPoint & p) const {
        return C2dPoint(x - p.x, y - p.y);
    };
    C2dPoint operator-() const {
        return C2dPoint(-x, -y);
    };

    double sum_square() const {
        return x * x + y*y;
    };

    double length() const {
        return sqrt(sum_square());
    };

    Eigen::Vector2d asVector() const { return Eigen::Vector2d(x, y); }

    void asVector(Eigen::Vector2d & vec) const;
    void asVector(Eigen::Vector3d & vec) const;

    void addNoise(double dSD);

    bool inBox(const double left, const double bottom, const double right, const double top) const {
        return x >= left && x < right && y >= bottom && y < top;
    }
        
    void normalise() //normalise to length 1
    {
        double dLength = length();
        if (fabs(dLength - 1) > 0.0000001)
            *this *= 1.0/dLength;
    }
};

inline C2dPoint operator*(double p, const C2dPoint & T) { return C2dPoint(p*T.getX(), p*T.getY()); }

std::ostream& operator<<(std::ostream& s, const CCamera& X);

class CCamera {
    friend CCamera C3dRotationMat::operator|(const C3dPoint & p) const;
    friend CCamera C3dRotationQuat::operator|(const C3dPoint & p) const;

    friend class C3dRotationQuat;
    friend class C3dPoint;

    friend C2dPoint operator*(const CCamera & P, const C3dPoint & p);
    friend std::ostream& operator<<(std::ostream& s, const CCamera& X);

    double adCam[12]; double dFocalLength; //Normally dFocalLength = 1 unless has been combined with a calibration matrix
    
    double at1(int r,int c) const { return at(r-1,c-1); }
public:
    double at(int r,int c) const { return rowData(r)[c]; }

    CCamera() : dFocalLength(1) {
        for (int i = 1; i < 12; i++)
            adCam[i] = 0;
        adCam[0] = adCam[5] = adCam[10] = 1;
    }

    const double * rowData(int nRow) const {
        return adCam + nRow * 4;
    }

    double * rowData(int nRow) {
        return adCam + nRow * 4;
    }

    C3dPoint translation() const {
        return C3dPoint(adCam[3], adCam[7], adCam[11]);
    }

    C3dRotationQuat rotation() const {
        return C3dRotation(*this);
    }
    void calibrate(const CCamCalibMatrix & K);

    void addOffset(const C3dPoint &dt) {
        adCam[3] += dt.getX();
        adCam[7] += dt.getY();
        adCam[11] += dt.getZ();
    }

    void makePureRot() {
        adCam[3] = adCam[7] = adCam[11] = 0;
    }
    
    void setFocalLength(double f) { dFocalLength = f; };
    double focalLength() const { return dFocalLength; };
    
    //Convert image coordinates to world
    C3dPoint imageToWorld(const C2dPoint & imPoint, const double dDepth) const;
};

typedef CDynArray<C2dPoint> TPointVec2d;

class CPointVec2d : public TPointVec2d {
public:

    CPointVec2d(int n) : TPointVec2d(n) {
    }

    CPointVec2d() {
    }

    const T2dPoints & operator()() const {
        CHECK_SIZES_EQUAL_RT(CSimple2dPoint, C2dPoint, CPointVec2d);
        return *reinterpret_cast<const T2dPoints *> (this);
    }
};

//TODO: Move me to geom_Eigen.h. Not the declaration.
void Xmat(Eigen::Vector3d const & e, Eigen::Matrix3d & X_matrix);

std::ostream& operator<<(std::ostream& s, const CCamCalibMatrix & X);
std::ostream& operator<<(std::ostream& s, const C2dPoint& X);
std::ostream& operator<<(std::ostream& s, const C3dPoint& X);
//ostream& operator<<(ostream& s, const C3dRotationMat& X);
std::ostream& operator<<(std::ostream& s, const C3dRotationQuat& X);
std::ostream& operator<<(std::ostream& s, const CLocation& X);

void printMat2(const C3dPoint &G, const char * caption);
void printMat2(const C3dRotation &G, const char * caption);

int dist(CLocation l1, CLocation l2);

bool clip(CLocation loc, const int IM_WIDTH, const int IM_HEIGHT);
#define CLIP(loc) clip(loc, IM_WIDTH, IM_HEIGHT)

CLocation colVecToLoc(const C2dPoint & v);
//Unused I think? Causing compiler errors CLocation colVecToLoc(const C3dPoint & v);

C3dPoint reconstruct(const CCamera &P, const CCamera &Pp, const C2dPoint &p, const C2dPoint &pp) HOT;
bool testPair(const CCamera &P, const CCamera &Pp, const C2dPoint &p1, const C2dPoint &p2);
bool testPointInFront(const CCamera &P, const CCamera &Pp, const C3dPoint &Q);

template<class T> T sign(T x) {
    return x > 0 ? 1 : -1;
}

void getCamsFromE(const double * pdE, CCamera aPp[4]); //Still used in autotests
void getCamsFromE(const Eigen::Matrix3d & E, CCamera aPp[4]);
/*int chooseCamMat(const CCamera * aPp, CPointVec2d & ap1, CPointVec2d & ap2);*/

bool chooseCamFromE(const double * pdE, CPointVec2d & ap1, CPointVec2d & ap2, CCamera & Pp, const bool bAllowPureRotation = true);
bool chooseCamFromE(const Eigen::Matrix3d & E, CPointVec2d & ap1, CPointVec2d & ap2, CCamera & Pp, const bool bAllowPureRotation = true);
bool chooseCamFromRT(const C3dRotation & R, const C3dPoint & t, CPointVec2d & ap1, CPointVec2d & ap2, CCamera & Pp, const bool bAllowPureRotation);

double camReprojectionErr(const CCamera & P, const CCamera & Pp, const CPointVec2d & p1, const CPointVec2d & p2, CMask * pMask);
void camReprojectionErr(const CCamera & P, const CCamera & Pp, const CPointVec2d & p1, const CPointVec2d & p2, Eigen::VectorXd & vResiduals, CMask * pMask);
//void camReprojectionErr(const CCamera & P, const CCamera & Pp, const CPointVec2d & p1, const CPointVec2d & p2, Eigen::VectorXd & vResiduals);

bool isRotMat(const Eigen::Matrix3d & R);

class C3dPose {
public:
    C3dRotation R; //Todo protect
    C3dPoint t;

    C3dPose(const C3dRotation R, const C3dPoint t) : R(R), t(t) {
    }

    C3dPose() {
    }

    C3dPose operator+(const C3dPose & p2) const {
        return C3dPose(p2.R * R, t + R.t() * p2.t);
        //Works worse return C3dPose(p2.R * R, t + R * p2.t);
    }

    void write(std::ostream & outFile, const bool b2d) const {
        outFile << t.getX() << ' ';
        if (!b2d)
            outFile << t.getY() << ' ';
        outFile << t.getZ() << ' ';

        double phi, theta, psi;
        R.t().toPhiThetaPsi2(phi, theta, psi); //Transpose because we compose rotations in different order
        //They are read in by TORO as phi, theta, psi, I think that makes them (http://en.wikipedia.org/wiki/Yaw,_pitch_and_roll):
        if (b2d) {
            outFile << theta << ' ';
        } else
            outFile << phi << ' ' << theta << ' ' << psi << ' ';
    }
};

class C3dNormalisedPose : public C3dPose {
public:

    C3dNormalisedPose(const C3dRotation R, const C3dPoint t) : C3dPose(R, t) {
        if(IS_DEBUG) CHECK(!zero(t.sum_square() - 1) && t.sum_square() != 0, "C3dNormalisedPose: t not normalised");
    }

    C3dPose scale(double dScale) const {
        return C3dPose(R, dScale * t);
    }

    C3dNormalisedPose reverse() const //Todo: test on simn data
    {
        //C3dNormalisedPose pose(R.t(), -t);
        C3dNormalisedPose pose(R.t(), -(R * t)); //Experimental change 8-2-10, appears to work better now
        return pose;
    }
};

class CRelPoseSD //Uncertainty in a relative pose with *unknown length*
{
    double dRelOrientationSD, dCameraMotionAngleSD, dBaselineFromPointDepthMean, dBaselineFromPointDepthVar;
public:
    //Errors in reconstructing cameras (Lev-Mar for Sampson error) from these calibrated image points.
    CRelPoseSD(const TPointVec2d & m0, const TPointVec2d & m1, const CCamCalibMatrix & K, const double dFeatureLocalisationAccuracy, double dPointDepthMean = -1, double dPointDepthSD = -1);

    CRelPoseSD() : dRelOrientationSD(-HUGE), dCameraMotionAngleSD(-HUGE), dBaselineFromPointDepthMean(-1), dBaselineFromPointDepthVar(-1) {
    }

    void makeUninformative() {
        dRelOrientationSD = M_PI_2, dCameraMotionAngleSD = M_PI_2, /* Means we don't trust the estimate */ dBaselineFromPointDepthMean = -1, dBaselineFromPointDepthVar = -1;
    }

    double relOrientationSD() const {
        if(IS_DEBUG) CHECK(!init(), "relOrientationSD: Uninit orientation")
        return dRelOrientationSD;
    }

    double cameraMotionAngleSD() const {
        if(IS_DEBUG) CHECK(!init(), "relOrientationSD: Uninit orientation")
        return dCameraMotionAngleSD;
    }

    bool init() const {
        return dRelOrientationSD > 0 && dCameraMotionAngleSD > 0;
    } //0 not allowed

    double baselineMean() const {
        return dBaselineFromPointDepthMean;
    }

    double baselineVar() const {
        return dBaselineFromPointDepthVar;
    }

    bool usePointDepths() const {
        return dBaselineFromPointDepthVar >= -1;
    }

    bool havePointDepth() const {
        return dBaselineFromPointDepthVar > 0;
    }

};

class C3dNormalisedPoseWithSD : public C3dNormalisedPose {
public:
    const CRelPoseSD SD;

    C3dNormalisedPoseWithSD(const C3dNormalisedPose & P, const CRelPoseSD & SD) : C3dNormalisedPose(P), SD(SD) {
        if(IS_DEBUG) CHECK(!SD.init(), "C3dNormalisedPoseWithSD: Uninit SD")
    }

    C3dNormalisedPoseWithSD(const C3dRotation R, const C3dPoint t, const CRelPoseSD & SD) : C3dNormalisedPose(R, t), SD(SD) {
        if(IS_DEBUG) CHECK(!SD.init(), "C3dNormalisedPoseWithSD: Uninit SD")
    }

    C3dNormalisedPoseWithSD reverse() const {
        return C3dNormalisedPoseWithSD(C3dNormalisedPose::reverse(), SD);
    }
};

//#define N_CAM_MATS 4
#define CAM_SELECT_ERROR -2

static const C3dPoint Origin(0, 0, 0);
static const C3dRotation R0;


//Select vector elements based on a mask
template<class T>
void maskVector(const CMask & mask, const CDynArray<T> & vIn1, CDynArray<T> & vOut1) {
    size_t size = (size_t)mask.size();
    if(IS_DEBUG) CHECK(vIn1.size() != size, "Vector sizes don't match");
    if(IS_DEBUG) CHECK(vOut1.size() != 0, "Vector not empty");

    CMask::const_iterator pMask = mask.begin();
    typename CDynArray<T>::const_iterator pIn1 = vIn1.begin();

    for (; pIn1 != vIn1.end(); pIn1++, pMask++) {
        if (*pMask) {
            vOut1.push_back(*pIn1);
        }
    }
    if(IS_DEBUG) CHECK(vOut1.size() != mask.countInliers(), "Vector wrong size");
}

template<class T>
void mask2Vectors(const CMask & mask, const CDynArray<T> & vIn1, const CDynArray<T> & vIn2, CDynArray<T> & vOut1, CDynArray<T> & vOut2) {
    int size = mask.size();
    if(IS_DEBUG) CHECK(vIn1.size() != size, "Vector sizes don't match");
    if(IS_DEBUG) CHECK(vIn2.size() != size, "Vector sizes don't match");
    if(IS_DEBUG) CHECK(vOut1.size() != 0, "Vector not empty");
    if(IS_DEBUG) CHECK(vOut2.size() != 0, "Vector not empty");
    vOut1.reserve(size);
    vOut2.reserve(size);

    CMask::const_iterator pMask = mask.begin();
    typename CDynArray<T>::const_iterator pIn1 = vIn1.begin();
    typename CDynArray<T>::const_iterator pIn2 = vIn2.begin();

    for (; pIn1 != vIn1.end(); pIn1++, pIn2++, pMask++) {
        if (*pMask) {
            vOut1.push_back(*pIn1);
            vOut2.push_back(*pIn2);
        }
    }
    if(IS_DEBUG) CHECK(vOut1.size() != mask.countInliers(), "Vector wrong size");
    if(IS_DEBUG) CHECK(vOut2.size() != mask.countInliers(), "Vector wrong size");
}

void calibrateBoWCorr(const CCamCalibMatrix & K_inv, const CBoWCorrespondences * pCorr, CPointVec2d & aCalibratedPoints1, CPointVec2d & aCalibratedPoints2, CInlierProbs & adArrLikelihood, CPointIdentifiers & pointIds, const bool CORRECT_RD);
bool decomposeHomography(const Eigen::Matrix3d & H_in, C3dRotation & rotation, C3dPoint & planeNormal, C3dPoint & camMotion);

std::ostream& operator<<(std::ostream& s, const C3dPose & X);

#endif
