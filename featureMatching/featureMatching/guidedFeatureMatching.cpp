/* 
 * File:   guidedFeatureMatching.cpp
 * Author: tom
 * 
 * Created on 9 March 2012, 10:24 AM
 */

#include <boost/smart_ptr/scoped_array.hpp>
#include <boost/smart_ptr/scoped_ptr.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/hypergeometric.hpp>

#include "guidedFeatureMatching.h"
#include "geom/mvNormalSampler.h" // featureDescription library needs to have ../cameraGeom on the include path
#include "geom/geom_eigen.h"
#include "ransac/refineEOnRTManifold.h"
#include "opencv2/opencv.hpp"
#include <Eigen/Dense>


//#include "../pruningSoftware/src/geometry.h"
double closestDistance(const Eigen::Vector2d & A, const Eigen::Vector2d & B, const Eigen::Vector2d & point) {
	Eigen::Vector2d AP = point - A;
	Eigen::Vector2d AB = B - A;
	double ab2 = AB.squaredNorm();
	double ap_ab = AP.dot(AB); // AP.x*AB.x + AP.y*AB.y;
	double t = ap_ab / ab2;
	const bool segmentClamp = true;
	if (segmentClamp) {
		if (t < 0.0) t = 0.0;
		else if (t > 1.0) t = 1.0;
	}
	Eigen::Vector2d closestPointOnLine = A + AB * t;

	return (point - closestPointOnLine).norm();
}

template<bool bSquared>
double sampsonsErr(const Eigen::Vector3d & x, const Eigen::Vector3d & xp, Eigen::Matrix3d & E) {
    Eigen::Vector3d xpE = xp.transpose() * E;
    Eigen::Vector3d Ex = E * x;

    double numerator = xpE.dot(x);
    double denom = xpE.segment(0, 2).squaredNorm() + Ex.segment(0, 2).squaredNorm();

    if(denom==0)
    {
        cout << "E: " << E << endl;
        cout << "x: " << x.transpose() << endl;
        cout << "xp: " << xp.transpose() << endl;
        THROW("Bad essential matrix or points");
    }
    
    if (bSquared)
        return sqr(numerator) / denom;
    else
        return numerator / sqrt(denom);
}

Eigen::Matrix3d phiThetaPsiToRotationMat(const Eigen::Vector3d & phiThetaPsi)
{
    Eigen::Matrix3d Rx, Rz1, Rz2;
    
    double phi=phiThetaPsi(0);
    double c=cos(phi), s=sin(phi);
    Rz1 << c, -s,0,
           s, c, 0,
           0, 0, 1;

    double theta=phiThetaPsi(1);
    c=cos(theta), s=sin(theta);
    Rx << 1, 0, 0,
          0, c, -s,
          0, s, c;
    
    double psi=phiThetaPsi(2);
    c=cos(psi), s=sin(psi);
    Rz2 << c, -s,0,
           s, c, 0,
           0, 0, 1;    
    
 
    Eigen::Matrix3d R = Rz1*Rx*Rz2; //TODO: what axes, what order? these came from http://en.wikipedia.org/wiki/Euler_angles#Euler_angles_as_composition_of_extrinsic_rotations

    return R;
}

C3dRotation phiThetaPsiToQuat(const Eigen::Vector3d & phiThetaPsi)
{
    Eigen::Matrix3d R = phiThetaPsiToRotationMat(phiThetaPsi);
    C3dRotation q(R);
    return q;
}




void calibrateAll(const CMatchableDescriptors * pDS, const CCamCalibMatrix & K, T2dPoints & aCalibratedPoints1)
{
    aCalibratedPoints1.clear();
    aCalibratedPoints1.reserve(pDS->Count());
    for (int i=0; i<pDS->Count(); i++)
    {
        const CLocation loc1 = pDS->get_const(i)->location();

        CSimple2dPoint p(loc1);

        p.calibrate(K, true);

        aCalibratedPoints1.push_back(p);
    }
}

class CLikelihoodFn
{
public:
    virtual void setupMatch(const Eigen::Vector3d & loc1) = 0;
    virtual bool canMatch(const Eigen::Vector3d & loc1, const Eigen::Vector3d & loc2, const double p) const = 0;
        
    virtual void visualise() 
    {
        const Eigen::Vector3d loc1(0,0,1);
        
        setupMatch(loc1);
        
        double dRad = IS_DEBUG ? 50 : 250;
        
        cv::Mat M(dRad*2, dRad*2, CV_8UC3, cv::Scalar());
        
        for(int r=0;r<2*dRad;r++)
            for(int c=0;c<2*dRad;c++)
            {
                const Eigen::Vector3d loc2((r-dRad)/dRad,(c-dRad)/dRad,1);
                if(canMatch(loc1, loc2, 0.01))
                {
                    M.at<cv::Vec3b>(r,c) = cv::Vec3b(255,0,0);
                    if(canMatch(loc1, loc2, 0.05))
                    {
                        M.at<cv::Vec3b>(r,c) = cv::Vec3b(255,255,0);
                        if(canMatch(loc1, loc2, 0.1))
                        {
                            M.at<cv::Vec3b>(r,c) = cv::Vec3b(255,0,255);
                        }
                    }
                }
            }
        
        cv::imshow("matchableRegion", M);
        cv::waitKey(0);
    }
//    virtual void visualise() {}

};

class CMatchAll : public CLikelihoodFn
{
public:
    virtual void setupMatch(const Eigen::Vector3d & loc1) {}
    virtual bool canMatch(const Eigen::Vector3d & loc1, const Eigen::Vector3d & loc2, const double p) const {return true; }
};

// _hat denotes means of distributions
class CLinearisedLikelihoodFn : public CLikelihoodFn
{
    typedef Eigen::Matrix<double, 9, 1> TEVec;
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> E_hat, DE[6], A_xyz[3]; //Borrow A_x, A_y, A_z notation from Wikipedia: http://en.wikipedia.org/wiki/Rotation_matrix#Lie_algebra
    TEVec E_hat_vec; //=E_hat as a vector
    
    typedef Eigen::Matrix<double, 9, 9> TJ_C_JT;
    typedef Eigen::Matrix<double, 9, 6> TJ_cov;
    typedef Eigen::Matrix<double, 3, 9> TX_cov;
    TJ_C_JT J_C_JT;
    Eigen::Vector3d l_hat,l_norm;
    Eigen::Matrix3d l_cov;
public:
    CLinearisedLikelihoodFn(const Eigen::Vector3d & t, const Eigen::Vector3d & phiThetaPsi, const CGuidedFeatureMatching::TCovMat66 & C)
    {
        Eigen::Matrix3d R_hat = phiThetaPsiToRotationMat(phiThetaPsi);
        Eigen::Matrix3d T_cross;
        const Eigen::Vector3d translationSample = t / t.norm();
        Xmat(translationSample, T_cross);
        //makeE(R_hat, translationSample, E_hat);
        E_hat = T_cross * R_hat;
        
        E_hat_vec = TEVec(E_hat.data());
        
        A_xyz[0] << 0,  0, 0,
                    0,  0, -1,
                    0,  1, 0;
        
        A_xyz[1] << 0, 0,  1,
                    0, 0,  0,
                    -1, 0,  0;
        
        A_xyz[2] << 0, -1, 0,
                    1,  0, 0,
                    0,  0, 0;
                
        for(int i=0;i<3;i++)
        {
            DE[i] = A_xyz[i] * R_hat;
            DE[i+3] = T_cross * A_xyz[i];//Only if phiThetaPsi are Rx,Ry and Rz
        }
        TJ_cov J;
        for(int i=0;i<6;i++)
            J.col(i) = TEVec(DE[i].data());
        
        /*epsilon has covariance C, mean 0
         * E = E_hat + J*epsilon
         * E ~ N(E_hat, J*C*J^T
        */
        
        J_C_JT = J*C*J.transpose();
    }
    /*To restore this must change rowMajor->colMajor
    virtual void setupMatch(const Eigen::Vector3d & loc2)
    {
        TX_cov Xp = TX_cov::Zero();
        Xp.block<1,3>(0,0) = loc2.transpose();
        Xp.block<1,3>(1,3) = loc2.transpose();
        Xp.block<1,3>(2,6) = loc2.transpose();
        l_hat = Xp*E_hat_vec; 
        l_cov = Xp * J_C_JT *Xp.transpose();
        
        Eigen::Matrix3d P = Eigen::Matrix3d::Identity();
        l_norm=l_hat/l_hat.norm();
        
        for(int i=0;i<3;i++)
            P.row(i) -= l_norm(i)*l_norm.transpose(); //Todo: check me
        
        Eigen::Matrix3d l_norm_cov = P*l_cov*P.transpose();
        if(IS_DEBUG) CHECK((P*l_hat).squaredNorm() > 1e-5, "P setup failed");
        if(IS_DEBUG) CHECK((P*l_norm).squaredNorm() > 1e-5, "P setup failed");
    }*/
    
    virtual void setupMatch(const Eigen::Vector3d & loc1)
    {
        TX_cov X = TX_cov::Zero();
        X.block<1,3>(0,0) = loc1.transpose();
        X.block<1,3>(1,3) = loc1.transpose();
        X.block<1,3>(2,6) = loc1.transpose();
        l_hat = X*E_hat_vec; 
        l_cov = X * J_C_JT *X.transpose();
        
        Eigen::Matrix3d P = Eigen::Matrix3d::Identity();
        l_norm=l_hat/l_hat.norm();
        
        for(int i=0;i<3;i++)
            P.row(i) -= l_norm(i)*l_norm.transpose(); //Todo: check me
        
        Eigen::Matrix3d l_norm_cov = P*l_cov*P.transpose();
        if(IS_DEBUG) CHECK((P*l_hat).squaredNorm() > 1e-5, "P setup failed");
        if(IS_DEBUG) CHECK((P*l_norm).squaredNorm() > 1e-5, "P setup failed");
    }
    
    virtual bool canMatch(const Eigen::Vector3d & loc1, const Eigen::Vector3d & loc2, const double p) const
    {
        double dTestStat = l_norm.dot(loc1);
        double dVar = loc1.transpose()*l_cov*loc1; //TODO: this is wrong because dTestStat=mean of this distn.
        double dNormal01TestStat = dTestStat / sqrt(dVar);
        
        static const boost::math::normal_distribution<double> ND(0,1);
        static const double P_THRESH = boost::math::cdf(ND, 1-p);
        
        return dNormal01TestStat < fabs(P_THRESH);
    }
};

// _hat denotes means of distributions
class CLinearisedLikelihoodFn_H : public CLikelihoodFn
{
    typedef Eigen::Matrix<double, 9, 1> THVec;
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> H_hat, DH[6], A_xyz[3]; //Borrow A_x, A_y, A_z notation from Wikipedia: http://en.wikipedia.org/wiki/Rotation_matrix#Lie_algebra
    THVec H_hat_vec; //=E_hat as a vector
    
    typedef Eigen::Matrix<double, 9, 9> TJ_C_JT;
    typedef Eigen::Matrix<double, 9, 9> TJ_cov;
    typedef Eigen::Matrix<double, 3, 9> TX_cov;
    TJ_C_JT J_C_JT;
    Eigen::Vector3d l_hat,l_homo;
    Eigen::Matrix3d l_homo_cov, l_homo_info;
public:
    CLinearisedLikelihoodFn_H(const Eigen::Vector3d & t, const Eigen::Vector3d & phiThetaPsi,const Eigen::Vector3d & n, const CGuidedFeatureMatching::TCovMat99 & C)
    {
        Eigen::Matrix3d R_hat = phiThetaPsiToRotationMat(phiThetaPsi);

        H_hat = makeH(t, phiThetaPsi, n);
        
        H_hat_vec = THVec(H_hat.data());
        
        A_xyz[0] << 0,  0, 0,
                    0,  0, -1,
                    0,  1, 0;
        
        A_xyz[1] << 0, 0,  1,
                    0, 0,  0,
                    -1, 0,  0;
        
        A_xyz[2] << 0, -1, 0,
                    1,  0, 0,
                    0,  0, 0;
                
        for(int i=0;i<3;i++)
        {
            //dH/dt1 ...
            DH[i].setZero(); 
            DH[i].row(i) = n.transpose();
            
            //dH/dphi dH/dtheta dH/dpsi
            DH[i+3] = A_xyz[i];

            //dH/dt1 ...
            DH[i].setZero();
            DH[i].col(i) = t;
        }
        
        TJ_cov J;
        for(int i=0;i<6;i++)
            J.col(i) = THVec(DH[i].data());
        
        J_C_JT = J*C*J.transpose();
    }
    
    //H * loc1 = loc2
    virtual void setupMatch(const Eigen::Vector3d & loc1)
    {
        TX_cov X = TX_cov::Zero();
        X.block<1,3>(0,0) = loc1.transpose();
        X.block<1,3>(1,3) = loc1.transpose();
        X.block<1,3>(2,6) = loc1.transpose();
        l_hat = X*H_hat_vec; 
        //cout << "l_hat: " << l_hat << endl;
        //cout << "H*x: " << H_hat*loc1 << endl;
        Eigen::Matrix3d l_cov = X * J_C_JT *X.transpose();
        
        l_homo=l_hat/l_hat.z();
        
        l_homo_cov = l_cov/sqr(l_hat.z());
        l_homo_info = l_homo_cov.inverse();
    }
    
    virtual bool canMatch(const Eigen::Vector3d & loc1, const Eigen::Vector3d & loc2, const double p) const
    {
        Eigen::Vector3d dev = loc2 - l_homo;
        const double dDeviationSD_sq = dev.transpose() * l_homo_info * dev;
        
        static const boost::math::normal_distribution<double> ND(0,1);
        //TODO: Don't need to recompute this every time
        const double P_THRESH = fabs(boost::math::quantile(ND, 1-p)); // quantile = icdf
        
        return dDeviationSD_sq < P_THRESH;
    }
};

class CKernalDensityLikelihoodFn : public CLikelihoodFn
{
protected:
    const int NUM_SAMPLES;
    boost::scoped_array<Eigen::Matrix3d> aEssentialMatrices; //can't use a vector as must be aligned
    boost::scoped_array<CCamera> aCameras; //can't use a vector as must be aligned
    boost::scoped_array< Eigen::Vector2d > aEpilineStart, aEpilineEnd;
    boost::scoped_array<bool> abCanMatch;

    double dBandwidth;
    const bool bConstrainDepths;
public:
    CKernalDensityLikelihoodFn(const Eigen::Vector3d & t, const Eigen::Vector3d & phiThetaPsi, const CGuidedFeatureMatching::TCovMat66 & fullCovariance, const int NUM_SAMPLES, const double dBandwidth, const bool bConstrainDepths_in) 
    : NUM_SAMPLES(NUM_SAMPLES), aEssentialMatrices(new Eigen::Matrix3d[NUM_SAMPLES]), aCameras(!bConstrainDepths_in ? 0 : (new CCamera[NUM_SAMPLES])), aEpilineStart(new Eigen::Vector2d[NUM_SAMPLES]), aEpilineEnd(new Eigen::Vector2d[NUM_SAMPLES]), abCanMatch(new bool[NUM_SAMPLES]), dBandwidth(dBandwidth), bConstrainDepths(bConstrainDepths_in)
    {
        //First make a set of essential matrices
        CMVNormalSampler<CGuidedFeatureMatching::TCovMat66>::VecType fullState;
        fullState.head(3)= t; fullState.tail(3) = phiThetaPsi;

        CMVNormalSampler<CGuidedFeatureMatching::TCovMat66> normalSampler(fullState, fullCovariance);

        for (int nSample = 0; nSample < NUM_SAMPLES; nSample++) {
            CMVNormalSampler<CGuidedFeatureMatching::TCovMat66>::VecType fullSample = normalSampler.sample();

            Eigen::Vector3d phiThetaPsiSample = fullSample.tail(3), translationSample=fullSample.head(3);

            Eigen::Matrix3d Rsample = phiThetaPsiToRotationMat(phiThetaPsiSample);
            C3dRotation q_sample(Rsample);

            translationSample /= translationSample.norm();

            makeE(q_sample, translationSample, aEssentialMatrices[nSample]);
            
            if(bConstrainDepths)
                aCameras[nSample] = q_sample | translationSample;
                
            //cout << "E: " << aEssentialMatrices[nSample] << endl;
        }
    }

    virtual void setupMatch(const Eigen::Vector3d & loc1) 
    {
        if(bConstrainDepths)
        {
            const bool bCheckAndResetBandwidth = true;
            CDynArray<double> aNearMagnitude, aFarMagnitude;
            
            if(bCheckAndResetBandwidth)
            {
                aNearMagnitude.reserve(NUM_SAMPLES);
                aFarMagnitude.reserve(NUM_SAMPLES);
            }
            
            for (int nSample = 0; nSample < NUM_SAMPLES; nSample++) {
                const CCamera & Pp = aCameras[nSample];
                double dMaxDepth = 1000; //near infinity
                double dMinDepth = 0.5; //relative to camera motion
                
                const Eigen::Vector3d X_far_vec = loc1*dMaxDepth;
                const Eigen::Vector3d X_near_vec = loc1*dMinDepth;
                C3dPoint X_far(X_far_vec);
                C3dPoint X_near(X_near_vec);
                C2dPoint xp_far = X_far.photo(Pp);
                C2dPoint xp_near = X_near.photo(Pp);
                const Eigen::Vector2d xp_far_vec = xp_far.asVector();
                const Eigen::Vector2d xp_near_vec = xp_near.asVector();
                //const Eigen::Vector2d xp = loc2.head(2);
                
                // Find X where points pass from in front of to behind the 2nd camera
                // R.row(z) * x * lambda + t_z =0
                Eigen::Vector3d R_row_z(Pp.at(2,0), Pp.at(2,1), Pp.at(2,2));
                const double lambda = -Pp.at(2,3)/loc1.dot(R_row_z);
                //cout << "Depth where point passes behind cam: " << lambda << endl;
                Eigen::Vector2d xp_mid_vec_lo, xp_mid_vec_hi;
                if(lambda>dMinDepth)
                {
                    const Eigen::Vector3d X_mid_vec_lo = 0.99*lambda*loc1;
                    C3dPoint X_mid_lo(X_mid_vec_lo);
                    C2dPoint xp_mid_lo = X_mid_lo.photo(Pp);
                    xp_mid_vec_lo = xp_mid_lo.asVector();

                    const Eigen::Vector3d X_mid_vec_hi = 1.01*lambda*loc1;
                    C3dPoint X_mid_hi(X_mid_vec_hi);
                    C2dPoint xp_mid_hi = X_mid_hi.photo(Pp);
                    xp_mid_vec_hi = xp_mid_hi.asVector();
                }

                const bool bNearInFront = X_near.testInFront(Pp);
                const bool bFarInFront = X_far.testInFront(Pp);
                
                abCanMatch[nSample] = true;

                if(bFarInFront)
                {
                    aEpilineStart[nSample] = bNearInFront ? xp_near_vec : xp_mid_vec_hi;
                    aEpilineEnd[nSample] = xp_far_vec;
                }
                else if(bNearInFront)
                {
                    aEpilineStart[nSample] = xp_near_vec;
                    aEpilineEnd[nSample] = xp_mid_vec_lo;
                }
                else
                {
                    cout << "Neither point in front of 2nd camera" << endl;
                    abCanMatch[nSample] = false; //TODO: should we only count the ones in front for the CDF?
                }
                
                if(bCheckAndResetBandwidth)
                {
                    aNearMagnitude.push_back(aEpilineStart[nSample].norm());
                    aFarMagnitude.push_back(aEpilineEnd[nSample].norm());
                }
            }

            if(bCheckAndResetBandwidth)
            {
                double dMean, dNearVar, dFarVar;
                aNearMagnitude.getMeanVar(dMean, dNearVar);
                aFarMagnitude.getMeanVar(dMean, dFarVar);
                
                const double dSDMax = sqrt(std::max<double>(dNearVar, dFarVar));
                
                cout << "Bandwidth: " << dBandwidth << " SD max: " << dSDMax << endl;
                dBandwidth = 10*dSDMax/NUM_SAMPLES;
                cout << "Setting bandwidth to " << dBandwidth << endl;
            }
            
        }
    }    
    
    virtual bool canMatch(const Eigen::Vector3d & loc1, const Eigen::Vector3d & loc2, const double p) const
    {
        double dCDF = 0;
        const double dBandwidth_inv = 1.0/(dBandwidth);

        boost::math::normal_distribution<double> ND(0, dBandwidth);

        for (int nSample = 0; nSample < NUM_SAMPLES; nSample++) {
            double dSE = HUGE;
            if(bConstrainDepths)
            {
                const Eigen::Vector2d xp = loc2.head(2);
                
                if(abCanMatch[nSample])
                    dSE = closestDistance(aEpilineStart[nSample], aEpilineEnd[nSample], xp);
            }
            else
            {
                dSE = sampsonsErr<false>(loc1, loc2, aEssentialMatrices[nSample]);
            }

            const double d1CDF = boost::math::cdf(ND, dSE);
            dCDF += d1CDF;

        }
        dCDF /= NUM_SAMPLES;

        if(IS_DEBUG) CHECK(dCDF > 1.000001, "CDF OOB");

        return dCDF > p && dCDF < (1-p);
    }
    
    /*virtual void visualise() 
    {
        const Eigen::Vector3d loc1(0,0,1);
        
        setupMatch(loc1);
        
        double dRad = IS_DEBUG ? 50 : 250;
        
        cv::Mat M(dRad*2, dRad*2, CV_8UC3, cv::Scalar());
        
        for(int r=0;r<2*dRad;r++)
            for(int c=0;c<2*dRad;c++)
            {
                const Eigen::Vector3d loc2((r-dRad)/dRad,(c-dRad)/dRad,1);
                if(canMatch(loc1, loc2, 0.01))
                {
                    M.at<cv::Vec3b>(r,c) = cv::Vec3b(255,0,0);
                    if(canMatch(loc1, loc2, 0.05))
                    {
                        M.at<cv::Vec3b>(r,c) = cv::Vec3b(255,255,0);
                        if(canMatch(loc1, loc2, 0.1))
                        {
                            M.at<cv::Vec3b>(r,c) = cv::Vec3b(255,0,255);
                        }
                    }
                }
            }
        
        cv::imshow("matchableRegion", M);
        cv::waitKey(0);
    }*/
};

Eigen::Matrix3d makeH(const Eigen::Vector3d & translationSample, const Eigen::Vector3d & phiThetaPsiSample, const Eigen::Vector3d & normalSample)
{
    Eigen::Matrix3d R = phiThetaPsiToRotationMat(phiThetaPsiSample);
    return R-translationSample*normalSample.transpose();
}

class CKernalDensityLikelihoodFn_H : public CLikelihoodFn
{
protected:
    const int NUM_SAMPLES;
    boost::scoped_array<Eigen::Matrix3d> aHMatrices; //can't use a vector as must be aligned
    boost::scoped_array<Eigen::Vector2d> aHxVectors;
    boost::scoped_array<CCamera> aCameras; //can't use a vector as must be aligned

    double dBandwidth;
public:
    CKernalDensityLikelihoodFn_H(const Eigen::Vector3d & t, const Eigen::Vector3d & phiThetaPsi, const Eigen::Vector3d & n, const CGuidedFeatureMatching::TCovMat99 & fullCovariance, const int NUM_SAMPLES, const double dBandwidth) 
    : NUM_SAMPLES(NUM_SAMPLES), aHMatrices(new Eigen::Matrix3d[NUM_SAMPLES]), aHxVectors(new Eigen::Vector2d[NUM_SAMPLES]), dBandwidth(dBandwidth)
    {
        //First make a set of essential matrices
        CMVNormalSampler<CGuidedFeatureMatching::TCovMat99>::VecType fullState;
        fullState.head(3)= t; fullState.segment<3>(3) = phiThetaPsi; fullState.tail(3) = n;

        CMVNormalSampler<CGuidedFeatureMatching::TCovMat99> normalSampler(fullState, fullCovariance);

        for (int nSample = 0; nSample < NUM_SAMPLES; nSample++) {
            CMVNormalSampler<CGuidedFeatureMatching::TCovMat99>::VecType fullSample = normalSampler.sample();

            Eigen::Vector3d phiThetaPsiSample = fullSample.segment<3>(3), translationSample=fullSample.head(3), normalSample=fullSample.tail(3);

            aHMatrices[nSample] = makeH(translationSample, phiThetaPsiSample, normalSample);
                
            REPEAT(3,cout << "H: " << aHMatrices[nSample] << endl);
        }
    }
    
    void setBandwidth() 
    {
        Eigen::Vector2d total = Eigen::Vector2d::Zero(), totalSq = Eigen::Vector2d::Zero();
        
        for (int nSample = 0; nSample < NUM_SAMPLES; nSample++) {
            total += aHxVectors[nSample];
            totalSq += aHxVectors[nSample].array().square().matrix();
        }
        
        Eigen::Vector2d var = totalSq/NUM_SAMPLES - (total/NUM_SAMPLES).array().square().matrix();
        
        const double dArea = sqrt(var.x()*var.y());
        
        dBandwidth = 50*sqrt(dArea/NUM_SAMPLES);
        
        cout << "Set bandwidth = " << dBandwidth << endl;
    }    

    virtual void setupMatch(const Eigen::Vector3d & loc1) 
    {
        for (int nSample = 0; nSample < NUM_SAMPLES; nSample++) {
            aHxVectors[nSample] = homog(aHMatrices[nSample] * loc1);
        }
        
        setBandwidth(); //TODO: don't need to call this every time
    }
    
    virtual bool canMatch(const Eigen::Vector3d & loc1, const Eigen::Vector3d & loc2, const double p) const
    {
        double dCDF = 0;
        const double dBandwidth_inv = 1.0/(dBandwidth);

        boost::math::normal_distribution<double> ND(0, dBandwidth);

        for (int nSample = 0; nSample < NUM_SAMPLES; nSample++) {

            const double dDist = (aHxVectors[nSample] - loc2.head<2>()).norm();
            
            REPEAT(5,
            cout << "aHxVectors[nSample]" << aHxVectors[nSample].transpose() << endl;
            cout << loc1.transpose() << endl;
            cout << loc2.transpose() << endl);
            
            const double d1CDF = boost::math::cdf(ND, dDist);
            dCDF += d1CDF;

        }
        dCDF /= (NUM_SAMPLES * cube(dBandwidth)); //NOT A CDF --todo--compute confidence bound properly

        //if(IS_DEBUG) CHECK(dCDF > 1.000001, "CDF OOB");

        return dCDF > p/* && dCDF < (1-p)*/;
    }
};

const CBoWCorrespondences * CGuidedFeatureMatching::guidedFeatureMatch_Homography(const CMatchableDescriptors * pDS1, const CMatchableDescriptors * pDS2, const CCamCalibMatrix & K, const CMatchableDescriptors::CMatchSettings & MS, const Eigen::Vector3d & t, const Eigen::Vector3d & phiThetaPsi, const Eigen::Vector3d & n, const CGuidedFeatureMatching::TCovMat99 & fullCovariance)
{
    const double dBandwidth = 0.01; //radians
    boost::scoped_ptr<CLikelihoodFn> pFn;
    //pFn.reset( new CKernalDensityLikelihoodFn_H( t,  phiThetaPsi, n, fullCovariance, 2500, dBandwidth) );
    pFn.reset(new CLinearisedLikelihoodFn_H( t,  phiThetaPsi, n, fullCovariance) );
    return guidedFeatureMatch_int(pFn.get(), pDS1, pDS2, K, MS);
}

const CBoWCorrespondences * CGuidedFeatureMatching::guidedFeatureMatch(const CMatchableDescriptors * pDS1, const CMatchableDescriptors * pDS2, const CCamCalibMatrix & K, const CMatchableDescriptors::CMatchSettings & MS, const Eigen::Vector3d & t, const Eigen::Vector3d & phiThetaPsi, const TCovMat66 & fullCovariance) 
{
    const double dBandwidth = 0.005; //radians
    boost::scoped_ptr<CLikelihoodFn> pFn;
    pFn.reset( new CKernalDensityLikelihoodFn( t,  phiThetaPsi,  fullCovariance, 250, dBandwidth, true) );
    //pFn.reset( new CMatchAll );
    return guidedFeatureMatch_int(pFn.get(), pDS1, pDS2, K, MS);
}
    
const CBoWCorrespondences * CGuidedFeatureMatching::guidedFeatureMatch_int(CLikelihoodFn * pFn, const CMatchableDescriptors * pDS1, const CMatchableDescriptors * pDS2,  const CCamCalibMatrix & K, const CMatchableDescriptors::CMatchSettings & MS)     
{
    T2dPoints aCalibratedPoints1, aCalibratedPoints2;
    calibrateAll(pDS1, K, aCalibratedPoints1);
    calibrateAll(pDS2, K, aCalibratedPoints2);

    //pFn->visualise();
   
    //Now the brute force match:
    const double p=0.01;
    
    CBoWCorrespondences * pCorr = new CBoWCorrespondences();
    const int nCount1 = pDS1->Count();
    const int nCount2 = pDS2->Count();
    ARRAYZ(int, anMatchesLeft, nCount1);
    ARRAYZ(int, anMatchesRight, nCount2);
    
    ARRAY(int, anRightIdxMatchingLeft, nCount1);
    setConstant(anRightIdxMatchingLeft, -1, nCount1);
    
    int nNumPairs=0, nNumMatchablePairs=0;
 
    for (int i = 0; i < nCount1; i++) {

        Eigen::Vector3d point1homog(aCalibratedPoints1[i].getX(), aCalibratedPoints1[i].getY(), 1);
        const CDescriptor * pDesc1 = pDS1->get_const(i);
        CLocation loc1 = pDesc1->location();
        
        CDescriptor::TDist closestDist = MAX_ALLOWED_DIST;
        int nClosestMatch = -1;
        //CLocation bestLoc2;
        
        pFn->setupMatch(point1homog);
        for (int j = 0; j < nCount2; j++) {

            Eigen::Vector3d point2homog(aCalibratedPoints2[j].getX(), aCalibratedPoints2[j].getY(), 1);
            const CDescriptor * pDesc2 = pDS2->get_const(j);
            const CLocation loc2 = pDesc2->location();            
        
            nNumPairs++;
            
            if(pFn->canMatch(point1homog, point2homog, p))
            {
                nNumMatchablePairs++;
                CDescriptor::TDist dist = pDesc1->distance(pDesc2);
            
                if(dist < closestDist)
                {
                    closestDist = dist;
                    nClosestMatch = j;
                    //bestLoc2 = loc2;
                }
            }            
        }
        
        if(nClosestMatch >=0 && anMatchesRight[nClosestMatch] < NN_MAX) //anMatchesLeft[nClosestMatch] < NN_MAX is just a safety check that sholdn't actually happen often
        {
            anMatchesRight[nClosestMatch]++;
            anMatchesLeft[i]++;
            anRightIdxMatchingLeft[i] = nClosestMatch;
        }
    }    

    //Now convert matches to matches with prior probabilities
    for (int i = 0; i < nCount1; i++) {
        const CDescriptor * pDesc1 = pDS1->get_const(i);
        const CLocation loc1 = pDesc1->location();
        
        if(anRightIdxMatchingLeft[i] >=0)
        {
            const CDescriptor * pDesc2 = pDS2->get_const(anRightIdxMatchingLeft[i]);
            const CLocation loc2 = pDesc2->location();

            const int nMaxNN = std::max<int> (anMatchesRight[anRightIdxMatchingLeft[i]], anMatchesLeft[i]); //= nMatchesLeft[anLeftIdxMatchingRight[j]] in practice

            double dPriorProb = MS.PP / nMaxNN;

            pCorr->push_back(CCorrespondence(loc1, loc2, dPriorProb));
        }
    }

    cout << pCorr->size() << " features matched, " << nNumMatchablePairs << "/" << nNumPairs << "=" << ((double)nNumMatchablePairs / nNumPairs) << endl;
    return pCorr;
}