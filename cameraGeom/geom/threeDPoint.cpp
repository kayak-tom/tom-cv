/*
 * threeDPoint.cpp
 *
 *  Created on: 25 Nov 2010
 *      Author: tom
 */

#include "vectorPoints.h"
#include "lines.h"

#include <Eigen/LU>
#include <geom/mvNormalSampler.h>
#include <boost/math/distributions/chi_squared.hpp>
#include <opencv2/core/core.hpp>
#include <fstream>
#include <boost/filesystem.hpp>
#include <image/drawEllipse.h>

#define INSTANTIATE_CAMERAS_IPP //also consistency checks
#include "cameras.ipp"

using namespace std;
using namespace Eigen;

template<typename TM> void EigenNormalise(TM & M) {
    M /= M.stableNorm();
}

/*void C3dCovariance::checkCovarianceInit() const {
    //det is product of eigenvectors. 1e-18 is a sd of 1e-3 (1mm)
    if (squaredNorm() != 0 && (determinant() < 1e-24 || determinant() > 0.1)) {
        cout << "Covariance:\n" << *this << endl;
        THROW("Determinant of covariance outside reasonable range");
    }
    
    const bool bCheckAlignedWithMotion = true;
    if(bCheckAlignedWithMotion)
    {
        const double dPll = sys().TRACTOR_DIRECTION().transpose() * *this * sys().TRACTOR_DIRECTION();
        const double dPerp = sys().VINES_NORMAL().transpose() * *this * sys().VINES_NORMAL();
        
        CHECK_P(dPll < 0.99*dPerp, *this, "Covariance is greater up-down than side-to-side (this should be unlikely as main uncertainty is how far the tractor has moved)");
    }
}

void C2dCovariance::checkCovarianceInit() const {
    //0.01 is about 0.3px accuracy. Determinants can be pretty high (from e.g. predicted endpoint location covariances)
    if (squaredNorm() != 0)
    {
        const double dPx = covToSDPixels(*this);
        if(dPx < 1 || dPx > 800) { //TB: Changed this to 800 to cope with large time gaps between frames arriving
            ALWAYS_VERBOSE;
            cout << "Covariance:\n" << *this << endl;
            COUT(covToSDPixels(*this));
            cout << "Warning: Determinant of covariance outside reasonable range" << endl;
            cout << "This can be caused by a large SUBSAMPLE_FRAMES!" << endl;
            THROW("Determinant of covariance outside reasonable range for pixel location (check SUBSAMPLE_FRAMES value is reasonable)");
        }
    }
}

/ * I'm not quite sure what's happening here--it seems that projected covariances normally have higher uncertainty in the Y direction than X, even though the 3D point had higher uncertainty in X
 * * /
void C2dCovariance::checkXYCompatability() const 
{
    const bool bVerbose = true;

    static int s_nGood=0, s_nBad=0;

    PERIODIC(1000, if(s_nGood<0.5*s_nBad) { COUT(s_nGood); COUT(s_nBad); COUT(s_nGood/(double)s_nBad); / *THROW("Too many covariances are aligned the wrong way (despite portrait images)");* / });
    
    const double dPll = (*this)(0,0);
    const double dPerp = (*this)(1,1);
    
    if(dPll > 0.9*dPerp)
        s_nGood++;
    else
        s_nBad++;
        
    const bool bCheckAlignedWithMotion = false;
    if(bCheckAlignedWithMotion)
    {
        
        CHECK_P(dPll < 0.75*dPerp, *this, "Covariance is greater up-down than side-to-side (this should be unlikely as main uncertainty is how far the tractor has moved, but happens at image top and bottom because of projection)");
    }    
}

void C3dCovariance::checkPosDef() const {
    checkIsPosDef<Base > (*this);
}

void C2dCovariance::checkPosDef() const {
    checkIsPosDef<Base > (*this);
}

CPointLocation::CPointLocation(const C3dWorldPoint & loc, const C3dCovariance & cov, const CTime & time) : location(loc), covariance(cov), time(time) {

}

CPointLocation::~CPointLocation() {
    // TODO Auto-generated destructor stub
    VALGRIND_CHECK_INIT(covariance.sum());
}

//Assimilates a new measurement/covariance.

void CPointLocation::updateLocation(const C3dWorldPoint & loc2, const C3dCovariance & cov2) {
    cout << "Old location2: " << *this;

    const C3dCovariance cov1 = getCovariance();

    //Eigen::Vector3d locVec1, locVec2, locVecUpdated;
    //location.asVector(locVec1);
    //loc2.asVector(locVec2);

    //    cout << "Det 1: " << cov1.determinant() << endl;
    //    cout << "Det 2: " << cov2.determinant() << endl;
    //
    {
        Eigen::Vector3d resid = (loc2 - location);
        double dInnovation = resid.transpose() * (cov1 + cov2).inverse() * resid;
        if (dInnovation > chi2(3, 0.95)) {
            cout << "UPDATE FAILED (actually will try anyway) " << dInnovation << ">" << chi2(3, 0.95) << endl;
            //return;
        }
    }

    Eigen::Matrix3d PRECON;
    PRECON.setZero();
    Eigen::Matrix3d PRECON_INV;
    PRECON_INV.setZero();
    for (int i = 0; i < 3; i++) {
        double dDiagEl = sqrt(cov1(i, i) * cov2(i, i));
        PRECON(i, i) = 1.0 / dDiagEl;
        PRECON_INV(i, i) = dDiagEl;
    }
    checkIsPosDef(PRECON);
    checkIsPosDef(PRECON_INV);

    const Eigen::Matrix3d & info1_precon = ((cov1 * PRECON).inverse());
    const Eigen::Matrix3d & info2_precon = ((cov2 * PRECON).inverse());

    //Will not be pos def, as not symmetric
    //checkIsPosDef<C3dCovariance::Base>(info1_precon);
    //checkIsPosDef<C3dCovariance::Base>(info2_precon);

    //if(IS_DEBUG) CHECK(!zero((info1.inverse() - covariance).squaredNorm()), "Matrix inversion failed (prob ill-conditioned)");

    const Eigen::Matrix3d & covariance_precon = ((info1_precon + info2_precon).inverse());
    //checkIsPosDef<C3dCovariance::Base>(covariance_precon);

    covariance = covariance_precon*PRECON_INV;

    checkIsPosDef<C3dCovariance::Base > (covariance);

    C3dWorldPoint locVecUpdated = (covariance_precon * (info1_precon * location + info2_precon * loc2)).eval(); //This is correct MLE. Must be ill-conditioned????

    location = C3dWorldPoint(locVecUpdated);


    Eigen::Vector3d resid1 = (locVecUpdated - location);
    Eigen::Vector3d resid2 = (locVecUpdated - loc2);

    //Check this really is an MLE
    const Eigen::MatrixXd & info1(cov1.inverse());
    const Eigen::MatrixXd & info2(cov2.inverse());
    double dLL1_best = resid1.transpose() * info1 * resid1;
    double dLL2_best = resid2.transpose() * info2 * resid2;
    double dLL_best = dLL1_best + dLL2_best;

    for (double eps = -0.01; eps <= 0.011; eps += 0.02)
        for (int i = 0; i < 3; i++) {
            resid1(i) += eps;
            resid2(i) += eps;

            double dLL1 = resid1.transpose() * info1 * resid1;
            double dLL2 = resid2.transpose() * info2 * resid2;
            double dLL = dLL1 + dLL2;

            if (dLL <= dLL_best)
                cout << "Point is not an MLE";
            if(IS_DEBUG) CHECK(dLL <= dLL_best, "Point is not an MLE")

            resid1(i) -= eps;
            resid2(i) -= eps;
        }

    cout << "New location: " << *this;

}

void CPointLocation::updateLocation(const CPointLocation & fullloc) {

    //updateLocToNewTime(newCamPos, false); //TODO: Orientation
    if(IS_DEBUG) CHECK(fullloc.getTime() != time, "Updates must be in the same frame (same time)");

    updateLocation(fullloc.location, fullloc.getCovariance());
}*/

/*double C1dMeasurement::getLogLikelihood(double d) const
{
    return -0.5*sqr(d-dMean)/dVariance;
}

void C1dMeasurement::update(C1dMeasurement & diam) {
    double dInf1 = 1.0 / getVar();
    double dInf2 = 1.0 / diam.getVar();
    double dNewVar = 1.0 / (dInf1 + dInf2);
    double dNewDiam = (getMean() * dInf1 + diam.getMean() * dInf2) * dNewVar;
    dMean = dNewDiam;
    dVariance = dNewVar;
    checkOk();
}
    
double CPointLocation::getLogLikelihood(const C3dWorldPoint & p) const
{
    C3dWorldPoint diff = p - location;
    return -0.5*diff.transpose() * covariance.inverse() * diff;
}

optional<const C2dPointLocation> CPointLocation::projectIntoImage_MC(const CWorldCamera & P) const {
    optional<const C2dImagePointPx> ploc2d = P.projectToPx(location, eValidateIO, eSTAny);
    
    if(!ploc2d) //behind camera
        return optional<const C2dPointLocation>();
        
    const C2dImagePointPx loc2d = *ploc2d;

    optional<const double> pdDepth = P.depth(location);
    const double dDepth_z = sys().depth_z(location);

    if(!pdDepth / * behind cam * / || dDepth_z < 0.15 / * too close to plane of cams * /) 
    {
        const bool bVerbose =true;
        COUT("Can't project covariance--point is too close");
        COUT(dDepth_z);
        COUT(location);
        COUT(loc2d);
        return optional<const C2dPointLocation>();
    }
    
    const cv::Size size = sys().getCameraData(P.getFrameId().getCamId()).getImSize();
    
    const int MARGIN = -500; //should be enough to catch any post ends which are out of view
    if(loc2d.x() < MARGIN || loc2d.y() < MARGIN || loc2d.x() > size.width - MARGIN || loc2d.y() > size.height - MARGIN)
    {
        const bool bVerbose =true;
        COUT("Can't project covariance--point is well outside image");
        COUT(location);
        COUT(loc2d);
        COUT(MARGIN);
        return optional<const C2dPointLocation>();
    }
        
    if(getCovariance().determinant() < 1e-16) //essentially 0
        return optional<const C2dPointLocation>(C2dPointLocation(loc2d, C2dCovariance(0.0)));
        
    const bool bVerbose = false;
    
    MAT(outputCov);

    / *const int MARGIN = -5;
    if (!loc2d.inBox(MARGIN, MARGIN, P.getWidth() + MARGIN, P.getHeight() + MARGIN)) {
        cout << "Point outside image--should not happen during initialisation...\n";
        C2dCovariance cov2d;
        cov2d.setIdentity();
        cov2d *= 1000; //Uninformative, outside image anyway so can't track

        return C2dPointLocation(loc2d, cov2d);
    }* /

    // Simulate points, project into image, compute covariance mat
    const int NUM_POINTS = 25;

    const bool MEAN_FROM_MC = true; //Only statistically correct if TRUE
    C2dImagePointPx meanFromSim;

    T2dImPointPxVector aSimPoints2d(NUM_POINTS);

    const double COV_SCALE = 0.01; // Scale down distribution before projecting (appears to be needed still)
    C3dCovariance scaledCov3d(COV_SCALE * getCovariance());
    CMVNormalSampler<C3dCovariance::Base> normalSampler(location, scaledCov3d);

    int nNumPointsTruncated = 0;

    const double MIN_DEPTH = 0.5 * sys().depth_z(location); // TODO: Heuristic

    for (int i = 0; i < NUM_POINTS; i++) {
        C3dWorldPoint simPoint3d;

        bool bTooClose = false;
        do {
            simPoint3d = normalSampler.sample();
            const optional<const double> pdDepth = P.depth(simPoint3d);
            const double dDepth_z = sys().depth_z(simPoint3d);
            bTooClose = !pdDepth || (dDepth_z < MIN_DEPTH);
            if (bTooClose)
                nNumPointsTruncated++;
        } while (bTooClose);

        optional<const C2dImagePointPx> px =  P.projectToPx(simPoint3d, eCheckIO, eSTAny);
        
        if(!px) //will throw...
        {
            bool bVerbose = true;
            COUT(simPoint3d);
            COUT(P.depth(simPoint3d));
            COUT(loc2d);
        }
        aSimPoints2d[i] = *px;//...here.
        
        if(bVerbose) cv::circle(outputCov, aSimPoints2d[i], 2, colour("SampledPoint"));

        //cout << aSimPoints2d[i].x() << '\t' << aSimPoints2d[i].y() << '\t' << simPoint3d.depth(P.cam()) << endl;

        if (MEAN_FROM_MC)
            meanFromSim += aSimPoints2d[i];
    }

    meanFromSim *= 1.0 / NUM_POINTS;
    
    if(bVerbose)
    {
        cv::circle(outputCov, meanFromSim, 3, colour("SampleMean"));
        cv::circle(outputCov, loc2d, 3, colour("PointLocation"));
        SHOW(outputCov);
    }
    
    if (nNumPointsTruncated > NUM_POINTS / 50)
        cout << "Warning: truncating " << nNumPointsTruncated << " of " << NUM_POINTS << " points\n";

    C2dCovariance cov2d = C2dCovariance::Zero();

    for (int i = 0; i < NUM_POINTS; i++) {
        C2dImagePointPx simPoint2d = aSimPoints2d[i];

        if (MEAN_FROM_MC)
            simPoint2d -= meanFromSim;
        else
            simPoint2d -= loc2d;

        cov2d += simPoint2d * simPoint2d.transpose();
    }
    cov2d /= NUM_POINTS; //or NUM_POINTS-1

    if ((meanFromSim - loc2d).norm() > 20) {
        cout << "\n\n3d point: " << *this << endl;
        cout << "Mean: " << meanFromSim << ", projected mean: " << loc2d << endl;
        cout << cov2d << endl;
    }

    const bool CHECK_KURTOSIS = true;
    if (CHECK_KURTOSIS) {//For calc see http://en.wikipedia.org/wiki/Multivariate_normal_distribution#Multivariate_normality_tests
        // alternative formulae: http://www.technion.ac.il/docs/sas/stat/chap19/sect35.htm

        double kurt2d = 0;

        const Eigen::Matrix2d inf(cov2d.inverse());

        for (int i = 0; i < NUM_POINTS; i++) {
            C2dImagePointPx simPoint2d = aSimPoints2d[i];

            if (MEAN_FROM_MC)
                simPoint2d -= meanFromSim;
            else
                simPoint2d -= loc2d;

            const double dev = simPoint2d.transpose() * inf * simPoint2d;
            kurt2d += sqr(dev);
        }

        const double p = (double)cov2d.rows();

        kurt2d = (1.0 / NUM_POINTS) * kurt2d - (p * (p + 2));
        kurt2d *= sqrt(NUM_POINTS / (8 * p * (p + 2)));

        if (fabs(kurt2d) > 4) 
        {
            bool bVerbose = true;
            cout << "Mardia's kurtosis score (should be from N(0,1)): " << kurt2d << endl;
            COUT(meanFromSim);
            COUT(loc2d);
        }
        if(IS_DEBUG) CHECK(fabs(kurt2d) > 100, "Mardia's kurtosis score way too high; distn in image is not normal");
    }
    cov2d /= COV_SCALE;

    return optional<const C2dPointLocation>(C2dPointLocation(loc2d, cov2d));
    //return C2dPointLocation((MEAN_FROM_MC ? meanFromSim : loc2d), cov2d);
}*/

double covToSDMetres(const Eigen::Matrix3d & S)
{
    return pow(S.determinant(), 1.0/6.0); //Approx. sigma in metres
}

double covToSDPixels(const Eigen::Matrix2d & S)
{
    const double dProductOfEigenvectors = S.determinant();
    if(IS_DEBUG) CHECK(dProductOfEigenvectors < 0, "cov matrix not positive definite");
    return pow(dProductOfEigenvectors, 0.25); //Approx. sigma in pixels
}

/**
 * @brief Project a covariance at a 3D point into the image by linearising the homogeneous normalisation
 * @param P
 * @return 
 *
optional<const C2dPointLocation> CPointLocation::projectIntoImage_linearise(const CWorldCamera & P) const
{
    //First check the point is somewhere in the FOV
    const double dTooClose = 0.3;//TODO: too close parameter
    optional<const double> pdDepth_cam = P.depth(location);
    if(!pdDepth_cam || *pdDepth_cam < dTooClose) 
        return optional<const C2dPointLocation>();

    const double dDepth_z = sys().depth_z(location);
    if(dDepth_z < dTooClose) 
        return optional<const C2dPointLocation>();
    
    optional<const C2dImagePointPx> ploc2d = P.projectToPx(location, eValidateIO, eSTAny);

    if(!ploc2d)
        return optional<const C2dPointLocation>();

    const C2dImagePointPx loc2d = *ploc2d;

    const cv::Size size = sys().getCameraData(P.getFrameId().getCamId()).getImSize();
    
    const int MARGIN = -1000; //should be enough to catch any post ends which are out of view, and allow some autotests. Also wires are projected in their entirety into the image and only the sucessfully-projecting part is used
    if(loc2d.x() < MARGIN || loc2d.y() < MARGIN || loc2d.x() > size.width - MARGIN || loc2d.y() > size.height - MARGIN)
    {
        ALWAYS_VERBOSE;
        COUT2("TODO: point projects OOB, but we shouldn't even be trying to compute a covariance for this point? yes for some autotests", loc2d);
        //breakInCpp();
        return optional<const C2dPointLocation>();
    }
    
    / *
     * This C3dPointLocation describes the distribution of a random 3D point X with mean mX, with cov. S
     * 
     * P*X_homog is a 3D point with distribution N(P*mX, P*S*P^T)  (affine transformation of a normal distribution -- http://en.wikipedia.org/wiki/Multivariate_normal_distribution#Affine_transformation)
     * 
     * Write P*mX=(mx,my,mz). We want the distribution of the R.V. project(P, X) = (x,y)/z where P*X=(x,y,z)
     * 
     * Assume (x,y) and z independent, so that:
     * (x,y) has distribution N((mx,my), (P*S*P^T)_block12)
     * z has distribution N(mz, (P*S*P^T)_33)
     * => z has distribution mz*N(1, sz^2 ) where sz^2 = (P*S*P^T)_33 / mz^2
     * 
     * Write A ~ N(1, sz^2 ) so that z ~ mz*A
     * Linearising (assuming sz is small)
     * 
     * 1/A ~ N(1, sz^2 )
     * 
     * 1/z = 1/(mz*A) ~ N(1/mz, (sz/mz)^2)
     * 
     * Now want product of 2 rv's, (x,y) * 1/z. Linearise again:
     * 
     * N((mx,my), (P*S*P^T)_block12) * N(1/mz, (sz/mz)^2)
     * 
     * = 1/mz * N((mx,my), (P*S*P^T)_block12) + (mx,my)*N(0, (sz/mz)^2) + H.O.T.         (H.O.T.=product of 2 zero-mean normal distns)
     * 
     * ~ N((mx,my)/mz, (P*S*P^T)_block12 / mz^2   +   diag((mx*sz/mz)^2, (my*sz/mz)^2))
     * 
     * 
     * *********************
     * Problem with this is that the linearisation introduces too significant an error. 
     * 
     * Try again with P0 = K^-1 P (turn on bUncalibrateP)
     * -- solution is a distribution in normalised image coordinates.
     * 
     * /
     
    const bool bUncalibrateP = true;
    
    const bool bVerbose = false; 
    //PERIODIC(10000, bVerbose = true);

    optional<const C2dImagePointPx> mean2d = P.projectToPx(location, eValidateIO, eSTAny);
    if(!mean2d)
        return optional<const C2dPointLocation>();
    
    COUT(P);
    
    const Eigen::Matrix3d & K_inv = sys().getCameraData(P.getFrameId().getCamId()).getCalibration().calibrationMat().mat_inv();
    TCamMat Pmat;

    if(bUncalibrateP)
    {
        Pmat = K_inv * P.getPMat();
        COUT(Pmat);
    }
    else 
    {
        Pmat = P.getPMat();
    }
        
    const Eigen::Vector3d mean_xyz = Pmat * location.homog();
    const Eigen::Vector2d mean_xy = mean_xyz.head<2>();
    
    COUT(mean_xy);
    
    const double mean_z = mean_xyz.z();
    
    if(mean_z<=0)
        return optional<const C2dPointLocation>();

    C2dImagePoint mean_im;
    if(bUncalibrateP)
        mean_im=mean_xy/mean_z; //homogeneous normalisation

    Eigen::Matrix4d covariance_ext = Eigen::Matrix4d::Zero();
    covariance_ext.block<3,3>(0,0) = covariance;
    const Eigen::Matrix3d cov = Pmat * covariance_ext * Pmat.transpose();
    //if(IS_DEBUG) checkIsPosDef(cov);// This is critical loop so don't want too many pos-def checks
    
    COUT(cov);
    
    const Eigen::Matrix2d cov_xy = cov.block<2,2>(0,0);
    
    COUT(cov_xy);
    
    //if(IS_DEBUG) checkIsPosDef(cov_xy);
    
    const double sz_sq = cov(2,2)/sqr(mean_z);
    
    const Eigen::Matrix2d cov_xy_normalised = cov_xy/sqr(mean_z);
    / *if(P.getFrameId().getCamId() == CImagingSystemData::eCamTop)* /
    if(!bUncalibrateP) {
        const C2dCovariance covForTesting = cov_xy_normalised; 
        covForTesting.checkXYCompatability();
    }
    Eigen::Matrix2d cov_xy_fromZDistn = Eigen::Matrix2d::Identity() * sz_sq / sqr(mean_z);
    
    //if(IS_DEBUG) checkIsPosDef(cov_xy_fromZDistn);
    cov_xy_fromZDistn(0,0) *= sqr(mean_xy.x());
    cov_xy_fromZDistn(1,1) *= sqr(mean_xy.y());
    
    Eigen::Matrix2d cov_xy_total = cov_xy_normalised + cov_xy_fromZDistn;
    
    if(IS_DEBUG) checkIsPosDef<Eigen::Matrix2d>(cov_xy_total);
    
    if(bUncalibrateP)
    {
        / *
        x_im ~ N(mean_im, cov_xy_total)
        
        //Actually want distribution for pixel coordinates
        x_px = K*x_im.homog()
        
        x_px ~ N(K*mean_im.homog(), K*cov_xy_total_extended*K.transpose())
         
         K*mean_im.homog().z==1 so we can just clip: K*cov_xy_total_extended*K.transpose() -> sqr(K.focalLength())*cov_xy_total
        * /
        const Eigen::Matrix2d Ksub = sys().getCameraData(P.getFrameId().getCamId()).getCalibration().calibrationMat().mat().block<2,2>(0,0);
        
        cov_xy_total = (Ksub*cov_xy_total*Ksub.transpose()).eval();
    }
    
    //bVerbose |= (covToSDPixels(cov_xy_total) > 100);
    if (bVerbose)
    {
        COUT("#######");
        COUT(location);
        COUT(covariance);
        COUT(P.translationToCam());
        //COUT(sz);
        COUT(cov_xy_normalised);
        COUT(cov_xy_fromZDistn);
        COUT(cov_xy_total);
        
        COUT(covToSDPixels(cov_xy_normalised));
        COUT(covToSDPixels(cov_xy_fromZDistn));
        COUT(covToSDPixels(cov_xy_total));

        COUT(covToSDMetres(covariance));
        
        MAT2(projCovariance, P.getFrameId());
        optional<const C2dPointLocation> locMC = projectIntoImage_MC(P);
        if(locMC)
        {
            COUT(locMC->getCovariance());
            COUT(covToSDPixels(locMC->getCovariance()));
            if(covToSDPixels(cov_xy_total) > 0) if(IS_DEBUG) CHECK(fabs(log(covToSDPixels(locMC->getCovariance())) - log(covToSDPixels(cov_xy_total))) > 2, "Linearised covariance projection doens't agree with MC covariance projection (probably MC is inaccurate)");
            
            drawEllipse(projCovariance, locMC->getCovariance(), locMC->getLocation(), 3, colour("CovEllipse_3SigmaMC"), 2);
        }

        drawEllipse(projCovariance, cov_xy_total, *mean2d, 3, colour("CovEllipse_3Sigma"), 2);
        SHOW(projCovariance);
        
        MAT2(projCovarianceCompXY, P.getFrameId());
        drawEllipse(projCovarianceCompXY, cov_xy_normalised, *mean2d, 3, colour("CovEllipse_3Sigma"), 2);
        SHOW(projCovarianceCompXY);
        
        MAT2(projCovarianceCompZ, P.getFrameId());
        drawEllipse(projCovarianceCompZ, cov_xy_fromZDistn, *mean2d, 3, colour("CovEllipse_3Sigma"), 2);
        SHOW(projCovarianceCompZ);
    }
    
    //if(fabs(location.y()) < 0.1 && P.getFrameId().getCamId() == CImagingSystemData::eCamTop)
    C2dCovariance(cov_xy_total).checkXYCompatability();
    
    return optional<const C2dPointLocation>(C2dPointLocation(*mean2d, cov_xy_total));
}

std::ostream & operator<<(std::ostream & out, const C1dMeasurement & diam) {
    return out << "Measurement mean " << diam.getMean() << ", var " << diam.getVar() << ", sd " << sqrt(diam.getVar()) << endl;
}

std::ostream & operator<<(std::ostream & out, const CPointLocation & diam) {
    Eigen::Array3d var = diam.getCovariance().diagonal();
    return out << "3D location mean " << diam.getLocation() << ", time: " << diam.getTime() << ", covariance: " << var.sqrt().transpose() << endl;
}

void CPointLocation::checkCovariance() const {
    if(IS_DEBUG) CHECK(covariance.determinant() < 0, "Uninit/bad covariance");
    checkIsPosDef<C3dCovariance::Base > (covariance);
}

const C2dCovariance makeCov2d(const Eigen::Vector2d & pll, const double dVarPll, const double dVarPerp) {
    if(IS_DEBUG) CHECK(!zero(pll.squaredNorm() - 1), "pll not normalised");
    const Eigen::Vector2d perp(pll(1), -pll(0));
    Eigen::Matrix2d Q;
    Q << pll, perp;

    if(IS_DEBUG) CHECK(std::isnan(Q.sum()), "nan making Q-matrix")

    Eigen::Matrix2d LAMBDA;
    LAMBDA.setZero();
    LAMBDA(0, 0) = dVarPll;
    LAMBDA(1, 1) = dVarPerp;

    if(IS_DEBUG) CHECK(std::isnan(LAMBDA.sum()), "nan making LAMBDA-matrix");

    //Build covariance mat
    const C2dCovariance C(Q * LAMBDA * Q.transpose());

    checkIsPosDef<C2dCovariance::Base > (C);

    return C;
}*/

double chi2(const int nDOFs, const double p) {
    boost::math::chi_squared chi2dist(nDOFs);
    double quantile = boost::math::quantile(chi2dist, p);
    return quantile;
}

/*CCam::~CCam()
{
}
    
CCam::CCam(const CImagingSystemData::eCams camId, const cv::Size imSize, const double dFocalLength, const TEigen3dRotation & R, const C3dWorldPoint & t, const double dActiveMargin)
: camId(camId), P(R, t, CNormalisedCameraMat::eRTFromCalibration), newK(imSize, dFocalLength), width(imSize.width), height(imSize.height), dActiveMargin(dActiveMargin), KP(newK.mat(), P)
{
}
 
CCam::CCam(const CImagingSystemData::eCams camId, const std::string strCalibFileName, const TEigen3dRotation & R, const C3dWorldPoint & t, const double dActiveMargin)
: camId(camId), P(R, t, CNormalisedCameraMat::eRTFromCalibration), newK(strCalibFileName), width(-1), height(-1), dActiveMargin(dActiveMargin), KP(newK.mat(), P) {

}

void CCam::pp() const
{
    ALWAYS_VERBOSE;
    COUT(camId);
    COUT(KP.getPMat());
}

C2dImagePointPx CCam::imCentre() const { 
    return C2dImagePointPx(height * 0.5, width * 0.5);
}

bool CCam::inBox(const C2dImagePointPx & imLoc, const double margin) const {
    return imLoc.x() >= margin && imLoc.y() >= margin && imLoc.x() < getWidth() - margin && imLoc.y() < getHeight() - margin;
    //return imLoc.inBox(margin, margin, getWidth() - margin, getHeight() - margin);
}

const C3dWorldPoint CCam::pixelToLocalCoords(const C2dImagePointPx & x, const double depth) const {
    return KP.pxToWorld_depth(x, depth);
}

/ * 
 * P maps local coordinates to camera coordinates
 * P*localPoint = camPoint
 * 
 * R*localPoint + t = camPoint
 * 
 * localPoint = R^T * (camPoint-t)
 * 
 * * /
const C3dWorldPoint CCam::camToLocalCoords(const C3dWorldPoint & camPoint) const
{
    const C3dWorldPoint localPoint = P.inverseTransform(camPoint);
    //const C3dWorldPoint localPoint = P.transform(camPoint);
    return localPoint;
}
*/
C2dImagePointPx CNewCamCalibMatrix::toPxCoords(const C2dImagePoint & imPoint) const 
{
    //TODO: RD correction (also in CPixelCameraMat)
    return C2dImagePointPx((K*imPoint.homog()).eval());
}

C2dImagePoint CNewCamCalibMatrix::toImCoords(const C2dImagePointPx & imPointPx) const 
{
    //TODO: RD uncorrection (also in CPixelCameraMat)
    return C2dImagePoint((K_inv*imPointPx.homog()).eval());
}

CNewCamCalibMatrix::CNewCamCalibMatrix() : K(Eigen::Matrix3d::Zero()), K_inv(Eigen::Matrix3d::Zero())
{
    
}

CNewCamCalibMatrix::CNewCamCalibMatrix(const std::string strCalibFileName) : K(Eigen::Matrix3d::Zero()), K_inv(Eigen::Matrix3d::Zero())
{
    if (boost::filesystem::exists(strCalibFileName)) {
        std::ifstream calibFile(strCalibFileName.c_str());
        if(IS_DEBUG) CHECK(!calibFile.is_open() || calibFile.eof(), "Failed to open intrinsic calibration file");
        
        for(int r=0;r<3;r++)
        {
            for(int c=0;c<3;c++)
            {
                if(IS_DEBUG) CHECK(calibFile.eof(), "EOF reading calibration");
                calibFile >> K(r,c);
            }
        }
        K_inv = K.inverse();
        
        int nRDCoeffs=-1;
        calibFile >> nRDCoeffs;
        aRDCoeffs = Eigen::VectorXd(nRDCoeffs);
        
        for(int nRDCoeff=0;nRDCoeff<nRDCoeffs;nRDCoeff++)
        {
            if(IS_DEBUG) CHECK(calibFile.eof(), "EOF reading calibration");
            calibFile >> aRDCoeffs(nRDCoeff);
            
            THROW("RD correction not completely implemented (NB we multiply this K into camera matrices so need to correct when we use these too)");
        }
        
    } else
        cout << "No intrinsic calibration at " << strCalibFileName << "\n***Calibration only***\n";

        
    if(IS_DEBUG) CHECK(std::isnan(K_inv.sum()), "K_inv is nan (K is singular)");
}

/*For simulating calibration matrices*/
CNewCamCalibMatrix::CNewCamCalibMatrix(const cv::Size & imSize, const double dFocalLength) 
{
    K = Eigen::Matrix3d::Identity();
    K(0,0) = K(1,1) = dFocalLength;

    K(0,2) = imSize.width*0.5;
    K(1,2) = imSize.height*0.5;

    K_inv = K.inverse();
}

void CWorldCamera::setCamMatForAnimation(const TCamMat & Pnew)
{
    P = Pnew;

    checkInit();
}

/*double CWorldCamera::depth_z(const C3dWorldPoint & X) const
{
    return getCamPose().depth_z(X); //A bit of a convoluted reference back to the parent of this camera...
}

double CWorldCamera::relativeY(const C3dWorldPoint & X) const
{
    return getCamPose().relativeY(X); //A bit of a convoluted reference back to the parent of this camera...
}

const CCameraPose & CWorldCamera::getCamPose() const
{
    return *scene().getCameraPose(frameId.getTime());
}*/

/*double C2dPointLocation::sigmaFromResidual(const C2dImagePointPx & diff) const
{
    const double dSigmas_sq = diff.transpose() * covariance.inverse() * diff;
    return sqrt(dSigmas_sq);
}

double C2dPointLocation::sigma(const C2dImagePointPx & measuredLocation) const
{
    TEigen2dPoint diff = measuredLocation-location;
    return sigmaFromResidual(diff);
}

void C2dPointLocation::drawCovEllipse(cv::Mat & im) const
{
    CHECK_MAT_INIT(im);
    drawEllipse(im, covariance, location, 3, colour("CovEllipse_3Sigma"));
}

//Area of 1-sigma covariance ellipse, in pixels
double C2dPointLocation::covAreaPx() const 
{
    return sqrt(covariance.determinant()); 
}*/


void CPixelCameraMat::setPxOffset(const int nRows, const int nCols)
{
    const bool bVerbose = false;
    //Pre-multiply by a matrix which will apply the correct offset
    Eigen::Matrix3d Koffset = TEigen3dRotation::Identity();
    Koffset(0,2) = -nCols;
    Koffset(1,2) = -nRows;
    COUT(P);
    P = Koffset * P;
    COUT(P);
}

std::ostream& operator<<(std::ostream& s, const CCameraMat_base & P) { return s << P.P; }

C2dImagePoint::C2dImagePoint(const double x, const double y) : TEigen2dPoint(x,y) 
{
    if(squaredNorm() > sqr(4))
    {
        REPEAT(100, cout << "Warning: these look more like pixel coordinates than image coordinates: C2dImagePoint=" << transpose() << endl);
    }
}

C2dImagePointPx::operator cv::Point () const { 
    
    if(squaredNorm() < sqr(100000))
    {
        cv::Point p(cvRound((*this)(0)), cvRound((*this)(1)));  
        if(IS_DEBUG) CHECK(fabs(p.x - (*this)(0)) > 1 || fabs(p.y - (*this)(1)) > 1, "operator cv::Point () failed");
        if(fabs((double)p.x) > 30000 || fabs((double)p.y) > 30000)
        {
            REPEAT(10, 
            cout << "Very large pixel coordinate: " << p << " ";
            cout << this->transpose() << endl;);
        }
        return p;
    }
    else
    {
        REPEAT(10, cout << "Converting massive px coordinates to cv::Point failed--returning cv::Point(-10000,-10000)" << endl);
        return cv::Point(-10000,-10000);
    }        
} 

