/* 
 * Testing new matching code. See main.cpp for an example of how to use my code.
 */

#include <cstdlib>
#include "imageSource/imageSourceFromDir.h"
#include <boost/smart_ptr.hpp>
#include "featureExtract/featureExtractor.h"
#include "ransac/ransacParams.h"
#include "ransac/ransac.h"
#include "description/vectorDescriptor.h"
#include "geom/geom.h"
#include "ransac/refineEOnRTManifold.h"
#include "featureMatching/guidedFeatureMatching.h"
#include <fstream> 
#include <util/stats.h>
#include <Eigen/Dense>

using namespace std;

void loadRelPoses(vector<C3dPoint> & a3dPoints, vector<C3dRotation> & a3dRotation);

/*
 * 
 */
int main_testMatching() {
    vector<C3dPoint> a3dPoints; 
    vector<C3dRotation> a3dRotation;
    loadRelPoses(a3dPoints, a3dRotation);
    ofstream results("resultsRT");
    CRTErrorStats /* statsRANSAC(results, "RANSAC"),*/ statsRefine(results, "Refine");

    //Setup image loader from directory
    CImParams IM_PARAMS(0,0);
    IM_PARAMS.ImageDir.IMAGE_DIR = "/home/data/data/data/SLAMDUNK/bmvc_mono/images"; 
    boost::scoped_ptr<CImageSource> pImageSource(CImageSourceFromDir::makeImageSource(IM_PARAMS));
    int nImageId = 0;

    IplImage * pLastFrame = pImageSource->createImage();
    IplImage * pFrame = pImageSource->createImage();

    //Make a feature extractor
    CCornerParams CORNERPARAMS(0,0);
    CPatchDescriptorParams PATCHDESCRIPTORPARAMS(0,0);
    CDescriptorSetClusteringParams DSCPARAMS(0,0);

    boost::scoped_ptr<CFeatureExtractor> pFeatureExtractor ( CFeatureExtractor::makeFeatureExtractor(IM_PARAMS, CORNERPARAMS, PATCHDESCRIPTORPARAMS, DSCPARAMS) );
    
    //Camera calibration settings
    T2dPoints aCalibratedPoints1, aCalibratedPoints2;
    CInlierProbs adArrLikelihood;
    CPointIdentifiers pointIds;

    //Settings for feature matching
    const int NN = 4; // N-M feature matches
    CMatchableDescriptors::CMatchSettings MS(0.8, 0.6, NN);
    
    //Settings for RANSAC
    CRANSACParams RANSAC_PARAMS(0,0);
    RANSAC_PARAMS.HypothesiseAlg = CRANSACParams::e5Pt_GradientDesc;
    RANSAC_PARAMS.TOPDOWN_ITERS = 0; //This turns off the topdown solutions refinement, which doesn't work very well.
    
    CDescriptorSet * pLastDescriptors = 0;
    const int nSkip = 10;
    while (pImageSource->loadImage(nImageId, pFrame)) {
        CDescriptorSet * pDescriptors = pFeatureExtractor->getDescriptors(pFrame);

        if (pLastDescriptors) //we're past the first frame
        {
            //Find matches between features
            //boost::scoped_ptr<const CBoWCorrespondences> pCorr ( pLastDescriptors->getBruteForceCorrespondenceSet(pDescriptors,MS,0) );
            C3dPoint T_GT = a3dPoints[nImageId] - a3dPoints[nImageId-nSkip];
            T_GT.rotate(a3dRotation[nImageId-nSkip].t());
            T_GT.normalise();
            C3dRotation R_GT =  a3dRotation[nImageId] * a3dRotation[nImageId-nSkip].t();
            cout << "GT orientation: " << R_GT << endl;
            cout << "GT direction: " << T_GT << endl;
            
            Eigen::Vector3d phiThetaPsi(0,0,0), t(CRandom::Uniform(-1.5,1.5)*0,0,1); // forward moving with no relative orientation

            CGuidedFeatureMatching::TCovMat66 fullCov;
            fullCov.setIdentity();
            fullCov *= 0.1;
            
            boost::scoped_ptr<const CBoWCorrespondences> pCorr ( CGuidedFeatureMatching::guidedFeatureMatch(pDescriptors, pLastDescriptors, IM_PARAMS.getCamCalibrationMat(), MS, t, phiThetaPsi, fullCov) );

            //Convert points to normalised image coordinates:
            pCorr->calibrate(IM_PARAMS.getCamCalibrationMat(), aCalibratedPoints1, aCalibratedPoints2, adArrLikelihood, pointIds, true);
            
            //Do RANSAC
            C3x3MatModel E; //Essential matrix
            CMask inlierMask(aCalibratedPoints1.size());
            int nInliers = getE(aCalibratedPoints1, aCalibratedPoints2, adArrLikelihood, pointIds, RANSAC_PARAMS, E, inlierMask, -1, IM_PARAMS.getCamCalibrationMat().focalLength(), 1);
            
            cout << nInliers << "/" << inlierMask.size() << " inliers after RANSAC" << endl;

            C3dPoint T;
            C3dRotation R;
            bool bSuccess = RTFromE(inlierMask, E, aCalibratedPoints1, aCalibratedPoints2, pointIds, R, T, RANSAC_PARAMS.E_INLIER_THRESH_PX / IM_PARAMS.getCamCalibrationMat().focalLength());
            if(!bSuccess)
            {
                cout << "Pure rotation detected, setting T=(0,0,1)" << endl;
                T=C3dPoint(0,0,1);
            }
            
            CRefineEOnRTManifold::refineRobustOnMask(aCalibratedPoints1, aCalibratedPoints2, R, T, inlierMask);
            
            cout << "Relative orientation: " << R << endl;
            cout << "Translation direction: " << T << endl; 

            double dErrR = diff(R_GT, R);
            double dErrT = angle(T_GT, T);
            cout << dErrR << " = error in R, " << dErrT << "=error in T" << endl;
            statsRefine.addRT(dErrR, dErrT);
        }

        delete pLastDescriptors; pLastDescriptors=0;
        std::swap(pDescriptors, pLastDescriptors);
        std::swap(pLastFrame, pFrame);
        nImageId += nSkip; //Skip a few frames, otherwise have essentially pure rotation
    }
    return 0;
}

// Simulated descriptor (just an integer, with abs distances between descriptors). Used for simulated data.
class CSymDescriptor : public CDescriptor
{
	int id_;
	CLocation loc_;
public:
	CSymDescriptor(CLocation loc, int id) : id_(id), loc_(loc) {}
	virtual TDist distance(const CDescriptor * pd) const { return abs(id_ - CAST<const CSymDescriptor *>(pd)->id_) ; }
	virtual CLocation location() const { return loc_; }
	virtual int size() const { return sizeof(CSymDescriptor); }
	virtual int length() const { return 1; };
};

bool OOB(const CLocation p, const int nWidth, const int nHeight)
{
    return p.x() < 0 || p.x() >= nWidth || p.y() < 0 || p.y() >= nHeight;
}

int main_simData(int argc, char** argv) {

    const bool bRefineH = true; //Estimate a homography instead
    
    ofstream results("resultsSim");
    CRTErrorStats /* statsRANSAC(results, "RANSAC"),*/ statsRefine(results, "Refine");

    //Setup image loader from directory (just for calibration)
    CImParams IM_PARAMS(0,0);
    IM_PARAMS.ImageDir.IMAGE_DIR = "/home/data/data/data/SLAMDUNK/bmvc_mono/images"; 
    boost::scoped_ptr<CImageSource> pImageSource(CImageSourceFromDir::makeImageSource(IM_PARAMS));
    int nImageId = 0;

    //Settings for feature matching
    const int NN = 4; // N-M feature matches
    CMatchableDescriptors::CMatchSettings MS(0.8, 0.6, NN);
    CDescriptorSetClusteringParams DSCPARAMS(0,0);
    
    //Settings for RANSAC
    CRANSACParams RANSAC_PARAMS(0,0);
    RANSAC_PARAMS.HypothesiseAlg = CRANSACParams::e5Pt_GradientDesc;
    RANSAC_PARAMS.TOPDOWN_ITERS = 0; //This turns off the topdown solutions refinement, which doesn't work very well.
    RANSAC_PARAMS.VERBOSE = true;
    RANSAC_PARAMS.E_INLIER_THRESH_PX = 2;
    
    CRANSACHomographyParams RANSAC_HOMOG_PARAMS(0,0);
    RANSAC_HOMOG_PARAMS.VERBOSE = true;

    const int NUM_ITERS = 1000;
    for(int nIter=0; nIter<NUM_ITERS; nIter++)
    {
        T2dPoints aCalibratedPoints1, aCalibratedPoints2;
        CInlierProbs adArrLikelihood;
        CPointIdentifiers pointIds;
        
        //Make 2 cameras:
        C3dPoint t_gt;
        t_gt.setRandomNormal();
        //t_gt = C3dPoint(t_gt.getX(), t_gt.getY(), 0);
        t_gt.normalise();

        //if(bRefineH)
            //t_gt *= 0;
        
        Eigen::Vector3d eulerAngles_gt(0,0,0), n_gt(0,0,-1);
        eulerAngles_gt.setRandom();
        eulerAngles_gt *= 0.1;
        
        C3dRotation R_gt = phiThetaPsiToQuat(eulerAngles_gt);
        
        CCamera P, Pp = R_gt | t_gt;
        
        const int NUM_POINTS=50;
        
        //boost::scoped_ptr<CBoWCorrespondences> pCorr = new CBoWCorrespondences;
        CMetricSpaceDescriptorSet aDescriptorsLeft(DSCPARAMS, NUM_POINTS), aDescriptorsRight(DSCPARAMS, NUM_POINTS);
        
        std::set<CLocation> aPointsLeft, aPointsRight;
        const double PLANE_DEPTH=4;
        
        const double EFFECTIVE_DEPTH = PLANE_DEPTH;

        Eigen::Vector3d t_Homography = t_gt.asVector() / EFFECTIVE_DEPTH;
        
        Eigen::Matrix3d H_gt=makeH(t_Homography, eulerAngles_gt, n_gt);
        
        for(int nId=0; aDescriptorsLeft.Count() < NUM_POINTS; nId++)
        {
            C3dPoint p;
            if(bRefineH)
                p.setRandomPlanar(1);
            else    
                p.setRandomNormal();
            
            p*=PLANE_DEPTH;
            
            //if(p.testInFront(P, Pp))
            const double dMinDepth = 0.25;

            if(p.depth(P) > dMinDepth && p.depth(Pp) > dMinDepth)
            {
                C2dPoint x=p.photo(P);
                C2dPoint xp=p.photo(Pp);
                C2dPoint x_uncalib=x;
                C2dPoint xp_uncalib=xp;

                //Simulate 2 descriptors:
                x_uncalib.uncalibrate(IM_PARAMS.getCamCalibrationMat());
                xp_uncalib.uncalibrate(IM_PARAMS.getCamCalibrationMat());

                CLocation locx(x_uncalib.getX(), x_uncalib.getY());
                CLocation locxp(xp_uncalib.getX(), xp_uncalib.getY());
                
                if(!OOB(locx, 640, 480) && !OOB(locxp, 640, 480)) 
                {
                    if(aPointsLeft.find(locx) == aPointsLeft.end() && aPointsRight.find(locxp) == aPointsRight.end() )
                    {
                        aPointsLeft.insert(locx);
                        aPointsRight.insert(locxp);

                        aCalibratedPoints1.push_back(x);
                        aCalibratedPoints2.push_back(xp);
                        adArrLikelihood.push_back(0.5);
                        pointIds.push_back(CPointIds(nId, nId));

                        //pCorr->push_back();

                        CDescriptor * pDescLeft = new CSymDescriptor(locx, nId);
                        CDescriptor * pDescRight = new CSymDescriptor(locxp, nId);
                        aDescriptorsLeft.Push(pDescLeft);
                        aDescriptorsRight.Push(pDescRight);
                    }
                }

                if(bRefineH)
                {
                    Eigen::Vector3d vx; x.asVector(vx);
                    Eigen::Vector3d vxp; xp.asVector(vxp);

                    REPEAT(2,
                            cout << "H: " << H_gt << endl;
                    cout << "vx:" << vx.transpose() << endl;
                    cout << "vxp:" << vxp.transpose() << endl;
                    cout << "H vx:" << homog(H_gt*vx).transpose() << endl;
                    cout << "H vxp:" << homog(H_gt*vxp).transpose() << endl;
                    cout << "H-1 vx:" << homog(H_gt.inverse()*vx).transpose() << endl;
                    cout << "H-1 vxp:" << homog(H_gt.inverse()*vxp).transpose() << endl);

                    CHECK((homog(H_gt*vx) - homog(vxp)).squaredNorm() > 0.001, "Error generating data");
                }
            }
        }        
        
    
        Eigen::Vector3d phiThetaPsi = eulerAngles_gt, t = t_gt.asVector(); 

        CGuidedFeatureMatching::TCovMat66 fullCov;
        fullCov.setIdentity();
        fullCov *= 0.01;
        
        CGuidedFeatureMatching::TCovMat99 fullCov_H;
        fullCov_H.setIdentity();
        fullCov_H *= 0.01;
        
        boost::scoped_ptr<const CBoWCorrespondences> pCorr;

        if(bRefineH)
        {
            Eigen::Vector3d t_Homography_perturbed = t_Homography, phiThetaPsi_perturbed = phiThetaPsi;
            t_Homography_perturbed += Eigen::Vector3d::Random() * 0.01;
            phiThetaPsi_perturbed += Eigen::Vector3d::Random() * 0.01;
            
            pCorr.reset( CGuidedFeatureMatching::guidedFeatureMatch_Homography(&aDescriptorsLeft, &aDescriptorsRight, IM_PARAMS.getCamCalibrationMat(), MS, /*((nIter % 2 == 0) ? 1 : -1) */ t_Homography_perturbed, phiThetaPsi_perturbed, n_gt, fullCov_H) );
        }
        else
            pCorr.reset( CGuidedFeatureMatching::guidedFeatureMatch(&aDescriptorsLeft, &aDescriptorsRight, IM_PARAMS.getCamCalibrationMat(), MS, /*((nIter % 2 == 0) ? 1 : -1) */ t, phiThetaPsi, fullCov) );
            
        //boost::scoped_ptr<const CBoWCorrespondences> pCorr ( aDescriptorsLeft.getBruteForceCorrespondenceSet(&aDescriptorsRight,MS,0) );

        //Convert points to normalised image coordinates:
        for(int i=0;i<5;i++)
            cout << aCalibratedPoints1[i] << endl;
        
        pCorr->calibrate(IM_PARAMS.getCamCalibrationMat(), aCalibratedPoints1, aCalibratedPoints2, adArrLikelihood, pointIds, true);

        if(pCorr->size() >= 5)
            for(int i=0;i<5;i++)
                cout << aCalibratedPoints1[i] << endl;

        //Do RANSAC
        C3x3MatModel E; //Essential matrix
        CMask inlierMask(aCalibratedPoints1.size());
        int nInliers = 0;

        C3dPoint T;
        C3dRotation R;
        
        if(bRefineH)
        {
            nInliers = getH(aCalibratedPoints1, aCalibratedPoints2, adArrLikelihood, pointIds, RANSAC_HOMOG_PARAMS, E, inlierMask, IM_PARAMS.getCamCalibrationMat().focalLength(), 1);
            C3dPoint n;
            Eigen::Matrix3d H(E.asDouble9());
            H.transposeInPlace();
            decomposeHomography(H,R,n,T);
            cout << "Plane normal: " << n << endl;
            //T.normalise();
        } else {
            nInliers = getE(aCalibratedPoints1, aCalibratedPoints2, adArrLikelihood, pointIds, RANSAC_PARAMS, E, inlierMask, -1, IM_PARAMS.getCamCalibrationMat().focalLength(), 1);
            bool bSuccess = (nInliers>0) ? RTFromE(inlierMask, E, aCalibratedPoints1, aCalibratedPoints2, pointIds, R, T, RANSAC_PARAMS.E_INLIER_THRESH_PX / IM_PARAMS.getCamCalibrationMat().focalLength()) : false;
            if(!bSuccess)
            {
                cout << "Pure rotation detected, setting T=(0,0,1)" << endl;
                T=C3dPoint(0,0,1);
            }

            if(nInliers > 5)
                CRefineEOnRTManifold::refineRobustOnMask(aCalibratedPoints1, aCalibratedPoints2, R, T, inlierMask);
        }

        cout << nInliers << "/" << inlierMask.size() << " inliers after RANSAC" << endl;

        cout << "Relative orientation: " << R << endl;
        cout << "Translation direction: " << T << endl; 

        double dErrR = diff(R_gt, R);
        double dErrT = bRefineH ? (C3dPoint(t_Homography)-T).length() : angle(t_gt, T);
        cout << dErrR << " = error in R, " << dErrT << "=error in T" << endl;
        statsRefine.addRT(dErrR, dErrT);
    }

    return 0;
}

int main(int argc, char** argv) {
    return main_simData(argc,argv);
    return main_testMatching();
}
