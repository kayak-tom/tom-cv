/* 
 * Evaluate RANSAC on data with GT
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
#include <fstream>
#include <util/stats.h>

using namespace std;

void loadRelPoses(vector<C3dPoint> & a3dPoints, vector<C3dRotation> & a3dRotation)
{
    std::ifstream slamdunkGT("/home/data/data/data/SLAMDUNK/repository/trajectories/circle_handmade.xml");
    
    while(!slamdunkGT.eof())
    {
        double dTime,qx,qy,qz,qw,x,y,z;
        slamdunkGT >> dTime;
        slamdunkGT >> qw;
        slamdunkGT >> qx;
        slamdunkGT >> qy;
        slamdunkGT >> qz;
        slamdunkGT >> x;
        slamdunkGT >> y;
        slamdunkGT >> z;
        
                
        a3dPoints.push_back(C3dPoint(x,y,z));
        a3dRotation.push_back(C3dRotation(qx,qy,qz,qw));
    }
}

/*
 * 
 */
int main2(int argc, char** argv) {

    vector<C3dPoint> a3dPoints; 
    vector<C3dRotation> a3dRotation;
    loadRelPoses(a3dPoints, a3dRotation);
    
    //Setup image loader from directory
    CImParams IM_PARAMS(0,0);
    IM_PARAMS.ImageDir.IMAGE_DIR = "/home/data/data/data/SLAMDUNK/bmvc_mono/images"; 
    boost::scoped_ptr<CImageSource> pImageSource(CImageSourceFromDir::makeImageSource(IM_PARAMS));

    IplImage * pLastFrame = pImageSource->createImage();
    IplImage * pFrame = pImageSource->createImage();

    //Make a feature extractor
    CCornerParams CORNERPARAMS(0,0);
    CPatchDescriptorParams PATCHDESCRIPTORPARAMS(0,0);
    CDescriptorSetClusteringParams DSCPARAMS(0,0);

    //CORNERPARAMS.MAX_FEATURES = 600;
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
    
    int nSkip = 20;
    
    double dThresh = RANSAC_PARAMS.E_INLIER_THRESH_PX / IM_PARAMS.getCamCalibrationMat().focalLength();
    
    ofstream results("resultsRT");
    CRTErrorStats statsRANSAC(results, "RANSAC"), statsRefine(results, "Refine");
    
    for(int nOffset=0; nOffset < nSkip; nOffset++)
    {
        int nImageId = nOffset;
        CDescriptorSet * pLastDescriptors = 0;
        while (pImageSource->loadImage(nImageId, pFrame)) {
            CDescriptorSet * pDescriptors = pFeatureExtractor->getDescriptors(pFrame);

            if (pLastDescriptors) //we're past the first frame
            {
                //Find matches between features
                boost::scoped_ptr<const CBoWCorrespondences> pCorr ( pDescriptors->getBruteForceCorrespondenceSet(pLastDescriptors,MS,0) );

                //Convert points to normalised image coordinates:
                pCorr->calibrate(IM_PARAMS.getCamCalibrationMat(), aCalibratedPoints1, aCalibratedPoints2, adArrLikelihood, pointIds, true);

                //Do RANSAC
                C3x3MatModel E; //Essential matrix
                CMask inlierMask(aCalibratedPoints1.size());
                int nInliers = getE(aCalibratedPoints1, aCalibratedPoints2, adArrLikelihood, pointIds, RANSAC_PARAMS, E, inlierMask, -1, IM_PARAMS.getCamCalibrationMat().focalLength(), 1);

                cout << nInliers << "/" << inlierMask.size() << " inliers after RANSAC" << endl;

                C3dPoint T_GT = a3dPoints[nImageId] - a3dPoints[nImageId-nSkip];
                T_GT.rotate(a3dRotation[nImageId-nSkip].t());
                T_GT.normalise();
                C3dRotation R_GT =  a3dRotation[nImageId] * a3dRotation[nImageId-nSkip].t();
                cout << "GT orientation: " << R_GT << endl;
                cout << "GT direction: " << T_GT << endl;

                C3dPoint T;
                C3dRotation R;
                bool bSuccess = RTFromE(inlierMask, E, aCalibratedPoints1, aCalibratedPoints2, pointIds, R, T, sqr(dThresh));
                if(!bSuccess)
                {
                    cout << "Pure rotation detected, setting T=(0,0,1)" << endl;
                    T=C3dPoint(0,0,1);
                    continue;
                }

                double dErrR = diff(R_GT, R);
                double dErrT = angle(T_GT, T);
                statsRANSAC.addRT(dErrR, dErrT);

                const int nIters = 2;
                for(int i=0; i<nIters;i++)
                {
                    //CRefineEOnRTManifold::refineRobustOnMask(aCalibratedPoints1, aCalibratedPoints2, R, T, inlierMask);
                    CRefineEOnRTManifold::refineLSOnMask(aCalibratedPoints1, aCalibratedPoints2, R, T, inlierMask);
                    RTFromRT(inlierMask, aCalibratedPoints1, aCalibratedPoints2, pointIds, R, T, sqr(dThresh));
                }
                if(nIters == 0)
                {
                    Eigen::Matrix3d E_fromLinearLS;
                
                    refineWeightedLinearLSOnMask(aCalibratedPoints1, aCalibratedPoints2, E_fromLinearLS, inlierMask, false);
                    if (!RTFromE(inlierMask, E_fromLinearLS, aCalibratedPoints1, aCalibratedPoints2, pointIds, R, T, sqr(dThresh))) {
                        cout << "Error recovering R,T after linear weighted LS" << endl;
                        statsRefine.addFail();
                        continue;
                    }
                }

                cout << "Relative orientation: " << R << endl;
                cout << "Translation direction: " << T << endl;


                dErrR = diff(R_GT, R);
                dErrT = angle(T_GT, T);
                cout << dErrR << " = error in R, " << dErrT << "=error in T" << endl;
                statsRefine.addRT(dErrR, dErrT);
            }

            delete pLastDescriptors; pLastDescriptors=0;
            std::swap(pDescriptors, pLastDescriptors);
            std::swap(pLastFrame, pFrame);

            nImageId += nSkip; //Skip a few frames, otherwise have essentially pure rotation
        }
    }
    return 0;
}


