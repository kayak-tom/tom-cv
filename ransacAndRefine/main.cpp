/* 
 * Example of how to use my code to:
 *  - Load images
 *  - Extract patch descriptors
 *  - Match features (by brute-force comparison)
 *  - Use RANSAC/BaySAC/my fivepoint solver to compute the relative pose
 *  - Use a nonlinear optimiser to refine the pose
 * 
 * Must be linked to: libparams libransac libcamerageom boost_system boost_filesystem boost_thread opencv_core opencv_highgui opencv_imgproc 
 * libfeatureextract libfeaturedescription libimage libimagesource libtiming libutil opencv_calib3d 
 * 
 * Include path must include: -I../util -I../ransac -I../params -I../cameraGeom -I../image -I../imageSource -I../featureDescription -I../featureExtract
 * plus Boost, Eigen, and OpenCV
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

using namespace std;

/*
 * Example program to demonstrate feature matching, RANSAC for E, and relative pose refinement
 */
int main(int argc, char** argv) {

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
    while (pImageSource->loadImage(nImageId, pFrame)) {
        CDescriptorSet * pDescriptors = pFeatureExtractor->getDescriptors(pFrame);

        if (pLastDescriptors) //we're past the first frame
        {
            //Find matches between features
            boost::scoped_ptr<const CBoWCorrespondences> pCorr ( pLastDescriptors->getBruteForceCorrespondenceSet(pDescriptors,MS,0) );

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
        }

        delete pLastDescriptors; pLastDescriptors=0;
        std::swap(pDescriptors, pLastDescriptors);
        std::swap(pLastFrame, pFrame);
        nImageId += 10; //Skip a few frames, otherwise have essentially pure rotation
    }
    return 0;
}

