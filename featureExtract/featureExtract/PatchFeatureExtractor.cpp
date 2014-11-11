/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * PatchFeatureExtractor.cpp
 *
 *  Created on: 29/06/2009
 *      Author: tom
 */

#include "PatchFeatureExtractor.h"
#include "description/patchFactory.h"
#include "opencv2/opencv.hpp"
#include "time/SpeedTest.h"
#include "image/convert_OpenCV.h"

using namespace std;

CPatchFeatureExtractor::CPatchFeatureExtractor(const int MAX_FEATURES, const CImParams & IM_PARAMS, CCornerDetector ** ppCornerDetector, const CPatchDescriptorParams & PATCH_PARAMS, const CDescriptorSetClusteringParams & DSC_PARAMS) : CFeatureExtractor(MAX_FEATURES, IM_PARAMS, DSC_PARAMS), pCornerDetector(*ppCornerDetector), aCorners(new CLocation[MAX_FEATURES]), PATCH_PARAMS(PATCH_PARAMS) {
    //	cout << PATCH_PARAMS.Patch.PATCH_SCALE << (int)PATCH_PARAMS.Patch.PATCH_SCALE << endl;
    if(IS_DEBUG) CHECK(!pCornerDetector, "No corner detector");
    //	bColour = (3 == IM_PARAMS.IM_CHANNELS && !bMonoDescriptor);
    *ppCornerDetector = 0;
}

CPatchFeatureExtractor::~CPatchFeatureExtractor() {
    delete pCornerDetector;
    delete [] aCorners;
}

CDescriptorSet * CPatchFeatureExtractor::getDescriptors_int(const IplImage * pImage) {
    int nCorners = MAX_FEATURES;
    //const bool bColour = pImage->nChannels == 3 && !bMonoDescriptor;

    CStopWatch s;
    s.startTimer();

    pCornerDetector->getCorners(pGreyImg, nCorners, aCorners);

    s.stopTimer();
    cout << "Extract corners took " << s.getElapsedTime() << " seconds" << endl;

    /*if (PATCH_PARAMS.Patch.PATCH_BLUR > 0) {
//        cvGaussianBlur(!PATCH_PARAMS.Patch.MONO_DESCRIPTOR ? pImage : pGreyImg, !PATCH_PARAMS.Patch.MONO_DESCRIPTOR ? (void*) pImage : pGreyImg, PATCH_PARAMS.Patch.PATCH_BLUR * 2 + 1, PATCH_PARAMS.Patch.PATCH_BLUR * 2 + 1);
//        static bool bWarned = false;
//        if (!bWarned)
//            cout << "Warning: dodgy smooth\n", bWarned = true;
        THROW("PATCH_PARAMS.Patch.PATCH_BLUR is deprecated, set it to 0")
    }*/

    CDescriptorSet * pDS = newDescriptorSet(nCorners);

    for (int i = 0; i < nCorners; i++) {
        CHECK(PATCH_PARAMS.Patch.MONO_DESCRIPTOR && !pGreyImg, "Mono descriptor selected but no grey image here");
        CDescriptor * pDescriptor = TDescriptorFactory::makeDescriptor(PATCH_PARAMS.Patch.MONO_DESCRIPTOR ? pGreyImg : pImage, aCorners[i], 1.0, PATCH_PARAMS);

        pDS->Push(pDescriptor);
        pDescriptor = 0;
    }
    //cvShowImage("Test", pGreyImg);
    //cvWaitKey(0);

    return pDS;
}

CDSPatchFeatureExtractor::CDSPatchFeatureExtractor(const int MAX_FEATURES, const CImParams & IM_PARAMS, CCornerDetector ** ppCornerDetector, CPatchDescriptorParams & PATCH_PARAMS, const CDescriptorSetClusteringParams & DSC_PARAMS) : CPatchFeatureExtractor(MAX_FEATURES, IM_PARAMS, ppCornerDetector, PATCH_PARAMS, DSC_PARAMS), pDownsampledImage(0), nScale(PATCH_PARAMS.Patch.PATCH_SCALE) {
    //if(nScale > 1)
    {
        CvSize DSSize = cvSize(IM_PARAMS.IM_WIDTH / nScale, IM_PARAMS.IM_HEIGHT / nScale);
        pDownsampledImage = cvCreateImage(DSSize, IPL_DEPTH_32S, PATCH_PARAMS.Patch.MONO_DESCRIPTOR ? 1 : 3);
    }

    THROW("TODO: Need the next line:PATCH_PARAMS.PATCH_SCALE = 1")
            //PATCH_PARAMS.PATCH_SCALE = 1; //Set globally
}

CDSPatchFeatureExtractor::~CDSPatchFeatureExtractor() {
    if (pDownsampledImage)
        cvReleaseImage(&pDownsampledImage);
}

CDescriptorSet * CDSPatchFeatureExtractor::getDescriptors_int(const IplImage * pImage) {
    int nCorners = MAX_FEATURES;

    pCornerDetector->getCorners(pGreyImg, nCorners, aCorners);

    const double dScaleInv = 1.0 / nScale;
    /*if(nScale > 1)
    {*/
    if (!PATCH_PARAMS.Patch.MONO_DESCRIPTOR)
        doDownSample(nScale, pImage, pDownsampledImage);
    else
        doDownSample(nScale, pGreyImg, pDownsampledImage);
    /*}
    else
            pDownsampledImage = bColour ? pImage : pGreyImg;*/

    CStopWatch s;
    s.startTimer();

    CDescriptorSet * pDS = newDescriptorSet(nCorners);

    for (int i = 0; i < nCorners; i++) {
        CHECK(!pDownsampledImage, "No downsampled image");
        CDescriptor * pDescriptor = TDescriptorFactory::makeDescriptor(pDownsampledImage, aCorners[i], dScaleInv, PATCH_PARAMS);
        pDS->Push(pDescriptor);
        pDescriptor = 0;
    }
    s.stopTimer();
    cout << "Extract " << nCorners << " patches took " << s.getElapsedTime() << " seconds\n";

    return pDS;
}

