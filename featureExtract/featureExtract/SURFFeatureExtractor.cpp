/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * CSURFFeatureExtractor.cpp
 *
 *  Created on: 29/06/2009
 *      Author: tom
 */

#include "util/opencv.h"
#include "util/exception.h"
#include "SURFFeatureExtractor.h"
#include "description/vectorDescriptorFactory.h"
#include "util/smallHashTable.h"
#include "cornerDetect/goodFeaturesFast.h"
#include "time/SpeedTest.h"

#ifndef USE_OLD_OPENCV
#include "opencv2/features2d/features2d.hpp"
#endif
#include "opencv/cv.h" //CvSurfParams has been moved here in the latest OpenCV

/*bool CSURFFeatureExtractor::bInitialised = false;
CvSURFParams CSURFFeatureExtractor::SURFParams;*/

using namespace std;
CSURFFeatureExtractor::CSURFFeatureExtractor(const int MAX_FEATURES, const CImParams & IM_PARAMS, const CDescriptorSetClusteringParams & DSC_PARAMS, const CCornerParams::CSURFParams & SURFParams) : CFeatureExtractor(MAX_FEATURES, IM_PARAMS, DSC_PARAMS), SURFParams(SURFParams) {
    pCvMemStorage = cvCreateMemStorage(0);
    //	if(IM_CHANNELS == 3)
    //		pGreyImg = cvCreateImage(cvSize(IM_WIDTH, IM_HEIGHT), IPL_DEPTH_8U, 1);
}
CSURFFeatureExtractor::~CSURFFeatureExtractor() {
    cvReleaseMemStorage(&pCvMemStorage);
    //	if(IM_CHANNELS == 3)
    //		cvReleaseImage(&pGreyImg);
}
void CCornerParams::CSURFParams::getSurfParams(CvSURFParams & surfParams) const {
    surfParams.extended = EXTENDED;
    surfParams.hessianThreshold = HESSIAN_THRESH;
    surfParams.nOctaves = OCTAVES;
    surfParams.nOctaveLayers = OCTAVE_LAYERS;
    surfParams.upright = 1;
}

typedef pair<CvSURFPoint *, float *> CSURFPair;
bool betterThan(const CSURFPair & p1, const CSURFPair & p2) {
    return p1.first->hessian > p2.first->hessian;
}
CDescriptorSet * CSURFFeatureExtractor::getDescriptors_int(const IplImage * pImage) {
    if(IS_DEBUG) CHECK(!pGreyImg, "SURF: no grey image here");

    CvSeq * pKeypoints = 0, * pSURFDescs = 0;
    cvClearMemStorage(pCvMemStorage);

    CvSURFParams surfParams;
    SURFParams.getSurfParams(surfParams);

    CStopWatch s;
    s.startTimer();
    cvExtractSURF(pGreyImg, 0,
            &pKeypoints, &pSURFDescs,
            pCvMemStorage, surfParams);
    s.stopTimer();
    /*double tAll = s.getElapsedTime();
        s.startTimer();

    cvExtractSURF( pGreyImg, 0,
                &pKeypoints, &pSURFDescs,
                        pCvMemStorage, surfParams, 1);
    s.stopTimer();
    cout << "Extract SURF took " << tAll << " seconds, blobs only took " << tAll - s.getElapsedTime() << " seconds" << endl;*/

    int nCorners = min<int>(MAX_FEATURES, pKeypoints->total);
    ARRAY(CSURFPair, aSortedKeypoints, pKeypoints->total);

    CvSeqReader keyPointReader, SURFDescReader;
    cvStartReadSeq(pKeypoints, &keyPointReader, 0);
    cvStartReadSeq(pSURFDescs, &SURFDescReader, 0);
    for (int i = 0; i < pKeypoints->total; i++) {
        CvSURFPoint * pKeyPoint = (CvSURFPoint*) keyPointReader.ptr;
        float * pfDesc = (float*) SURFDescReader.ptr;

        aSortedKeypoints[i] = CSURFPair(pKeyPoint, pfDesc);

        CV_NEXT_SEQ_ELEM(pKeypoints->elem_size, keyPointReader);
        CV_NEXT_SEQ_ELEM(pSURFDescs->elem_size, SURFDescReader);
    }
    cout << pKeypoints->total << " SURF blobs found\n";

    //if(nCorners < pKeypoints->total) //Choose the best only
    std::sort(PTR(aSortedKeypoints), PTR(aSortedKeypoints) + pKeypoints->total, betterThan);

    CPointBin < 32, 32, 100 > pointBin(pImage->width, pImage->height);
    pointBin.reset();

    CDescriptorSet * pDS = newDescriptorSet(nCorners);

    //CSmallHashTable<int, 1024, 1, intHash<1023> > detectDuplicates;
    int nRad = SURFParams.CORNER_MIN_DIST;
    for (int i = 0; i < nCorners; i++) {
        CvSURFPoint * pKeyPoint = aSortedKeypoints[i].first;
        double x = pKeyPoint->pt.x;
        double y = pKeyPoint->pt.y;
        //if(x > MARGIN && y > MARGIN && x < pGreyImg->width-2*MARGIN && y < pGreyImg->height-2*MARGIN )
        CLocation loc(x, y);
        //if(!detectDuplicates.insert(loc.id()))
        if (!pointBin.isTooClose((int) x, (int) y, nRad)) {
            //TODO: use laplacian (prob won't actually make much difference)
            pDS->Push(CVectorDescriptorFactory::makeDescriptor(aSortedKeypoints[i].second, SURFParams.SURF_LENGTH(), pKeyPoint->dir, loc));
        }
    }


    return pDS;
}
