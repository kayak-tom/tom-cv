/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * openCVCornerDetector.cpp
 *
 *  Created on: 19/05/2009
 *      Author: tom
 */

#include "subSampledCornerDetector.h"
#include "description/descriptor.h"
#include "time/SpeedTest.h"
#include <iomanip>

using namespace std;

CSubSampledCornerDetector::CSubSampledCornerDetector(const CImParams & IM_PARAMS, int MARGIN, const CCornerParams & FEATURE_PARAMS) : CCornerDetector(IM_PARAMS), MARGIN(MARGIN), MARGIN_INT(MARGIN),
	eigImg(cvCreateImage(cvSize(IM_PARAMS.IM_WIDTH/2 - MARGIN, IM_PARAMS.IM_HEIGHT/2 - MARGIN), IPL_DEPTH_32F, 1)),
	tempImg(cvCreateImage(cvSize(IM_PARAMS.IM_WIDTH/2 - MARGIN, IM_PARAMS.IM_HEIGHT/2 - MARGIN), IPL_DEPTH_32F, 1)),
    aCvPointCorners(new CvPoint2D64f[FEATURE_PARAMS.MAX_FEATURES]),
    pointBin(IM_PARAMS.IM_WIDTH/2 - MARGIN, IM_PARAMS.IM_HEIGHT/2 - MARGIN), CORNER_PARAMS(FEATURE_PARAMS.CornerDetector)
{
	pGreySSImage = cvCreateImage(cvSize(IM_PARAMS.IM_WIDTH/2-MARGIN, IM_PARAMS.IM_HEIGHT/2-MARGIN), IPL_DEPTH_8U, 1);
}

CSubSampledCornerDetector::~CSubSampledCornerDetector()
{
	delete [] aCvPointCorners;
	cvReleaseImage(&pGreySSImage);
}

#define PRINTWH(im) cout << #im " width=" << (im)->width << " height=" << (im)->height << endl;
void CSubSampledCornerDetector::getCorners(IplImage * pGreyImgUse, int & nCorners, CLocation * aCorners)
{
	CHECK( pGreyImgUse->nChannels != 1, "Grey image required");
	//if(IS_DEBUG) CHECK(pGreySSImage->width != pGreySSImage->widthStep, "Width assumption failed");

	CStopWatch s; s.startTimer();
	char * pcRowStart = pGreyImgUse->imageData + MARGIN_INT * pGreyImgUse->widthStep + MARGIN_INT;
	char * pcDest = pGreySSImage->imageData;
	for(int r = 0; r < IM_PARAMS.IM_HEIGHT/2-MARGIN_INT; r++)
	{
		char * pcRow = pcRowStart + r * 2 * pGreyImgUse->widthStep;
		char * pcDestRow = pcDest + r * pGreySSImage->widthStep;
		for(int c = IM_PARAMS.IM_WIDTH/2-MARGIN_INT; c>0; c--)
		{
			*pcDestRow = *pcRow; //Sampling 'top left' of 4--need to change below +0.5 as well if interpolate.
			pcRow += 2;
			pcDestRow++;
		}
	}
	int nPatchSize = CORNER_PARAMS.CORNER_BLOCK_RAD > 1 ? (CORNER_PARAMS.CORNER_BLOCK_RAD/2)*2+1 : 3;
	
    fastCvGoodFeaturesToTrack(pGreySSImage, eigImg, tempImg, aCvPointCorners, &nCorners, CORNER_PARAMS.CORNER_QUAL, CORNER_PARAMS.CORNER_MIN_DIST/2, 0, nPatchSize,  CORNER_PARAMS.USE_HARRIS, CORNER_PARAMS.CORNER_HARRIS_K, pointBin);

    CvPoint2D64f * pCvPointFloat = aCvPointCorners;
    for (int n=nCorners; n>0; n--)
    {
		pCvPointFloat->x = pCvPointFloat->x*2 + MARGIN + 0.5;
		pCvPointFloat->y = pCvPointFloat->y*2 + MARGIN + 0.5;
		pCvPointFloat++;
	}
	s.stopTimer();
    cout << setprecision(6) << "SubSampled Corner detect took " << s.getElapsedTime() << " secs, " << nCorners << " detected\n";

	if(IS_DEBUG) CHECK(!CLocation::SUBPIX_SUPPORT(), "Not enough precision for subpix");
	const int SUBPIX_SIZE=5, ZZ_SIZE=-1;

	s.startTimer();
	fastCvFindCornerSubPix( pGreyImgUse, aCvPointCorners,
							 nCorners, cvSize(SUBPIX_SIZE, SUBPIX_SIZE), cvSize(ZZ_SIZE, ZZ_SIZE),
							 cvTermCriteria(CV_TERMCRIT_ITER | CV_TERMCRIT_EPS, 8, 0.15) );
	s.stopTimer();
    cout << "Subpix took " << s.getElapsedTime()/nCorners << " secs per corner (" << s.getElapsedTime() << " total)\n";

    pCvPointFloat = aCvPointCorners;
    CLocation * pLoc = aCorners;

    for (int n=nCorners; n>0; n--)
    {
    	if((pCvPointFloat->x < MARGIN || pCvPointFloat->y < MARGIN || pCvPointFloat->x >= IM_PARAMS.IM_WIDTH-MARGIN || pCvPointFloat->y >= IM_PARAMS.IM_HEIGHT-MARGIN))
    	{
    		pCvPointFloat++;
    		continue;
    	}

        double x = pCvPointFloat->x;
        double y = pCvPointFloat->y;
        *pLoc = CLocation(x, y);
    	if(IS_DEBUG) CHECK(pLoc->x() == 0 || pLoc->y() == 0, "Found OOB corner");

        pCvPointFloat++;
        pLoc++;
    }
    nCorners = pLoc - aCorners;
    for (int n=0; n<nCorners; n++)
    {
    	if(IS_DEBUG) CHECK(aCorners[n].x() == 0 || aCorners[n].y() == 0, "Found OOB corner");
    }
}
