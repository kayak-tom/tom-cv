/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * openCVCornerDetector.cpp
 *
 *  Created on: 19/05/2009
 *      Author: tom
 */

#include "openCVCornerDetector.h"
#include "description/descriptor.h"
#include "time/SpeedTest.h"
#include <iomanip>

using namespace std;

COpenCVCornerDetector::COpenCVCornerDetector(const CImParams & IM_PARAMS, int MARGIN, const CCornerParams & FEATURE_PARAMS) : CCornerDetector(IM_PARAMS), MARGIN(MARGIN),
	eigImg(cvCreateImage(cvSize(IM_PARAMS.IM_WIDTH - 2*MARGIN, IM_PARAMS.IM_HEIGHT - 2*MARGIN), IPL_DEPTH_32F, 1)),
	tempImg(cvCreateImage(cvSize(IM_PARAMS.IM_WIDTH - 2*MARGIN, IM_PARAMS.IM_HEIGHT - 2*MARGIN), IPL_DEPTH_32F, 1)),
    aCvPointCorners(new CvPoint2D64f[FEATURE_PARAMS.MAX_FEATURES]),
    pointBin(IM_PARAMS.IM_WIDTH - 2*MARGIN, IM_PARAMS.IM_HEIGHT - 2*MARGIN), CORNER_PARAMS(FEATURE_PARAMS.CornerDetector)
{}

COpenCVCornerDetector::~COpenCVCornerDetector()
{
	delete [] aCvPointCorners;
}

#define PRINTWH(im) cout << #im " width=" << (im)->width << " height=" << (im)->height << endl;
void COpenCVCornerDetector::getCorners(IplImage * pGreyImgUse, int & nCorners, CLocation * aCorners)
{
	CHECK( pGreyImgUse->nChannels != 1, "Grey image required");

	IplImage greySubImage = *pGreyImgUse;
	CIplPx<uchar>::cropImage(greySubImage, (int)MARGIN);

	CStopWatch s; s.startTimer();
    fastCvGoodFeaturesToTrack(&greySubImage, eigImg, tempImg, aCvPointCorners, &nCorners, CORNER_PARAMS.CORNER_QUAL, CORNER_PARAMS.CORNER_MIN_DIST, 0, CORNER_PARAMS.CORNER_BLOCK_RAD*2+1, CORNER_PARAMS.USE_HARRIS, CORNER_PARAMS.CORNER_HARRIS_K, pointBin);
	s.stopTimer();
    cout << setprecision(6) << "Corner detect took " << s.getElapsedTime() << " secs, " << nCorners << " detected\n";


    if(CORNER_PARAMS.USE_SUBPIX)
    {
    	if(IS_DEBUG) CHECK(!CLocation::SUBPIX_SUPPORT(), "Not enough precision for subpix");
		const int SUBPIX_SIZE=3, ZZ_SIZE=-1;
    	s.startTimer();
    	fastCvFindCornerSubPix( &greySubImage, aCvPointCorners,
								 nCorners, cvSize(SUBPIX_SIZE, SUBPIX_SIZE), cvSize(ZZ_SIZE, ZZ_SIZE),
								 cvTermCriteria(CV_TERMCRIT_ITER | CV_TERMCRIT_EPS, 8, 0.15) );
    	s.stopTimer();
        cout << "Subpix took " << s.getElapsedTime()/nCorners << " secs per corner (" << s.getElapsedTime() << " total)\n";
    }

    CvPoint2D64f * pCvPointFloat = aCvPointCorners;
    CLocation * pLoc = aCorners;

    for (int n=nCorners; n>0; n--)
    {
    	if(CORNER_PARAMS.USE_SUBPIX && (pCvPointFloat->x < 0 || pCvPointFloat->y < 0 || pCvPointFloat->x >= IM_PARAMS.IM_WIDTH-2*MARGIN || pCvPointFloat->y >= IM_PARAMS.IM_HEIGHT-2*MARGIN))
    	{
    		pCvPointFloat++;
    		continue;
    	}

        double x = pCvPointFloat->x+MARGIN;
        double y = pCvPointFloat->y+MARGIN;
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
