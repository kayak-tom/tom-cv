/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * FASTCornerDetector.cpp
 *
 *  Created on: 17/06/2009
 *      Author: tom
 */

#include "FASTCornerDetector.h"
#include "fast.h"
#include "description/descriptor.h"
#include <boost/thread/mutex.hpp>
#include "../cornerDetect/goodFeaturesFast.h"
#include <boost/smart_ptr.hpp>

using namespace std;

class CFastCornerSorter
{
	static int* scores;
	static boost::mutex mxScoresStatic;

	static bool betterThan(const int i, const int j) { return scores[i] > scores[j]; };
public:
	static void sort(int * pnBegin, int * pnEnd, int * scores_in)
	{
		boost::mutex::scoped_lock lock(mxScoresStatic);
		scores=scores_in;
		std::sort(pnBegin, pnEnd, CFastCornerSorter::betterThan);
		scores=0;
	}
};

int * CFastCornerSorter::scores = 0;
boost::mutex CFastCornerSorter::mxScoresStatic;

xy * fast9_detect_nonmax_binned(const byte* im, int xsize, int ysize, int stride, int FASTCORNER_T, const CCornerParams::CCornerDetectorParams & CORNER_PARAMS, int* pnTargetFeaturesInOut)
{
	xy * corners = 0;
	xy * nonmaxCorners = 0;
	int num_corners = -1;
	int * scores = 0;

	do
	{
		free(corners);
		free(nonmaxCorners);
		corners = fast9_detect(im, xsize, ysize, stride, FASTCORNER_T, &num_corners);
		scores = fast9_score(im, stride, corners, num_corners, FASTCORNER_T);
		nonmaxCorners = nonmax_suppression(corners, scores, num_corners, &num_corners);
		free(scores);
		FASTCORNER_T /= 2;
	} while (FASTCORNER_T > 0 && num_corners < *pnTargetFeaturesInOut/2);

	scores = fast9_score(im, stride, nonmaxCorners, num_corners, FASTCORNER_T);

	ARRAY(int, anIndices, num_corners);
	for(int i=0; i<num_corners; i++)
		anIndices[i] = i;

	CFastCornerSorter::sort(PTR(anIndices), PTR(anIndices)+num_corners, scores);

	CPointBin<32, 32, 100> pointBin(xsize, ysize);
	pointBin.reset();

	//FIRST PASS: Go for SPATIAL SPREAD
	int nSeperation = CORNER_PARAMS.FASTCORNER_TWO_PASS_BINNING ? 64 : CORNER_PARAMS.FASTCORNER_MIN_SEPERATION;
	int nNumGoodCorners=0;
	xy* aNewCorners = corners, * pNewCorners = aNewCorners;
	for(int i=0; i<num_corners; i++)
	{
		int nCornerIdx = anIndices[i]; //should now be in decreasing order of score
		if(nCornerIdx>=0)
		{
			if(!pointBin.isTooClose(nonmaxCorners[nCornerIdx].x, nonmaxCorners[nCornerIdx].y, nSeperation))
			{
				*pNewCorners = nonmaxCorners[nCornerIdx];
				pNewCorners++;
				nNumGoodCorners++;
				anIndices[i] = -1; //We're using this corner, don't re-add later
				if(nNumGoodCorners==*pnTargetFeaturesInOut)
					break;
			}
			if(i==num_corners-1 && nSeperation != CORNER_PARAMS.FASTCORNER_MIN_SEPERATION)
			{
				//We've got to the end--drop min seperation and start again with the corners that were too close:
				i=-1;
				nSeperation = CORNER_PARAMS.FASTCORNER_MIN_SEPERATION;
			}
		}
	}
	//SECOND PASS: Go for BEST CORNERS to make up number


	*pnTargetFeaturesInOut = nNumGoodCorners;

	free(nonmaxCorners);
	free(scores);
	return aNewCorners;
}

xy* getFastCorners(const IplImage & im, const int FASTCORNER_T, const CCornerParams::CCornerDetectorParams & CORNER_PARAMS, int* pnTargetFeaturesInOut)
{
	const int w = im.width, h = im.height, stride=im.widthStep;
	const byte * data = (const byte *)im.imageData;
	if(CORNER_PARAMS.FASTCORNER_MIN_SEPERATION<=0)
	{
		return fast9_detect_nonmax(data, w, h, stride, FASTCORNER_T, pnTargetFeaturesInOut); //Could switch 9/10/11/12...
	}
	else
	{
		return fast9_detect_nonmax_binned(data, w, h, stride,  FASTCORNER_T, CORNER_PARAMS, pnTargetFeaturesInOut); //Could switch 9/10/11/12...
	}
}

CFASTCornerDetector::CFASTCornerDetector(const CImParams& IM_PARAMS, int MARGIN, const CCornerParams::CCornerDetectorParams & CORNER_PARAMS) : CCornerDetector(IM_PARAMS), CORNER_PARAMS(CORNER_PARAMS), MARGIN(MARGIN)
{
}

void CFASTCornerDetector::getCorners(IplImage * pImage, int & nCorners, CLocation * aCorners)
{
	CHECK(pImage->nChannels != 1, "Grey image expected");

	IplImage greySubImage = *pImage;
	CIplPx<uchar>::cropImage(greySubImage, MARGIN);

	if(IM_PARAMS.SAVE_FRAMES != CImParams::eDontSave)
	{
		cout << "Cut corners to save mem...\n";
		nCorners /= 20;
	}

	int nFastCorners = nCorners;

	if(CORNER_PARAMS.FASTCORNER_MIN_SEPERATION <= 1)
	{
		xy * axy = getFastCorners(greySubImage, CORNER_PARAMS.FASTCORNER_T, CORNER_PARAMS, &nFastCorners);
		if (nFastCorners < nCorners / 2)
		{
			free(axy);
			axy = getFastCorners(greySubImage, CORNER_PARAMS.FASTCORNER_T / 2, CORNER_PARAMS, &nFastCorners);
		}
		else if (nFastCorners > nCorners * 2)
		{
			free(axy);
			axy = getFastCorners(greySubImage, CORNER_PARAMS.FASTCORNER_T * 2, CORNER_PARAMS, &nFastCorners);
		}

		if(CORNER_PARAMS.USE_SUBPIX)
		{
			cout << "Warning--FAST and SUBPIX not supported together...";
			if(IS_DEBUG) CHECK(!CLocation::SUBPIX_SUPPORT(), "Not enough precision for subpix");
			/*const int SUBPIX_SIZE=3, ZZ_SIZE=-1;
			cvFindCornerSubPix( &greySubImage, aCvPointCorners,
									 nCorners, cvSize(SUBPIX_SIZE, SUBPIX_SIZE), cvSize(ZZ_SIZE, ZZ_SIZE),
									 cvTermCriteria(CV_TERMCRIT_ITER | CV_TERMCRIT_EPS, 10, 0.01) );*/
		}

		xy * pxy = axy;
		CLocation * pLoc = aCorners;

		for (int n=nCorners; n>0; n--)
		{
			if(!(CORNER_PARAMS.USE_SUBPIX && (pxy->x < 0 || pxy->y < 0 || pxy->x >= IM_PARAMS.IM_WIDTH-2*MARGIN || pxy->y >= IM_PARAMS.IM_HEIGHT-2*MARGIN)))
			{
				double x = pxy->x+MARGIN;
				double y = pxy->y+MARGIN;
				*pLoc = CLocation(x, y);

				pLoc++;
			}
			pxy++;
		}
		nCorners = pLoc - aCorners;
		free(axy);
	}
	else
	{
		xy * axy = getFastCorners(greySubImage, CORNER_PARAMS.FASTCORNER_T, CORNER_PARAMS, &nFastCorners);

		nCorners = min<int>(nCorners, nFastCorners);

		xy * pxy = axy;
		CLocation * pLoc = aCorners;

		for (int i=0; i<nCorners; )
		{
			if(pxy->x >= 0 || pxy->y >= 0 )
			{
				double x = pxy->x+MARGIN;
				double y = pxy->y+MARGIN;
				*pLoc = CLocation(x, y);
				pLoc++;
				i++;
			}
			pxy++;
		}
		free(axy);
	}
}

CFASTCornerDetector::~CFASTCornerDetector()
{
	// TODO Auto-generated destructor stub
}
