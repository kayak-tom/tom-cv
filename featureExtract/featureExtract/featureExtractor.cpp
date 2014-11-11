/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * featureExtractor.cpp
 *
 *  Created on: 19/05/2009
 *      Author: tom
 */

#include "featureExtractor.h"
#include "util/exception.h"
#include "description/vectorDescriptor.h"
#include "image/convert_OpenCV.h"
#include "time/SpeedTest.h"

using namespace std;

CDescriptorSet * CFeatureExtractor::newDescriptorSet(int nDescriptorEstimate)
{
	CDescriptorSet * pDS = new CMetricSpaceDescriptorSet(DSC_PARAMS, nDescriptorEstimate);

	return pDS;
}

CDescriptorSet * CFeatureExtractor::getDescriptors(const IplImage * pImage)
{
	CStopWatch s; s.startTimer();

	const int nChannels = pImage->nChannels;

	if(nChannels == 1)
	{
		pGreyImg = const_cast<IplImage *>(pImage);
	}
	else
	{
		CHECK(!pGreyImg, "Haven't created a grey image for corner detection");
		//cvCvtColor(pImage, pGreyImg, CV_RGB2GRAY);
		greyScaler.greyScale(pImage, pGreyImg);
	}

	CDescriptorSet * pDS = getDescriptors_int(pImage);

	s.stopTimer();
	REPEAT(20, cout << "Extract corners and describe features took " << s.getElapsedTime() << " seconds\n");

	return pDS;
}

CFeatureExtractor::CFeatureExtractor(const int MAX_FEATURES, const CImParams & IM_PARAMS_IN, const CDescriptorSetClusteringParams & DSC_PARAMS) :
	greyScaler(IM_PARAMS_IN.IM_CHANNELS>1, IM_PARAMS_IN.Greyscale.R, IM_PARAMS_IN.Greyscale.G, IM_PARAMS_IN.Greyscale.B, IM_PARAMS_IN.Greyscale.GAMMA), DSC_PARAMS(DSC_PARAMS), MAX_FEATURES(MAX_FEATURES), IM_PARAMS(IM_PARAMS_IN), pGreyImg(0)
{
	if(IM_PARAMS.IM_CHANNELS == 3)
		pGreyImg = cvCreateImage(cvSize(IM_PARAMS.IM_WIDTH, IM_PARAMS.IM_HEIGHT), IPL_DEPTH_8U, 1);
}

CFeatureExtractor::~CFeatureExtractor()
{
	if(IM_PARAMS.IM_CHANNELS == 3)
		cvReleaseImage(&pGreyImg);
}

void markDescriptors(IplImage * pImage, const CDescriptorSet * pDesc)
{
	CHECK(!pImage || !pDesc, "markDescriptors: Null parameter--has the descriptor set been given to the BoW?");

	for (int nPoint = 0; nPoint < pDesc->Count(); nPoint++)
	{
		CvPoint point = locToCvPoint(pDesc->get_const(nPoint)->location());
		cvCircle( pImage, point, 3, CV_RGB( 255, 0, 0 ), -1 );
	}
}

