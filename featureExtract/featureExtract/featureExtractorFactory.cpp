/*
 * featureExtractorFactory.cpp
 *
 *  Created on: 9/12/2009
 *      Author: tom
 */
#include "featureExtractor.h"
#include "util/exception.h"
#include "description/vectorDescriptor.h"

#include "SURFFeatureExtractor.h"
#include "PatchFeatureExtractor.h"
#include "cornerDetect/openCVCornerDetector.h"
#include "cornerDetect/orfeoCornerDetector.h"
#include "cornerDetect/subSampledCornerDetector.h"
#include "fast/FASTCornerDetector.h"

using namespace std;

CFeatureExtractor * CFeatureExtractor::makeFeatureExtractor(const CImParams & IMPARAMS, const CCornerParams & CORNERPARAMS,
			const CPatchDescriptorParams & PATCHDESCRIPTORPARAMS, const CDescriptorSetClusteringParams & DSCPARAMS)
{
	CCornerDetector * pCornerDetector = 0;
	CFeatureExtractor * pFeatureExtractor = 0;

	switch(CORNERPARAMS.SALIENT_FEATURE_TYPE)
	{
	/*case CCornerParams::eShiTomasiCorners:
		switch(CORNERPARAMS.CornerDetector.CORNER_MODE)
		{
		case CCornerParams::CCornerDetectorParams::eFasterOpenCVGoodFeatures:
			//pCornerDetector = new COpenCVCornerDetector(IMPARAMS, PATCHDESCRIPTORPARAMS.margin(), CORNERPARAMS);
			//TODO: Faster corners...
			//break;
		case CCornerParams::CCornerDetectorParams::eSubSampledCorners:
			//pCornerDetector = new CSubSampledCornerDetector(IMPARAMS, PATCHDESCRIPTORPARAMS.margin(), CORNERPARAMS);
			//break;
			cout << "TODO: Restore faster corner detection code...\n";
		case CCornerParams::CCornerDetectorParams::eOpenCVGoodFeatures:
			//pCornerDetector = (new COpenCVGoodFeaturesCornerDetector(IMPARAMS, PATCHDESCRIPTORPARAMS.margin(), CORNERPARAMS));
			//break;
		default:
			THROW("Unhandled corner detection mode");
		}
		break;*/
	case CCornerParams::eFastCorners:
		pCornerDetector = (new CFASTCornerDetector(IMPARAMS, PATCHDESCRIPTORPARAMS.margin(), CORNERPARAMS.CornerDetector));
		break;
	case CCornerParams::eSURFBlobs:
		pFeatureExtractor = (new CSURFFeatureExtractor(CORNERPARAMS.MAX_FEATURES, IMPARAMS, DSCPARAMS, CORNERPARAMS.SURF));
		break;
	default:
		THROW("Unhandled corner mode");
	}

	if(pCornerDetector)
	{
		pFeatureExtractor = (new CPatchFeatureExtractor(CORNERPARAMS.MAX_FEATURES, IMPARAMS, &pCornerDetector, PATCHDESCRIPTORPARAMS, DSCPARAMS));
	}

	CHECK(!pFeatureExtractor, "makeFeatureExtractor: No feature extractor made");
	return pFeatureExtractor;
}
