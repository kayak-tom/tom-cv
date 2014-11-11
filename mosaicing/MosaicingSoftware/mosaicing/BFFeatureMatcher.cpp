/*
 * BFFeatureMatcher.cpp
 *
 *  Created on: 12 Apr 2010
 *      Author: tom
 */

#include "BFFeatureMatcher.h"
#include "description/descriptor.h"
#include "bow/bagOfWordsParam.h"

namespace grc
{

BFFeatureMatcher::BFFeatureMatcher( ImageSourceSimple & imSource, CFeatureExtractor & featureExtractor, const CBOWMatchingParams & PARAMS)
 : FeatureMatcher2(imSource, featureExtractor), PARAMS(PARAMS)
{
}

BFFeatureMatcher::~BFFeatureMatcher() {
}

const CBoWCorrespondences * BFFeatureMatcher::getCorrespondences(const CDescriptorSet * pDS1, const CDescriptorSet * pDS2) const
{
        const CMatchableDescriptors::CMatchSettings matchSettings(0.6, PARAMS.BF_CORNER_CONDITION, PARAMS.MATCH_NN);
	const CBoWCorrespondences * pCorr = pDS1->getBruteForceCorrespondenceSet(pDS2, matchSettings, 0);
	std::cout << pCorr->size() << " correspondences found from " << pDS1->Count() << 'x' << pDS2->Count() << " descriptors\n";
	return pCorr;
}

}
