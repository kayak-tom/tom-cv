/*
 * BoWFeatureMatcher.cpp
 *
 *  Created on: 12 Apr 2010
 *      Author: tom
 */

#include "BoWFeatureMatcher.h"
#include "description/descriptor.h"
#include "bow/bagOfWordsParam.h"
#include "bow/bagOfWords.h"
#include "boost/smart_ptr.hpp"

namespace grc
{

BoWFeatureMatcher::BoWFeatureMatcher( ImageSourceSimple & imSource, CFeatureExtractor & featureExtractor, CBoW & bow, const CBOWMatchingParams & PARAMS)
 : FeatureMatcher2(imSource, featureExtractor), bow(bow), PARAMS(PARAMS)
{
}

BoWFeatureMatcher::~BoWFeatureMatcher() {
}

const CBoWCorrespondences * BoWFeatureMatcher::getCorrespondences(const CDescriptorSet * pDS1, const CDescriptorSet * pDS2) const
{
	THROW("BoWFeatureMatcher::getCorrespondences: Do not call: DSs should be owned by BoW")
	/*const CBoWCorrespondences * pCorr = pDS1->getBruteForceCorrespondenceSet(pDS2, PARAMS.BoW_CORNER_CONDITION, 0.6, PARAMS.MATCH_NN, 0);
	std::cout << pCorr->size() << " correspondences found from " << pDS1->Count() << 'x' << pDS2->Count() << " descriptors\n";
	return pCorr;*/
}

const CBoWCorrespondences * BoWFeatureMatcher::getCorrespondenceSet_int(IdPair & idPair)
{
	const int imageId1 = idPair.im1Id();
	const int imageId2 = idPair.im2Id();

	//Make sure they're both in BoW DB
	for(int i = 0; i<2; i++)
	{
		const int id = i ? imageId2 : imageId1;

		if(!bow.contains(id))
		{
			const IplImage * pIm1 = imSource.getImage(id);
			CDescriptorSet *features1 = featureExtractor_.getDescriptors(pIm1);
			if(features1->Count() < MIN_CORRESPONDENCES)
			{
				cout << "Failed to enough features in " << imageId1 << " " << endl;
				CDescriptorSet::deleteDS(&features1);
				return 0;
			}
			else
				bow.addImage(&features1, id);
		}
	}

	const CBoWCorrespondences * pCorrs = 0;

	// check there's enough to bother with, return 0 and don't recompute if not
	//CStopWatch s;
	//s.startTimer();
	pCorrs = bow.getCorrespondences(imageId2, imageId1, PARAMS);
	//s.stopTimer();
	//double numMatches = features1->size() * features2->size();
	//std::cout << "Match took " << s.getElapsedTime()/numMatches << " secs\n";
	//cout << "Found " << pCorrs->size() << " correspondences between " << imageId1 << " and " << imageId2 << endl;

	if(!pCorrs || pCorrs->size() < MIN_CORRESPONDENCES)
	{
		cout << "Failed to enough correspondences between " << imageId1 << " and " << imageId2 << endl;
		delete pCorrs; pCorrs=0;
	}
	return pCorrs;
}

void BoWFeatureMatcher::getSimilarIds(const int nImageId, const int nTopN, TCloseIdSet & ids) const
{
	int nLastId = nImageId;
	if(ids.size() == 1)
		nLastId = ids[0];
	else if (ids.size() > 1)
		THROW("Expected 0 or 1 id here");

	boost::scoped_ptr<TBoWMatchVector> pMatches ( bow.getMatches(nImageId, nTopN+5));
	if(pMatches)
	{
		for(TBoWMatchVector::const_iterator pMatch = pMatches->begin(); pMatch != pMatches->end() && (int)ids.size() < nTopN; pMatch++)
		{
			if(pMatch->id() < nLastId)
			{
				ids.push_back(pMatch->id());
			}
		}
	}
}

}
