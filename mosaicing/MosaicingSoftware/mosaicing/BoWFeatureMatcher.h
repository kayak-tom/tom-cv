/*
 * BoWFeatureMatcher.h
 *
 *  Created on: 12 Apr 2010
 *      Author: tom
 */

#ifndef BoWFEATUREMATCHER_H_
#define BoWFEATUREMATCHER_H_

#include "FeatureMatcher.h"
class CBoWCorrespondences;
class CDescriptorSet;
class CBOWMatchingParams;
class CBoW;

namespace grc
{

class BoWFeatureMatcher: public grc::FeatureMatcher2 {
	CBoW & bow;
	const CBOWMatchingParams & PARAMS;
protected:
	const CBoWCorrespondences * getCorrespondenceSet_int(IdPair & idPair);
	const CBoWCorrespondences * getCorrespondences(const CDescriptorSet * pDS1, const CDescriptorSet * pDS2) const;
public:
	BoWFeatureMatcher(ImageSourceSimple & imSource, CFeatureExtractor & featureExtractor, CBoW & bow, const CBOWMatchingParams & PARAMS);
	virtual ~BoWFeatureMatcher();
	virtual void getSimilarIds(const int nImageId, const int nTopN, TCloseIdSet & ids) const;
};

}

#endif /* BoWFEATUREMATCHER_H_ */
