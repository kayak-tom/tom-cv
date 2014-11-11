/*
 * BFFeatureMatcher.h
 *
 *  Created on: 12 Apr 2010
 *      Author: tom
 */

#ifndef BFFEATUREMATCHER_H_
#define BFFEATUREMATCHER_H_

#include "FeatureMatcher.h"
class CBoWCorrespondences;
class CDescriptorSet;
class CBOWMatchingParams;

namespace grc
{
class BFFeatureMatcher: public grc::FeatureMatcher2 {
	const CBOWMatchingParams & PARAMS;
protected:
    virtual const CBoWCorrespondences * getCorrespondences(const CDescriptorSet * pDS1, const CDescriptorSet * pDS2) const;
public:
	BFFeatureMatcher(ImageSourceSimple & imSource, CFeatureExtractor & featureExtractor, const CBOWMatchingParams & PARAMS);
	virtual ~BFFeatureMatcher();
};
}

#endif /* BFFEATUREMATCHER_H_ */
