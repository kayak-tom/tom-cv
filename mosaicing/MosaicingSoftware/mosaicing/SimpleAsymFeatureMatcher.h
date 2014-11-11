#pragma once
#ifndef _SIMPLE_ASYM_FEATURE_MATCHER
#define _SIMPLE_ASYM_FEATURE_MATCHER

#include "FeatureMatcher.h"

namespace grc {

//! Implements simple algorithm to find correspondences between 2 frames by brute-force comparisons. 
class SimpleAsymFeatureMatcher : public FeatureMatcher
{
    //! Simple algorithm to find correspondences between 2 frames by brute-force comparisons. Accurate and not especially slow with <200 features per image.
    virtual CorrespondenceSet * getCorrespondences(const DescriptorSet * pDS1, const DescriptorSet * pDS2) const;
public:
    SimpleAsymFeatureMatcher(FeatureExtractor & featureExtractor) : FeatureMatcher(featureExtractor) {};
};

}


#endif // _SIMPLE_ASYM_FEATURE_MATCHER
