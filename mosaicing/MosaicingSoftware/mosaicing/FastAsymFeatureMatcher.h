#pragma once
#ifndef _FAST_ASYM_FEATURE_MATCHER
#define _FAST_ASYM_FEATURE_MATCHER

#include "FeatureMatcher.h"

namespace grc {

//! Does an initial approximate match between 2 features to speed up matching (by only looking at 1 in 4 pixel vals on the first pass)
class FastAsymFeatureMatcher : public FeatureMatcher
{
    double conditionNum_; //!< Only potential matches with the next best match less than conditionNum*best match strength will be accepted. Values of about 0.75 are suitable.

    virtual CorrespondenceSet * getCorrespondences(const DescriptorSet * pDS1, const DescriptorSet * pDS2) const;
public:
    /*! \param conditionNum Only potential matches with the next best match less than conditionNum*best match strength will be accepted. Values of about 0.75 are suitable. */
    FastAsymFeatureMatcher(FeatureExtractor & featureExtractor, double conditionNum) : FeatureMatcher(featureExtractor), conditionNum_(conditionNum) {};
};

}


#endif // _FAST_ASYM_FEATURE_MATCHER
