#pragma once
#ifndef _GRC_DUMB_TRANSFORM_ESTIMATOR_
#define _GRC_DUMB_TRANSFORM_ESTIMATOR_

#include "TransformEstimator.h"
namespace grc {
	class DumbTransformEstimator : public TransformEstimator
    {
    public:
        DumbTransformEstimator(FeatureMatcher & featureMatcher) : TransformEstimator(featureMatcher) { };
		const Transform * getTransform(size_t imageId1, size_t imageId2);
	};
}

#endif
