#pragma once
#ifndef _GRC_RANSAC_PERSPECTIVE_TRANSFORM_ESTIMATOR_
#define _GRC_RANSAC_PERSPECTIVE_TRANSFORM_ESTIMATOR_

#include "TransformEstimator.h"

namespace grc {
	
    //! Use RANSAC to estimate transformation from a correspondence set.
	class RansacTransformEstimator: public TransformEstimator 
    {
        //! Implements RANSAC (Least-squares estimation implemented by transformation).
		virtual Transform * getTransform(const CorrespondenceSet * correspondences) const;
		int maxRansacIters_;
    public:

        RansacTransformEstimator(FeatureMatcher & featureMatcher, int maxRansacIters) : TransformEstimator(featureMatcher), maxRansacIters_(maxRansacIters) { }
	};

}

#endif