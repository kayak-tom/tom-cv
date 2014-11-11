#pragma once
#ifndef _GRC_BaySAC_PERSPECTIVE_TRANSFORM_ESTIMATOR_
#define _GRC_BaySAC_PERSPECTIVE_TRANSFORM_ESTIMATOR_

#include "TransformEstimator.h"
#include "ransac/ransacParams.h"
#include "imageSource/imageSource.h"

namespace grc {
	
    //! Use BaySAC to estimate transformation from a correspondence set.
	class BaySACTransformEstimator: public TransformEstimator2
    {
        //! Implements BaySAC (Least-squares estimation implemented by transformation).
		virtual Transform * getTransform(const CBoWCorrespondences * correspondences, size_t imageId1, size_t imageId2, const Transform  * lastTrans, double dLastTransScale) const;
		const CRANSACHomographyParams & PARAMS;
		const CCamCalibMatrix & K;
		CImageSource * imSource;
		int nRansacSuccess , nRansacFail , nTotalInliers;
    public:
		~BaySACTransformEstimator();
        BaySACTransformEstimator(FeatureMatcher2 & featureMatcher, const CRANSACHomographyParams & PARAMS, const CCamCalibMatrix & K, CImageSource * imSource) : TransformEstimator2(featureMatcher), PARAMS(PARAMS), K(K), imSource(imSource), nRansacSuccess(0), nRansacFail(0), nTotalInliers(0) { }
	};

}

#endif
