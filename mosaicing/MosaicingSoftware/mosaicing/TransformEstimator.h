#pragma once
#ifndef _TRANSFORM_ESTIMATOR
#define _TRANSFORM_ESTIMATOR

#include "util/set2.h"
#include "FeatureMatcher.h"
#include "TransformSet.h"
#include "TransformInfo.h"
#include "util/opencv.h"

namespace grc {


	//! Derived classes estimate transforms between pairs of images. Implements cacheing of estimated transforms.
	class TransformEstimator
	{
	protected:
		
		FeatureMatcher & featureMatcher_;

        typedef std::map<IdPair, TransformInfo *> CachedTransformations; 
        CachedTransformations transformCache_; //!< Cache transforms.

        //! Returns transform between a pair of images
		/*! 
		 */
		virtual Transform * getTransform(const CorrespondenceSet * correspondences) const = 0;
        
        //! Test if a transform is good or not (transforms corrupted by outliers are likely to be highly skewed and have large scale changes)
        bool transIsGood(const TransformInfo * trans);
	public:
		TransformEstimator(FeatureMatcher & featureMatcher) : featureMatcher_(featureMatcher) {};
		~TransformEstimator();

        //! Returns transform between 2 images. Also implements cacheing.
		const TransformInfo * getTransform(size_t imageId1, size_t imageId2); 
	};

    //typedef std::set<size_t, std::greater<size_t> > TIntSetDesc;

	class TransformEstimator2
	{
	protected:

		FeatureMatcher2 & featureMatcher_;

        typedef std::map<IdPair, TransformInfo2 *> CachedTransformations;
        CachedTransformations transformCache_; //!< Cache transforms.

        //! Returns transform between a pair of images
		/*!
		 */
		virtual Transform * getTransform(const CBoWCorrespondences * correspondences, size_t imageId1, size_t imageId2, const Transform  * lastTrans, double dLastTransScale) const = 0;

        //! Test if a transform is good or not (transforms corrupted by outliers are likely to be highly skewed and have large scale changes)
        bool transIsGood(const TransformInfo2 * trans);
	public:
		TransformEstimator2(FeatureMatcher2 & featureMatcher) : featureMatcher_(featureMatcher) {};
		~TransformEstimator2();

        //! Returns transform between 2 images. Also implements cacheing.
		const TransformInfo2 * getTransform(size_t imageId1, size_t imageId2);

		void getSimilarIds(const int nImageId, const int nTopN, TCloseIdSet & ids) const;
	};

}

#endif //_TRANSFORM_ESTIMATOR
