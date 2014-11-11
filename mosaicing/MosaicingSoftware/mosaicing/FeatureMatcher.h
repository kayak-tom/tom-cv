#pragma once
#ifndef _FEATURE_MATCHER
#define _FEATURE_MATCHER

#include "CorrespondenceSet.h"
#include "FeatureExtractor.h"
#include "TransformSet.h"

#include "featureExtract/cornerDetector.h"
#include "featureExtract/featureExtractor.h"

namespace grc {

    //! Classes derived from this find sets of correspondences between image pairs. Cacheing implemented here.
    class FeatureMatcher
    {
        FeatureExtractor & featureExtractor_; //!< Source for features/descriptors.

        typedef std::map<IdPair, CorrespondenceSet *> CachedCorrespondences; //!< Cache correspondence sets
        CachedCorrespondences correspondenceCache_; //!<Cache correspondence sets

        boost::mutex lockMatcherCache_; //!< Lock while matching features: only one thread can get corresps. at a time at the moment (likewise the feature xtor)
    protected:
        //! Returns correspoondences between 2 images. 
        virtual CorrespondenceSet * getCorrespondences(const DescriptorSet * pDS1, const DescriptorSet * pDS2) const = 0;
    public:
        FeatureMatcher(FeatureExtractor & featureExtractor) : featureExtractor_(featureExtractor) {};
        virtual ~FeatureMatcher();

        //! Returns correspoondences between 2 images, implements caching.
        const CorrespondenceSet * getCorrespondenceSet(size_t imageId1, size_t imageId2);
    };

    class FeatureMatcher2
    {
    protected:
        static const int MIN_CORRESPONDENCES = 4;
    	ImageSourceSimple & imSource;
    	CFeatureExtractor & featureExtractor_; //!< Source for features/descriptors.

        typedef std::map<IdPair, const CBoWCorrespondences *> CachedCorrespondences; //!< Cache correspondence sets
        CachedCorrespondences correspondenceCache_; //!<Cache correspondence sets

        boost::mutex lockMatcherCache_; //!< Lock while matching features: only one thread can get corresps. at a time at the moment (likewise the feature xtor)
        //! Returns correspoondences between 2 images.
        virtual const CBoWCorrespondences * getCorrespondences(const CDescriptorSet * pDS1, const CDescriptorSet * pDS2) const = 0;

        virtual const CBoWCorrespondences * getCorrespondenceSet_int(IdPair & idPair);
    public:
        FeatureMatcher2(ImageSourceSimple & imSource, CFeatureExtractor & featureExtractor) : imSource(imSource), featureExtractor_(featureExtractor) {};
        virtual ~FeatureMatcher2();

        //! Returns correspoondences between 2 images, implements caching.
        const CBoWCorrespondences * getCorrespondenceSet(size_t imageId1, size_t imageId2);

        virtual void getSimilarIds(const int nImageId, const int nTopN, TCloseIdSet & ids) const = 0;
    };
}

#endif // _FEATURE_MATCHER
