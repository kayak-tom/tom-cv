#pragma once
#ifndef _FEATURE_EXTRACTOR_
#define _FEATURE_EXTRACTOR_

#include "Descriptor.h"
#include "boost/thread.hpp"
//#include "ImageSource.h" will cause circular refs
#include "ImageSourceSimple.h"
#include "FeatureDetector.h"

namespace grc {

//! Implements feature extractors and returns sets of extracted features (and their locations)
/*! TFeatureDetector should be derived from FeatureDetector
 *  TDescriptor should be derived from Descriptor
 */
class FeatureExtractor
{
public:
    //! Returns a set of descriptors from an image.
    virtual const DescriptorSet * getFeatures(size_t imageId) = 0;
};

//! Parent of feature extractor classes: detects corners, extracts descriptors, caches extracted features.
class FeatureExtractorImpl : public FeatureExtractor
{
	ImageSourceSimple & imageSource_; //! Used to retreive frames

    FeatureDetector * featureDetector_; //! Any feature/corner/blob detector can be set here.

    //! Returns set of features (descriptors and locations)
    /*! Gets image from imageSource_, extracts corners, descriptors, returns feature set.
     */
    const DescriptorSet * getFeaturesFromImage(size_t imageId);

    boost::mutex featureCacheLock_; //!< Allow concurrent users
    typedef std::map<size_t, const DescriptorSet *> TCachedFeatureSets;
    TCachedFeatureSets cachedDescriptorSets_; //!< Cache extracted descriptor sets
public:
    FeatureExtractorImpl(ImageSourceSimple & imageSource, FeatureDetector * featureDetector) : imageSource_(imageSource), featureDetector_(featureDetector) {};
    ~FeatureExtractorImpl();

    //! Returns set of features (descriptors and locations) and manages cacheing.
    const DescriptorSet * getFeatures(size_t imageId);
};

}

#endif //_FEATURE_EXTRACTOR_
