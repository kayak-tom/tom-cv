#include "util/exception.h"
//! Implementation: Implements feature extractors and returns sets of extracted features (and their locations)
pragma_warning(push)
pragma_warning(disable:4996)  // warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)

#include "FeatureExtractor.h"
#include "GRCException.h"

namespace grc {


FeatureExtractorImpl::~FeatureExtractorImpl()
{
    boost::mutex::scoped_lock featureCacheLock(featureCacheLock_);

    for(TCachedFeatureSets::iterator ppFS = cachedDescriptorSets_.begin(); ppFS != cachedDescriptorSets_.end(); ppFS++)
        delete ppFS->second;
}

//! Returns set of features (descriptors and locations), caches those that are extracted
/*!
 */
const DescriptorSet * FeatureExtractorImpl::getFeatures(size_t imageId)
{
    boost::mutex::scoped_lock featureCacheLock(featureCacheLock_);

    if(cachedDescriptorSets_.find(imageId) == cachedDescriptorSets_.end())
    {
        cachedDescriptorSets_[imageId] = getFeaturesFromImage(imageId);
    }

    return cachedDescriptorSets_[imageId];
}

//! Returns set of features (descriptors and locations)
/*! Gets image from imageSource_, extracts corners, descriptors, returns feature set.
 * TODO: Shorten function
 */
const DescriptorSet * FeatureExtractorImpl::getFeaturesFromImage(size_t imageId)
{
    const IplImage * pRGBImage = imageSource_.getImage(imageId);
    if(!pRGBImage || pRGBImage->nChannels < 3) throw new GRCException("FeatureExtractor::getFeatures: No RGB image returned from image source");

    //Todo: reuse one grey image: (Threadsafe)
    IplImage * pGreyImage = cvCreateImage(cvGetSize(pRGBImage), pRGBImage->depth, 1);

    cvCvtColor(pRGBImage, pGreyImage, pRGBImage->nChannels == 4 ? CV_BGRA2GRAY : CV_BGR2GRAY);

    if(!pGreyImage || pGreyImage->nChannels != 1 /*|| !CV_ARE_SIZES_EQ(pGreyImage, pRGBImage)*/)
        throw new GRCException("FeatureExtractor::getFeatures: Grey image returned from image source is the wrong size or depth");

    FeatureVector * cornerVector = featureDetector_->findFeatures(pGreyImage, cDescriptor::getMargin());

    if(!cornerVector) throw new GRCException("FeatureExtractor::getFeatures: No corner vector");

	// This used to raise an exception, but changed to a warning since it was causing problems with DTA's camera
    //if(!cornerVector->size()) throw new GRCException("FeatureExtractor::getFeatures: 0 corners returned");
	if (!cornerVector->size()) {
		cout << "WARNING: FeatureExtractor::getFeatures: 0 corners returned" << endl;
	}

    //Create new empty descriptor set
    DescriptorSet * descriptorSet = new DescriptorSet;

    //Compute a descriptor for each corner and add to descriptor set
    for(FeatureVector::iterator pCorner = cornerVector->begin(); pCorner < cornerVector->end(); pCorner++)
    {
        descriptorSet->push_back(cDescriptor::newDescriptor(pRGBImage, *pCorner));
    }

    if(cornerVector->size() != descriptorSet->size()) throw new GRCException("FeatureExtractor::getFeatures: Number of corners does no match number of features");

    cvReleaseImage(&pGreyImage);
    delete cornerVector;

    return descriptorSet;
}

}


pragma_warning(pop) // warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)
