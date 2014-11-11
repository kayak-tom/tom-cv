#include "util/exception.h"
//! Implements transform engine: Most logic is done here. Requests and organises sets of image+transformation (from 0) pairs that the rendering engine will draw the latest of.
pragma_warning(push)
pragma_warning(disable:4996)  // warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)

#include "FeatureMatcher.h"
#include "GRCException.h"
#include "SpeedTest.h"
#include "ImageSource.h"

namespace grc {

FeatureMatcher::~FeatureMatcher()
{
	for(CachedCorrespondences::iterator pCorresp = correspondenceCache_.begin(); pCorresp != correspondenceCache_.end(); pCorresp++)
		delete pCorresp->second;
};
FeatureMatcher2::~FeatureMatcher2()
{
	for(CachedCorrespondences::iterator pCorresp = correspondenceCache_.begin(); pCorresp != correspondenceCache_.end(); pCorresp++)
		delete pCorresp->second;
};

//! FeatureMatcher Cacheing. Maybe we don't need to cache correspondences--the transform stuff will do this for us and they will only be needed once
/*! Will return 0 if no/insufficient (<4) correspondences are found. Need 4 for all kinds of transform I think. */
const CorrespondenceSet * FeatureMatcher::getCorrespondenceSet(size_t imageId1, size_t imageId2)
{
    boost::mutex::scoped_lock lock(lockMatcherCache_);//lock while matching features: only one thread can get corresps. at a time at the moment (likewise the feature xtor so would need to change this too to get much advantage)

    const int MIN_CORRESPONDENCES = 4;

    IdPair idPair(imageId1, imageId2);

    if(correspondenceCache_.find(idPair) == correspondenceCache_.end())
    {
        //No cached correspondences
	    const DescriptorSet *features1 = featureExtractor_.getFeatures(imageId1);
	    const DescriptorSet *features2 = featureExtractor_.getFeatures(imageId2);

        CorrespondenceSet * pCorrs = 0;

        if((int)features1->size() < MIN_CORRESPONDENCES || (int)features2->size() < MIN_CORRESPONDENCES)
        {
            cout << "Failed to enough features in " << imageId1 << " or " << imageId2 << endl;
        }
        else
        {
            // check there's enough to bother with, return 0 and don't recompute if not
            //CStopWatch s;
            //s.startTimer();
            pCorrs = getCorrespondences(features1, features2);
            //s.stopTimer();
            //double numMatches = features1->size() * features2->size();
            //std::cout << "Match took " << s.getElapsedTime()/numMatches << " secs\n";
            //cout << "Found " << pCorrs->size() << " correspondences between " << imageId1 << " and " << imageId2 << endl;

            if((int)pCorrs->size() < MIN_CORRESPONDENCES)
            {
                cout << "Failed to enough correspondences between " << imageId1 << " and " << imageId2 << endl;
                delete pCorrs; pCorrs=0;
            }

        }
        correspondenceCache_[idPair] = pCorrs; // Now it exists in map we won't attempt to compute them again
    }

    return correspondenceCache_[idPair]; //May be 0 if we couldn't find enough
}

const CBoWCorrespondences * FeatureMatcher2::getCorrespondenceSet_int(IdPair & idPair)
{
	const int imageId1 = idPair.im1Id();
	const int imageId2 = idPair.im2Id();

	//No cached correspondences
	const IplImage * pIm1 = imSource.getImage(imageId1);
	const IplImage * pIm2 = imSource.getImage(imageId2);
	const CDescriptorSet *features1 = featureExtractor_.getDescriptors(pIm1);
	const CDescriptorSet *features2 = featureExtractor_.getDescriptors(pIm2);

	const CBoWCorrespondences * pCorrs = 0;

	if(features1->Count() < MIN_CORRESPONDENCES || features2->Count() < MIN_CORRESPONDENCES)
	{
		cout << "Failed to enough features in " << imageId1 << " or " << imageId2 << endl;
	}
	else
	{
		// check there's enough to bother with, return 0 and don't recompute if not
		//CStopWatch s;
		//s.startTimer();
		pCorrs = getCorrespondences(features1, features2);
		//s.stopTimer();
		//double numMatches = features1->size() * features2->size();
		//std::cout << "Match took " << s.getElapsedTime()/numMatches << " secs\n";
		//cout << "Found " << pCorrs->size() << " correspondences between " << imageId1 << " and " << imageId2 << endl;

		if(pCorrs->size() < MIN_CORRESPONDENCES)
		{
			cout << "Failed to enough correspondences between " << imageId1 << " and " << imageId2 << endl;
			delete pCorrs; pCorrs=0;
		}
	}
	return pCorrs;
}

const CBoWCorrespondences * FeatureMatcher2::getCorrespondenceSet(size_t imageId1, size_t imageId2)
{
    boost::mutex::scoped_lock lock(lockMatcherCache_);//lock while matching features: only one thread can get corresps. at a time at the moment (likewise the feature xtor so would need to change this too to get much advantage)

    IdPair idPair(imageId1, imageId2);

    if(correspondenceCache_.find(idPair) == correspondenceCache_.end())
    {
    	const CBoWCorrespondences * pCorrs = getCorrespondenceSet_int(idPair);

        correspondenceCache_[idPair] = pCorrs; // Now it exists in map we won't attempt to compute them again

        return pCorrs;
    }

    return correspondenceCache_[idPair]; //May be 0 if we couldn't find enough
}

}
pragma_warning(pop) // warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)
