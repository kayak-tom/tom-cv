#include "util/exception.h"
//! Implementation: Estimates transforms between pairs of images (possibly between larger sets of images in the future??)
pragma_warning(push)
pragma_warning(disable:4996)  // warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)

#include "TransformEstimator.h"
#include "GRCException.h"

namespace grc {

//! Checks whether a transform is sensible (RANSAC may have found a transform that fits bad points, either in a degenerate configuration or from a set that should be outliers)
bool TransformEstimator::transIsGood(const TransformInfo * trans)
{
    const Transform * transform = trans->transform();

    const PerspectiveTransform * PT = dynamic_cast<const PerspectiveTransform *>(transform);
    bool bSequential = (abs((int)trans->id1() - (int)trans->id2()) <= 2);

    if(PT)
    {
        const double MAX_ALLOWED_DET = bSequential ? 1.3 : 2;
        const double MIN_ALLOWED_DET = 1.0/MAX_ALLOWED_DET;
        double det = cvDet(*PT);
        if(det > MAX_ALLOWED_DET || det < MIN_ALLOWED_DET)
            return false;

		const double MAX_ALLOWED_DIAG=MAX_ALLOWED_DET, MIN_ALLOWED_DIAG=MIN_ALLOWED_DET; 
		const double MAX_ALLOWED_SKEW = 0.15;

		for(int i=0;i<3;i++)
		{
			for(int j=0;j<2;j++)
			{
				double Hij = cvmGet(*PT, i, j);
				if(i==j)
				{
					if(Hij > MAX_ALLOWED_DIAG || Hij < MIN_ALLOWED_DIAG)
						cout << "Bad diagonal: " << Hij << "\n";
						return false;
				}
				else
				{
					if(fabs(Hij) > MAX_ALLOWED_SKEW)
					{
						cout << "Too much skew: " << Hij << "\n";
						return false;
					}
				}
			}
		}

    }

    //Todo: constraint on rotation angle, skew, maybe difference from I

    return true;
}

bool TransformEstimator2::transIsGood(const TransformInfo2 * trans)
{
    const Transform * transform = trans->transform();

    const PerspectiveTransform * PT = dynamic_cast<const PerspectiveTransform *>(transform);
    bool bSequential = (abs((int)trans->id1() - (int)trans->id2()) <= 2);

    if(PT)
    {
        const double MAX_ALLOWED_DET = bSequential ? 1.5 : 2;
        const double MIN_ALLOWED_DET = 1.0/MAX_ALLOWED_DET;
        double det = cvDet(*PT);
        if(det > MAX_ALLOWED_DET || det < MIN_ALLOWED_DET)
        {
			cout << "Transform rejected: Bad scale: " << det << "\n";
            return false;
        }

        if(bSequential)
        {
			for(int j=0;j<2;j++)
			{
				double Hjj = cvmGet(*PT, j, j);
				if(Hjj<0.4)
				{
					cout << "Transform rejected: Too much rotation: " << Hjj << "\n";
					return false;
				}
			}
		}

		/* This includes rotation
		 * const double MAX_ALLOWED_DIAG=MAX_ALLOWED_DET, MIN_ALLOWED_DIAG=MIN_ALLOWED_DET;
		const double MAX_ALLOWED_SKEW = 0.15;

		for(int i=0;i<3;i++)
		{
			for(int j=0;j<2;j++)
			{
				double Hij = cvmGet(*PT, i, j);
				if(i==j)
				{
					if(Hij > MAX_ALLOWED_DIAG || Hij < MIN_ALLOWED_DIAG)
					{
						cout << "Transform rejected: Bad diagonal: " << Hij << "\n";
						return false;
					}
				}
				else
				{
					if(fabs(Hij) > MAX_ALLOWED_SKEW)
					{
						cout << "Transform rejected: Too much skew: " << Hij << "\n";
						return false;
					}
				}
			}
		}*/
		CvPoint2D64f trans = PT->translation();
		if(fabs(trans.x) > 650 || fabs(trans.y) > 500)
		{
			cout << "Transform rejected: Too much translation: " << trans.x << ',' << trans.y << "\n";
			return false;
		}

    }

    //Todo: constraint on rotation angle, skew, maybe difference from I

    return true;
}

//! Returns transform between a pair of images, 0 on failure.
/*! Todo: not corrently threadsafe
 */
const TransformInfo * TransformEstimator::getTransform(size_t imageId1, size_t imageId2)
{
    //todo lock. Lock a map of id pairs that are currently being created--then we can be running multi threads here.
    IdPair ids(imageId1, imageId2);

    if(transformCache_.find(ids) == transformCache_.end())
    {
        //No cached transformation
        const CorrespondenceSet *corresp = featureMatcher_.getCorrespondenceSet(imageId1, imageId2);

        //check we actually found a transform, make sure we don't recompute if not (by setting map to 0)
        Transform * trans = 0;
        TransformInfo * transInfo = 0;

        if(corresp)
        {
            trans = getTransform(corresp);

            if(trans)
            {
                transInfo = new TransformInfo(trans, corresp, imageId1, imageId2); //corresp copied, trans given away
                if(!transIsGood(transInfo))
                {
                    delete transInfo; transInfo=0;
                }
            }
        }

        if(!transInfo)
        {
            cout << "Failed to find a transform between " << imageId1 << " and " << imageId2 << endl;
        }
        transformCache_[ids] = transInfo;
    }

    return transformCache_[ids]; //May be 0, if we've (previously) failed to find a transform
}

TransformEstimator::~TransformEstimator()
{
    for(CachedTransformations::iterator ppTrans = transformCache_.begin(); ppTrans != transformCache_.end(); ppTrans++)
        delete ppTrans->second;

};


//! Checks whether a transform is sensible (RANSAC may have found a transform that fits bad points, either in a degenerate configuration or from a set that should be outliers)
/*bool TransformEstimator2::transIsGood(const TransformInfo2 * trans)
{
    const Transform * transform = trans->transform();

    const PerspectiveTransform * PT = dynamic_cast<const PerspectiveTransform *>(transform);
    bool bSequential = (abs((int)trans->id1() - (int)trans->id2()) <= 2);

    if(PT)
    {
        const double MAX_ALLOWED_DET = bSequential ? 1.5 : 2;
        const double MIN_ALLOWED_DET = 1.0/MAX_ALLOWED_DET;
        double det = cvDet(*PT);
        if(det > MAX_ALLOWED_DET || det < MIN_ALLOWED_DET)
        {
        	cout << "Invalid transformation, det=" << det << endl;
            return false;
        }
    }

    //Todo: constraint on rotation angle, skew, maybe difference from I

    return true;
}*/

void TransformEstimator2::getSimilarIds(const int nImageId, const int nTopN, TCloseIdSet & ids) const
{
	if(IS_DEBUG) CHECK(ids.size() != 0, "Expected empty TCloseIdSet");
	featureMatcher_.getSimilarIds(nImageId, nTopN, ids);
}


//! Returns transform between a pair of images, 0 on failure.
/*! Todo: not corrently threadsafe
 */
const TransformInfo2 * TransformEstimator2::getTransform(size_t imageId1, size_t imageId2)
{
    //todo lock. Lock a map of id pairs that are currently being created--then we can be running multi threads here.
    IdPair ids(imageId1, imageId2);

    if(transformCache_.find(ids) == transformCache_.end())
    {
        //No cached transformation
        const CBoWCorrespondences *corresp = featureMatcher_.getCorrespondenceSet(imageId1, imageId2);

        //check we actually found a transform, make sure we don't recompute if not (by setting map to 0)
        Transform * trans = 0;
        TransformInfo2 * transInfo = 0;

        if(corresp)
        {
        	Transform const * lastTrans = 0;
        	double dLastTransScale = -1;


        	/*IdPair lastIds(imageId1-1, imageId1); //Do we have one to extrap from?
        	CachedTransformations::const_iterator pLastTrans = transformCache_.find(lastIds);
        	if(pLastTrans != transformCache_.end())
        	{
        		lastTrans = pLastTrans->second->transform();
        		dLastTransScale =
        	}*/

        	if(transformCache_.size() > 0)
        	{
        		CachedTransformations::const_iterator pLastTrans = transformCache_.end();
        		pLastTrans--;
        		const IdPair & lastIds = pLastTrans->first;
        		//cout << lastIds.im1Id() << ' ' << lastIds.im2Id() << endl;
        		if(pLastTrans->second && (int)lastIds.im2Id() > ((int)imageId1)-4 && (int)lastIds.im1Id() > ((int)imageId1)-8)
        		{
        			double dOldDiff = lastIds.im2Id() - lastIds.im1Id();
        			double dNewDiff = imageId2-imageId1;
        			dLastTransScale = dNewDiff/dOldDiff; //Todo a bit dodge when lastIds.im2Id() < imageId1
        			lastTrans = pLastTrans->second->transform();
        		}

        	}


            trans = getTransform(corresp, imageId1, imageId2, lastTrans, dLastTransScale);
            if(trans)
            {
                transInfo = new TransformInfo2(trans, corresp, imageId1, imageId2); //corresp copied, trans given away
                if(!transIsGood(transInfo))
                {
                    delete transInfo; transInfo=0;
                }
            }
        }

        if(!transInfo)
        {
            cout << "Failed to find a transform between " << imageId1 << " and " << imageId2 << endl;
        }
        transformCache_[ids] = transInfo;
    }

    return transformCache_[ids]; //May be 0, if we've (previously) failed to find a transform
}

TransformEstimator2::~TransformEstimator2()
{
    for(CachedTransformations::iterator ppTrans = transformCache_.begin(); ppTrans != transformCache_.end(); ppTrans++)
        delete ppTrans->second;

};


}
pragma_warning(pop)  // warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)

