#include "util/exception.h"
//!
pragma_warning(push)
pragma_warning(disable:4996)  // warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)

#include "SimpleAsymFeatureMatcher.h"
#include "GRCException.h"
#include "cvGeom.h"
//#include "Convert.h"

namespace grc {

//!
CorrespondenceSet * SimpleAsymFeatureMatcher::getCorrespondences(const DescriptorSet * features1, const DescriptorSet * features2) const
{
	CorrespondenceSet * corresp = new CorrespondenceSet();

	for (DescriptorSet::const_iterator iter1 = features1->begin(); iter1 != features1->end(); iter1++) {

	    cDescriptor::TDist bestDistance = MAX_INT;
	    cDescriptor::TDist nextDistance = MAX_INT;
        CvPoint2D32f bestMatch2_32;

        for(DescriptorSet::const_iterator iter2 = features2->begin(); iter2 != features2->end(); iter2++)
        {
			cDescriptor::TDist thisDistance = (*iter1)->distance(*iter2);
			if (thisDistance < bestDistance) {
				nextDistance = bestDistance;
				bestDistance = thisDistance;
				bestMatch2_32 = (*iter2)->location();
			}
		}

		if (/*bestDistance < 0.7*nextDistance -- equiv. but faster:*/ 10*bestDistance < 7*nextDistance)
        {
			CvPoint2D64f bestMatch2_64 = cvPoint2D64f(bestMatch2_32.x, bestMatch2_32.y);
            if(!corresp->point2Exists(bestMatch2_64))//Check best match isn't already matched to a feature in the first image
            {
                CvPoint2D32f bestMatch1_32 = (*iter1)->location();
			    corresp->insertCorresp(cvPoint2D64f(bestMatch1_32.x, bestMatch1_32.y), bestMatch2_64);
            }
		}
	}

    return corresp;
}


}
pragma_warning(pop) // warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)
