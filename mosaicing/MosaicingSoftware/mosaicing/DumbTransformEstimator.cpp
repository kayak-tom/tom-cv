//pragma_warning(push)
//pragma_warning(disable:4996)  // warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)

#include "DumbTransformEstimator.h"

namespace grc {
	const Transform * DumbTransformEstimator::getTransform(size_t imageId1, size_t imageId2) {
		int diff = imageId2 - imageId1;
		Transform *result = new PerspectiveTransform();
		cvmSet(*result, 0, 2, -diff*10);
		cvmSet(*result, 1, 2, -diff*10);
		return result;
	}
}

//pragma_warning(pop) // warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)
