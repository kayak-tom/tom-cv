#include "util/exception.h"
pragma_warning(push) // Disable warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)
pragma_warning(disable:4996)

#include "SumSquaredErrorEvaluationFunction.h"

namespace grc {

	double SumSquaredErrorEvaluationFunction::evaluate(IplImage *mosaic, TransformSet *TS, ImageSource &imageSource) {
		if (TS->size() == 0) {
			return 0;
		}
		IplImage *warp = cvCreateImage(cvGetSize(mosaic), IPL_DEPTH_8U, mosaic->nChannels);
		IplImage *mask = cvCreateImage(cvGetSize(mosaic), IPL_DEPTH_8U, 1);
		IplImage *error = cvCreateImage(cvGetSize(mosaic), IPL_DEPTH_8U, mosaic->nChannels);

		const IplImage *frame = imageSource.getImage(TS->begin()->first.im1Id());

		IplImage *white = cvCreateImage(cvGetSize(frame), IPL_DEPTH_8U, 1);
		cvSet(white, cvScalarAll(255));

		double totalError = 0;
		double totalSize = 0;

		// Iterate over each image in the mosaic
		for (TransformSet::iterator iter = TS->begin(); iter != TS->end(); iter++) {
			// Get the image from the Transform set and warp it to the mosaic
			frame = imageSource.getImage(iter->first.im1Id());
			iter->second->transform()->applyToImage(frame, warp);
			// Warp the mask to the mosaic also
			iter->second->transform()->applyToImage(white, mask);
			// Sum up the squared differences
			CvScalar size = cvSum(mask);
			cvSubRS(mask, cvScalarAll(255), mask);
			cvAbsDiff(warp, mosaic, error);
			cvSet(error, cvScalarAll(0), mask);
			CvScalar sum = cvSum(error);
			// Accumulate the result
			totalSize += size.val[0] * frame->nChannels;
			for (int c = 0; c < frame->nChannels; c++) {
				totalError += sum.val[0];
			}
		}

		//Release storage
		cvReleaseImage(&white);
		cvReleaseImage(&mask);
		cvReleaseImage(&warp);
		cvReleaseImage(&error);

		// Scale the result to the range [0,1]
		return 1.0 - totalError/totalSize;
		
	}
}
pragma_warning(pop) // Disable warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)

