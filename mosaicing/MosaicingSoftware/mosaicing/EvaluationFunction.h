#ifndef _GRC_EVALUATION_FUNCTION_
#define _GRC_EVALUATION_FUNCTION_

#include "TransformSet.h"
#include "ImageSource.h"
#include "util/opencv.h"

namespace grc {

    //! Parent class for functions evalkuating the quality of a mosaic.
	class EvaluationFunction {
	public:
        //! Evaluate the quality of a mosaic w.r.t. a transform set.
		virtual double evaluate(IplImage *mosaic, TransformSet *TS, ImageSource &imageSource) = 0;
	};

}

#endif
