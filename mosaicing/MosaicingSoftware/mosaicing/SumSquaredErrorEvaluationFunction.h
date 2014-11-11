#ifndef _GRC_SUM_SQUARED_ERROR_EVALUATION_FUNCTION_
#define _GRC_SUM_SQUARED_ERROR_EVALUATION_FUNCTION_

#include "EvaluationFunction.h"

namespace grc {
	
	/*!
	 * \brief Evaluation of a mosaic by sum of squared differences between the mosaic and individual frames.
	 *
	 * This class implements an EvaluationFunction that is based around the idea that the mosaic should
	 * be locally similar to each of the input images. This is not necessarily a good indicator of the
	 * visual quality of the mosaic, however, and so should be considered with caution.
	 */
	class SumSquaredErrorEvaluationFunction : public EvaluationFunction {
	public:
		/*!
		 * \brief Function to evaluate a mosaic via sum of squared differences metric
		 * 
		 * For each component image of the mosaic, the differences between the mosaic and
		 * the image warped according to the appropriate transform is found. These differences
		 * are squared, summed, and averaged over all the images. This value is then scaled
		 * to the range [0,1].
		 *
		 * \param mosaic The mosaic to evaluate
		 * \param TS The set of transforms used to form the mosaic
		 * \param imageSource The source of the frames used to construct the mosaic
		 * \return An estimate of the quality of the mosaic from 0 (terrible) to 1 (perfect). 
		 */
		double evaluate(IplImage *mosaic, TransformSet *TS, ImageSource &imageSource);
	};

}

#endif