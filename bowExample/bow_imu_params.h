/*
 * bowslamParams.h
 *
 *  Created on: 19/10/2009
 *      Author: tom
 */

#ifndef BOWSLAMPARAMS_H_
#define BOWSLAMPARAMS_H_

#include "params/param.h"
#include "util/convert.h"
#include "ransac/ransacParams.h"
#include "bow/bagOfWordsParam.h"
#include "image/imageAccess.h"
#include "description/vectorDescriptor.h"
#include "description/patchParams.h"
#include "featureExtract/cornerDetector.h"

WRAPPERPARAMCLASS(BoWIMU)
	PARAM(MIN_INLIERS_NEARBY, 8, 50, 20, "Unused at the moment")
	PARAM(MIN_INLIERS_DISTANT, 8, 100, 28, "Min number of inlier correspondences needed before a model will be used")

    PARAMB(CORRECT_RD, true, "Turn on/off radial distortion correction")

    CHILDCLASS(BOW, "Bag-of-Words")
	CHILDCLASS(BOWMatching, "Params for correspondences from Bag-of-Words")
	CHILDCLASS(PatchDescriptor, "Image patch selection params")
	CHILDCLASS(RANSAC, "RANSAC for E params")
	CHILDCLASS(RANSACHomography, "RANSAC for H params")
	CHILDCLASS(Im, "Image sizes, channels, etc.")
	CHILDCLASS(Corner, "Corner detection params")
	CHILDCLASS(DescriptorSetClustering, "BoW dictionary generation params")
	{}

	CNumParam<int> MIN_INLIERS_NEARBY, MIN_INLIERS_DISTANT;
	CNumParam<bool> CORRECT_RD;

	//BoWSLAM requires these parameterisation so instantiate inside
	MAKECHILDCLASS(BOW)
	MAKECHILDCLASS(BOWMatching)
	MAKECHILDCLASS(PatchDescriptor)
	MAKECHILDCLASS(RANSAC)
	MAKECHILDCLASS(RANSACHomography)
	MAKECHILDCLASS(Im)
	MAKECHILDCLASS(Corner)
	MAKECHILDCLASS(DescriptorSetClustering)
};

#endif /* BOWSLAMPARAMS_H_ */
