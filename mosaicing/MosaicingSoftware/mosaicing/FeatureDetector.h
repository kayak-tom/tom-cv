#pragma once
#ifndef _GRC_FEATURE_DETECTOR_
#define _GRC_FEATURE_DETECTOR_

#include "util/opencv.h"
#include <vector>

namespace grc {

    typedef std::vector<CvPoint2D32f> FeatureVector;

    //!Parent of all salient feature (corner/blob) detection classes.
	class FeatureDetector {
	public:
        //! Find some salient features in an image, no closer than margin to the edge.
		virtual FeatureVector * findFeatures(IplImage * pGreyImage, int margin) = 0;
	};	
}
#endif