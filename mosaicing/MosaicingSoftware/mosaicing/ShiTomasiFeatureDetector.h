#pragma once
#ifndef _GRC_SHI_TOMASI_FEATURE_DETECTOR_
#define _GRC_SHI_TOMASI_FEATURE_DETECTOR_

#include "FeatureDetector.h"
#include "boost/thread.hpp"

namespace grc {
	
	class ShiTomasiFeatureDetector : public FeatureDetector {
	public:
		ShiTomasiFeatureDetector(int maxCorners, bool subPix);
		~ShiTomasiFeatureDetector();

		FeatureVector * findFeatures(IplImage *pGreyImage, int margin);

	private:
		bool isInitialised_;
		IplImage *eigImage_;
		IplImage *tempImage_;
		int maxCorners_;
		double qualityLevel_;
		double minDistance_;
		CvPoint2D32f *cornerStore_;
        boost::mutex findFeatures_;
		IplImage *mask_;
        const bool subPix_; //Detect features with subpixel accuracy. Makes a small difference.
	};
}

#endif
