#include "util/exception.h"
pragma_warning(push) // Disable warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)
pragma_warning(disable:4996)  

#include "ShiTomasiFeatureDetector.h"
#include "GRCException.h"

namespace grc {

    ShiTomasiFeatureDetector::ShiTomasiFeatureDetector(int maxCorners, bool subPix) : isInitialised_(false), eigImage_(0), tempImage_(0), maxCorners_(maxCorners), subPix_(subPix) {
        
        if(maxCorners<=0) throw new GRCException("ShiTomasiFeatureDetector::ShiTomasiFeatureDetector: Bad parameter");

        cornerStore_ = new CvPoint2D32f[maxCorners_];
	}

	ShiTomasiFeatureDetector::~ShiTomasiFeatureDetector() {
		if (isInitialised_) {
			//cvReleaseImage(&greyImage_);
			cvReleaseImage(&eigImage_);
			cvReleaseImage(&tempImage_);
			delete [] cornerStore_;
		}
	}

	FeatureVector * ShiTomasiFeatureDetector::findFeatures(IplImage * pGreyImage, int margin) {
		
        if(!pGreyImage || margin < 0) throw new GRCException("ShiTomasiFeatureDetector::findFeatures: Bad parameters");

        boost::mutex::scoped_lock findFeatures(findFeatures_); //Lock here: 1) So only initialise once, 2) because eigImage_ etc. are shared

	    if (!isInitialised_) { 
		    // Allocate some storage
		    //greyImage_ = cvCreateImage(cvGetSize(pGreyImage), img->depth, 1);
		    eigImage_  = cvCreateImage(cvGetSize(pGreyImage), IPL_DEPTH_32F, 1);
		    tempImage_ = cvCreateImage(cvGetSize(pGreyImage), IPL_DEPTH_32F, 1);
		    qualityLevel_ = 0.1;
		    minDistance_ = 10;

            //Enforce margin:

			mask_ = cvCreateImage(cvGetSize(pGreyImage), IPL_DEPTH_8U, 1);
			cvSetZero(mask_);
			cvDrawRect(mask_, cvPoint(margin, margin), 
				       cvPoint(pGreyImage->width-(2*margin+1), 
					           pGreyImage->height-(2*margin+1)), 
					   CV_RGB(255,255,255), -1);

            isInitialised_ = true;
		}

        //if(!CV_ARE_SIZES_EQ(pGreyImage, eigImage_) ) throw new GRCException("ShiTomasiFeatureDetector::findFeatures: Image size has changed");
		
		int cornerCount = maxCorners_;
		cvGoodFeaturesToTrack(pGreyImage, eigImage_, tempImage_,
			cornerStore_, &cornerCount, qualityLevel_, minDistance_, mask_);

        if(subPix_)
            cvFindCornerSubPix(pGreyImage, cornerStore_, cornerCount, cvSize(3,3), cvSize(-1,-1), cvTermCriteria(CV_TERMCRIT_ITER+CV_TERMCRIT_EPS,5,.1));

        FeatureVector * featureVector = new FeatureVector(cornerCount);

		for (int c = 0; c < cornerCount; c++) {//Todo: iterator/std::copy?
			(*featureVector)[c] = (cornerStore_[c]);
		}

        return featureVector;
	}
}

pragma_warning(pop) // Disable warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)
