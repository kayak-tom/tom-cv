#include "util/exception.h"
pragma_warning(push) // Disable warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)
pragma_warning(disable:4996)  

#include "DoGBlobDetector.h"
#include "GRCException.h"

#include "util/opencv_highgui.h"

namespace grc {

    DoGBlobDetector::DoGBlobDetector(int maxCorners, int levels) 
        : levels_(levels + 2) /* the +2 is so we have a 'top' and 'bottom' blurred image--and levels with one on either side*/,
          isInitialised_(false), blurredImageArray_(new IplImage * [levels_]), maxCorners_(maxCorners)
    {
        if(maxCorners<=0 || levels_ < 3) throw new GRCException("DoGBlobDetector::DoGBlobDetector: Bad parameter");
	}

	DoGBlobDetector::~DoGBlobDetector() {
		if (isInitialised_) {
            for(int i=0; i<levels_; i++)
			    cvReleaseImage(blurredImageArray_ + i);
		}

        delete [] blurredImageArray_;
	}

	FeatureVector * DoGBlobDetector::findFeatures(IplImage * pGreyImage, int margin) {
		
        if(!pGreyImage || margin < 0) throw new GRCException("DoGBlobDetector::findFeatures: Bad parameters");

        boost::mutex::scoped_lock findFeatures(findFeatures_); //Lock here: 1) So only initialise once, 2) because eigImage_ etc. are shared

	    if (!isInitialised_) { 
            for(int i=0; i<levels_; i++)
            {
                blurredImageArray_[i] = cvCreateImage(cvGetSize(pGreyImage), pGreyImage->depth, 1);

            }

            isInitialised_ = true;
		}

        //if(!CV_ARE_SIZES_EQ(pGreyImage, blurredImageArray_[0]) ) throw new GRCException("DoGBlobDetector::findFeatures: Image size has changed");
		
		//int cornerCount = maxCorners_;
        
        double sigma = 1.0;
        int kernelSize = 4;
        for(int i=0; i<levels_; i++)
        {
            cvSmooth(pGreyImage, blurredImageArray_[i], CV_GAUSSIAN, kernelSize+1, kernelSize+1, sigma);
            
            kernelSize = (kernelSize-1)*2;
            sigma *= 2.0;
        }

        //For every pixel, see if it is max (or min) of 8 neighbours, and 2*9 neighbours in other images
        FeatureVector * featureVector = new FeatureVector;

        for(int i=levels_-2; i>0; i--)
        {
            for(int x=margin; x < pGreyImage->width-margin; x++)
            {        
                for(int y=margin; y < pGreyImage->height-margin; y++)
                {   
                    if(localMaximum<std::less<double> >(i, x, y))
                    {
                        featureVector->push_back(cvPoint2D32f(x, y));
                        //cvCircle(pGreyImage, cvPoint(x,y), 3, CV_RGB(0,0,0));
                    }

                    else if(localMaximum<std::greater<double> >(i, x, y))
                    {
                        featureVector->push_back(cvPoint2D32f(x, y));
                        //cvCircle(pGreyImage, cvPoint(x,y), 3, CV_RGB(255,255,255));
                    }
                }
            }
            
            //Stop descending when we have enough corners:
            if(maxCorners_ < (int)featureVector->size()) break;
        }
        //cvSaveImage("Blobs.bmp", pGreyImage);

        return featureVector;
	}
}

pragma_warning(pop) // Disable warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)
