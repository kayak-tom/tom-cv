#include "util/exception.h"
pragma_warning(push) // Disable warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)
pragma_warning(disable:4996)

#include "MultiScaleFeatheredRenderer.h"
#include "myGetSet2D.h"

//#ifndef USE_OLD_OPENCV
//#include "opencv2/imgproc/imgproc_c.h"
//#endif


namespace grc {

	MultiScaleFeatheredRenderer::~MultiScaleFeatheredRenderer() {
		if (initialised_) {
			cvReleaseImage(&mosaic_);
			cvReleaseImage(&mask_);
			cvReleaseImage(&warpedMask_);
			cvReleaseImage(&totalMask_);
			cvReleaseImage(&tempMask_);
			cvReleaseImage(&warpedHi_);
			cvReleaseImage(&warpedLo_);
			cvReleaseImage(&frameLo_);
			cvReleaseImage(&mosaicHi_);
			cvReleaseImage(&mosaicLo_);
		}
	}

	IplImage * MultiScaleFeatheredRenderer::renderImagesToMosaic(grc::TransformSet *TS) {
		if (!TS)
			throw new GRCException("MultiScaleFeatheredRenderer::renderImagesToMosaic: Rendering 0-transform -- should have died");

		if (!initialised_) {
			// Allocate storage
			const IplImage *frame = imageSource_.getImage(TS->begin()->first.im1Id());
			mosaic_ = cvCreateImage(mosaicSize_, frame->depth, frame->nChannels);
			cvSetZero(mosaic_);

			mask_ = cvCreateImage(cvGetSize(frame), IPL_DEPTH_8U, 1);
			warpedMask_ = cvCreateImage(mosaicSize_, IPL_DEPTH_8U, 1);
			totalMask_ = cvCreateImage(mosaicSize_, IPL_DEPTH_8U, 1);
			tempMask_  = cvCreateImage(mosaicSize_, IPL_DEPTH_8U, 1);

			warpedHi_ = cvCreateImage(mosaicSize_, mosaic_->depth, mosaic_->nChannels);
			warpedLo_ = cvCreateImage(mosaicSize_, mosaic_->depth, mosaic_->nChannels);

			frameLo_ = cvCreateImage(cvGetSize(frame), frame->depth, frame->nChannels);

			mosaicHi_ = cvCreateImage(mosaicSize_, mosaic_->depth, mosaic_->nChannels);
			mosaicLo_ = cvCreateImage(mosaicSize_, mosaic_->depth, mosaic_->nChannels);
			
			if ( (featherRadius_ > mask_->width/2) || (featherRadius_ > mask_->height/2) ) {
				throw new GRCException("FeatheredRenderer::renderImagesToMosaic: Feather radius greater than half the image size");	
			}

			cvSet(mask_, cvScalarAll(255));

			// Set up a feathering mask
			for (int y = 0; y < mask_->height; y++) {
				for (int x = 0; x < (featherRadius_-1); x++) {
					mySet2D(mask_, x, y, (x+1)*255/featherRadius_);
					mySet2D(mask_, mask_->width-(x+1), y, (x+1)*255/featherRadius_);
				}
			}

			for (int x = 0; x < mask_->width; x++) {
				for (int y = 0; y < (featherRadius_-1); y++) {
					if ( (y < x) && (y < mask_->width-x) ) {
						mySet2D(mask_, x, y, (y+1)*255/featherRadius_);
						mySet2D(mask_, x, mask_->height-(y+1), (y+1)*255/featherRadius_);
					}
				}
			}
			initialised_ = true;
		}


		IplImage *result = cvCreateImage(mosaicSize_, mosaic_->depth, mosaic_->nChannels);
		cvSetZero(result);
		
		// There's a lot of code duplication between the if and else below. Could be moved to a separte function
		if (TS->onlyRenderLatestTransform()) {
			TS->mosaicTransform()->applyToImage(mosaicLo_, result);
			cvCopy(result, mosaicLo_);
			TS->mosaicTransform()->applyToImage(mosaic_, result);
			TS->mosaicTransform()->applyToImage(totalMask_, tempMask_);
			cvCopy(tempMask_, totalMask_);
			const IplImage *frameHi = imageSource_.getImage(TS->latestTransform()->id1());
			const Transform *T = TS->latestTransform()->transform();
			
			// Warp the image and separate it into low (smoothed) and high (original-smoothed) frequency components
			cvSetZero(warpedHi_);
			T->applyToImage(frameHi, warpedHi_);
			cvSmooth(frameHi, frameLo_, CV_BLUR, 3);
			cvSetZero(warpedLo_);
			T->applyToImage(frameLo_, warpedLo_);
			cvSetZero(warpedMask_);
			T->applyToImage(mask_, warpedMask_);
			// Feather the low-frequency components and add back the high-frequency elements
			for (int y = 0; y < result->height; y++) {
				for (int x = 0; x < result->width; x++) {
					if (myGet2D(warpedMask_, x, y) > 0) {
						if (myGet2D(totalMask_, x, y) > 0) {
							for (int c = 0; c < result->nChannels; c++) {
								int v = myGet2D(warpedMask_, x, y);
								v = v*myGet2D(warpedLo_, x, y, c) + (255-v)*myGet2D(mosaicLo_, x, y, c);
								v = v / 255;
								mySet2D(mosaicLo_, x, y, c, v);
								v = v + myGet2D(warpedHi_, x, y, c) - myGet2D(warpedLo_, x, y, c);
								if (v > 255) v = 255;
								if (v < 0) v = 0;
								mySet2D(result, x, y, c, v);
							}
						} else {
							for (int c = 0; c < result->nChannels; c++) {
								mySet2D(mosaicLo_, x, y, c, myGet2D(warpedLo_, x, y, c));
								mySet2D(result, x, y, c, myGet2D(warpedHi_, x, y, c));
							}
							mySet2D(totalMask_, x, y, 255);
						}
					}
				}
			}

		} else {
			cvSetZero(totalMask_);
			cvSetZero(mosaicLo_);
			const Transform *T;
			// For each image making up the mosaic...
			for (TransformSet::iterator iter = TS->begin(); iter != TS->end(); iter++) {
				T = iter->second->transform();
				const IplImage *frameHi = imageSource_.getImage(iter->first.im1Id());
				cvSmooth(frameHi, frameLo_, CV_BLUR, 3);

				cvSetZero(warpedLo_);
				T->applyToImage(frameLo_, warpedLo_);
				cvSetZero(warpedHi_);
				T->applyToImage(frameHi, warpedHi_);

				cvSetZero(warpedMask_);
				T->applyToImage(mask_, warpedMask_);
				for (int y = 0; y < result->height; y++) {
					for (int x = 0; x < result->width; x++) {
						if (myGet2D(warpedMask_, x, y) > 0) {
							if (myGet2D(totalMask_, x, y) > 0) {
								for (int c = 0; c < result->nChannels; c++) {
									int v = myGet2D(warpedMask_, x, y);
									v = v*myGet2D(warpedLo_, x, y, c) + (255-v)*myGet2D(mosaicLo_, x, y, c);
									v = v / 255;
									mySet2D(mosaicLo_, x, y, c, v);
									v = v + myGet2D(warpedHi_, x, y, c) - myGet2D(warpedLo_, x, y, c);
									if (v > 255) v = 255;
									if (v < 0) v = 0;
									mySet2D(result, x, y, c, v);
								}
							} else {
								for (int c = 0; c < result->nChannels; c++) {
									mySet2D(mosaicLo_, x, y, c, myGet2D(warpedLo_, x, y, c));
									mySet2D(result, x, y, c, myGet2D(warpedHi_, x, y, c));
								}
								mySet2D(totalMask_, x, y, 255);
							}
						}
					}
				}
			}
		}
		cvCopy(result, mosaic_);

		return result;

		
	}



}

pragma_warning(pop) // Disable warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)
