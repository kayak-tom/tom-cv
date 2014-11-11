#include "util/exception.h"
pragma_warning(push) // Disable warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)
pragma_warning(disable:4996)

#include "FeatheredRenderer.h"
#include "myGetSet2D.h"

namespace grc {

	FeatheredRenderer::~FeatheredRenderer() {
		if (initialised_) {
			cvReleaseImage(&mosaic_);
			cvReleaseImage(&warped_);
			cvReleaseImage(&mask_);
			cvReleaseImage(&warpedMask_);
		}
	}

	IplImage * FeatheredRenderer::renderImagesToMosaic(grc::TransformSet *TS) {
		if (!TS)
			throw new GRCException("FeatheredRenderer::renderImagesToMosaic: Rendering 0-transform -- should have died");
		
		if (!initialised_) {
			// Allocate storage
			const IplImage *frame = imageSource_.getImage(TS->begin()->first.im1Id());
			mosaic_ = cvCreateImage(mosaicSize_, frame->depth, frame->nChannels);
			cvSetZero(mosaic_);
			warped_ = cvCreateImage(mosaicSize_, frame->depth, frame->nChannels);
			mask_ = cvCreateImage(cvGetSize(frame), IPL_DEPTH_8U, 1);
			warpedMask_   = cvCreateImage(mosaicSize_, IPL_DEPTH_8U, 1);
			initialised_ = true;

			// Initialise the feathering mask
			if ( (featherRadius_ > mask_->width/2) || (featherRadius_ > mask_->height/2) ) {
				throw new GRCException("FeatheredRenderer::renderImagesToMosaic: Feather radius greater than half the image size");	
			}

			cvSet(mask_, cvScalarAll(255));

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
		}

		IplImage *result = cvCreateImage(mosaicSize_, mosaic_->depth, mosaic_->nChannels);
		cvSetZero(result);

		// There's a lot of code duplication between the if and else below, could separate this out
		if (TS->onlyRenderLatestTransform()) {
			// Update the mosaic
			TS->mosaicTransform()->applyToImage(mosaic_, result);
			cvSetZero(warped_);
			// Warp the latest frame
			TS->latestTransform()->transform()->applyToImage(imageSource_.getImage(TS->latestTransform()->id1()), warped_);
			cvSetZero(warpedMask_);
			TS->latestTransform()->transform()->applyToImage(mask_, warpedMask_);	
			// Feather it with the previous mosaic
			for (int y = 0; y < result->height; y++) {
				for (int x = 0; x < result->width; x++) {
					for (int c = 0; c < result->nChannels; c++) {
						int v = myGet2D(warpedMask_, x, y);
						if (v > 0) {
								if (myGet2D(result, x, y, c) > 0) {
									v = v*myGet2D(warped_, x, y, c) + (255-v)*myGet2D(result, x, y, c);
									v = v / 255;
									if (v <1 ) v = 1;
									mySet2D(result, x, y, c, v);
								} else {
									mySet2D(result, x, y, c, myGet2D(warped_, x, y, c));
								}
						}
					}
				}
			}
		} else {
			// For each frame in the mosaic
			for (TransformSet::iterator iter = TS->begin(); iter != TS->end(); iter++) {
				// Warp this frame
				cvSetZero(warped_);
				iter->second->transform()->applyToImage(imageSource_.getImage(iter->first.im1Id()), warped_);
				cvSetZero(warpedMask_);
				iter->second->transform()->applyToImage(mask_, warpedMask_);
				// Feather it into the mosaic
				for (int y = 0; y < result->height; y++) {
					for (int x = 0; x < result->width; x++) {
						for (int c = 0; c < result->nChannels; c++) {
							int v = myGet2D(warpedMask_, x, y);
							if (v > 0) {
								if (myGet2D(result, x, y, c) > 0) {
									v = v*myGet2D(warped_, x, y, c) + (255-v)*myGet2D(result, x, y, c);
									v = v / 255;
									if (v <1 ) v = 1;
									mySet2D(result, x, y, c, v);
								} else {
									mySet2D(result, x, y, c, myGet2D(warped_, x, y, c));
								}
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
