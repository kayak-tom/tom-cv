#include "FeatheredRenderer2.h"
#include "myGetSet2D.h"

namespace grc {

	IplImage * FeatheredRenderer2::renderImagesToMosaic(grc::TransformSet *TS) {
		if (!TS)
			throw new GRCException("FeatheredRenderer::renderImagesToMosaic: Rendering 0-transform -- should have died");

		IplImage *result = cvCreateImage(mosaicSize_, IPL_DEPTH_8U, 3);
		cvSetZero(result);

		CvSize frameSize = cvGetSize(imageSource_.getImage(TS->begin()->first.im1Id()));
		IplImage *warped = cvCreateImage(mosaicSize_, IPL_DEPTH_8U, 3);
		IplImage *weight = cvCreateImage(frameSize, IPL_DEPTH_8U, 1);
		IplImage *mask = cvCreateImage(mosaicSize_, IPL_DEPTH_8U, 1);
		cvSet(weight, cvScalarAll(255), 0);

		int r = 16;

		for (int y = 0; y < weight->height; y++) {
			for (int x = 0; x < (r-1); x++) {
				mySet2D(weight, x, y, (x+1)*255/r);
				mySet2D(weight, weight->width-(x+1), y, (x+1)*255/r);
			}
		}

		for (int x = 0; x < weight->width; x++) {
			for (int y = 0; y < (r-1); y++) {
				if ( (y < x) && (y < weight->width-x) ) {
					mySet2D(weight, x, y, (y+1)*255/r);
					mySet2D(weight, x, weight->height-(y+1), (y+1)*255/r);
				}
			}
		}

		cvSaveImage("weight.bmp", weight);

		const Transform *T;
		
		for (TransformSet::iterator iter = TS->begin(); iter != TS->end(); iter++) {
			T = iter->second->transform();
			cvSetZero(warped);
			T->applyToImage(imageSource_.getImage(iter->first.im1Id()), warped);
			cvSetZero(mask);
			T->applyToImage(weight, mask);
			for (int y = 0; y < result->height; y++) {
				for (int x = 0; x < result->width; x++) {
					for (int c = 0; c < result->nChannels; c++) {
						int v = myGet2D(mask, x, y);
						v = v*myGet2D(warped, x, y, c) + (255-v)*myGet2D(result, x, y, c);
						v = v / 255;
						mySet2D(result, x, y, c, v);
					}
				}
			}
		}

		cvReleaseImage(&weight);
		cvReleaseImage(&warped);
		cvReleaseImage(&mask);

		return result;
	}

}