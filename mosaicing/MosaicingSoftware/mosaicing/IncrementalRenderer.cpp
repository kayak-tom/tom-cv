#include "IncrementalRenderer.h"

namespace grc {

	IplImage * IncrementalRenderer::renderImagesToMosaic(TransformSet *TS) {

		IplImage *result = cvCreateImage(mosaicSize_, IPL_DEPTH_8U, 3);
		cvSetZero(result);
		if (TS->onlyRenderLatestTransform()) {
			cout << "RENDERING LATEST ONLY" << endl;
			TS->mosaicTransform()->applyToImage(mosaic_, result);
			TS->latestTransform()->transform()->applyToImage(imageSource_.getImage(TS->latestTransform()->id1()), result);
		} else {
			cout << "RENDERING ALL" << endl;
			for (TransformSet::iterator iter = TS->begin(); iter != TS->end(); iter++) {
				iter->second->transform()->applyToImage(imageSource_.getImage(iter->first.im1Id()), result);
			}
		}

		cvCopy(result, mosaic_);
		return result;
	}

}
