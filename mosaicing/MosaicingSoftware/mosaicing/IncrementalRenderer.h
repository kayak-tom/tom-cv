#ifndef _GRC_INCREMENTAL_RENDERER_
#define _GRC_INCREMENTAL_RENDERER_

#include "Renderer.h"

namespace grc {
	
	class IncrementalRenderer : public Renderer {

		IplImage * renderImagesToMosaic(TransformSet *TS);

		IplImage *mosaic_;

	public:

		IncrementalRenderer(ImageSource & imageSource, TransformEngine & transformEngine, CvSize mosaicSize) :
		  Renderer(imageSource, transformEngine, mosaicSize) {
			mosaic_ = cvCreateImage(mosaicSize_, IPL_DEPTH_8U, 3);
			cvSetZero(mosaic_);
		}

		~IncrementalRenderer() {
			cvReleaseImage(&mosaic_);
		}

	};
}

#endif