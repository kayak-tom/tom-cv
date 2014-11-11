#ifndef _GRC_FEATHERED_RENDERER_2_
#define _GRC_FEATHERED_RENDERER_2_

#include "Renderer.h"

namespace grc {

	class FeatheredRenderer2 : public Renderer {

		IplImage * renderImagesToMosaic(TransformSet *TS);

	public:

		FeatheredRenderer2(ImageSource & imageSource, TransformEngine & transformEngine, CvSize mosaicSize) :
		  Renderer(imageSource, transformEngine, mosaicSize) {};

	};

}

#endif
