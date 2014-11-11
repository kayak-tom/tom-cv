#ifndef _GRC_MULTI_SCALE_FEATHERED_RENDERER_
#define _GRC_MULTI_SCALE_FEATHERED_RENDERER_

#include "Renderer.h"

namespace grc {

	/*!
	 * \brief Reads the latest TransformSet from the TransformEngine and renders an image from it using multi-scale feathering.
	 *
	 * Renderers are responsible for creating a mosaic from a set of images and Transforms that align 
	 * them to a common co-ordinate frame. The MultiScaleFeatheredRenderer reduces the appearance of the 
	 * borders between images by blurring across image boundaries. In order to try and preserve fine detail
	 * only low-frequency image components are blurred, high-frequency components are not.
	 */	
	class MultiScaleFeatheredRenderer : public Renderer {

		/*!
		 * \brief Creates a mosaic from a set of images.
		 *
		 * This method creates a mosaic from a set of images. The TransformSet parameter specifies which images to
		 * render and the transforms between them. The MultiScaleFeatheredRenderer reduces the appearance of the 
		 * borders between images by blurring across image boundaries. In order to try and preserve fine detail
		 * only low-frequency image components are blurred, high-frequency components are not.
		 *
		 * \param TS a pointer to a TransformSet to be rendered
		 * \return The rendered mosaic
		 */
		IplImage * renderImagesToMosaic(TransformSet *TS);
		int featherRadius_;    //!< Width of border across which to smooth images
		bool initialised_;     //!< Flags whether or not internal storage has been allocated 
		IplImage * mosaic_;    //!< Internal storage for rendering mosaic
		IplImage * mask_;      //!< Internal storage for feathering mask
		IplImage * warpedMask_;//!< Internal storage for warping the feathering mask
		IplImage * totalMask_; //!< Internal storage for the overall mosaic mask
		IplImage * tempMask_;  //!< Temporary storage for in-place warping
		IplImage * warpedHi_;  //!< Internal storage for the warped high-frequency component of each frame
		IplImage * warpedLo_;  //!< Internal storage for the warped low-frequency component of each frame
		IplImage * frameLo_;   //!< Internal storage for the low-frequency component of each frame
		IplImage * mosaicHi_;  //!< Internal storage for the high-frequency component of the mosaic
		IplImage * mosaicLo_;  //!< Internal storage for the low-frequency component of the mosaic

	public:

		/*!
		 * \brief Constructor, initialises the DijkstraCutRenderer and creates a thread for it to execute in
		 *
		 * As with all Renderers the constructor also sets doRendering in action within renderingThread_.
		 * The MultiScaleFeatheredRenderer also takes a parameter which specifies how much to blur the boundaries 
		 * between images. A larger value means more blurring but a smoother transition.
		 * 
		 * \param imageSource The source of images for the renderer
		 * \param transformEngine The source of TransformSets to be rendered
		 * \param mosaicSize The size of the mosaic to render
		 * \param featherRadius Width of border across which to blur the images
		 * \param evaluationFunction A function with which to evaluate the resulting mosaics. If NULL, no evaluation is made
		 */
		MultiScaleFeatheredRenderer(ImageSource & imageSource, TransformEngine & transformEngine, CvSize mosaicSize, int featherRadius, EvaluationFunction *evaluationFunction) :
		  Renderer(imageSource, transformEngine, mosaicSize, evaluationFunction), featherRadius_(featherRadius), initialised_(false),
		  mosaic_(0), mask_(0), warpedMask_(0), totalMask_(0), warpedHi_(0), warpedLo_(0), frameLo_(0), mosaicHi_(0), mosaicLo_(0) {};

		//! Destructor, waits for the renderingThread_ to terminate then releases resources
		~MultiScaleFeatheredRenderer();

	};

}

#endif
