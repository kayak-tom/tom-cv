#pragma once
#ifndef _GRC_RENDERER_
#define _GRC_RENDERER_

#include "ImageSource.h"
#include "TransformEngine.h"
#include "EvaluationFunction.h"
#include "boost/thread.hpp"

namespace grc {

/*!
 * \brief Reads the latest TransformSet from the TransformEngine and renders an image from it.
 *
 * The Renderer is responsible for creating a mosaic from a set of images and Transforms that align 
 * them to a common co-ordinate frame. The basic Renderer implemented here simply overlays the 
 * images on one another, but sub-classes of Renderer perform more complicated compositions.
 */
class Renderer
{
protected:
    ImageSource & imageSource_; //!< Source of images for rendering
    TransformEngine & transformEngine_; //!< Source of transforms for rendering
    IplImage * imToRender_; //!< Storage for rendered mosaic
	EvaluationFunction *evaluationFunction_; //!< Function to evaluate the rendered mosaic

    CvSize mosaicSize_; //!< Size of the mosaic to render

	IplImage *mosaic_; //!< Internal workspace for rendering the mosaic
	bool mosaicAllocated_; //!< Flags if storage has been allocated or not

    boost::thread * renderingThread_; //!< Thread in which to conduct the rendering

	/*!
	 * \brief Reads a TransformSet from transformEngine_, renders it, and passes it off to imageSource_ for display.
	 * 
	 * Extracts TransformSets from transformEngine_ and passes them off to renderImagesToMosaic in order to generate
	 * mosaics. The resulting mosaics are passed to evaluationFunction_ (if one is provided), and then back
	 * to imageSource_ for display. This process continues until there are no more TransformSets available to 
	 * render. Typically, doRendering is launched as a thread by the constructor.
	 *
	 */
    void doRendering();
    
	/*!
	 * \brief Creates a mosaic from a set of images.
	 *
	 * This method creates a mosaic from a set of images. The TransformSet parameter specifies which images to
	 * render and the transforms between them. The implementation in Renderer simply overlays each image on to
	 * the mosaic, but this function may be over-ridden in order to create more complex mosaic Renderers.
	 *
	 * \param TS a pointer to a TransformSet to be rendered
	 * \return The rendered mosaic
	 */
    virtual IplImage * renderImagesToMosaic(TransformSet * TS);

    boost::mutex renderLock_; //!< Lock to prevent multiple render calls conflicting
    
    static const int MAX_FRAMES_AHEAD = 3; //!< Initial value for semWaitForFrame_
    boost::interprocess::interprocess_semaphore semWaitForFrame_; //!< Semaphore to ensure that the process waits for engine to have rendered previous frames, and block it until a new one is ready

	//static const char *MOSAIC_DIR; //!< Directory to save rendered mosaics to

public:

	/*!
	 * \brief Constructor, initialises the Renderer and creates a thread for it to execute in
	 *
	 * The constructor for a Renderer also sets doRendering in action within renderingThread_.
	 * 
	 * \param imageSource The source of images for the renderer
	 * \param transformEngine The source of TransformSets to be rendered
	 * \param mosaicSize The size of the mosaic to render
	 * \param evaluationFunction A function with which to evaluate the resulting mosaics. If NULL, no evaluation is made
	 */
    Renderer(ImageSource & imageSource, TransformEngine & transformEngine, CvSize mosaicSize, EvaluationFunction *evaluationFunction = NULL);

	//! Destructor, waits for the renderingThread_ to terminate then releases resources
    ~Renderer();

    /*!
	 * \brief Get the next mosaic to display.
	 *
	 * Renders the new mosaic, if available. Depending on the value of block this may wait until
	 * there is a new mosaic to render, or may return a null image if no new frames are available.
	 *
	 * \param block Flag determining whether or not to wait for a new frame.
	 * \return The next rendered image, or NULL if not blocking and there are no new frames.
	 */
	IplImage * getMosaicToDisplay(bool block); //If block may wait for rendering to finish (this blocks image display/capture thread)


}; // class Renderer

} // namespace grc

#endif //_GRC_RENDERER_
