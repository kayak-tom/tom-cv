#pragma once
#ifndef _TRANSFORM_ENGINE
#define _TRANSFORM_ENGINE

#include "TransformEstimator.h"
#include "TransformSet.h"
#include "boost/thread.hpp"
#include <boost/interprocess/sync/interprocess_semaphore.hpp>
#include <set>
#include "Enums.h"

//! Namespace enclosing GRC code
namespace grc {

//! Drives the application: chooses sets of images to render, constructs and refines a set of transforms between them, converts to transforms into a mosaic, and passes to renderer
class TransformEngine
{
    TransformEstimator2 & transformEstimator_; //!< Called to estimate pairwise transforms between images
    size_t latestFrameId_; //!< Keep track of the latest frame--updated by ImageSource as frames arrive
    static const size_t NO_FRAME = 2147483647; //!< For initialising frame ids that are not yet known
    bool finished_; //! We've been notified that video has finished
    static const bool renderAllFrames_ = true; //! Wait for renderer to render every frame
    
    const int maxFrames_; //!< Max number of frames to render per mosiac
    const bool doFastRender_; //!< Controls whether incremental rendering is used
    const int maxLMIterations_; //!< Maximum number of Levenberg-Marquardt iterations
    const CvSize mosaicSize_; //!< Dimensions of mosaic (for positioning images)
    CvSize imageSize_;//!< Image size (set dynamically by image source before the first frame has arrived)
    const Enums::ePositionMethod positionMethod_;//!< Controls method used to position mosaic
    const Enums::eFrameChoiceMethod frameChoiceMethod_;//!< Controls whether to render all recent frames, or to selectively render a subset in order to cover a larger area.
    TransformSet * latestTS_;//!< The complete computed transform set that will (hopefully) be passed to the renderer.
	const int fullFrameUpdateFreq_; //!< Frequency at which to update all frames in incremental mosaicing
	const int maxSearchForTransform_; //!< Maximum number of frames to search when trying to link in a new frame

    boost::mutex notifyLock_; //!< Lock while being notified that new images have arrived
    boost::mutex giveAwayTransformLock_; //!< Lock while creating new/giving latest transform to renderer
    boost::interprocess::interprocess_semaphore semWaitForFrame_; //!< Wait for notification of a new frame
    boost::interprocess::interprocess_semaphore semWaitForTransformSet_; //!< Make renderer block until transforms are ready
    boost::interprocess::interprocess_semaphore semWaitForRenderer_; //!< Block until renderer has grabbed next frame (used if renderAllFrames_)
    boost::interprocess::interprocess_semaphore semRendererFinished_; //!< Once image source is finished, wait until renderer is finished

    boost::thread * mainEngineThread_; //!< Main engine thread--runs mainTransformLoop

    typedef set2_NF<int, std::greater<int> > TIntSetDesc;

    //typedef CFixedArray<int, MAX_BF> TCloseIdSet;

    typedef map2_NF<int, TCloseIdSet > TCloseIdMap;
    TCloseIdMap aBestBoWMatches;

    //! Main loop computing a transform every time an image arrives, and passing it to renderer
    void mainTransformLoop();

    //! Set latest transform to be passed to renderer
    void setCurrentTransform(TransformSet * TS);

    //! Calculate transform set (image+transformations) to be rendered, including images up to latestImageId.
    TransformSet * computeTransforms(int latestImageId);

    //! Delete each transformation x->y with y != latestImageId
    void cleanupTransformSet(TransformSet * TS, int alignToImageId);

    //! For each transformation x->pivotId add x->latestImageId and recurse with x as pivot
    void alignTransformToLastImage(TransformSet * TS, int pivotId, int alignToImageId, bool bFirstTime);

    //! Translate all transforms so latest image is visible and as many previous images as possible are visible
    /*! Cull images that may be out of view */
    void positionMosaic(TransformSet * TS, int alignToId, int latestImageId);

    //! Get 4 corners of source image that can be transformed to give position of source in dest image.
    //CvMat * getSourceCorners() const;

    //! Apply transform to 4 corners of source image. If any are further away than current max bounds then adjust max bounds
    void findDistantCorners(const Transform * pTrans, double & globalMaxX, double & globalMaxY, double & globalMinX, double & globalMinY) const;

    //!Calculate shift so that shift+localMin and shift+localMax fall within mosaic, and hopefully global ones do too
    static double position1d(double localMin, double localMax, double globalMin, double globalMax, double width);

    /*! Convert all transformations to align to one image */
    /*! First add all transformations x->alignToImageId for all images
     *  Then delete all transformations not pointing to alignToImageId
     *  alignToImageId is normally the latest image
     */
    void alignTransformsToOneImage(TransformSet * TS, int alignToImageId);

    //! Use Levenberg-Marquardt to refine transform set
    void refineTranformSet(TransformSet * TS, int alignTo, int latestId);

    //! Choose pairs of images that we're likely to get pairs of transforms between, and are likely to fall within the rendered window
    void chooseTransforms(TransformSet * TS, int latestImageId, std::set<size_t> & permittedAlignToIds);

    //!Choose an id to align to that will be linked to the latest frame, and is fairly central in the sequence of frames that will be rendered
    int getAlignToId(int latestImageId, const std::set<size_t> & permittedAlignToIds);
    
    //!Choose an id to align to that will be linked to the latest frame and will (normally) persist for several frames, allowing incremental rendering
    int getFastAlignToId(int latestImageId, const std::set<size_t> & permittedAlignToIds);

    //! Remove any transforms that are totally disconnected from the aligned-to frame (we won't be able to draw them)
    void removeDisconnectedTransforms(TransformSet * TS, int alignToImageId);

    //! Get a local copy of the transform between 2 ids (0 if can't find one)
    TransformInfo2 * tryGetTransformCopy(int id, int lastId);

    //! Setup global mosaic transform for fast rendering in transform set
    void setFastRender(TransformSet * TS, int alignToId, int latestId) const;

    //! Discard transforms that are heavily overlapping neighbours
    /*!These should be transforms to alignToId for all images, but we haven't refined yet so we want to maintain links between transforms */
    void thinTransforms(TransformSet * TS, int alignToId, int latestId) const;

    //! Called by the main engine thread when it has finished.
    void notifyEngineFinished();

    void addAllIds(TransformSet * TS, std::set<size_t> & permittedAlignToIds, const int latestImageId, TIntSetDesc & candidateAlignToIds, const int DEPTH);
public:
    //! Will spawn a thread that generates a new transformation set every time we are notified that a new image is available
    TransformEngine(TransformEstimator2 & transformEstimator, int maxFrames, bool useIncrementalRendering, int LMIterations, CvSize mosaicSize, Enums::eFrameChoiceMethod frameChoiceMethod, int fullFrameUpdateFreq, int maxSearchForTransform);

    ~TransformEngine();

    //!Called when new images are available to render
    void notifyImageAvailable(size_t id);
    
    //!Called (by image source) to end everything. Then blocks until the engine has ensured everything else (the renderer) is finished.
    void notifyFinished();

    //! Called by the renderer after it has finished.
    void notifyRendererFinished();
    
    //! Set image size once a frame has arrived--we don't know image sizes until capture has started
    bool notifySize(CvSize s); 

    //! Returns the latest set of image+transformations.
    /*! The rendering engine calls this. It will block if we haven't 
     *  finished computing a new one. Will always return the latest transform: it 
     *  will discard a transformation if a new one is calculated before the previous one has been geto by the renderer.
     *  this class has generated a new transform before the last one has been rendered.
     */
    grc::TransformSet * getLatestTransforms(); 
};

}

#endif // _TRANSFORM_ENGINE
