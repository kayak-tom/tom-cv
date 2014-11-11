#pragma once

#ifndef _IMAGE_SOURCE
#define _IMAGE_SOURCE

#include "util/opencv_highgui.h"
//#ifndef USE_OLD_OPENCV
//#include "opencv2/opencv.hpp"
//#endif

#include <vector>
#include <string>
#include "TransformEngine.h"
pragma_warning(push)
pragma_warning(disable:4996)
#include "boost/thread.hpp"
pragma_warning(pop)

#include "ImageSourceSimple.h"
#include "Enums.h"

namespace grc {

class Renderer; //!< Forward declaration avoids circular refs in header files

//! Reads images from video device/file/directory containing images. 
/*! Capturing from camera is done in doCaptureImages (blocking). Quitting is initiated from here.
 *  Recent frames are cached in memory, frames are all saved to disk in case they are needed later.*/
class ImageSource : public ImageSourceSimple
{
    TransformEngine * engine_; //!< The ImageSource notifys the engine when new frames are available
    Renderer * mosaicRenderer_; //!< Display of rendered images is done here in the main thread, to avoid OpenCV threading issues.

    //! Check if a rendered image is available for display, display it if so
    void tryDisplayRenderedMosaic(bool keepUp);

    //! Check for a new frame
    void getImage(bool & finished);

    //! Get new frame from capture (file or cam) object, return 0 when no more are available
    IplImage * tryGetImageFromCamOrAVI();
    
    //! Save and cache this frame.
    void storeImage(IplImage * pFrame);

    //! Remove frames that have not ben accessed recently from memory.
    void doCleanup();

    CvCapture * captureDevice_; //!< OpenCV structure for capturing frames from an attached camera, or video file.
	vector<string> fileList_; //!< List of image filenames in directory.

    typedef std::vector<IplImage *> TFrameVector; 
    typedef std::vector<size_t> TFrameAccessTimeVector; 
    TFrameVector frameVector_;//!< Frames currently available in memory
    TFrameAccessTimeVector frameAccessTimeVector_;//!< Record frame access-times (so we can remove frames not accessed recently)
    boost::mutex lockFrameVector_; //!< Lock for frame cache

    const grc::Enums::eDisplayMode saveMosaic_; //!< Controls whether mosaics are saved.
    const char* mosaicSaveDir_;//!< Folder in which mosaics are saved.

    static const char * MOSAIC_WIN_NAME;//!< Mosaic window title
    static const char * CACHE_DIR;//!< Folder where frames are cached.

    Enums::eVideoSource eVidSource_;//!< Controls whether images come from .
public:
    ImageSource(Enums::eVideoSource eVidSource, const char* imDirName, const char* vidFilename, grc::Enums::eDisplayMode saveMosaic, const char* mosaicSaveDir );
    ~ImageSource();

    //! Read images from video device/file/whatever.
    /*! Calling this from main ensures images are captured by the main thread, otherwise OpenCv can't cope.
     *  Only returns once queue is empty (to ensure tidy cleanup, this will be implemented by blocking within the queue).
     *  Notifies engine once images have been captures.
     */
    void doCaptureImages(TransformEngine * engine, grc::Renderer * mosaicRenderer);

    //! Returns images for either rendering of for feature extraction.
    /*! Images will always be requested after we have notified engine they exist, so should always have been captured already.
     *  Eventually returns 0 to indicate it has finished.
     */
    const IplImage * getImage(size_t id);
};

}

#endif //_IMAGE_SOURCE
