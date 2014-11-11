#include "util/exception.h"
pragma_warning(push) // Disable warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)
pragma_warning(disable:4996)

#include "ImageSource.h"
#include "Renderer.h"
#include "GRCException.h"
#include "util/opencv_highgui.h"
#include <boost/filesystem.hpp>
#include "SpeedTest.h"

namespace grc {

const char * ImageSource::MOSAIC_WIN_NAME = "GRCMosaic";
const char * ImageSource::CACHE_DIR = "frameCache";

//! Reads images from video device/file/whatever. Capturing is done in doCaptureImages (blocking). We are all finished once this returns.

//! Read images from video device/file/whatever.
/*! Calling this from main ensures images are captured by the main thread, otherwise OpenCv can't cope.
 *  Only returns once queue is empty (to ensure tidy cleanup, this will be implemented by blocking within the queue).
 *  Notifies engine once images have been captures.
 */
void ImageSource::doCaptureImages(TransformEngine * engine, Renderer * mosaicRenderer)
{
    if(engine_ != 0 || mosaicRenderer_ != 0) throw new GRCException("ImageSource::doCaptureImages: Already initialised");
    if(engine == 0 || mosaicRenderer == 0) throw new GRCException("ImageSource::doCaptureImages: Null parameter");

    engine_ = engine;
    mosaicRenderer_ = mosaicRenderer;

	int delay = 1; // Number of milliseconds to wait for a keypress;
	if (eVidSource_ == Enums::eImageDirectory) {
		// For debugging it can be useful to set the delay to be long for a test image set
		delay = 5;
	}

    bool keepUp = (eVidSource_ != Enums::eAttachedCam);

    //Capture frame/render mosaic loop
    for(bool finished = false; !finished;) // Todo: keep on displaying rendered images until done--we may be behind (?)
    {
        //Render
        tryDisplayRenderedMosaic(keepUp);

        getImage(finished);

        int keyPressed = cvWaitKey(delay); //Need to wait: either for rendered image to display, or to avoid looping more than we need
        if(keyPressed == 'Q' || keyPressed == 'q')
            finished = true;
    }

    std::cout << "Image source finished\n";

    engine_->notifyFinished();
}

//! Returns images for either rendering of for feature extraction.
/*! Images will always be requested after we have notified engine they exist, so should always have been captured already.
 *  Eventually returns 0 to indicate it has finished.
 */
const IplImage * ImageSource::getImage(size_t id)
{
    boost::mutex::scoped_lock lockFrameVec(lockFrameVector_);

	if(id >= frameVector_.size()) {
		cout << id << endl;
        throw new GRCException("ImageSource::getImage: Frame not yet available");
	}

    //Log image access times so we know which frames to delete
    const int currentTime = frameVector_.size(); //Use current frame count as frame access time
    frameAccessTimeVector_[id] = currentTime;

    if(frameVector_[id] == 0)
    {
        //This frame has been dumped to disk--reload into memory
        char filename[50];
#ifdef GCC
        sprintf(filename, "%s/%05d.jpg", CACHE_DIR, id);
#else
        sprintf_s(filename, 50, "%s/%05d.jpg", CACHE_DIR, (int)id);
#endif
        CvPtr<IplImage> frame(cvLoadImage(filename));
		IplImage * frameBGRA = cvCreateImage(cvSize(frame->width, frame->height), IPL_DEPTH_8U, 4);
		cvCvtColor(frame, frameBGRA, CV_RGB2RGBA);
	    //cout << "Creating image with " << frameBGRA->nChannels << " channels\n";

        frameVector_[id] = frameBGRA;
    }

    if(!frameVector_[id]) throw new GRCException("ImageSource::getImage: Missing image");

    //cout << "Returning image with " << frameVector_[id]->nChannels << " channels\n";
    return frameVector_[id];
}

//! Free's memory used by images that haven't been accessed for a while.
/*! All images are saved on arrival, so will be reloaded if necessary
 */
void ImageSource::doCleanup()
{
    boost::mutex::scoped_lock lockFrameVec(lockFrameVector_);

    const size_t MEM_CACHE_TIME = 30;
    //Remove all images not accessed in the last MEM_CACHE_TIME frames
    const size_t currentTime = frameVector_.size();
    size_t timeCutoff = currentTime - MEM_CACHE_TIME;

    for(size_t i = 0; i < frameAccessTimeVector_.size(); i++)
    {
        if(frameAccessTimeVector_[i] < timeCutoff)
        {
            cvReleaseImage(&frameVector_[i]);
        }
    }
}

void removeBorder(IplImage * pFrame)
{
	const int WIDTH=pFrame->width;
	const int HEIGHT=pFrame->height;

	for(int x=0;x<WIDTH;x++)
	{
		CIplPx<uchar>::copyPx(pFrame, x, 1, pFrame, x, 0);
		CIplPx<uchar>::copyPx(pFrame, x, HEIGHT-2, pFrame, x, HEIGHT-1);
	}
	for(int y=0;y<HEIGHT;y++)
	{
		CIplPx<uchar>::copyPx(pFrame, 1, y, pFrame, 0, y);
		CIplPx<uchar>::copyPx(pFrame, WIDTH-2, y, pFrame, WIDTH-1, y);
	}
}

void ImageSource::storeImage(IplImage * pFrame)
{
    if(!pFrame || pFrame->height <= 0 || pFrame->width <= 0)
        throw new GRCException("ImageSource::storeImage: Bad image supplied");

    static const int CLEANUP_INTERVAL = 50;
    if(frameVector_.size() % CLEANUP_INTERVAL == 0) //Before we lock: cleanup (doesn't matter if size changes)
        doCleanup();

    boost::mutex::scoped_lock lockFrameVec(lockFrameVector_);

    if(!pFrame || pFrame->height <= 0 || pFrame->width <= 0)
        throw new GRCException("ImageSource::storeImage: Cleanup has cleaned up this frame while it is being processed");

    frameVector_.push_back(pFrame);

    size_t id = frameVector_.size()-1;
    frameAccessTimeVector_.push_back(id); //Set last access time to now

    if(frameAccessTimeVector_.size() != frameVector_.size()) throw new GRCException("ImageSource::storeImage: Vector lengths out-of-sync");

    //Setup getSize method (BEFORE notifying engine)
    static const bool setupSizeOnceOnly = engine_->notifySize(cvSize(pFrame->width, pFrame->height));

    //Notify engine that image has arrived:
    engine_->notifyImageAvailable(id);

    if(!pFrame || pFrame->height <= 0 || pFrame->width <= 0)
        throw new GRCException("ImageSource::storeImage: Image corrupted during notification");

    //Save a copy
    char filename[50];

#ifdef GCC
    sprintf(filename, "%s/%05d.jpg", CACHE_DIR, id);
#else
    sprintf_s(filename, 50, "%s/%05d.jpg", CACHE_DIR, (int)id);
#endif

    //if(eVidSource_ != Enums::eImageDirectory)
    	//cvSaveImage(filename, pFrame);
    REPEAT(1, cout << "NOT saving frames (for speed)\n");
}

//! Check if a rendered image is available for display, display it if so
void ImageSource::tryDisplayRenderedMosaic(bool keepUp)
{
    IplImage * pRenderedMosiac = mosaicRenderer_->getMosaicToDisplay(keepUp);
    if(pRenderedMosiac)
    {
    	//pRenderedMosiac->nChannels = 3;
    	//IplImage * pIm = cvLoadImage("/media/Ubuntu2/data/albert-B-laser-vision-dataset/albertB_1127373486-408821.jpg");
        //cvShowImage(MOSAIC_WIN_NAME, pIm);
    	//cv::showImage(MOSAIC_WIN_NAME, pRenderedMosiac);
        //cvShowImage(MOSAIC_WIN_NAME, pRenderedMosiac);
    	//cv::Mat m(pIm);
    	//cv::imshow("TEST", m);
    	cvShowImage(MOSAIC_WIN_NAME, pRenderedMosiac);
        cvWaitKey(1);

        if(saveMosaic_ == grc::Enums::eSaveImages)
        {
        	if(!boost::filesystem::exists(mosaicSaveDir_))
        		boost::filesystem::create_directory(mosaicSaveDir_);

    		char filename[500];
            sprintf_s(filename, 500, "%s/Mosaic%d.jpg", mosaicSaveDir_, (int)(frameVector_.size()-1));
    		cvSaveImage(filename, pRenderedMosiac);
        }
		cvReleaseImage(&pRenderedMosiac);
    }
}

//Actually blocking I think?? Todo: test
IplImage * ImageSource::tryGetImageFromCamOrAVI()
{
	if(cvGrabFrame(captureDevice_))
    {
	    IplImage *tmp = cvRetrieveFrame(captureDevice_);
        if(!tmp) throw new GRCException("tryGetImageFromCamOrAVI: No frame retreived");
		//std::cout << "Frame grabbed" << std::endl;
	    return cvCloneImage(tmp);
    }
    else
        return 0;
}



//! Check for a new frame
void ImageSource::getImage(bool & finished)
{
    if(eVidSource_ == Enums::eAttachedCam || eVidSource_ == Enums::eVideoFile)
    {
        IplImage * pFrame = tryGetImageFromCamOrAVI();
        if(!pFrame)
        {
            finished = true;
        }
        else
        {
			if (pFrame->origin) {
				cvFlip(pFrame);
				pFrame->origin = 0;
			}

            storeImage(pFrame);
        }
    }

	if (eVidSource_ == Enums::eImageDirectory) {
		static size_t ix = 0;

		IplImage *pFrame = NULL;
		while (!pFrame && ix < fileList_.size()) {
			pFrame = cvLoadImage(fileList_[ix].c_str());
			cout << "Loading image " << fileList_[ix].c_str() << "...";
                        cout << (pFrame ? "success" : "fail") << endl;
			ix++;
		}
		if (pFrame) {
			removeBorder(pFrame);
			{
				IplImage * frameBGRA = cvCreateImage(cvSize(pFrame->width, pFrame->height), IPL_DEPTH_8U, 4);
				cvCvtColor(pFrame, frameBGRA, CV_RGB2RGBA);
			    //cout << "Creating image with " << frameBGRA->nChannels << " channels\n";
			    storeImage(frameBGRA);
			    cvReleaseImage(&pFrame);
			}
		} else {
			finished = true;
		}
	}
}

ImageSource::ImageSource(Enums::eVideoSource eVidSource, const char* imDirName, const char* vidFilename, grc::Enums::eDisplayMode saveMosaic, const char* mosaicSaveDir) : engine_(0), captureDevice_(0), eVidSource_(eVidSource), mosaicRenderer_(0), saveMosaic_(saveMosaic), mosaicSaveDir_(mosaicSaveDir)
{
    if(eVidSource == Enums::eVideoFile && !vidFilename) throw new GRCException("ImageSource::ImageSource: No video source filename given");
	if(eVidSource == Enums::eImageDirectory && !imDirName) throw new GRCException("ImageSource::ImageSource: No directory given for JPEG files");

	// Create a directory to cache the frames in
    boost::filesystem::path dirPath(CACHE_DIR);
    boost::filesystem::create_directory(dirPath);

    cvNamedWindow(MOSAIC_WIN_NAME);

	if (eVidSource == Enums::eAttachedCam) {
        captureDevice_ = cvCaptureFromCAM(-1); //choose any cam
	} else if (eVidSource == Enums::eVideoFile) {
        captureDevice_ = cvCaptureFromAVI(vidFilename);
	} else if (eVidSource == Enums::eImageDirectory) {
		boost::filesystem::path dir(imDirName);
		if (!exists(dir) || !is_directory(dir)) {
			throw new GRCException("ImageSource::ImageSource: Directory not found");
		}
		boost::filesystem::directory_iterator endIter;
		boost::filesystem::directory_iterator iter(dir);
		for (; iter != endIter; iter++) {
			fileList_.push_back(iter->path().string());
		}
		sort( fileList_.begin(), fileList_.end() );

	} else {
		throw new GRCException("ImageSource::ImageSource: unexpected video source type");
	}

};

ImageSource::~ImageSource()
{
    boost::mutex::scoped_lock lockFrameVec(lockFrameVector_);

    for(size_t i = 0; i < frameVector_.size(); i++)
    {
        cvReleaseImage(&frameVector_[i]);
    }

    if(captureDevice_)
        cvReleaseCapture(&captureDevice_);

    cvDestroyWindow(MOSAIC_WIN_NAME);
};

}

pragma_warning(pop) // Disable warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)
