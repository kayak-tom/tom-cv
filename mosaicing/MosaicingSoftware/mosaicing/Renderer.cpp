#include "util/exception.h"
//! Implementation: Reads latest transformation set from Transform engine and renders all of it.
pragma_warning(push) // Disable warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)
pragma_warning(disable:4996)

#include "Renderer.h"
#include "GRCException.h"
#include "boost/bind.hpp"
#include "boost/filesystem.hpp"
#include <typeinfo>
#include "cvGeom.h"

#include "SpeedTest.h"

namespace grc {

	Renderer::Renderer(ImageSource & imageSource, TransformEngine & transformEngine, CvSize mosaicSize, EvaluationFunction *evaluationFunction)
		: imageSource_(imageSource), transformEngine_(transformEngine), imToRender_(0), mosaicSize_(mosaicSize), evaluationFunction_(evaluationFunction),
		  mosaicAllocated_(false), mosaic_(0),
		  semWaitForFrame_(MAX_FRAMES_AHEAD), renderingThread_(0 /* renderingThread_(boost::bind(&grc::Renderer::doRendering, this)) gives a "this isn't constructed yet" warning */ )

	{
		/* Instantiate rendering thread last: */
		renderingThread_ = new boost::thread(boost::bind(&grc::Renderer::doRendering, this));
	};

	Renderer::~Renderer()
	{
        std::cout << "Rejoining rendering thread\n";
		renderingThread_->join();
		delete renderingThread_; renderingThread_=0;
        std::cout << "Done rejoining rendering thread, notifying engine...\n";
        transformEngine_.notifyRendererFinished();

		if(mosaicAllocated_) {
			cvReleaseImage(&mosaic_);
		}
	};

	IplImage * Renderer::renderImagesToMosaic(TransformSet * TS)
	{
		if(!TS)
			throw new GRCException("Renderer::renderImagesToMosaic: Rendering 0-transform--should have died");
		

		if (!mosaicAllocated_) {
			const IplImage *frame = imageSource_.getImage(TS->begin()->first.im1Id());
			/*CvPtr<IplImage> frame ( cvCreateImage(cvSize(frame_temp->width, frame_temp->height), IPL_DEPTH_8U, 4) );
			cvCvtColor(frame_temp, frame, CV_RGB2RGBA);*/

			mosaic_ = cvCreateImage(mosaicSize_, frame->depth, frame->nChannels);
			cvSetZero(mosaic_);
			mosaicAllocated_ = true;
		}	

		IplImage *result = cvCreateImage(mosaicSize_, mosaic_->depth, mosaic_->nChannels);
		cvSetZero(result);

        const int NONE = -1;
        int idAlignedTo = NONE; //For consistency check
		if (TS->onlyRenderLatestTransform()) {
			TS->mosaicTransform()->applyToImage(mosaic_, result);


			const IplImage *src = imageSource_.getImage(TS->latestTransform()->id1());
			/*CvPtr<IplImage> src ( cvCreateImage(cvSize(src_temp->width, src_temp->height), IPL_DEPTH_8U, 4) );
			cvCvtColor(src_temp, src, CV_RGB2RGBA);*/

			TS->latestTransform()->transform()->applyToImage(src, result);
		} else {
			for (TransformSet::const_iterator iter = TS->begin(); iter != TS->end(); iter++) {
				const IplImage *src = imageSource_.getImage(iter->first.im1Id());
				/*CvPtr<IplImage> src ( cvCreateImage(cvSize(src_temp->width, src_temp->height), IPL_DEPTH_8U, 4) );
				cvCvtColor(src_temp, src, CV_RGB2RGBA);*/


				if(!iter->second)
					throw new GRCException("Renderer::renderImagesToMosaic: TransformInfo structure missing");

				const Transform * transform = iter->second->transform();

				if(!transform)
					throw new GRCException("Renderer::renderImagesToMosaic: TransformInfo structure is missing a transform");

				//Check bookkeeping/ba has worked and everything is being warped to the same frame
				if(idAlignedTo == NONE) idAlignedTo=iter->second->id2();
				if(idAlignedTo != (int)iter->second->id2())
					throw new GRCException("Renderer::renderImagesToMosaic: Images in transform set not all warped to same frame");

				if(iter->second->id1() != iter->first.im1Id() || iter->second->id2() != iter->first.im2Id())
					throw new GRCException("Renderer::renderImagesToMosaic: Inconsistent image ids");

    			transform->applyToImage(src, result);
			}
		}
		cvCopy(result, mosaic_);
		return result;
	}

	void Renderer::doRendering()
	{
		//int frameNum = 0;
		try
		{
			for(;;)
			{
				TransformSet * TS  = transformEngine_.getLatestTransforms(); //will block until one is ready
				if(!TS) //die
				{
					break; //will return
				}

				IplImage * pImg = renderImagesToMosaic(TS);
				if (evaluationFunction_) {
					cout << "Evaluation: " << evaluationFunction_->evaluate(pImg, TS, imageSource_) << endl;

				}

				const bool bMarkLatest = true;
				if(bMarkLatest)
				{
					static const IplImage * temp = imageSource_.getImage(0);
					static CvSize imageSize_ = cvSize(temp->width, temp->height);
					static const CvMat * sourceCorners = getSourceCorners(imageSize_);
				    static CvMat * projectedCorners = getSourceCorners(imageSize_);

				    const TransformInfo2 * pLatestTrans = TS->latestTransform();
				    if(pLatestTrans)
				    {
						pLatestTrans->transform()->applyToPoints(sourceCorners, projectedCorners);

						for(int i=0; i<4; i++)
						{
							const CvScalar RED = CV_RGB(255,0,0);
							const int WIDTH=5, LENGTH=15;
							double x = cvmGet(projectedCorners, 0, i);
							double y = cvmGet(projectedCorners, 1, i);
							CvPoint corner = cvPoint((int)x,(int)y);
							CvPoint corner1 = corner;
							CvPoint corner2 = corner;
							if(i % 2 == 0)
								corner1.x -= LENGTH;
							else
								corner1.x += LENGTH;

							if(i < 2)
								corner2.y += LENGTH;
							else
								corner2.y -= LENGTH;

							cvLine(pImg, corner1, corner, RED, WIDTH);
							cvLine(pImg, corner2, corner, RED, WIDTH);

							//cvCircle(pImg, corner, 2, RED);

						}
				    }

				}

				delete TS;

				boost::mutex::scoped_lock renderDisplayLock( renderLock_ );

				if(imToRender_ != 0)
				{
					std::cout << "Warning: Renderer::doRendering: Display thread is not keeping up. Another rendered frame is ready--deleting previous rendered frame\n";
					cvReleaseImage(&imToRender_);
					//throw new GRCException("Renderer::doRendering: Haven't rendered last frame yet");
				}
                // imToRender_ should be 0 now
				imToRender_ = pImg;
                semWaitForFrame_.post(); //May be doing nothing
			}
		}
		catch(grc::GRCException * pEx)
		{
			std::cout << "ERROR: Unhandled exception in renderer main thread: " << *(pEx->GetErrorMessage()) << std::endl ;
			delete pEx;
		}

        transformEngine_.notifyRendererFinished();
	}

	IplImage * Renderer::getMosaicToDisplay(bool block)
	{
        if(block)
            //Make the engine wait for another frame before grabbing another:
            semWaitForFrame_.wait();

		boost::mutex::scoped_lock renderDisplayLock( renderLock_ );

		IplImage * mosaicToShow = imToRender_;
		imToRender_ = 0;
		return mosaicToShow; // may be 0 if we've already displayed the most recent rendered frame
	}


}

pragma_warning(pop) // Disable warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)
