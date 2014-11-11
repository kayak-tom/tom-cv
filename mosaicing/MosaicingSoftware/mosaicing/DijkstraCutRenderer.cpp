#include "util/exception.h"

pragma_warning(push)
pragma_warning(disable:4996)  // warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)

#include "DijkstraCutRenderer.h"
#include "SpeedTest.h"
#include <limits>
#include <vector>
#include <stdio.h>
#include "myGetSet2D.h"
#include "Warp.h"
//#ifndef USE_OLD_OPENCV
//#include "opencv2/opencv.hpp"
//#endif

template<typename T>
inline void bound(T & val, const T min, const T max)
{
	if(val < min)
		val=min;
	else if(val > max)
		val=max;
}

static const double infinity = 1E+12;// numeric_limits<double>::infinity();

namespace grc {

    /*template<unsigned int SCALE>
	void DijkstraCutRenderer::resizeAndThresh()
    {
        //The code below is roughly equivalent to, but faster than, the following two commented lines.
        //cvResize(totalMask_, largeMask_, CV_INTER_LINEAR);
        //cvThreshold(largeMask_, largeMask_, 128, 255, CV_THRESH_BINARY);
		const double scale_inv = 1.0/SCALE;

		for (int y = 0; y < totalMask_->height; y++) {
			for (int x = 0; x < totalMask_->width; x++) {
				double p00 = myGet2D(totalMask_, x, y);
				double p01 = myGet2D(totalMask_, x, y+1);
				double p10 = myGet2D(totalMask_, x+1, y);
				double p11 = myGet2D(totalMask_, x+1, y+1);
				double a=0;
//				uchar * pData =


				for (unsigned int ndx = 0; ndx < SCALE; ndx++, a += scale_inv) {
					//double a = ndx*scale_inv;
					double v1 = a*p10 + (1-a)*p00;
					double v2 = a*p11 + (1-a)*p01;
					double b=0;
					for (unsigned int ndy = 0; ndy < SCALE; ndy++, b += scale_inv) {
						//double b = ndy*scale_inv;
						double v = b*v2 + (1-b)*v1;
						const int nx=SCALE*x+ndx, ny = SCALE*y+ndy;

						unsigned char src_px = myGet2D(largeMosaic_, nx, ny, 0);

						//largeMask_->imageData[ny*largeMask_->widthStep + nx*largeMask_->nChannels] = value;


						if ( (v > 128) || !src_px ) {
							mySet2D(largeMask_, nx, ny, 255);
						} else {
							mySet2D(largeMask_, nx, ny, 0);
						}
					}
				}
			}
		}
    }*/

    template<unsigned int SCALE, unsigned int CHANNELS>
	void DijkstraCutRenderer::resizeAndThresh()
    {
        //The code below is roughly equivalent to, but faster than, the following two commented lines.
        //cvResize(totalMask_, largeMask_, CV_INTER_LINEAR);
        //cvThreshold(largeMask_, largeMask_, 128, 255, CV_THRESH_BINARY);
		//const double scale_inv = 1.0/SCALE;

		for (int y = 0; y < totalMask_->height; y++) {
			for (int x = 0; x < totalMask_->width; x++) {
				int p00 = myGet2D(totalMask_, x, y);
				int p01 = myGet2D(totalMask_, x, y+1);
				int p10 = myGet2D(totalMask_, x+1, y);
				int p11 = myGet2D(totalMask_, x+1, y+1); //These may be off the edge, but creates artifacts if we only go to totalMask_->width-1
				//int a=0; use ndx

				uchar * pOutDataStart = (uchar *)largeMask_->imageData + SCALE*y*largeMask_->widthStep + SCALE*x*largeMask_->nChannels;
				uchar * pInDataStart = (uchar *)largeMosaic_->imageData + SCALE*y*largeMosaic_->widthStep + SCALE*x*CHANNELS;

				int sum = p00 + p01 + p10 + p11;
				if(sum == 0)
				{
					for (unsigned int ndx = 0; ndx < SCALE;
							ndx++,
							pOutDataStart += largeMask_->widthStep,
							pInDataStart += largeMosaic_->widthStep
							) {

						uchar * pInData = pInDataStart;
						uchar * pOutData = pOutDataStart;

						for (unsigned int ndy = 0; ndy < SCALE;
								ndy++,
								pOutData++, pInData += CHANNELS) {

							/*if(CHANNELS == 4)
							{
								volatile unsigned int src_px = *reinterpret_cast<const unsigned int *>(pInData); //nonzero flagged channel (blue) indicates point in image
								*pOutData = ((v > sqr((int)SCALE)*128) || !src_px ) ? '\255' : '\0'; //2nd critical loop
							}
							else*/
							{
								volatile unsigned char src_px = *pInData;
								*pOutData = !src_px ? '\255' : '\0';
							}
						}
					}
				}
				else if (sum == 4*255)
				{
					for (unsigned int ndx = 0; ndx < SCALE;
							ndx++,
							pOutDataStart += largeMask_->widthStep,
							pInDataStart += largeMosaic_->widthStep
							) {

						uchar * pOutData = pOutDataStart;

						for (unsigned int ndy = 0; ndy < SCALE;
								ndy++,
								pOutData++) {

							/*if(CHANNELS == 4)
							{
								volatile unsigned int src_px = *reinterpret_cast<const unsigned int *>(pInData); //nonzero flagged channel (blue) indicates point in image
								*pOutData = ((v > sqr((int)SCALE)*128) || !src_px ) ? '\255' : '\0'; //2nd critical loop
							}
							else*/
							{
								*pOutData = '\255';
							}
						}
					}
				}
				else
				{

					int v1 = SCALE*p00;
					int v2 = SCALE*p01; //Are these 1's and 0's? Yes 255's and 0's

					int v1inc = p10-p00;
					int v2inc = p11-p01;

					for (unsigned int ndx = 0; ndx < SCALE;
							ndx++,
							//a += scale_inv,
							pOutDataStart += largeMask_->widthStep,
							pInDataStart += largeMosaic_->widthStep,
							v1 += v1inc,
							v2 += v2inc
							) {

						uchar * pInData = pInDataStart;
						uchar * pOutData = pOutDataStart;

						//int v1 = ndx*p10 + (SCALE-ndx)*p00;
						//int v2 = ndx*p11 + (SCALE-ndx)*p01; //Are these 1's and 0's? Yes 255's and 0's

						//double b=0;
						int v = SCALE*v1;
						int v_inc = v2 - v1;
						for (unsigned int ndy = 0; ndy < SCALE;
								ndy++,
								//b += scale_inv,
								pOutData++, pInData += CHANNELS, v += v_inc) {

							//volatile int v = ndy*v2 + (SCALE-ndy)*v1;
							//v in 0...SCALE^2 * 255
							//const int nx=SCALE*x+ndx, ny = SCALE*y+ndy;

							if(CHANNELS == 4)
							{
								volatile unsigned int src_px = *reinterpret_cast<const unsigned int *>(pInData); //nonzero flagged channel (blue) indicates point in image
								//cout << hex << src_px << ' ';
								*pOutData = ((v > sqr((int)SCALE)*128) || !src_px ) ? '\255' : '\0'; //2nd critical loop
							}
							else
							{
								volatile unsigned char src_px = *pInData; //myGet2D(largeMosaic_, nx, ny, 0); //flagged channel (blue) indicates point in image

								//largeMask_->imageData[ny*largeMask_->widthStep + nx*largeMask_->nChannels] = value;

								*pOutData = ((v > sqr((int)SCALE)*128) || !src_px ) ? '\255' : '\0';
							}
	/*						if ( (v > 128) || !src_px ) {
								mySet2D(largeMask_, nx, ny, 255);
							} else {
								mySet2D(largeMask_, nx, ny, 0);
							}*/
						}
					}
				}
			}
		}
    }

    template<unsigned int CHANNELS>
    void DijkstraCutRenderer::findAndApplyCut2(const IplImage *frame, const Transform *transform)
    {
    	switch(scale_)
    	{
    	case 1:
    	    findAndApplyCut_template<1, CHANNELS>(frame, transform);
    	    break;
    	case 2:
    	    findAndApplyCut_template<2, CHANNELS>(frame, transform);
    	    break;
    	case 3:
    	    findAndApplyCut_template<3, CHANNELS>(frame, transform);
    	    break;
    	case 4:
    	    findAndApplyCut_template<4, CHANNELS>(frame, transform);
    	    break;
    	case 5:
    	    findAndApplyCut_template<5, CHANNELS>(frame, transform);
    	    break;
    	case 6:
    	    findAndApplyCut_template<6, CHANNELS>(frame, transform);
    	    break;
    	case 7:
    	    findAndApplyCut_template<7, CHANNELS>(frame, transform);
    	    break;
    	case 8:
    	    findAndApplyCut_template<8, CHANNELS>(frame, transform);
    	    break;
    	default:
    		THROW("Scale too low or too high (maybe should be supported?)")
    	}
    }
    void DijkstraCutRenderer::findAndApplyCut(const IplImage *frame, const Transform *transform)
    {
    	switch(frame->nChannels)
    	{
    	case 1:
    		findAndApplyCut2<1>(frame, transform);
    	    break;
    	case 3:
    		findAndApplyCut2<3>(frame, transform);
    	    break;
    	case 4:
    		findAndApplyCut2<4>(frame, transform);
    	    break;
    	default:
    		THROW("Channels not 1 or 3 or 4 (maybe should be supported?)")
    	}
    }
    template<unsigned int SCALE, unsigned int CHANNELS>
    void DijkstraCutRenderer::findAndApplyCut_template(const IplImage *frame, const Transform *transform)
    {
    	const bool DEBUG_IMAGES = false;

    	CvSize dsize = cvGetSize(largeMosaic_);
        BB bbLatestFrame(getBB(((const CvMat *)(*transform))->data.db, dsize, frameSize_));


        vector<pair<int,int> > search;

        enum eStates { eInMosaic = 0, eOutsideImage = 1, eOutsideMosaic = 2, eFirstPx = 4 };
		int x_dest_old = SCALE;
		int y_dest_old = SCALE;

		for(int y=0; y<frameSize_.height; y+=frameSize_.height-1)
		{
	        eStates status = eFirstPx;
			for(int x=0; x<frameSize_.width; x++)
			{
				CvPoint2D64f ptSrc = cvPoint2D64f(x, y);
				CvPoint2D64f ptDest = transform->applyToPoint(ptSrc);
				int x_dest = doubleToInt(ptDest.x), y_dest = doubleToInt(ptDest.y);

				eStates newStatus = eOutsideMosaic;
				if(x_dest<(int)SCALE || y_dest<(int)SCALE || x_dest >= largeMosaic_->width-2*(int)SCALE || y_dest >= largeMosaic_->height-2*(int)SCALE)
					newStatus = eOutsideImage;
				else if(myGet2D(largeMosaic_, x_dest, y_dest, 0))
					newStatus = eInMosaic;

				if(newStatus != status)
				{
					if(status != eFirstPx && ((newStatus | status) != (eOutsideImage | eOutsideMosaic)))
					{
						int x_mask = (newStatus == eInMosaic) ? x_dest/SCALE : x_dest_old/SCALE ;
						int y_mask = (newStatus == eInMosaic) ? y_dest/SCALE : y_dest_old/SCALE ;
						bound(x_mask, 1, scaledMosaicSize_.width-2);
						bound(y_mask, 1, scaledMosaicSize_.height-2);
                        search.push_back(make_pair<int,int>(x_mask, y_mask));
					}
					status = newStatus;
				}
				x_dest_old = x_dest;
				y_dest_old = y_dest;
			}
		}
		for(int x=0; x<frameSize_.width; x+=frameSize_.width-1)
		{
	        eStates status = eFirstPx;
			for(int y=0; y<frameSize_.height; y++)
			{
				CvPoint2D64f ptSrc = cvPoint2D64f(x, y);
				CvPoint2D64f ptDest = transform->applyToPoint(ptSrc);
				int x_dest = doubleToInt(ptDest.x), y_dest = doubleToInt(ptDest.y);

				eStates newStatus = eOutsideMosaic;
				if(x_dest<(int)SCALE || y_dest<(int)SCALE || x_dest >= largeMosaic_->width-(int)SCALE || y_dest >= largeMosaic_->height-(int)SCALE)
					newStatus = eOutsideImage;
				else if(myGet2D(largeMosaic_, x_dest, y_dest, 0))
					newStatus = eInMosaic;

				if(newStatus != status)
				{
					if(status != eFirstPx && ((newStatus | status) != (eOutsideImage | eOutsideMosaic)))
					{
						int x_mask = (newStatus == eInMosaic) ? x_dest/SCALE : x_dest_old/SCALE ;
						int y_mask = (newStatus == eInMosaic) ? y_dest/SCALE : y_dest_old/SCALE ;
						bound(x_mask, 1, scaledMosaicSize_.width-2);
						bound(y_mask, 1, scaledMosaicSize_.height-2);
                        search.push_back(make_pair<int,int>(x_mask, y_mask));
					}
					status = newStatus;
				}
				x_dest_old = x_dest;
				y_dest_old = y_dest;
			}
		}

        if (search.size() < 2)
        {
			transform->applyToImage(frame, largeMosaic_); //Warp straight into mosaic
			if (DEBUG_IMAGES) cout << search.size() << " No cut found. Happens when one image obliterates another, first frame, sometimes when cut search fails\n";

		} else {
			if (DEBUG_IMAGES) cout << search.size() << " Cut found\n";

			cvSetZero(largeWarped_); //could prob just apply to BB...
			transform->applyToImage(frame, largeWarped_);

			if (DEBUG_IMAGES) cvSaveImage("largeWarped_before.jpg", largeWarped_);
			if (DEBUG_IMAGES) cvSaveImage("largeMosaic_before.jpg", largeMosaic_);

			cvSetZero(totalMask_); //Excessive as we're about to set the whole thing almost

			{// Extra braces stop MS C++ complaining about ANSI scoping of for
				// Compute the graph costs
				const int SUBPIX_STEP = (SCALE + 1)/2;

				for (int y = 1; y < scaledMosaicSize_.height-1; ++y) {
					//unsigned char *in1 = myGetRowStart(mosaicMask_, y)+1;
					//unsigned char *in2 = myGetRowStart(warpedMask_, y)+1;

					const int SCALE_y=SCALE*y;

					for (int x = 1; x < scaledMosaicSize_.width-1; ++x/*, ++in1, ++in2*/) {

						const int SCALE_x=SCALE*x;

						uchar * pLargeMosaicPxStart = (uchar *)(largeMosaic_->imageData + SCALE_y*largeMosaic_->widthStep + SCALE_x*CHANNELS);
						uchar * pLargeWarpedPxStart = (uchar *)(largeWarped_->imageData + SCALE_y*largeWarped_->widthStep + SCALE_x*CHANNELS);
						//if (*in1 && *in2)
						bool bInNewIm = (*pLargeWarpedPxStart) ? true : false;
						bool bInMosaic = (*pLargeMosaicPxStart) ? true : false;
						bool bInBoth = bInNewIm && bInMosaic;
						bool bInEither = bInNewIm || bInMosaic;

						mySet2D(warpedMask_, x, y, bInNewIm ? 255 : 0);
						mySet2D(totalMask_, x, y, bInEither ? 255 : 0);

						if (bInBoth)
						{
							cost_[x][y] = 0;
							int nCost = 0;

							for (unsigned int dy = 0; dy < SCALE; dy+=SUBPIX_STEP, pLargeMosaicPxStart += SUBPIX_STEP*largeMosaic_->widthStep, pLargeWarpedPxStart += SUBPIX_STEP*largeWarped_->widthStep)
							{
								uchar * pLargeMosaicPx = pLargeMosaicPxStart;
								uchar * pLargeWarpedPx = pLargeWarpedPxStart;
								for (unsigned int dx = 0; dx < SCALE; dx+=SUBPIX_STEP, pLargeMosaicPx+=SUBPIX_STEP*CHANNELS, pLargeWarpedPx+=SUBPIX_STEP*CHANNELS)
								{
									for (int c = 0; c < (int)CHANNELS; c++) //3rd critical loop
									{ //Todo: alpha channel is always 0 here...
										volatile const int val1 = pLargeMosaicPx[c];
										if(IS_DEBUG) CHECK(myGet2D(largeMosaic_, SCALE_x+dx, SCALE_y+dy, c) != val1, "Ptr conversion error");

										volatile const int val2 = pLargeWarpedPx[c];
										if(IS_DEBUG) CHECK(myGet2D(largeWarped_, SCALE_x+dx, SCALE_y+dy, c) != val2, "Ptr 2 conversion error");

										//int thisCost = int(myGet2D(largeMosaic_, SCALE_x+dx, SCALE_y+dy, c)) - int(myGet2D(largeWarped_, SCALE_x+dx, SCALE_y+dy, c));
										nCost += sqr(val2 - val1);
									}
								}
							}
							cost_[x][y] = (double)nCost;
						} else {
							cost_[x][y] = infinity;
						}
					}
				}
				//cvErode(totalMask_, totalMask_); // totalMask_ contains latest image AND MOSAIC. erode should stop flood-fill overflowing
			}
			if (DEBUG_IMAGES) cvSaveImage("totalMask_1.jpg", totalMask_);
			if (DEBUG_IMAGES) cvSaveImage("warpedMask_1.jpg", warpedMask_);


			CDynArray<bool> solved(search.size(), false);
			// For each crossing point, make sure there's a cut leading to it
			for (size_t i = 0; i < search.size(); i++) {
				if (!solved[i]) {
					// No cut yet, initialise the graph search

					{// Extra braces stop MS C++ complaining about ANSI scoping of for
						for (int x = 0; x < scaledMosaicSize_.width; x++) {
							for (int y = 0; y < scaledMosaicSize_.height; y++) {
								distance_[x][y] = infinity;
								expanded_[x][y] = false;
								dx_[x][y] = 0;
								dy_[x][y] = 0;
							}
						}
						for (int x = 0; x < scaledMosaicSize_.width; x++) {
							expanded_[x][0] = expanded_[x][scaledMosaicSize_.height-1] = true;
						}
						for (int y = 0; y < scaledMosaicSize_.height; y++) {
							expanded_[0][y] = expanded_[scaledMosaicSize_.width-1][y] = true;
						}

					}
				
					// Get the starting point for the search

					int x = search[i].first;
					int y = search[i].second;
					if(IS_DEBUG) CHECK(!(x > 0 && y > 0 && x < scaledMosaicSize_.width -1 && y < scaledMosaicSize_.height - 1), "xy OOB")

					distance_[x][y] = 0;
					solved[i] = true;

					//Make surrounding weights low in case we're just off the edge...
					for (int x2 = x-1; x2 <= x+1; x2++) {
						for (int y2 = y-1; y2 <= y+1; y2++)
						{
							if(x2 > 0 && y2 > 0 && x2 < scaledMosaicSize_.width -1 && y2 < scaledMosaicSize_.height)
								if (cost_[x2][y2] == infinity)
									cost_[x2][y2] = sqr(255);
						}
					}

					// Run Dijkstra's algorithm across the graph

					PriorityQueue Q;

					Q.push(x, y, infinity, 0);
					bool done = false;
					while (!done) {

						Q.pop(&x, &y);

						if(IS_DEBUG) CHECK(!(x > 0 && y > 0 && x < scaledMosaicSize_.width -1 && y < scaledMosaicSize_.height - 1), "xy OOB")

						expanded_[x][y] = true;
						for (int x2 = x-1; x2 <= x+1; x2++) {
							for (int y2 = y-1; y2 <= y+1; y2++) {
								if (!expanded_[x2][y2]) {
									double newDist = distance_[x][y] + cost_[x2][y2];
									if (newDist < distance_[x2][y2]) {
										Q.push(x2, y2, distance_[x2][y2], newDist);
										distance_[x2][y2] = newDist;
										dx_[x2][y2] = x-x2;
										dy_[x2][y2] = y-y2;
									}
								}
							}
						}

						for (size_t j = 0; j < search.size(); j++) {

							if ( (j != i) && (search[j].first == x) && (search[j].second == y) ) {
								// Reached a destination
								done = true;
								solved[j] = true;
								// Trace a path back to the start
								mySet2D(totalMask_, x, y, 200); //This draws a line back through the total mask so we can cut off the existing mosaic
								while ( (dx_[x][y] != 0) || (dy_[x][y] != 0) ) {
									int t = x;
									x += dx_[x][y];
									y += dy_[t][y];
									mySet2D(totalMask_, x, y, 200);
								}
							}
						}
						if (!done && Q.isEmpty()) {
                            std::cout << "WARNING: DijkstraCutRenderer::renderImagesToMosaic: Empty queue while searching for path\n";
                            done = true;
						}
					}
				} // For each unsolved crossing point
				CvPoint p = cvPoint(search[i].first, search[i].second);
				cvCircle(totalMask_, p, 1, cvScalarAll(200), -1);
			} // For each crossing point
			
			if (DEBUG_IMAGES) cvSaveImage("totalMask_.jpg", totalMask_);
			if (DEBUG_IMAGES) cvSaveImage("warpedMask_.jpg", warpedMask_);

			//cvErode(totalMask_, totalMask_); // totalMask_ contains latest image AND MOSAIC. helps flood-fill not sneak around edge. Or could put a littls square at each X point
			//Anything in totalmask AND NOT image (i.e. in old mosaic) set 0
			//int temp=0;
			for (int x = 0; x < scaledMosaicSize_.width; x++) {
				for (int y = 0; y < scaledMosaicSize_.height; y++) {
					if ( (myGet2D(totalMask_, x, y) == 255) && (myGet2D(warpedMask_, x, y) == 0) ) {
						// Point inside the mosaic, but outside the new image
						cvFloodFill(totalMask_, cvPoint(x,y), cvScalarAll(0));

						/*char fn[100];
						sprintf(fn, "FF%d.jpg", temp++);
						if (DEBUG_IMAGES) cvSaveImage(fn, totalMask_);*/
					}
				}
			}
			//cvErode(totalMask_, totalMask_); // totalMask_ contains latest image AND MOSAIC. erode gets rid of little artifacts

			if (DEBUG_IMAGES) cvSaveImage("totalMask_filled.jpg", totalMask_);

			//Vals of totalMask_ are used--could improve this...
			cvThreshold(totalMask_, totalMask_, 1, 255, CV_THRESH_BINARY);
			
			//cvOr(warpedMask_, totalMask_, mosaicMask_);
			if (DEBUG_IMAGES) cvSaveImage("mosaicMask_.jpg", mosaicMask_);

			//The code below is roughly equivalent to, but faster than, the following two commented lines.
			//cvResize(totalMask_, largeMask_, CV_INTER_LINEAR);
			//cvThreshold(largeMask_, largeMask_, 128, 255, CV_THRESH_BINARY);
			resizeAndThresh<SCALE, CHANNELS>(); //uses totalMask_ and largeMosaic_ to set largeMask_
			/*const double scale_inv = 1.0/SCALE;

			for (int y = 0; y < totalMask_->height; y++) {
				for (int x = 0; x < totalMask_->width; x++) {
					double p00 = myGet2D(totalMask_, x, y);
					double p01 = myGet2D(totalMask_, x, y+1);
					double p10 = myGet2D(totalMask_, x+1, y);
					double p11 = myGet2D(totalMask_, x+1, y+1);
					double a=0;
					for (unsigned int ndx = 0; ndx < SCALE; ndx++, a += scale_inv) {
						//double a = ndx*scale_inv;
						double v1 = a*p10 + (1-a)*p00;
						double v2 = a*p11 + (1-a)*p01;
						double b=0;
						for (unsigned int ndy = 0; ndy < SCALE; ndy++, b += scale_inv) {
							//double b = ndy*scale_inv;
							double v = b*v2 + (1-b)*v1;
							const int nx=SCALE*x+ndx, ny = SCALE*y+ndy;

							unsigned char src_px = myGet2D(largeMosaic_, nx, ny, 0);

							//img->imageData[ny*img->widthStep + nx*img->nChannels] = value;

							if ( (v > 128) || !src_px ) {
								mySet2D(largeMask_, nx, ny, 255);
							} else {
								mySet2D(largeMask_, nx, ny, 0);
							}
						}
					}
				}
			}*/
			
			//cvCopy(largeWarped_, largeMosaic_, largeMask_);
			// We replace the above line with a check to make sure that the 1st channel of the 
			// mosaic is only zero outside of the rendered region - this lets us use the mosaic to check for errors
			// around the borders of images above.

			//cvSetZero(checkMask_);
			//transform->applyToImage(checkFrame_, checkMask_); //TODO set this up earlier
			//cvAnd(checkMask_, largeMask_, checkMask_); //Could use ROI
			if (DEBUG_IMAGES) cvSaveImage("largeMask_.jpg", largeMask_);
			if (DEBUG_IMAGES) cvSaveImage("mosaicBefore.jpg", largeMosaic_);

			for (int y = bbLatestFrame.minY; y < bbLatestFrame.maxY; y++)
			{
				uchar * pSrcMask = (uchar*)(largeMask_->imageData + bbLatestFrame.minX + y*largeMask_->widthStep);
				uchar * pSrc = (uchar*)(largeWarped_->imageData + bbLatestFrame.minX*CHANNELS + y*largeWarped_->widthStep);
				uchar * pDst = (uchar*)(largeMosaic_->imageData + bbLatestFrame.minX*CHANNELS + y*largeMosaic_->widthStep);

				for (int x = bbLatestFrame.minX; x < bbLatestFrame.maxX; x++, pSrcMask++, pSrc += CHANNELS, pDst += CHANNELS) { //only need to iterate over BB
					if(IS_DEBUG) CHECK(*pSrcMask != myGet2D(largeMask_, x, y), "Bad conversion to ptr");
					if (*pSrcMask)// && myGet2D(largeMask_, x, y))
					{ //1st critical loop
						if(CHANNELS == 4)
						{
							volatile const int val = *(reinterpret_cast<int *>(pSrc));// *(myGet2D_int(largeWarped_, x, y));
							if(IS_DEBUG) CHECK(val != *(myGet2D_int(largeWarped_, x, y)), "Bad conversion to ptr");
							//val |= 0x00000001; //Could be DEBUGONLY?
							//*(myGet2D_int(largeMosaic_, x, y)) = val;
							volatile int * pnDst = reinterpret_cast<int *>(pDst);
							if(IS_DEBUG) CHECK(myGet2D_int(largeMosaic_, x, y) != pnDst, "Bad conversion to ptr");
							*pnDst = val;
						}
						else
						{
							for (int c = 0; c < (int)CHANNELS; c++) {
								unsigned char val = pSrc[c]; //myGet2D(largeWarped_, x, y, c);
								if (c == 0 && val == 0) val = 1;
								//mySet2D(largeMosaic_, x, y, c, val);
								pDst[c] = val;
							}
						}
					}
				}
			}
			if (DEBUG_IMAGES) cvSaveImage("mosaicAfter.jpg", largeMosaic_);
			if (DEBUG_IMAGES) cout << "";
		}
        //cvReleaseMat(&T);
	}

	void DijkstraCutRenderer::initialiseStorage(TransformSet *TS) {
		const IplImage *frame = imageSource_.getImage(TS->begin()->first.im1Id());

		frameSize_ = cvGetSize(frame);
		scaledFrameSize_ = cvSize(frameSize_.width/scale_, frameSize_.height/scale_);
		scaledMosaicSize_ = cvSize(mosaicSize_.width/scale_, mosaicSize_.height/scale_);

		// Masks, 8 bit single channel images
		mosaicMask_ = cvCreateImage(scaledMosaicSize_, IPL_DEPTH_8U, 1);
		warpedMask_ = cvCreateImage(scaledMosaicSize_, IPL_DEPTH_8U, 1);
		//mosaicEdge_ = cvCreateImage(scaledMosaicSize_, IPL_DEPTH_8U, 1);
		//warpedEdge_ = cvCreateImage(scaledMosaicSize_, IPL_DEPTH_8U, 1);
		//crossingMask_ = cvCreateImage(scaledMosaicSize_, IPL_DEPTH_8U, 1);
		totalMask_ = cvCreateImage(cvSize(scaledMosaicSize_.width, scaledMosaicSize_.height+1), IPL_DEPTH_8U, 1); //so access OOB is harmless
		largeMask_ = cvCreateImage(mosaicSize_, IPL_DEPTH_8U, 1);
		//checkMask_ = cvCreateImage(mosaicSize_, IPL_DEPTH_8U, 1);
		//checkFrame_ = cvCreateImage(frameSize_, IPL_DEPTH_8U, 1);
		//cvSet(checkFrame_, cvScalarAll(255));

		// Other images have same bit depth and channels as the first frame
		largeMosaic_ = cvCreateImage(mosaicSize_, frame->depth, frame->nChannels);
		largeWarped_ = cvCreateImage(mosaicSize_, frame->depth, frame->nChannels);

		// Various matrices of use
		scaleTransform_ = cvCreateMat(3, 3, CV_64F);
		cvSetIdentity(scaleTransform_, cvRealScalar(1.0/scale_));
		cvmSet(scaleTransform_, 2, 2, 1.0);

		corners_ = cvCreateMat(3, 4, CV_64F);
		cvmSet(corners_, 0, 0, 0);
		cvmSet(corners_, 1, 0, 0);
		cvmSet(corners_, 2, 0, 1);
		cvmSet(corners_, 0, 1, frameSize_.width-1);
		cvmSet(corners_, 1, 1, 0);
		cvmSet(corners_, 2, 1, 1);
		cvmSet(corners_, 0, 2, frameSize_.width-1);
		cvmSet(corners_, 1, 2, frameSize_.height-1);
		cvmSet(corners_, 2, 2, 1);
		cvmSet(corners_, 0, 3, 0);
		cvmSet(corners_, 1, 3, frameSize_.height-1);
		cvmSet(corners_, 2, 3, 1);

		warpedCorners_ = cvCreateMat(3, 4, CV_64F);

		// Things we need for the search
		distance_ = new double * [scaledMosaicSize_.width];
		cost_ = new double * [scaledMosaicSize_.width];
		expanded_ = new bool * [scaledMosaicSize_.width];
		dx_ = new signed char * [scaledMosaicSize_.width];
		dy_ = new signed char * [scaledMosaicSize_.width];
		for (int x = 0; x < scaledMosaicSize_.width; x++) {
			distance_[x] = new double[scaledMosaicSize_.height];
			cost_[x] = new double[scaledMosaicSize_.height];
			expanded_[x] = new bool[scaledMosaicSize_.height];
			dx_[x] = new signed char[scaledMosaicSize_.height];
			dy_[x] = new signed char[scaledMosaicSize_.height];
		}
	}

	IplImage * DijkstraCutRenderer::renderImagesToMosaic(TransformSet *TS) {
		
		if (!TS || TS->size() == 0) {
			throw new GRCException("DijkstraCutRenderer::renderImagesToMosaic: Cannot render an empty TransformSet"); 
		}

		if (!initialised_) {
			initialiseStorage(TS);
			initialised_ = true;
		}

		IplImage *result = cvCreateImage(mosaicSize_, largeMosaic_->depth, largeMosaic_->nChannels);
		cvSetZero(result);

		//const double infinity = numeric_limits<double>::infinity();
		static int nNumTransforms=0, nNumTimesRendered = 0;
		nNumTransforms += TS->size(); nNumTimesRendered++;

		cout << TS->size() << " transforms " << nNumTimesRendered << " frames rendered " << nNumTransforms << " nNumTransforms total " << nNumTransforms/(double)nNumTimesRendered << " images per mosaic\n";

		if (TS->onlyRenderLatestTransform()) {
			// Update mosaic
			TS->mosaicTransform()->applyToImage(largeMosaic_, result);
			// Find latest frame and transform then add them to the mosaic
			const Transform *transform = TS->latestTransform()->transform();
			const IplImage * frame = imageSource_.getImage(TS->latestTransform()->id1());
			findAndApplyCut(frame, transform);
			
		} else {

			// Reset the mosaic and related images
			cvSetZero(mosaicMask_);
			cvSetZero(warpedMask_);
			//cvSetZero(mosaicEdge_);
			//cvSetZero(warpedEdge_);
			//cvSetZero(crossingMask_);
			cvSetZero(totalMask_); totalMask_->height = scaledMosaicSize_.height-1;
			cvSetZero(largeMask_);
			cvSetZero(largeMosaic_);
			cvSetZero(largeWarped_);

			// For each image in the mosaic
			//cout << TS->size() << " transforms\n";
			for (TransformSet::const_iterator iter = TS->begin(); iter != TS->end(); iter++) {
				// Get the image and related transform, then render it into the mosaic
				const IplImage * frame = imageSource_.getImage(iter->first.im1Id());

				const Transform * transform = iter->second->transform();

				findAndApplyCut(frame, transform);

			} // For each entry in the TransformSet
		}


		cvCopy(largeMosaic_, result);

		return result;
	}


	DijkstraCutRenderer::DijkstraCutRenderer(ImageSource &imageSource, TransformEngine &transformEngine, CvSize mosaicSize, unsigned int scale, EvaluationFunction *evaluationFunction)  
		: Renderer(imageSource, transformEngine, mosaicSize, evaluationFunction), initialised_(false), scale_(scale), frameSize_(), scaledFrameSize_()  {}

	DijkstraCutRenderer::~DijkstraCutRenderer() {
		if (initialised_) {
			cvReleaseImage(&mosaicMask_);
			cvReleaseImage(&warpedMask_);
			//cvReleaseImage(&mosaicEdge_);
			//cvReleaseImage(&warpedEdge_);
			//cvReleaseImage(&crossingMask_);
			cvReleaseImage(&totalMask_);
			cvReleaseImage(&largeMask_);
			cvReleaseImage(&largeMosaic_);
			cvReleaseImage(&largeWarped_);
			//cvReleaseImage(&checkMask_);
			//cvReleaseImage(&checkFrame_);

			cvReleaseMat(&scaleTransform_);
			cvReleaseMat(&corners_);
			cvReleaseMat(&warpedCorners_);

			for (int x = 0; x < scaledMosaicSize_.width; x++) {
				delete [] distance_[x];
				delete [] cost_[x];
				delete [] expanded_[x];
				delete [] dx_[x];
				delete [] dy_[x];
			}
			delete [] distance_;
			delete [] cost_;
			delete [] expanded_;
			delete [] dx_;
			delete [] dy_;
		}
	}

	DijkstraCutRenderer::PriorityQueue::PriorityQueue() : Q_() { }

	DijkstraCutRenderer::PriorityQueue::~PriorityQueue() { }

	void DijkstraCutRenderer::PriorityQueue::push(int x, int y, double oldDistance, double newDistance)
	{
		if(oldDistance != infinity)
		{
			TQueue::iterator iter = Q_.find(oldDistance);
			while (iter != Q_.end() && (iter->first == oldDistance)) {
				if ( (iter->second.first == x) && (iter->second.second ==y) ) {
					Q_.erase(iter);
					Q_.insert( make_pair<double, pair<int, int> >(newDistance, make_pair<int, int>(x,y)) );
					return;
				}
				iter++; //TODO: I think this is invalid now...
			}
		}
	
		Q_.insert( make_pair<double, pair<int, int> >(newDistance, make_pair<int, int>(x,y)) );
	
	}

	bool DijkstraCutRenderer::PriorityQueue::isEmpty() const {
		return Q_.empty();
	}

	void DijkstraCutRenderer::PriorityQueue::pop(int *x, int *y) {
		if (isEmpty()) {
			throw new GRCException("ERROR: DijkstraCutRenderer::PriorityQueue::pop - Attempt to pop off an empty Queue");
		}

		TQueue::iterator iter = Q_.begin();
		*x = iter->second.first;
		*y = iter->second.second;
		Q_.erase(iter);
	}

	void DijkstraCutRenderer::PriorityQueue::clear() {
		Q_.clear();
	}
};

pragma_warning(pop) // warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)
