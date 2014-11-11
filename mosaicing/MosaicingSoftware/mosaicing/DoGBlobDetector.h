#pragma once
#ifndef _GRC_DOG_BLOB_DETECTOR_
#define _GRC_DOG_BLOB_DETECTOR_

#include "FeatureDetector.h"
#include "boost/thread.hpp"

namespace grc {
	
	//! Find Difference-of-Gaussian blobs in a greyscale image
    class DoGBlobDetector : public FeatureDetector {
	public:
        /*! \param maxCorners Target maximum number of corners (stop detecting them when we've found this many)
           \param levels Number of Gaussian-smoothed images to compute (2 more than the number of scales we detect blobs at)
         */
		DoGBlobDetector(int maxCorners, int levels);
		~DoGBlobDetector();

		FeatureVector * findFeatures(IplImage *pGreyImage, int margin);

	private:
		inline double getVal(int id, int x, int y) 
        {
			return (double)blurredImageArray_[id]->imageData[y*blurredImageArray_[id]->widthStep + x];
		};
		
        template<class lessThan>
        inline float localMaximum(int id, int x, int y) 
        {
            lessThan lessThanFn;
            double val = getVal(id, x, y);

            //Bit of duplication to avoid comparison to itself, and to search the same image first
            for(int x_neighbor = x-1; x_neighbor <= x+1; x_neighbor++)
                for(int y_neighbor = y-1; y_neighbor <= y+1; y_neighbor++)
                {
                    if(y_neighbor == y && x_neighbor == x) y_neighbor++; //Avoid comparing to itself
                    if(!lessThanFn(val, getVal(id, x_neighbor, y_neighbor))) return false;
                }
            for(int im_neighbor = id-1; im_neighbor <= id+1; im_neighbor+=2)
                for(int x_neighbor = x-1; x_neighbor <= x+1; x_neighbor++)
                    for(int y_neighbor = y-1; y_neighbor <= y+1; y_neighbor++)
                        if(!lessThanFn(val, getVal(im_neighbor, x_neighbor, y_neighbor))) return false;

            return true;
        };

		const int levels_; //!< Number of Gaussian-smoothed images to compute (2 more than the number of scales we detect blobs at)
		bool isInitialised_;  //!< Flag when images have been allocated
		IplImage ** blurredImageArray_;  //!< Allocate blurred images once only
		const int maxCorners_; //!< Target maximum number of corners

        boost::mutex findFeatures_; //!< Lock to allow multiple concurrent users
    };
}

#endif //_GRC_DOG_BLOB_DETECTOR_
