#pragma once
#ifndef _WARP
#define _WARP
#include "util/opencv.h"

//! Perspective warp implemented with faster version of various internal OpenCV warping functions
/*! Nearest-neighbour interpolation supported. A bounding box limits the area we iterate over. Division 
 *  is implemented with a binomial expansion. Floating-point maths done with doubles (faster than floats). Repeated arithmatic, 
 *  filling of outliers, and unnecessary logic removed. 
 * \param srcarray Source image--1 or 3 channel and 8-bit.
 * \param dstarray Destination image--1 or 3 channel and 8-bit.
 * \param matrix Perspective transformation matrix mapping srcarray into dstarray
 * \param flags Interpolation method: CV_INTER_NN or CV_INTER_LINEAR
*/
void WarpPerspective( const CvArr* srcarr, CvArr* dstarr,
                   const CvMat* matrix, int flags );

//! Affine warp implemented with faster version of various internal OpenCV warping functions
/*! Nearest-neighbour interpolation supported. A bounding box limits the area we iterate over. 
 *  Floating-point maths done with doubles (faster than floats). Repeated arithmatic, 
 *  filling of outliers, and unnecessary logic removed. 
 * \param srcarray Source image--1 or 3 channel and 8-bit.
 * \param dstarray Destination image--1 or 3 channel and 8-bit.
 * \param matrix 2x3 Affine transformation matrix mapping srcarray into dstarray
 * \param flags Interpolation method: CV_INTER_NN or CV_INTER_LINEAR
*/
void
WarpAffine( const CvArr* srcarr, CvArr* dstarr, const CvMat* matrix,
              int flags );

class BB
{
public:
	int minX, minY, maxX, maxY; //Todo protect...
	BB(const int minX, const int minY, const int maxX, const int maxY)
	  : minX(minX), minY(minY), maxX(maxX), maxY(maxY)
	{}

	BB() : minX(-1), minY(-1), maxX(-1), maxY(-1)
	{}

	void operator=(const BB & bb)
	{
		minX = bb.minX;
		maxX = bb.maxX;
		minY = bb.minY;
		maxY = bb.maxY;
	}

	void split(BB & firstBB, BB & secondBB) const // For threading
	{
        int midY = minY + (maxY - minY)/2;

        firstBB = *this; firstBB.maxY = midY;
		secondBB = *this; secondBB.minY = midY;
	}
};
const BB getBB(const double * inv_mat, const CvSize ssize, const CvSize dsize);

#endif
