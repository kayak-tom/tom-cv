/*! \file cvGeom.h \brief Simple vector addition/multiplication operators for CvPoint2D??f. */
#pragma once
#ifndef _CVGEOM
#define _CVGEOM

#include "util/opencv.h"

namespace grc {

CvPoint2D64f operator+(const CvPoint2D64f &p1, const CvPoint2D64f &p2);

CvPoint2D64f operator*(const CvPoint2D64f &p1, double dScale);

bool operator==(const CvPoint2D64f &p1, const CvPoint2D64f &p2);

double SSD(const CvPoint2D64f &p1, const CvPoint2D64f &p2);

CvPoint2D32f operator+(const CvPoint2D32f &p1, const CvPoint2D32f &p2);

CvPoint2D32f operator*(const CvPoint2D32f &p1, double dScale);

bool operator==(const CvPoint2D32f &p1, const CvPoint2D32f &p2);

//! Get 4 corners of source image that can be transformed to give position of source in dest image.
CvMat * getSourceCorners(const CvSize & imageSize);

}

#endif // _CVGEOM
