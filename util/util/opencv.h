#pragma once
/* 
 * File:   opencv.h
 * Author: tom
 *
 * Created on 27 January 2011, 13:31
 */

#ifndef OPENCV_H
#define	OPENCV_H

//#define USE_OLD_OPENCV

#ifdef USE_OLD_OPENCV
#include "cv.h"
#else //Version 2.2 or later
#include "opencv2/core/core.hpp"
#include "opencv2/core/core_c.h"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/imgproc/imgproc_c.h"
#endif

#endif	/* OPENCV_H */

