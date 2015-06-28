#pragma once
/* 
 * File:   opencv.h
 * Author: tom
 *
 * Created on 27 January 2011, 13:31
 */

#ifndef OPENCV_HIGHGUI_H
#define    OPENCV_HIGHGUI_H

#include "opencv.h"

#ifdef USE_OLD_OPENCV
#include "highgui.h"
#else //Version 2.2 or later
#include "opencv2/highgui/highgui_c.h"
#include "opencv2/highgui/highgui.hpp"
#endif

#endif    /* OPENCV_HIGHGUI_H */

