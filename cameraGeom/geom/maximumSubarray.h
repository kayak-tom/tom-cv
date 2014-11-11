#ifndef MAXIMUMSUBARRAY_H
#define MAXIMUMSUBARRAY_H

#include <opencv2/opencv.hpp>

void findMaximalSubarray2d(const cv::Mat & image_in, const int nMean, cv::Rect & BB_best);
int percentile(const cv::Mat & M, const int nPercentile);

#endif // MAXIMUMSUBARRAY_H
