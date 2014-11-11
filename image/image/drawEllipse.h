#ifndef DRAW_ELLIPSE_H
#define DRAW_ELLIPSE_H

#include "util/exception.h"
#include <Eigen/Core>
#include <opencv2/core/core.hpp>


/**
 * @brief Draws covariance ellipse with mean b, covariance A, at sigma s.d.'s on image
 * @param im
 * @param A
 * @param b
 * @param sigma
 * @param col
 */
void drawEllipse(cv::Mat & im, const Eigen::Matrix2d & A, const Eigen::Vector2d & b, const double sigma, const cv::Scalar col, const int thickness = 1, const bool bMarkCentre = true);

#endif //DRAW_ELLIPSE_H