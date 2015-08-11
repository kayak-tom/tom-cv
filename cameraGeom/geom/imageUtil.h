#ifndef IMAGEUTIL_H
#define    IMAGEUTIL_H

/* 
 * File:   imageUtil.h
 * Author: tom
 *
 * Created on 1 March 2012, 11:28 AM
 */
#include <opencv2/opencv.hpp>

//Output stats on the image
//void pp(const cv::Mat & M);
void dumpPixels(const cv::Mat & M); //Dump grey pixel values to cout
void dumpPixels(const std::string & filename);

class C2dPoint;
const cv::Point toCvPoint(const C2dPoint & p);

class C2dImagePointPx;

cv::Scalar getCol(int nId); //Returns a range of colours (mod 30)
/* For colour depth maps. dDepth should normally be in 0...1, but outside that is ok */
cv::Scalar getDepthCol(const double dDepth);

cv::Mat cropRect(const cv::Mat & image, const C2dImagePointPx & top1, const C2dImagePointPx & top2, const C2dImagePointPx & bottom2, const C2dImagePointPx & bottom1, const int nOutputRows, const int nOutputCols);
cv::Vec3b getSubPx_C3(const cv::Mat & M, const C2dImagePointPx & pt);
Eigen::Vector3d getVec(const cv::Mat & M, const int r, const int c);

cv::Scalar vecToScalar(const cv::Vec3b & v);
cv::Vec3b scalarToVec(const cv::Scalar & v);
cv::Scalar vector3dToScalar(const Eigen::Vector3d & v);
inline Eigen::Vector3d vec3bToVector(const cv::Vec3b & BGRb)
{
    return Eigen::Vector3d((double)BGRb[0],(double)BGRb[1],(double)BGRb[2]);
}

inline bool inIm(const cv::Point & p, const cv::Mat & im) { 
    return p.x>=0 && p.y >=0 && p.x < im.cols && p.y < im.rows;
}
inline bool inIm(const cv::Point & p, const cv::Mat & im, const int nMargin) { 
    return p.x>= nMargin && p.y >= nMargin && p.x < im.cols-nMargin && p.y < im.rows-nMargin;
}
/* 
 * Extra drawing functions--e.g. with safety checks or corrections for invalid arguments (these occur when projecting points very close to the camera sensor plane)
 * */

void drawText(cv::Mat & image, const cv::Point location, const C3dWorldPoint & toolTipend, const std::string label="");
void drawText(cv::Mat & image, const cv::Point location, const std::string label, const cv::Scalar col = CV_RGB(255,255,255), const int nSize=1, const bool bOutline=false);
void drawNumber(cv::Mat & image, const cv::Point location, const double dNum, const char * szUnit = "", const cv::Scalar col =  CV_RGB(255,255,255), const int nSize=1);
void drawLine(cv::Mat & image, const cv::Point & lastCentre, const cv::Point & pointOnBranch, const cv::Scalar & colToUse, int nThickness);
void drawCircle(cv::Mat & image, const cv::Point & centre, int nRadius, const cv::Scalar & colToUse, int nThickness=1);
void drawX(cv::Mat & image, const cv::Point & centre, const cv::Scalar & colToUse, const int nRadius=4, const int nThickness=2);
void drawPlus(cv::Mat & image, const cv::Point & centre, const cv::Scalar & colToUse, const int nRadius=4, const int nThickness=2);
void labelImagePointLocation(cv::Mat & image, const C2dImagePointPx & pointToDisplay, const std::string label);


void interpMax(const Eigen::VectorXd & aScales, double & dMax, double & dMagnitude, double & dVar);
bool quadraticMax(const double a, const double b, const double c, double & dMax, double & dMagnitude, double & dVar, const bool bTruncate);

cv::Mat getCircleStrel(const int erosion_size); //Strel with size cv::Size( 2*erosion_size + 1, 2*erosion_size+1 )

void addGaussianNoise(cv::Mat & M, const double dSigma);

//while keeping within image
void expandROI(cv::Rect & ROI, const int nExpansion, const cv::Mat & M);

int toCircleRad(const double dRad);
cv::Rect scaleRect(const cv::Rect & rect_in, const double dScale);

void imageToTSV(const std::string filename, const std::string label);
void imageToTSV_float(const cv::Mat & M, const std::string label);
void imageToTSV_uchar(const cv::Mat & M, const std::string label);

//For finding wire pixels
template<typename T> void nonMaxSuppression_vertical(cv::Mat & M, const int nStartRow, const int nEndRow)
{
    if(nStartRow==0) 
        M.row(0).setTo(cv::Scalar(0));
    
    if(nEndRow==M.rows)
        M.row(M.rows-1).setTo(cv::Scalar(0));
     
    const int nEnd = std::min<int>(nEndRow, M.rows-1);
    
    for(int r=std::max<int>(nStartRow, 1); r < nEnd; r++)
        for(int c=0; c < M.cols; c++)
        {
            T & val = M.at<T>(r,c);
            if(val <= M.at<T>(r+1,c) || val <= M.at<T>(r-1,c))
                val = 0;
        }
}

#endif    /* IMAGEUTIL_H */

