/*
 * File:   imageUtil.cpp
 * Author: tom
 *
 * Created on 1 March 2012, 11:28 AM
 */

#include "newCamera.h"
#include <opencv2/highgui/highgui.hpp>
#include <iomanip>
//#include "geometry.h"
//#include "colourLookup.h"
#include <opencv2/core/core_c.h>
//#include "pruningUI.h"
#include "lines.h"
//#include "probImage.h"
#include "imageUtil.h"
#include <util/random.h>
#include <util/pp.h>
#include <fstream>

Eigen::Vector3d getPx(const cv::Mat & M, const int x, const int y)
{
    cv::Vec3b BGR = M.at<cv::Vec3b>(y,x);
    return Eigen::Vector3d(BGR[0], BGR[1], BGR[2]);
}

Eigen::Vector3d getSubPx_intx(const cv::Mat & M, const int nx, const double y)
{
    const int ny_lo = intFloor(y);
    const int ny_hi = ny_lo + 1;

    const double y_lo = (double) ny_lo;

    if(IS_DEBUG) CHECK(y_lo > y, "floor failed");

    double dWeightHi = y - y_lo;
    double dWeightLo = 1 - dWeightHi;

    if(IS_DEBUG) CHECK(dWeightHi < 0 || dWeightHi > 1, "Weight OOB");
    if(IS_DEBUG) CHECK(dWeightLo < 0 || dWeightLo > 1, "Weight OOB");

    return dWeightLo * getPx(M, nx, ny_lo) + dWeightHi * getPx(M, nx, ny_hi);
}

Eigen::Vector3d getSubPx(const cv::Mat & M, const double x, const double y)
{
    const int nx_lo = intFloor(x);
    const int nx_hi = nx_lo + 1;

    const double x_lo = (double) nx_lo;

    double dWeightHi = x - x_lo;
    double dWeightLo = 1 - dWeightHi;

    if(IS_DEBUG) CHECK(dWeightHi < 0 || dWeightHi > 1, "Weight OOB");
    if(IS_DEBUG) CHECK(dWeightLo < 0 || dWeightLo > 1, "Weight OOB");

    return dWeightLo * getSubPx_intx(M, nx_lo, y) + dWeightHi * getSubPx_intx(M, nx_hi, y);
}

cv::Vec3b getSubPx_C3(const cv::Mat & M, const C2dImagePointPx & pt)
{
    Eigen::Vector3d BGR = getSubPx(M, pt.x(), pt.y());
    return cv::Vec3b(cv::saturate_cast<uchar>(BGR(0)),cv::saturate_cast<uchar>(BGR(1)),cv::saturate_cast<uchar>(BGR(2)));
}

cv::Mat cropRect(const cv::Mat & image, const C2dImagePointPx & top1, const C2dImagePointPx & top2, const C2dImagePointPx & bottom2, const C2dImagePointPx & bottom1, const int nOutputRows, const int nOutputCols)
{

    cv::Mat outMat(nOutputRows, nOutputCols, CV_8UC3);

    for(int r=0; r<nOutputRows; r++) {
        C2dImagePointPx rowStart = (bottom1 + (bottom2-bottom1)*(r/(double)nOutputRows)).eval();
        C2dImagePointPx rowEnd = (top1 + (top2-top1)*(r/(double)nOutputRows)).eval();
        for(int c=0; c<nOutputCols; c++) {
            C2dImagePointPx rowPos = (rowStart + (rowEnd-rowStart)*(c/(double)nOutputCols)).eval();
            outMat.at<cv::Vec3b>(r, c) = getSubPx_C3(image, rowPos);
        }
    }

    return outMat;
}

int val(int & nShift)
{
    nShift++;
    if (nShift == 5)
        nShift = 1;

    return 50 + (205 * nShift) / 5;
}

cv::Scalar getCol(int nId)
{
    int nShift = nId % 5;
    nId %= 6;
    nId++; //3 bits, 1 or 2 set
    int R = (nId & 1) ? val(nShift) : 0;
    int G = (nId & 2) ? val(nShift) : 0;
    int B = (nId & 4) ? val(nShift) : 0;
    int m = std::max(std::max(R, G), B);
    m = 255 - m;

    return CV_RGB(R + m, G + m, B + m);
}

void dumpPixels(const cv::Mat & M)
{
    for(int r=0; r<M.rows; r++) {
        for(int c=r % 3; c<M.cols; c+=3) { //1/3 pixels
            int val=(int)M.at<uchar>(r,c);
            if(val < 250)
                cout << val << '\n';
        }
    }
    cout << endl;
}

void dumpPixels(const std::string & filename)
{
    cv::Mat M = cv::imread(filename, cv::IMREAD_GRAYSCALE);
    dumpPixels(M);
}

Eigen::Vector3d getVec(const cv::Mat & M, const int r, const int c)
{
    switch(M.type()) {
    case CV_8UC3: {
        cv::Vec3b BGRb = M.at<cv::Vec3b>(r,c);
        return vec3bToVector(BGRb);
    }
    case CV_32FC3: {
        cv::Vec3f BGRf = M.at<cv::Vec3f>(r,c);
        return Eigen::Vector3d(BGRf[0],BGRf[1],BGRf[2]);
    }
    case CV_64FC3: {
        cv::Vec3d BGRd = M.at<cv::Vec3d>(r,c);
        return Eigen::Vector3d(BGRd[0],BGRd[1],BGRd[2]);
    }
    default: {
        THROW("Not a C3 colour image");
    }
    }

}

cv::Scalar vecToScalar(const cv::Vec3b & v) //TODO: these should be used everywhere... Search [0].*[1]
{
    return cv::Scalar( v[0], v[1], v[2] );
}

cv::Vec3b scalarToVec(const cv::Scalar & v)
{
    return cv::Vec3b( cv::saturate_cast<uchar>(v[0]), cv::saturate_cast<uchar>(v[1]), cv::saturate_cast<uchar>(v[2]) );
}

cv::Scalar vector3dToScalar(const Eigen::Vector3d & v)
{
    return cv::Scalar(v(0), v(1), v(2));
}


void drawText(cv::Mat & image, const cv::Point location, const std::string label, const cv::Scalar colour, const int nSize, const bool bOutline)
{
    if(bOutline)
        cv::putText(image, label, location, CV_FONT_HERSHEY_PLAIN, 1.2*nSize, cv::Scalar(), nSize+6);
    cv::putText(image, label, location, CV_FONT_HERSHEY_PLAIN, 1.2*nSize, colour, nSize);
}

void drawText(cv::Mat & image, const cv::Point location, const C3dWorldPoint & pointToDisplay, const std::string label)
{
    std::ostringstream ss;
    ss << std::setprecision(3) << label << " (" << pointToDisplay.x() << ", " << pointToDisplay.y() << ", " << pointToDisplay.z() << ")";
    drawText(image, location, ss.str(), CV_RGB(255,255,255), 2);
}

void labelImagePointLocation(cv::Mat & image, const C2dImagePointPx & pointToDisplay, const std::string label)
{
    std::ostringstream ss;
    ss << std::setprecision(3) << "  " << label << " (" << pointToDisplay.x() << ", " << pointToDisplay.y() << ")";
    drawText(image, pointToDisplay, ss.str(), CV_RGB(255,255,255), 3);
}

void drawNumber(cv::Mat & image, const cv::Point location, const double dNum, const char * szUnit, const cv::Scalar colour, const int nSize)
{
    std::ostringstream ss;
    ss << std::setprecision(2) << dNum << szUnit;
    drawText(image, location, ss.str(), colour, nSize);
}

void drawLine(cv::Mat & image, const cv::Point & lastCentre, const cv::Point & pointOnBranch, const cv::Scalar & colToUse, int nThickness)
{
    const C2dImagePointPx p1=lastCentre;
    const C2dImagePointPx p2=pointOnBranch;
    double dDistLineToOrigin = -1;
    if(p1 != p2) {
        const C2dBoundedLine lineBounded(p1, p2);
        dDistLineToOrigin = lineBounded.closestDistance(C2dImagePointPx::Zero());
    }

    const double dMaxSize=image.rows+image.cols;
    if(fabs(dDistLineToOrigin) > dMaxSize || p1.norm() > 10*dMaxSize || p2.norm() > 10*dMaxSize) {
        return; //Line doesn't pass near origin
    }

    if(nThickness > 255) {
        REPEAT(10, cout << "Warning, line too thick: " << nThickness << " pixels. Object too close to camera. Resetting to 3" << endl);
        nThickness = 3;
    } else if (nThickness<1)
        nThickness = 1;

    cv::line(image, lastCentre, pointOnBranch, colToUse, nThickness);
}
int toCircleRad(const double dRad)
{
    return clip<int>(cvRound(dRad), 1, 255);
}

void drawCircle(cv::Mat & image, const cv::Point & centre, int nRadius, const cv::Scalar & colToUse, int nThickness)
{
    const C2dImagePointPx p=centre;
    const double dLength = p.squaredNorm();
    if(dLength > sqr(10000))
        return; //Centre too far away

    if(nRadius < 0) {
        REPEAT(10, cout << "Negative circle radius: " << nRadius << endl;);
        nRadius = 1;
    }

    if((nThickness < 1 && nThickness != -1) || nThickness > 255) {
        REPEAT(10, cout << "Circle thickness too low or high: " << nThickness << endl;);
        nThickness = 3;
    }

    cv::circle(image, centre, nRadius, colToUse, nThickness);
}

void drawX(cv::Mat & image, const cv::Point & centre, const cv::Scalar & colToUse, const int nRadius, const int nThickness)
{
    cv::Point aPointOffsets[4];
    aPointOffsets[0] = cv::Point(-nRadius,-nRadius);
    aPointOffsets[1] = cv::Point(-nRadius,nRadius);
    aPointOffsets[2] = cv::Point(nRadius,-nRadius);
    aPointOffsets[3] = cv::Point(nRadius,nRadius);

    for(int i=0; i<4; i++) {
        cv::line(image, centre, centre+aPointOffsets[i], colToUse, nThickness);
    }
}

void drawPlus(cv::Mat & image, const cv::Point & centre, const cv::Scalar & colToUse, const int nRadius, const int nThickness)
{
    cv::Point aPointOffsets[4];
    aPointOffsets[0] = cv::Point(0,-nRadius);
    aPointOffsets[1] = cv::Point(0,nRadius);
    aPointOffsets[2] = cv::Point(nRadius,0);
    aPointOffsets[3] = cv::Point(-nRadius,0);

    for(int i=0; i<4; i++) {
        cv::line(image, centre, centre+aPointOffsets[i], colToUse, nThickness);
    }
}

/*Fit a quadratic to a,b,c and return the position of its maximum. If unbounded then return false. dVar is a heuristic (TODO)*/
bool quadraticMax(const double a, const double b, const double c, double & dMax, double & dMagnitude, double & dVar, const bool bTruncate)
{
    dVar = HUGE;

    //f(-1) = a, f(0)=b, f(1)=c, f(x) = A*x^2+B*x+C
    const double C=b;
    const double A=0.5*(a-2*b+c);
    const double B=A+b-a;
    if(IS_DEBUG) CHECK(!zero(A+B+C-c), "math error");

    if(A==0) { //Flat line (e.g. simulated images sometimes)
        dMagnitude = b;
        dMax = 0;

        return false;
    }
    if(A>0) {
        //Take an endpoint as a maximum
        dMagnitude = std::max<double>(a,c);
        dMax = (dMagnitude == a) ? -1 : 1;

        return false;
    }

    //Df(x) = 2Ax+B -> max at x=-B/2A
    dMax = -B/(2*A);

    //We can;t reasonably interpolate further than 1
    const double TRUNCATE_MAX = 2;
    if(bTruncate && fabs(dMax) > TRUNCATE_MAX)
        dMax = TRUNCATE_MAX*sign(dMax);

    dMagnitude = A*sqr(dMax) + B*dMax + C;

    if(IS_DEBUG) CHECK(dMagnitude > 2*max3(a,b,c)-min3(a,b,c), "Magnitude extrapolated too high");

    if(IS_DEBUG && fabs(dMax) < TRUNCATE_MAX) {
        const double dTemp = A*sqr(dMax+0.001) + B*(dMax+0.001) + C;
        CHECK(dTemp > dMagnitude, "Find extrema failed");
    }

    if(A<0 && dMagnitude > 0)
        dVar = -dMagnitude/A;

    return true;
}

/* Approximate about maximum with a quadratic, return subpixel max. dMax is the (subpixel) index into aScales */
void interpMax(const Eigen::VectorXd & aScales, double & dMax, double & dMagnitude, double & dVar)
{
    int nMax=-1;
    aScales.maxCoeff(&nMax);
    int nInterpMiddle = nMax;
    if(nMax==0)
        nInterpMiddle++;
    if(nMax==aScales.size()-1)
        nInterpMiddle--;

    quadraticMax(aScales(nInterpMiddle-1), aScales(nInterpMiddle), aScales(nInterpMiddle+1), dMax, dMagnitude, dVar, true);
    dMax += nInterpMiddle;

    //REPEAT(100000, if(rand() % 100 == 0) {cout << aScales.transpose() << " has max " << dMagnitude << ", var " << dVar << " at index " << dMax << endl;});
}

cv::Mat getCircleStrel(const int erosion_size)
{
    return cv::getStructuringElement( cv::MORPH_ELLIPSE,
                                      cv::Size( 2*erosion_size + 1, 2*erosion_size+1 ),
                                      cv::Point( erosion_size, erosion_size ) );
}

void addGaussianNoise(cv::Mat & M, const double dSigma)
{
    CRandom::fast_srand((unsigned int)cv::sum(M.row(0))[0]); //make it more deterministic

    for(int r=0;r<M.rows;r++)
        for(int c=0;c<M.cols;c++)
            switch(M.type())
            {
                case CV_8UC1:
                    M.at<uchar>(r,c) = cv::saturate_cast<uchar>((int)M.at<uchar>(r,c) + (int)CRandom::Normal(0, dSigma));
                    break;
                case CV_8UC3:
                    for(int i=0;i<3;i++)
                    {
                        M.at<cv::Vec3b>(r,c)[i] = cv::saturate_cast<uchar>((int)M.at<cv::Vec3b>(r,c)[i] + (int)CRandom::Normal(0, dSigma));
                    }
                    break;
                case CV_32FC1:
                    M.at<float>(r,c) += (float)CRandom::Normal(0, dSigma);
                    break;
                default:
                    THROW("Unhandled image type");
            }
}
void expandROI(cv::Rect & ROI, const int nExpansion, const cv::Mat & M)
{
    ROI.x = std::max<int>(0, ROI.x-nExpansion);
    ROI.y = std::max<int>(0, ROI.y-nExpansion);
    ROI.width = std::min<int>(M.cols - ROI.x, ROI.width+2*nExpansion);
    ROI.height = std::min<int>(M.rows - ROI.y, ROI.height+2*nExpansion);
}

cv::Rect scaleRect(const cv::Rect & rect_in, const double dScale)
{
    return cv::Rect(cvRound(rect_in.x*dScale), cvRound(rect_in.y*dScale), cvRound(rect_in.width*dScale), cvRound(rect_in.height*dScale));
}

void imageToTSV_float(const cv::Mat & M, const std::string label)
{
    const std::string filename = label + ".tsv";
    std::ofstream tsv(filename.c_str()); 
    
    for(int r=0;r<M.rows;r++)
    {
        for(int c=0;c<M.cols;c++)
        {
            tsv << M.at<float>(r,c) << '\t';
        }
        tsv << endl;
    }
}

void imageToTSV_8UC3(const cv::Mat & M, const std::string label)
{
    const bool bVerbose = true;

    const std::string filename = label + ".tsv";
    const std::string filename_b = "blue-" + filename;
    const std::string filename_g = "green-" + filename;
    const std::string filename_r = "red-" + filename;
    
    std::ofstream aChannels[3]; 
    aChannels[0].open(filename_b.c_str());
    aChannels[1].open(filename_g.c_str());
    aChannels[2].open(filename_r.c_str());
    
    for(int r=0;r<M.rows;r++)
    {
        for(int c=0;c<M.cols;c++)
        {
            cv::Vec3b px = M.at<cv::Vec3b>(r,c);
            for(int nChannel=0;nChannel<3;nChannel++)
            {
                aChannels[nChannel] << (int)px[nChannel] << '\t';
            }
        }
        for(int nChannel=0;nChannel<3;nChannel++)
        {
            aChannels[nChannel] << endl;
        }        
    }
    
    cv::Scalar mean = cv::mean(M);
    COUT(mean);
    
    M *= 100/mean[1]; //normalise
    cv::Mat Mat, med;
    cv::medianBlur(M, med, 3);
    
    cv::absdiff(M, med, M);
    COUT(cv::sum(M)/M.size().area());
}

void imageToTSV(const std::string filename, const std::string label)
{
    cv::Mat M = cv::imread(filename);
    imageToTSV_8UC3(M, label);
}

