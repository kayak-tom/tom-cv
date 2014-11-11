#include "cvGeom.h"

namespace grc {

CvPoint2D64f operator+(const CvPoint2D64f &p1, const CvPoint2D64f &p2)
{
    return cvPoint2D64f(p1.x+p2.x, p1.y+p2.y);
}

CvPoint2D64f operator*(const CvPoint2D64f &p1, double dScale)
{
    return cvPoint2D64f(p1.x*dScale, p1.y*dScale);
}

bool operator==(const CvPoint2D64f &p1, const CvPoint2D64f &p2)
{
    return p1.x==p2.x && p1.y==p2.y;
}

double SSD(const CvPoint2D64f &p1, const CvPoint2D64f &p2)
{
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return (dx*dx) + (dy*dy);
}

CvPoint2D32f operator+(const CvPoint2D32f &p1, const CvPoint2D32f &p2)
{
    return cvPoint2D32f(p1.x+p2.x, p1.y+p2.y);
}

CvPoint2D32f operator*(const CvPoint2D32f &p1, double dScale)
{
    return cvPoint2D32f(p1.x*dScale, p1.y*dScale);
}

bool operator==(const CvPoint2D32f &p1, const CvPoint2D32f &p2)
{
    return p1.x==p2.x && p1.y==p2.y;
}

//! Get 4 corners of source image that can be transformed to give position of source in dest image.
CvMat * getSourceCorners(const CvSize & imageSize_)
{
    CvMat * sourceCorners = cvCreateMat(3, 4, CV_64FC1);
    for(int i=0; i<4; i++)
    {
        cvmSet(sourceCorners, 0, i, imageSize_.width * (i%2==0));
        cvmSet(sourceCorners, 1, i, imageSize_.height * (i>1));
        cvmSet(sourceCorners, 2, i, 1);
    }

    return sourceCorners;
}

}
