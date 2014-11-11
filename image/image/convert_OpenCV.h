/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * convert_OpenCV.h
 *
 *  Created on: 12/10/2009
 *      Author: tom
 */

#ifndef CONVERT_OPENCV_H_
#define CONVERT_OPENCV_H_

#include "util/exception.h"
#include "image/HSI.h"
#include "util/opencv.h"
#include "util/location.h"
#include <boost/noncopyable.hpp>
#include <iostream>
#include "imageAccess.h"

struct CvPoint;
struct CvPoint2D32f;
struct CvPoint2D64f;

inline CvPoint operator+(const CvPoint &p1, int n)
{
    return cvPoint(p1.x+n, p1.y+n);
};

inline CvPoint locToCvPoint (CLocation const &Loc)
{
    return cvPoint(Loc.x(), Loc.y());
}

//Don't use templates cos then matches loads of other classes using same operators:
#define MAKE_OPENCV_OPERATORS( TCvPoint )\
    inline bool operator==(const TCvPoint & p1, const TCvPoint & p2) { return p1.x==p2.x && p1.y == p2.y; }\
    inline bool operator!=(const TCvPoint & p1, const TCvPoint & p2) { return !(p1 == p2); }\
    inline TCvPoint operator+(const TCvPoint & p, const TCvPoint & p2) { \
        TCvPoint out; \
        out.x=p.x+p2.x; \
        out.y=p.y+p2.y; \
        return out;}\
    inline TCvPoint operator-(const TCvPoint & p, const TCvPoint & p2) { \
        TCvPoint out; \
        out.x=p.x-p2.x; \
        out.y=p.y-p2.y; \
        return out;}\
    template<typename TNum> inline void operator*=(TCvPoint & p, TNum d) { p=p*d; }\
    inline void operator+=(TCvPoint & p, const TCvPoint & d) { p=p+d; }\
    inline double sum_square(const TCvPoint & p) { return sqr(p.x) + sqr(p.y); }\
    inline double length(const TCvPoint & p) { return sqrt(sum_square(p)); }\
    template<typename TNum> inline TNum SSD(const TCvPoint &p1, const TCvPoint &p2)\
    {    return (sqr(p1.x-p2.x) + sqr(p1.y-p2.y)); }

MAKE_OPENCV_OPERATORS(CvPoint)
MAKE_OPENCV_OPERATORS(CvPoint2D32f)
MAKE_OPENCV_OPERATORS(CvPoint2D64f)

inline CvPoint2D32f operator*(const CvPoint2D32f & p, const float d)
{
    CvPoint2D32f out;
    out.x=p.x*d;
    out.y=p.y*d;
    return out;
}

inline CvPoint2D64f operator*(const CvPoint2D64f & p, const double d)
{
    CvPoint2D64f out;
    out.x=p.x*d;
    out.y=p.y*d;
    return out;
}


template<typename TCvPoint>
bool inImage(const TCvPoint & p, const CvSize & s)
{
    return p.x>=0 && p.y >=0 && p.x <= s.width && p.y <= s.height;
}

template<typename TCvPoint>
inline void normalise(TCvPoint & p)
{
    double dLength_inv = (1.0/length(p));
    p.x *= dLength_inv;
    p.y *= dLength_inv;
}

inline CvPoint cvPointFrom64f(const CvPoint2D64f p)
{
    return cvPoint(doubleToInt(p.x), doubleToInt(p.y));
}

//This is now redundant--TODO: Upgrade to cv::Mat throughout
template<class T>
class CvPtr : boost::noncopyable //smart ptr
{
protected:
    T * ptr;
public:
    CvPtr() : ptr(0) {}
    CvPtr(T * p) : ptr(p) {}

    ~CvPtr(); //MUST PROVIDE SPECIALISATION FOR THIS

    void operator=(T * p) {
        if(IS_DEBUG) CHECK(!p, "No CvMat/image supplied");
        if(IS_DEBUG) CHECK(ptr, "CvMat/Image already supplied");
        ptr=p;
    }
    operator const T * () const {
        return ptr;
    }
    operator T * () {
        return ptr;
    }
    T * operator->() {
        return (ptr);
    }
    const T * operator->() const {
        return (ptr);
    }

    void setRGB(int x, int y, uchar R, uchar G, uchar B) HARD_INLINE {
        CIplPx<uchar>::setRGB(ptr, x, y, R, G, B);
    }
    void setRed(int x, int y, uchar val) HARD_INLINE {
        CIplPx<uchar>::setRed(ptr, x, y, val);
    }
    void setGreen(int x, int y, uchar val) HARD_INLINE {
        CIplPx<uchar>::setGreen(ptr, x, y, val);
    }
    void setBlue(int x, int y, uchar val) HARD_INLINE {
        CIplPx<uchar>::setBlue(ptr, x, y, val);
    }
    void setGrey(int x, int y, uchar val) HARD_INLINE {
        CIplPx<uchar>::setGrey(ptr, x, y, val);
    }

    uchar getRed(int x, int y) const PURE_FN HARD_INLINE {
        return CIplPx<uchar>::getRed(ptr, x, y);
    }

    uchar getGreen(int x, int y) const PURE_FN HARD_INLINE {
        return CIplPx<uchar>::getGreen(ptr, x, y);
    }

    uchar getBlue(int x, int y) const PURE_FN HARD_INLINE {
        return CIplPx<uchar>::getBlue(ptr, x, y);
    }

    void getRGB(int x, int y, CRedGreenBlue & rgb) const HARD_INLINE {
        CIplPx<uchar>::getRGB(ptr, x, y, rgb);
    }

    uchar getPx(int x, int y, int channel) const {
        return CIplPx<uchar>::getPx(ptr, x, y, channel);
    }

    template <typename TProbType>
    double getSubPx(int x, double y) const {
        return CIplPx<uchar>::getSubPx<TProbType>(ptr, x, y, 0);
    }

    template <typename TProbType>
    double getSubPx(double x, double y) const {
        return CIplPx<uchar>::getSubPx<TProbType>(ptr, x, y, 0);
    }
    uchar getGrey(int x, int y) const HARD_INLINE {
        return CIplPx<uchar>::getGrey(ptr, x, y);
    }
};

template<int nRows, int nCols>
class CvMatFixed : boost::noncopyable //smart ptr
{
    CvMat mat;
    double data[nRows*nCols];
public:
    CvMatFixed() {
        mat = cvMat(nRows, nCols, CV_64FC1, data);
    }

    operator CvMat * () {
        return &mat;
    }
    operator const CvMat * () const {
        return &mat;
    }
    const CvMat * operator->() const {
        return (&mat);
    }

    double operator()(int r, int c) const {
        return cvmGet(&mat, r, c);
    }
};

template<typename T>
class CIplImIt
{
    const IplImage * pIm;
    T * pData;
    int x,y;
public:
    int getX() const {
        return x;
    }
    int getY() const {
        return y;
    }
    CvPoint getCvPoint() {
        return cvPoint(x,y);
    }

    CIplImIt(const IplImage * pIm) : pIm(pIm), pData(reinterpret_cast<T *>(pIm->imageData)), x(0), y(0) {
    }

    void reset() {
        *this = CIplImIt(pIm);
    }

    T & operator*() {
        return *pData;
    }
    T & operator[](int n) {
        return pData[n];    //TODO Change to RGB
    }
    void setBGR(uchar b, uchar g, uchar r) {
        if(IS_DEBUG) CHECK(pIm->nChannels != 3 && pIm->nChannels != 4, "setBGR: Image not colour");
        pData[0]=b;
        pData[1]=g;
        pData[2]=r;
    }

    //Returns 'false' at end if input
    bool operator++() {
        pData += pIm->nChannels;
        x++;

        if(x==pIm->width) {
            x=0;
            y++;
            if(y==pIm->height) {
                pData=0;
                return false;
            }

            pData = reinterpret_cast<T *>(pIm->imageData + pIm->widthStep*y);
        }
        return true;
    }
    bool operator++(int) {
        CIplImIt<T> & me=*this;
        return ++me;
    }

    operator bool () const {
        return pData != 0;
    }
};

class CGreyscaler : boost::noncopyable
{
    static const int LOOKUP_SIZE=1000;
    int anLookup[LOOKUP_SIZE+1];
    const bool init;
    const int coeffR, coeffG, coeffB, TOTAL;
    //const double gamma;

public:
    CGreyscaler(bool init, int coeffR, int coeffG, int coeffB, double gamma) : init(init), coeffR(coeffR), coeffG(coeffG), coeffB(coeffB), TOTAL(coeffR + coeffG + coeffB) { /*, gamma(gamma)*/
        if(init) {
            if(IS_DEBUG) CHECK(gamma<=0 || gamma > 1000, "CGreyscaler: Bad parameters" );
            for(int i=0; i<= LOOKUP_SIZE; i++) {
                double dValIn01 = (double)i/(double)LOOKUP_SIZE;
                if(gamma != 1)
                    dValIn01 = pow(dValIn01, gamma);
                anLookup[i] = doubleToInt(255.0*dValIn01);
                if(i>0 && anLookup[i] > anLookup[i-1]+1) {
                    std::cout << "Warning: losing GS conversion precision...\n";
                    std::cout << "anLookup[" << i-1 << "]=" << anLookup[i-1] << std::endl;
                    std::cout << "anLookup[" << i << "]=" << anLookup[i] << std::endl;
                }
            }
            if(IS_DEBUG) CHECK(anLookup[0] != 0 || anLookup[LOOKUP_SIZE] != 255, "CGreyscaler: Failed");
        }
    }

    void greyScale2(const IplImage * pBGR, IplImage * pGrey) const {
        if(IS_DEBUG) CHECK(!init, "CGreyscaler not initialised--must initialise if using colour images");

        if(IS_DEBUG) CHECK(pBGR->nChannels != 3 && pBGR->nChannels != 4, "Colour input image required");
        if(IS_DEBUG) CHECK(pGrey->nChannels != 1, "Greyscale output image required");

        if(IS_DEBUG) CHECK(pBGR->height != pGrey->height || pBGR->width != pGrey->width, "Size mismatch");

        CIplImIt<const uchar> colourIt(pBGR);
        for(CIplImIt<uchar /* type reused in cast 6 lines down... */> greyIt(pGrey); greyIt; greyIt++, colourIt++) {
            const int b=colourIt[0];
            const int g=colourIt[1];
            const int r=colourIt[2];
            const int weightedSum = ((coeffB*b+coeffG*g+coeffR*r))/TOTAL;
            if(IS_DEBUG) CHECK(weightedSum >= 0 && weightedSum < 256, "Conversion failed");
            *greyIt = (uchar)weightedSum;//(uchar)anLookup[weightedSum];
        }
    }

    void greyScale(const IplImage * pBGR, IplImage * pGrey) const {
        if(IS_DEBUG) CHECK(!init, "CGreyscaler not initialised--must initialise if using colour images");

        if(IS_DEBUG) CHECK(pBGR->nChannels != 3 && pBGR->nChannels != 4, "Colour input image required");
        if(IS_DEBUG) CHECK(pGrey->nChannels != 1, "Greyscale output image required");

        if(IS_DEBUG) CHECK(pBGR->height != pGrey->height || pBGR->width != pGrey->width, "Size mismatch");

        uchar * pSrcRow = (uchar *)pBGR->imageData;
        uchar * pDstRow = (uchar *)pGrey->imageData;
        for(int y=pBGR->height; y>0; y--, pSrcRow += pBGR->widthStep,  pDstRow += pGrey->widthStep) {
            uchar * pSrc = pSrcRow;
            uchar * pDst = pDstRow;
            for(int x=pBGR->width; x>0; x--, pSrc += pBGR->nChannels, pDst++) {
                const int b=pSrc[0];
                const int g=pSrc[1];
                const int r=pSrc[2];
                const int weightedSum = ((coeffB*b+coeffG*g+coeffR*r) * LOOKUP_SIZE)/(TOTAL * 255);
                *pDst = (uchar)anLookup[weightedSum];
            }
        }
    }
};
void setZero(IplImage * pIm);

void doDownSample(const int nScale, const IplImage * pSrc, IplImage * pDest);

//Can pass in & of null ptr and image will be created
void markCorrespondences(const CBoWCorrespondences * pCorr, const CMask & mask, const IplImage * pIm1, const IplImage * pIm2, IplImage ** pImOut);

#define CV_RED CV_RGB(255,0,0)
#define CV_GREEN CV_RGB(0,255,0)
#define CV_BLUE CV_RGB(0,0,255)
#define CV_BLACK CV_RGB(0,0,0)
#define CV_WHITE CV_RGB(255,255,255)

#endif /* CONVERT_OPENCV_H_ */
