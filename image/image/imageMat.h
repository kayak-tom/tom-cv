/* 
 * File:   imageMat.h
 * Author: tom
 *
 * Created on 28 May 2011, 3:28 PM
 */

#ifndef IMAGEMAT_H
#define	IMAGEMAT_H

#include "util/exception.h"
#include "util/opencv.h"
#include <boost/noncopyable.hpp>
#include "util/optimisation_attributes.h"
#include "util/floor.h"

//1-channel image with output in range [0, ..., 1.0]
template<typename TImType=uchar /* or float */, typename TOutputType=double /* or float */> 
class CImageMat : public cv::Mat /*, boost::noncopyable*/ {
    inline void checkInit() const
    {
#ifdef _DEBUG
        if((rows==0||cols==0||!isInit()))
            THROW("CImageMat not initialised");
        //else
          //  if(IS_DEBUG) CHECK(cv::Mat::type() != CV_8UC1, "CImageMat is not CV_8UC1");
#endif

        checkType();
    } 
    
    inline void checkType() const
    {
        if(IS_DEBUG) CHECK(channels() != 1, "should only be used with 1-channel images");
        if(IS_DEBUG) CHECK(type() != getCVType(), "Probability image initialised with wrong type");
    }
public:

    inline bool isInit() const { return data ? true : false; }
    
    void operator=(const cv::Mat & M) {
        cv::Mat::operator =(M);
        checkType();
    }

    void operator=(const cv::MatExpr & M) {
        cv::Mat::operator =(M);
        checkType();
    }

    CImageMat() { data = 0; }
    
    CImageMat(const cv::Size & size) : cv::Mat(size, getCVType()) {
    }
    CImageMat(int rows, int cols) : cv::Mat(cv::Size(cols, rows), getCVType()) {
    }
    
    inline static const int getCVType()  
    {
        return isUchar() ? CV_8UC1 : CV_32FC1;
    }
    
    void init(const cv::Size & size) 
    {
        if(IS_DEBUG) CHECK(size.height == 0 || size.width == 0, "CImageMat already initialised");
        if(IS_DEBUG) CHECK(data != 0, "CImageMat already initialised");
        
        //cv::Mat temp(size, getCVType());
        //*this = temp;
        create(size.height, size.width, getCVType());
        
        checkInit(); //failing here indicates a problem...
    }
    
    void setZero() {
        checkInit(); 
        //cv::Mat::operator =(cv::Scalar(0));
        this->setTo(0);
        checkInit(); 
    }

    inline TOutputType getPx(int x, int y) const {
        if(IS_DEBUG) CHECK(channels() != 1, "getPx should only be used with 1-channel images");
        
        const TOutputType dProb = imToProb(at<TImType> (y, x));
        //if(IS_DEBUG) CHECK(dProb < 0 || dProb > 1, "Prob OOB");

        return dProb;
    }
    inline TOutputType getPxWithDefault(int x, int y, TOutputType dDefault) const {
        if(OOB(x,y))
            return dDefault;
        return getPx(x,y);
    }
    inline void setProb(int x, int y, TOutputType dProb) {
        if(IS_DEBUG) CHECK(channels() != 1, "getPx should only be used with 1-channel images");
        //if(IS_DEBUG) CHECK(dProb < 0 || dProb > 1, "Prob OOB");
        
        at<TImType>(y, x) = probToIm(dProb);
    }
    
protected:
    /*inline uchar & getPx(int x, int y) {
        if(IS_DEBUG) CHECK(channels() != 1, "getPx should only be used with 1-channel images");
        return at<uchar> (y, x);
    }*/
    
    inline static TOutputType imToProb(TImType prob);
    inline static TImType probToIm(TOutputType prob);
public:
    inline TOutputType getSubPxWithDefault(TOutputType x, TOutputType y, TOutputType dDefault) const
    {
        if(OOB_subpix(x,y))
            return dDefault;
        else
            return getSubPx(x, y);
    }
    
    inline TOutputType getSubPx(int x, TOutputType y) const HOT HARD_INLINE;
    inline TOutputType getSubPx(TOutputType x, TOutputType y) const HOT HARD_INLINE;
    
protected:
//    inline uchar * getData() { checkInit(); return (uchar *)data; }
//    inline const uchar * getData() const { checkInit(); return (const uchar *)data; }
public:    
    inline bool OOB(int x, int y) const
    {
        return x < 0 || y < 0 || x >= cols || y >= rows;
    }
    inline bool OOB_subpix(TOutputType x, TOutputType y) const
    {
        return x < 1 || y < 1 || x >= cols-1.0 || y >= rows-1.0;
    }
    
    inline static const bool isUchar() { return sizeof(TImType) == sizeof(uchar); }
};

//For original probability images with probs as uchars
/*template<typename TOutputType>
class CImageMat<uchar, TOutputType>
{
    protected:
    inline uchar probToIm(TOutputType prob)
    {
        static const TOutputType dProbToUchar = 255;
        return cv::saturate_cast<uchar>(dProbToUchar*prob);
    }
    inline TOutputType imToProb(uchar pxVal)
    {
        static const double dUcharToProb = 1.0/255.0;
        return dUcharToProb*(TOutputType)pxVal;
    }
};*/
//template<typename TOutputType>
//inline uchar CImageMat<uchar, TOutputType>::probToIm(TOutputType prob)
//{
//    static const TOutputType dProbToUchar = 255;
//    return cv::saturate_cast<uchar>(dProbToUchar*prob);
//}
//template<typename TOutputType>
//inline TOutputType CImageMat<uchar, TOutputType>::imToProb(uchar pxVal)
//{
//    static const double dUcharToProb = 1.0/255.0;
//    return dUcharToProb*(TOutputType)pxVal;
//}

//For other floating point types:
template<typename TImType, typename TOutputType>
inline TImType CImageMat<TImType, TOutputType>::probToIm(TOutputType prob)
{
    if(IS_DEBUG) CHECK(prob < 0 || prob > 1, "Prob OOB");
    if(isUchar())
        return cv::saturate_cast<uchar>(255.0*prob);
    else
        return (TImType)prob;
}
template<typename TImType, typename TOutputType>
inline TOutputType CImageMat<TImType, TOutputType>::imToProb(TImType pxVal)
{
    if(isUchar())
    {
        static const TOutputType dUcharToProb =  1.0/255.0;
        return (TOutputType)(dUcharToProb * pxVal);
    }
    else
    {
        //if(IS_DEBUG) CHECK(pxVal < 0, "Prob OOB");
        return pxVal;
    }
}

template<typename TImType, typename TOutputType>
inline TOutputType CImageMat<TImType, TOutputType>::getSubPx(TOutputType x, TOutputType y) const {
    const int nx_lo = intFloor(x);
    const int nx_hi = nx_lo + 1;

    const TOutputType x_lo = (TOutputType) nx_lo;

    TOutputType dWeightHi = x - x_lo;
    TOutputType dWeightLo = 1 - dWeightHi;

    if(IS_DEBUG) CHECK(dWeightHi < 0 || dWeightHi > 1, "Weight OOB");
    if(IS_DEBUG) CHECK(dWeightLo < 0 || dWeightLo > 1, "Weight OOB");

    return dWeightLo * getSubPx(nx_lo, y) + dWeightHi * getSubPx(nx_hi, y);
}

template<typename TImType, typename TOutputType>
inline TOutputType CImageMat<TImType, TOutputType>::getSubPx(int x, TOutputType y) const {
    const int ny_lo = intFloor(y);
    const int ny_hi = ny_lo + 1;

    const TOutputType y_lo = (TOutputType) ny_lo;

    if(IS_DEBUG) CHECK(y_lo > y, "floor failed");

    TOutputType dWeightHi = y - y_lo;
    TOutputType dWeightLo = 1 - dWeightHi;

    if(IS_DEBUG) CHECK(dWeightHi < 0 || dWeightHi > 1, "Weight OOB");
    if(IS_DEBUG) CHECK(dWeightLo < 0 || dWeightLo > 1, "Weight OOB");

    return dWeightLo * (TOutputType) getPx(x, ny_lo) + dWeightHi * (TOutputType) getPx(x, ny_hi);
}


/**
 * @brief If dest is not setup, copy image size and format from src to dest
 * @param src
 * @param dest
 */
inline void initImage(const cv::Mat & src, cv::Mat & dest, int type = -1)
{
    if(type == -1)
        type = src.type();
    dest.create(src.rows, src.cols, type);
}

/**
 * @brief If dest is not setup, copy image size and format from src to dest. Also set 0
 * @param src
 * @param dest
 */
inline void initImage0(const cv::Mat & src, cv::Mat & dest, int type = -1)
{
    if(type == -1)
        type = src.type();
    /*cout << src.rows << endl;
    cout << src.cols << endl;
    cout << type << endl;
    cout << CV_8UC1 << endl;*/
    dest.create(src.rows, src.cols, type);
    dest.setTo(cv::Scalar());
}

inline std::ostream & operator<<(std::ostream & out, const cv::Scalar & s)
{
    for (int i = 0; i < 4; i++)
        out << s[i] << " ";
    return out;
}

inline void pp(const cv::Mat & M) {
    cout << "Image rows " << M.rows << ", cols " << M.cols << ", channels " << M.channels() << ", depth " << M.depth() << " ";
    if (M.depth() == CV_32F)
        cout << "CV_32F";
    else if (M.depth() == CV_8U)
        cout << "CV_8U";
    else if (M.depth() == CV_64F)
        cout << "CV_64F";

	if(M.channels() <= 4)
	{
		cout << " means: ";
		cout << cv::mean(M);
	}

    if(M.channels() == 1)
    {
        double dMin,dMax;
        cv::minMaxLoc(M, &dMin, &dMax);
        cout << "Min: " << dMin << " Max: " << dMax << endl;
    }

    cout << endl;
}
inline void pp(const std::string name, const cv::Mat & M) {
    cout << name << " ";
    pp(M);
}

#define PP(M) pp(#M, M)

/* Gaussian blur in-place, with suitable size.
*/
inline void blur(cv::Mat & im, const double dSigma)
{
    const int nSize = 2*(cvRound(2*(dSigma))) + 1;
    cv::GaussianBlur(im, im, cv::Size(nSize, nSize), dSigma);
}
inline void blur(const cv::Mat & src, cv::Mat & dest, const double dSigma)
{
    const int nSize = 2*(cvRound(2*(dSigma))) + 1;
    cv::GaussianBlur(src, dest, cv::Size(nSize, nSize), dSigma);
}

/* Check rows, cols suitable
*/
inline void checkImInit(const cv::Mat & im)
{
    if(IS_DEBUG) CHECK(!(im.rows>0 && im.rows < 10000 && im.cols>0 && im.cols < 10000 && im.channels() > 0 && im.channels() < 5), "Image is uninitialised (or extremely large)");
}










#endif	/* IMAGEMAT_H */

