/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * imageAccess.h
 *
 *  Created on: 13/10/2009
 *      Author: tom
 */

#ifndef IMAGEACCESS_H_
#define IMAGEACCESS_H_

#include "util/exception.h"
#include "util/opencv.h"
#include "util/location.h"
#include "params/param.h"
#include "util/calibration.h"
#include "HSI.h"
#include "util/convert.h"

#ifdef _DEBUG
#ifdef __OPTIMIZE__
//#error "Debug checks in optimised build will really slow down image access"
//Actually these are off. #  pragma message "WARNING: Debug checks in optimised build will really slow down image access"
#endif
#endif

#define CHECK2(x, message)

//Extends cv::Mat functionality (subpixel access, cleaner interface)

template<typename TImData>
class CIplPx {
    inline static TImData & index(IplImage * m_img, int x, int y, int channels, int c) HOT HARD_INLINE;
    inline static TImData val(const IplImage * m_img, int x, int y, int channels, int c) HOT HARD_INLINE;
public:
    //Todo: possibly inefficient

    inline static void setRGB(IplImage * pIm, int x, int y, TImData R, TImData G, TImData B) HARD_INLINE {
        CHECK2(!pIm || pIm->nChannels != 3, "CIplPx::setRGB: Wrong num of channels");
        setBlue(pIm, x, y, B);
        setGreen(pIm, x, y, G);
        setRed(pIm, x, y, R);
    }

    inline static void setRed(IplImage * pIm, int x, int y, TImData val) HARD_INLINE {
        CHECK2(!pIm || pIm->nChannels != 3, "CIplPx::setRed: Wrong num of channels");
        index(pIm, x, y, 3, 2) = val;
    }

    inline static void setGreen(IplImage * pIm, int x, int y, TImData val) HARD_INLINE {
        CHECK2(!pIm || pIm->nChannels != 3, "CIplPx::setGreen: Wrong num of channels");
        index(pIm, x, y, 3, 1) = val;
    }

    inline static void setBlue(IplImage * pIm, int x, int y, TImData val) HARD_INLINE {
        CHECK2(!pIm || pIm->nChannels != 3, "CIplPx::setBlue: Wrong num of channels");
        index(pIm, x, y, 3, 0) = val;
    }
    inline static void setGrey(IplImage * pIm, int x, int y, TImData val) HOT HARD_INLINE;

    inline static TImData getRed(const IplImage * pIm, int x, int y) PURE_FN HARD_INLINE {
        CHECK2(!pIm || pIm->nChannels != 3, "CIplPx::getRed: Wrong num of channels");
        return val(pIm, x, y, 3, 2);
    }

    inline static TImData getGreen(const IplImage * pIm, int x, int y) PURE_FN HARD_INLINE {
        CHECK2(!pIm || pIm->nChannels != 3, "CIplPx::getGreen: Wrong num of channels");
        return val(pIm, x, y, 3, 1);
    }

    inline static TImData getBlue(const IplImage * pIm, int x, int y) PURE_FN HARD_INLINE {
        CHECK2(!pIm || pIm->nChannels != 3, "CIplPx::getBlue: Wrong num of channels");
        return val(pIm, x, y, 3, 0);
    }

    inline static void getRGB(const IplImage * pIm, int x, int y, CRedGreenBlue & rgb) HARD_INLINE {
        CHECK2(!pIm || pIm->nChannels != 3, "CIplPx::getRGB: Wrong num of channels");
        rgb.B = getBlue(pIm, x, y);
        rgb.G = getGreen(pIm, x, y);
        rgb.R = getRed(pIm, x, y);
    }

    inline static TImData getPx(const IplImage * pIm, int x, int y, int channel) {
        return val(pIm, x, y, pIm->nChannels, channel);
    }

    template <typename TProbType>
    inline static TProbType getSubPx(const IplImage * pIm, int x, double y, int channel) {
        //int y_lo = (int)y;
        /*int y_lo = doubleToInt(y);*/
        //int y_hi = y_lo+1;

        const int ny_lo = intFloor(y);
        const int ny_hi = ny_lo + 1;

        const double y_lo = (double) ny_lo;

        CHECK2(y_lo > y, "floor failed");

        TProbType dWeightHi = y - y_lo;
        TProbType dWeightLo = 1 - dWeightHi;

        CHECK2(dWeightHi < 0 || dWeightHi > 1, "Weight OOB");
        CHECK2(dWeightLo < 0 || dWeightLo > 1, "Weight OOB");

        return dWeightLo * (TProbType) val(pIm, x, ny_lo, pIm->nChannels, channel) + dWeightHi * (TProbType) val(pIm, x, ny_hi, pIm->nChannels, channel);
    }

    template <typename TProbType>
    inline static double getSubPx(const IplImage * pIm, double x, double y, int channel) {
        //int x_lo = (int)x;

        //const int x_lo = doubleToInt(x); Doesn't work, rounds rather than truncates
        //const int x_hi = x_lo+1;

        const int nx_lo = intFloor(x);
        const int nx_hi = nx_lo + 1;

        const double x_lo = (double) nx_lo;

        TProbType dWeightHi = x - x_lo;
        TProbType dWeightLo = 1 - dWeightHi;

        CHECK2(dWeightHi < 0 || dWeightHi > 1, "Weight OOB");
        CHECK2(dWeightLo < 0 || dWeightLo > 1, "Weight OOB");

        return dWeightLo * getSubPx<TProbType > (pIm, nx_lo, y, channel) + dWeightHi * getSubPx<TProbType > (pIm, nx_hi, y, channel);
    }

    inline static TImData getGrey(const IplImage * pIm, int x, int y) HOT PURE_FN HARD_INLINE;

    inline static void copyPx(const IplImage * pIm, int xs, int ys, IplImage * pImDest, int xd, int yd) HARD_INLINE {
        CHECK2(!pIm || !pImDest, "CIplPx::setRGB: No im");
        if (pIm->nChannels >= 3) {
            setBlue(pImDest, xd, yd, getBlue(pIm, xs, ys));
            setGreen(pImDest, xd, yd, getGreen(pIm, xs, ys));
            setRed(pImDest, xd, yd, getRed(pIm, xs, ys));
        } else {
            setGrey(pImDest, xd, yd, getGrey(pIm, xs, ys));
        }
    }

    /*
        Access a subarray of an IplImage. Use: IplImage subIm=*pMyWholeImage; CIplPx<uchar>::cropImageToRect(subIm, cvRect(20,20,200,300));
     */
    static void cropImageToRect(IplImage & subImage, const CvRect & rect) {
        CHECK2(rect.x + rect.width > subImage.width, "Sub-imge width too big");
        CHECK2(rect.y + rect.height > subImage.height, "Sub-imge height too big");
        CHECK2(rect.x < 0 || rect.width < 0 || rect.y < 0 || rect.height < 0, "Rect contains negative values");

        subImage.imageData += rect.x * subImage.nChannels + rect.y * subImage.widthStep;
        subImage.width = rect.width;
        subImage.height = rect.height;
    }

    static void cropImage(IplImage & subImage, int nMargin) {
        subImage.imageData += nMargin * subImage.nChannels + nMargin * subImage.widthStep;
        subImage.width -= 2 * nMargin;
        subImage.height -= 2 * nMargin;
    }
};

inline bool OOB(const cv::Mat& img, const cv::Point pt1) {
    if (pt1.x < 0 || pt1.y < 0 || pt1.x >= img.cols || pt1.y >= img.rows)
        return true;
    return false;
}


template<typename TImData>
inline TImData & CIplPx<TImData>::index(IplImage * m_img, int x, int y, int channels, int c) {
    CHECK2(sizeof (TImData)*8 != m_img->depth, "Incorrect image depth(datatype) specified");
    CHECK2(!(x >= 0 && y >= 0 && x < m_img->width && y < m_img->height) || m_img->nChannels != channels, "CIplPx: Access OOB");
    return ((TImData *) (void *) (m_img->imageData + m_img->widthStep * y))[x * channels + c];
}

template<typename TImData>
inline TImData CIplPx<TImData>::val(const IplImage * m_img, int x, int y, int channels, int c) {
    CHECK2(sizeof (TImData)*8 != m_img->depth, "Incorrect image depth(datatype) specified");
    CHECK2(!(x >= 0 && y >= 0 && x < m_img->width && y < m_img->height) || m_img->nChannels != channels, "CIplPx: Access OOB");
    return ((TImData *) (void *) (m_img->imageData + m_img->widthStep * y))[x * channels + c];
}

template<typename TImData>
inline void CIplPx<TImData>::setGrey(IplImage * pIm, int x, int y, TImData val) {
    CHECK2(!pIm || pIm->nChannels != 1, "CIplPx::setGrey: Wrong num of channels");
    index(pIm, x, y, 1, 0) = val;
}

template<typename TImData>
inline TImData CIplPx<TImData>::getGrey(const IplImage * pIm, int x, int y) {
    CHECK2(!pIm || pIm->nChannels != 1, "CIplPx::getGrey: Wrong num of channels");
    return val(pIm, x, y, 1, 0);
}

PARAMCLASS(Im)
PARAM(IM_WIDTH, 32, 65536 / SUBPIX_RES, 640, "Set automatically by image loader. NB this has to fit in a short (allowing fast detection of duplicates), limiting max val. #define SUBPIX_RES to something lower to work with large images (see PX_SCALE_INT=SUBPIX_RES in location.h).")
PARAM(IM_HEIGHT, 32, 65536 / SUBPIX_RES, 480, "Set automatically by image loader. NB this has to fit in a short (allowing fast detection of duplicates), limiting max val. #define SUBPIX_RES to something lower to work with large images (see PX_SCALE_INT=SUBPIX_RES in location.h).")
PARAM(IM_CHANNELS, 1, 3, 3, "Set automatically by image loader. Num channels in images being loaded. 1 (grey) and 3 (RGB) supported in Linux, possibly only 3 in Windows.")
PARAME(IM_SOURCE, ImageDir, "Frame source")
PARAMB(CALIBRATION_INIT, false, "Flag whether calibration matrix is initialised (set automatically, not always needed)")
PARAM(SCALE_DOWN, 1, 10, 1, "Scale down images in image loader by this factor")
PARAME(ILLUMINATION_CORRECTION, NoIlluminationCorrection, "Correction--can also set SAVE_FRAMES to save correced frames for next run")
PARAMB(GREYSCALE_ON_LOAD, false, "Greyscale images in image loader (not good if colour frames for video needed). Turned on if ILLUMINATION_CORRECTION is set.")
PARAME(SAVE_FRAMES, DontSave, "Dumps all frames used to 'frames' folder, *includes illumination corrections, resizeing, etc.")
PARAM(RD2, -10, 10, 0, "RD correction coeffs. For tuning. Don't use, set in config file instead.")
PARAM(RD4, -10, 10, 0, "RD correction coeffs. For tuning. Don't use, set in config file instead.")
CHILDCLASS(Greyscale, "Params for greyscaling colour images (used for corner detection and usually for descriptors)")
CHILDCLASS(ImageDir, "Params for loading frames from images from a directory")
CHILDCLASS(VideoFile, "Params for loading frames from a video file (.avi usually works)")
CHILDCLASS(Camera, "Params for loading frames from an attached camera (untested)") {
}

CNumParamDerived<int> IM_WIDTH, IM_HEIGHT, IM_CHANNELS;

MAKEENUMPARAM5(IM_SOURCE, ImageSim, ImageDir, SimMap, VideoFile, Camera)

CvSize SIZE(int nMargin = 0) const {
    return cvSize(((int) IM_WIDTH) - 2 * nMargin, ((int) IM_HEIGHT) - 2 * nMargin);
}

private:
CNumParamDerived<bool> CALIBRATION_INIT;
CCamCalibMatrix K;
public:
CNumParam<int> SCALE_DOWN;
MAKEENUMPARAM3(ILLUMINATION_CORRECTION, NoIlluminationCorrection, /*EqualiseHist,*/ EqualiseMean, EqualiseMeanSD)
CNumParam<bool> GREYSCALE_ON_LOAD;
MAKEENUMPARAM4(SAVE_FRAMES, DontSave, SaveBMP, SavePNG, SaveJPG)
CNumParam<double> RD2, RD4;

const CCamCalibMatrix & getCamCalibrationMat() const {
    if(IS_DEBUG) CHECK(!CALIBRATION_INIT, "No calibration matrix setup, should have been found and initialised by image loader")
    return K;
}

/*const CInvCamCalibMatrix & getInvCamCalibrationMat() const
{
    if(IS_DEBUG) CHECK(!CALIBRATION_INIT, "No calibration matrix setup, should have been found and initialised by image loader")
    return K_inv; 
}*/

double getFocalLength() const {
    if(IS_DEBUG) CHECK(!CALIBRATION_INIT, "No calibration matrix setup, should have been found and initialised by image loader")
    return K.focalLength();
}

void initCalibration(const char * szFolderName) {
    std::cout << "Initialising calibration matrix from " << szFolderName << "...\n";
    if(IS_DEBUG) CHECK(!szFolderName, "Calibration matrix folder is null")

    try {
        K.init(szFolderName, SCALE_DOWN);
        CALIBRATION_INIT = true;
        K.testK(IM_WIDTH, IM_HEIGHT);
    } catch (...) {
        std::cout << "ERROR: Failed to find calibration matrix. Continuing without (may not be needed)...\n";
    }
}

PARAMCLASS(Greyscale)
PARAM(R, 0, 1024, 250, "Relative importance of red channel, default values work well experimentally")
PARAM(G, 0, 1024, 150, "Relative importance of green channel, default values work well experimentally")
PARAM(B, 0, 1024, 10, "Relative importance of blue channel, default values work well experimentally") //These default values work well experimentally
PARAM(GAMMA, 0, 5, 1, "Gamma correction, 1 turns it off (works best)") {
}

CNumParam<int> R, G, B;
CNumParam<double> GAMMA;
};

PARAMCLASS(ImageDir)
PARAMSTR(IMAGE_DIR, "Path to folder containing images (jpg, png, bmp, pgm, tiff if HAVE_TIFF defined)")
PARAM(LOAD_SUBSET, 1, 1000, 1, "Only load a subset of images, e.g. 2=every other one") {
}

CStringParam IMAGE_DIR;
CNumParam<int> LOAD_SUBSET;
};

PARAMCLASS(VideoFile)
PARAMSTR(FILENAME, "Filename and path to video file. Make sure calib.txt containing calibration is in the same folder.")
PARAM(LOAD_SUBSET, 1, 1000, 1, "Only load a subset of images, e.g. 2=every other one") {
}

CStringParam FILENAME;
CNumParam<int> LOAD_SUBSET;
};

PARAMCLASS(Camera)
PARAM(DEVICE_ID, 0, 1000000, 0, "Select one of several attached cameras")
PARAM(LOAD_SUBSET, 1, 1000, 1, "Only load a subset of images, e.g. 2=every other one (UNTESTED)") {
}

CNumParam<int> DEVICE_ID, LOAD_SUBSET;
};

MAKECHILDCLASS(Greyscale);
MAKECHILDCLASS(ImageDir);
MAKECHILDCLASS(VideoFile);
MAKECHILDCLASS(Camera);
};


#endif /* IMAGEACCESS_H_ */
