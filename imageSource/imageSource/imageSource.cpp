/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 *
 *
 * imageSource.cpp
 */

#include "util/opencv.h"
#include "util/opencv_highgui.h"

#include "imageSourceFromDir.h"
#include "imageSourceFromVideo.h"
#include "image/convert_OpenCV.h"

#include <boost/filesystem.hpp>
#include <string>
#include <boost/lexical_cast.hpp>
#include "stdio.h"

/*#ifndef USE_OLD_OPENCV
#include "opencv2/imgproc/imgproc_c.h"
#else
#include "highgui.h"
#endif*/

using namespace std;

CImageSource * CImageSource::makeImageSource(CImParams & IM_PARAMS)
{
    if(IM_PARAMS.IM_SOURCE == CImParams::eImageDir)
        return new CImageSourceFromDir(IM_PARAMS);
    else if(IM_PARAMS.IM_SOURCE == CImParams::eCamera || IM_PARAMS.IM_SOURCE == CImParams::eVideoFile)
        return new CImageSourceFromVideo(IM_PARAMS);
    else
    {
        cout << "No image source required--simulated images. Creating one just to setup calibration and image size params...\n";
        CImageSourceFromDir imSourceToSetParams(IM_PARAMS);
        return 0;
    }
    //THROW("Cannot create this kind of image source");
}
void CImageSource::doDS(IplImage * pDest, const bool bDoGreyscale) const
{
    if(IM_PARAMS.SCALE_DOWN > 1)
        doDownSample(IM_PARAMS.SCALE_DOWN, pImInternal, bDoGreyscale ? pImInternalColour : pDest);
}

void CImageSource::doCorrection(IplImage * pFrame, const bool bDoGreyscale) //Eventually greyscale here too, but want colour images for making video frames...
{
    CHECK(bDoGreyscale != (0 != pImInternalColour), "Asked to greyscale but no internal colour image available");

    doDS(pFrame, bDoGreyscale);

    CHECK(IM_PARAMS.IM_CHANNELS != 1 && IM_PARAMS.ILLUMINATION_CORRECTION != CImParams::eNoIlluminationCorrection, "Illumination correction only works on monochrome images")

    if(bDoGreyscale)
    {
        if(!pGreyscaler)
            pGreyscaler = new CGreyscaler(true, IM_PARAMS.Greyscale.R, IM_PARAMS.Greyscale.G, IM_PARAMS.Greyscale.B, IM_PARAMS.Greyscale.GAMMA);

        pGreyscaler->greyScale(pImInternalColour, pFrame);
        REPEAT(3, cout << "Greyscale in image loader applied\n");
    }

    /*cvSaveImage("pImInternalColour.png", pImInternalColour);
    cvSaveImage("pImInternal.png", pImInternal);
    cvSaveImage("pFrame.png", pFrame);*/
    //removeBorder(pFrame);

    /* deprecatedif(IM_PARAMS.ILLUMINATION_CORRECTION == CImParams::eEqualiseHist)
    {
        cvEqualizeHist(pFrame, pFrame);
    }
    else*/ if(IM_PARAMS.ILLUMINATION_CORRECTION == CImParams::eEqualiseMeanSD)
    {
        CvScalar adMean, adSD;
        cvAvgSdv(pFrame, &adMean, &adSD);
        const double TARGET_SD = 80, TARGET_MEAN=128;
        double dMean = adMean.val[0];
        double dSD = adSD.val[0];

        //1) Scale image so it has the same SD
        double dAdjustScale = TARGET_SD/dSD;
        double dMeanAfterScaling = dAdjustScale*dMean;

        double dAdjustMean = TARGET_MEAN-dMeanAfterScaling;
        cvConvertScale(pFrame, pFrame, dAdjustScale, dAdjustMean);
    }
    else if(IM_PARAMS.ILLUMINATION_CORRECTION == CImParams::eEqualiseMean)
    {
        CvScalar adMean = cvAvg(pFrame);
        const double TARGET_MEAN=128;
        double dMean = adMean.val[0];

        double dAdjustMean = TARGET_MEAN-dMean;
        cvConvertScale(pFrame, pFrame, 1, dAdjustMean);
    }

    if(IM_PARAMS.SAVE_FRAMES != CImParams::eDontSave)
    {
        //if(boost::filesystem::exists("frames"))
        remove("frames");
        boost::filesystem::create_directory("frames");
        string frameName("frames/frame");
        char num[10];
        static int s_nFrame = 0;
        sprintf_s(num, 10, "%07d", s_nFrame++);
        //frameName += boost::lexical_cast<std::string>(s_nFrame++);
        frameName += num;
        frameName += ".";
        switch(IM_PARAMS.SAVE_FRAMES)
        {
        case CImParams::eSaveBMP:
            frameName += "bmp";
            break;
        case CImParams::eSaveJPG:
            frameName += "jpg";
            break;
        case CImParams::eSavePNG:
            frameName += "png";
            break;
        default:
            THROW("Enum val not handled")
        }
        cvSaveImage(frameName.c_str(), pFrame);
    }
}

IplImage * CImageSource::createImage() const
{
    return cvCreateImage(IM_PARAMS.SIZE(), IPL_DEPTH_8U, IM_PARAMS.IM_CHANNELS);
}

bool CImageSource::loadImage_choose(int & nId, IplImage * pIm, bool bNext)
{
    boost::mutex::scoped_lock lock(mxImLoader); //Multiple threads may try to load images

    bool bGreyscaleHere = SOURCE_CHANNELS != IM_PARAMS.IM_CHANNELS;

    IplImage * pImToLoadInto = pIm;
    if(bGreyscaleHere)
    {
        if(!pImInternalColour)
            pImInternalColour = cvCreateImage(IM_PARAMS.SIZE(), IPL_DEPTH_8U, SOURCE_CHANNELS);

        pImToLoadInto = pImInternalColour;
    }

    bool res = false;
    if(IM_PARAMS.SCALE_DOWN>1)
    {
        if(!pImInternal)
            pImInternal = cvCreateImage(cvSize(SOURCE_WIDTH, SOURCE_HEIGHT), IPL_DEPTH_8U, SOURCE_CHANNELS);

        pImToLoadInto = pImInternal;
    }

    if(!bNext)
        res = loadImage_int(nId, pImToLoadInto);
    else
        res = loadNextImage_int(nId, pImToLoadInto);

    if(!res)
        return false;

    doCorrection(pIm, bGreyscaleHere); //Includes downsample
    return true;
}
