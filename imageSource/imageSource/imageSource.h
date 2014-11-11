/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * imageSource.h
 *
 *  Created on: 5/06/2009
 *      Author: tom
 */

#ifndef IMAGESOURCE_H_
#define IMAGESOURCE_H_

#include "image/imageAccess.h"
#include "boost/thread/mutex.hpp"

class CGreyscaler;

class CImageSource : boost::noncopyable
{
    IplImage * pImInternal, * pImInternalColour;
    CGreyscaler * pGreyscaler;
    boost::mutex mxImLoader;
protected:
    const CImParams & IM_PARAMS;
    int SOURCE_WIDTH, SOURCE_HEIGHT, SOURCE_CHANNELS;

    virtual bool loadImage_int(int nId, IplImage * pIm) = 0;
    virtual bool loadNextImage_int(int & nId, IplImage * pIm) = 0;
    virtual void setGlobalDims(int & x, int & y, int & channels) = 0;

    void doDS(IplImage * pDest, const bool bDoGreyscale) const;
    void doCorrection(IplImage * pFrame, const bool bDoGreyscale); //Eventually greyscale here too, but want colour images for making video frames...
    bool loadImage_choose(int & nId, IplImage * pIm, bool bNext);

    //Could make public...
    void setImageParams(CImParams & IM_PARAMS)
    {
        CHECK(SOURCE_CHANNELS != -1, "setting dims when already initialised");

        int nScaledown = IM_PARAMS.SCALE_DOWN;
        setGlobalDims(SOURCE_WIDTH, SOURCE_HEIGHT, SOURCE_CHANNELS);
        if(!IM_PARAMS.IM_WIDTH.isInit())
        {
            IM_PARAMS.IM_WIDTH = SOURCE_WIDTH/nScaledown;
            IM_PARAMS.IM_HEIGHT = SOURCE_HEIGHT/nScaledown;
            if(IM_PARAMS.ILLUMINATION_CORRECTION != CImParams::eNoIlluminationCorrection && !IM_PARAMS.GREYSCALE_ON_LOAD)
            {
                cout << "WARNING: Greyscale images needed for illumination correction, setting ";
                IM_PARAMS.GREYSCALE_ON_LOAD = true;
            }

            IM_PARAMS.IM_CHANNELS = IM_PARAMS.GREYSCALE_ON_LOAD ? 1 : SOURCE_CHANNELS;
        }
    }
public:

    bool loadImage(int nId, IplImage * pIm)
    {
        return loadImage_choose(nId, pIm, false);
    }

    bool loadNextImage(int & nId, IplImage * pIm)
    {
        return loadImage_choose(nId, pIm, true);
    }

    virtual ~CImageSource()
    {
        if(pImInternal)
            cvReleaseImage(&pImInternal);
        if(pImInternalColour)
            cvReleaseImage(&pImInternalColour);
    };

    CImageSource(const CImParams & IM_PARAMS) : pImInternal(0), pImInternalColour(0), pGreyscaler(0), IM_PARAMS(IM_PARAMS), SOURCE_WIDTH(-1), SOURCE_HEIGHT(-1), SOURCE_CHANNELS(-1)
    {
        /*if(IM_PARAMS.IM_WIDTH.isInit())
        {
            int nScaledown = IM_PARAMS.SCALE_DOWN;
            SOURCE_WIDTH = IM_PARAMS.IM_WIDTH * nScaledown;
            SOURCE_HEIGHT = IM_PARAMS.IM_HEIGHT * nScaledown;

            if(!IM_PARAMS.GREYSCALE_ON_LOAD)
                SOURCE_CHANNELS = IM_PARAMS.IM_CHANNELS;
            else
            {
                cout << "Warning: assuming images have 3 channels, as aren't greyscaling on load\n";
                SOURCE_CHANNELS = 3;
            }
        }*/
    };

    static CImageSource * makeImageSource(CImParams & IM_PARAMS);

    IplImage * createImage() const;
};

#endif /* IMAGESOURCE_H_ */
