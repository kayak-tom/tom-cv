/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * imageSourceFromVideo.h
 *
 */

#ifndef IMAGESOURCEFROMVIDEO_H_
#define IMAGESOURCEFROMVIDEO_H_

namespace cv
{
    class VideoCapture;
}

#include "imageSource.h"
#include <vector>
#include <string>
#include <boost/thread/mutex.hpp>

class CImageSourceFromVideo: public CImageSource
{
    boost::mutex imLoadMutex;
    const int nLoadSubset;
    int nLastId, nLastImageId;
    std::vector<int> aJumps;

    void init(CImParams & IM_PARAMS_in);

    cv::VideoCapture * pVideoCap;
    cv::Mat m;
public:
    CImageSourceFromVideo(CImParams & IM_PARAMS_in);
    virtual ~CImageSourceFromVideo();
    /*virtual void setGlobalDims(CImParams & IM_PARAMS);
    virtual bool loadImage(int nId, IplImage * pIm);
    virtual bool loadNextImage(int & nId, IplImage * pIm);*/
    virtual bool loadImage_int(int nId, IplImage * pIm);
    virtual bool loadNextImage_int(int & nId, IplImage * pIm);
    virtual void setGlobalDims(int & x, int & y, int & channels);

};

#endif /* IMAGESOURCEFROMVIDEO_H_ */
