/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 *
 * imageSourceFromVideo.cpp
 *
 */

#include "imageSourceFromVideo.h"
#include <boost/filesystem.hpp>
#include "util/opencv.h"
#include "util/opencv_highgui.h"
#include <boost/thread/mutex.hpp>
#include <fstream>
#include "image/convert_OpenCV.h"

using namespace std;
using namespace boost::filesystem;
using namespace cv;

CImageSourceFromVideo::CImageSourceFromVideo(CImParams & IM_PARAMS_in) : CImageSource(IM_PARAMS_in), nLoadSubset(IM_PARAMS_in.ImageDir.LOAD_SUBSET), pVideoCap(0)
{
    init(IM_PARAMS_in);
    CHECK(!pVideoCap, "Error initialising video capture");
    CHECK(!pVideoCap->isOpened(), "Error opening video capture");
}

void CImageSourceFromVideo::init(CImParams & IM_PARAMS_in)
{
    nLastId = -1; nLastImageId = -1;
    if(IM_PARAMS.IM_SOURCE == CImParams::eVideoFile)
    {
        CHECK(!IM_PARAMS.VideoFile.FILENAME.isInit() || ((const string &)IM_PARAMS.VideoFile.FILENAME).length()==0, "Video filename is required, specify Im.VideoFile.FILENAME=\"/path/to/video/file.avi.mov\" in config file");

        if(!boost::filesystem::exists((const string &)IM_PARAMS.VideoFile.FILENAME))
        {
            cout << "Video file \"" << IM_PARAMS.VideoFile.FILENAME.asSz() << "\" does not exist\n";
            THROW("Video file does not exist")
        }
        if(!boost::filesystem::is_regular_file((const string &)IM_PARAMS.VideoFile.FILENAME) && !boost::filesystem::is_symlink((const string &)IM_PARAMS.VideoFile.FILENAME))
        {
            cout << "Video file \"" << IM_PARAMS.VideoFile.FILENAME.asSz() << "\" is not a file\n";
            THROW("Video file not a file or link")
        }

        pVideoCap = new VideoCapture((const std::string &)IM_PARAMS.VideoFile.FILENAME);
        
        CHECK(!pVideoCap, "Error initialising video capture");

    }
    else if(IM_PARAMS.IM_SOURCE == CImParams::eCamera)
    {
        pVideoCap = new VideoCapture(IM_PARAMS.Camera.DEVICE_ID);
        
        CHECK(!pVideoCap, "Error initialising video capture from camera");
    }
    else
        THROW("Unhandled video source");

    CHECK(!pVideoCap->isOpened(), "Error opening video capture");

    if(!IM_PARAMS.IM_WIDTH.isInit())
    {
        cout << "Dynamically initialising image height and width and colour channels...\n";
        setImageParams(IM_PARAMS_in);

        boost::filesystem::path calibFolderName((const string &)IM_PARAMS.VideoFile.FILENAME);
        
        IM_PARAMS_in.initCalibration(calibFolderName.string().c_str());
    }
    else
        cout << "Image height and width and colour channels initialised via config file\n";
}

CImageSourceFromVideo::~CImageSourceFromVideo()
{
    delete pVideoCap;
}

void CImageSourceFromVideo::setGlobalDims(int & x, int & y, int & channels)
{
    *pVideoCap >> m;

    x = m.cols;
    y = m.rows;
    channels = m.channels();
}

//Also set id/
bool CImageSourceFromVideo::loadNextImage_int(int & nId, IplImage * pIm)
{
    if(nLastId < 0)
    {
        nLastId = 0;
        nLastImageId = 0;
    }
    else
    {
        nLastId++;
        nLastImageId++;
    }
    for(std::vector<int>::const_iterator pJump = aJumps.begin(); pJump != aJumps.end(); pJump++)
        if(nLastImageId == *pJump)
            nLastId+=1000;

    nId = nLastId;

    cout << "Returning image " << nLastImageId << " with id " << nId << endl;
    return loadImage_int(nLastImageId, pIm);
}

bool CImageSourceFromVideo::loadImage_int(int nId, IplImage * pIm)
{
    const uchar * pDataEnd = m.data+m.step*m.rows;
    for(;;)
    {
        for(int i=nLoadSubset; (i>0) && m.data; i--)
            *pVideoCap >> m; //Don't know any better way of doing this..?

        if(!m.data)
        {
            cout << "End of video...\n";
            return false;
        }

        uchar data = *(m.data);
        bool bGoodFrame = false;
        for(uchar * pData = m.data; pData < pDataEnd; pData++)
            if(data != *pData)
            {
                bGoodFrame = true;
                break;
            }

        if(!bGoodFrame)
        {
            cout << "Dropping bad video frame...\n";
            cvWaitKey(10);
        }
        else
            break;
    }
    for(uchar * pData = m.data, * pImData = (uchar *)pIm->imageData; pData < pDataEnd; pData++, pImData++)
        *pImData = *pData;

    return true;
}
