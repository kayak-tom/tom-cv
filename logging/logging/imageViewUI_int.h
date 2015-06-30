#ifndef PRUNINGUI_INT_H_
#define PRUNINGUI_INT_H_

#include "imageViewUI.h"
#include <util/exception.h>
#include <opencv2/core/core.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <string>
//#include <image/convert_OpenCV.h>
//#include "newCamera.h"
#include <boost/thread.hpp> 
#include <boost/circular_buffer.hpp> 


class redirectCout;

namespace cv
{
    class VideoWriter;
}

/**
 * @class CPruningUI
 * @brief Tiles 3 frames horizontally, then shows them
 */
class CFrameJoiner
{
    cv::Mat joinedFrame;
    int nNumFramesAdded;
public:
    CFrameJoiner() : nNumFramesAdded(0) {}
    void addFrame(const cv::Mat & frame, const std::string & title, const int id);
private:    
    void show(const std::string & title) const;
};

//TODO: Way too much state here--tidy me up 
class CImageViewUI : public CImageViewUI_base {
    std::string strFolderName, strLogfileName, strVidFrameDir;
    std::string WINDOW_TITLE;

    //Details on mouse handlers: http://www.cs.iit.edu/~agam/cs512/lect-notes/opencv-intro/index.html
    static void mouseHandler(int event, int x, int y, int flags, void* param);
    static void mouseHandlerMat(int event, int x, int y, int flags, void* param);
    redirectCout * pRedir;
    int nSave;
    
    eDisplayMode displayMode; //how long to show images for, output video, etc
    
    boost::scoped_ptr<cv::VideoWriter> pVideoWriter;
    boost::thread::id mainThreadId;
    
    typedef boost::circular_buffer<cv::Mat> TImHistoryBuffer;
    TImHistoryBuffer aImageHistory;
    static const int HISTORY_SIZE=20;
    cv::Size videoSize;
    
    cv::Size maxSize = cv::Size(640,480);
    
public:    
    static const int LONG_DISPLAY = 1000;
    static const int SHORT_DISPLAY = 10; //Display times in ms (10 insufficient for large images)
    
private:

    std::set<std::string> skipList; //set of image titles which won't be displayed
    std::set<std::string> fastList; //set of image titles which will only be displayed briefly
    
    CFrameJoiner frameJoiner;
        
    char showOneImage_int(const cv::Mat & image, const std::string TITLE, const int nWaitMS, const cv::Mat & actualIm);
public:
    CImageViewUI(); //Constructed by static variable initialisation
    virtual ~CImageViewUI();
    
    virtual void setDisplayMode(const eDisplayMode newDispMode);
    virtual void setDebugging();

    virtual char showOneImage(const cv::Mat & image, const std::string TITLE, const int nWaitMS=0);
    virtual char mergeAndShow(const cv::Mat & M1, const cv::Mat & M2, const cv::Mat & M3_optional, const int nWaitMS = 0);
    
    virtual char showText(const std::string text);
    
    //saves image in current output dir
    virtual void saveImage(const cv::Mat & image, const std::string & filename);
    
    virtual std::string getLogPath() { return strFolderName; }
    
    virtual bool isMT() const; //Is the caller in a thread other than the 'main' thread? (in which case we don't display anything as OpenCV's display fn doesn't work well with threads)
    
    virtual bool makingVideo() const; 

    virtual void joinFrame(const cv::Mat & frame, const std::string & title, const int cam) { frameJoiner.addFrame(frame, title, cam); }

    virtual void setNumBytesPrinted(const size_t nMaxNumBytesPrinted_in);
    virtual const cv::Size & getMaxSize() const { return maxSize;}
    
private:
    bool fastDisplay() const;
    //CFrameJoiner & getFrameJoiner();
    
    void tryStartGedit() const;
    void tryStartNautilus() const; 
    
    void setupWindow();
    
    const char * getVidFramePath() const;
};



#endif 
