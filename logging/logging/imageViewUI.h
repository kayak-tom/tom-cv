#ifndef PRUNINGUI_H_
#define PRUNINGUI_H_

#include <util/pp.h>
#include <boost/scoped_ptr.hpp>

#include <opencv2/core/core.hpp>

enum eDisplayMode { eFullDisplay, eSlowDisplay, eFastDisplay, eNoDisplay, eMakeVideo, eMakeFastVideo, eMakeSplitVideo };
enum eDisplaySize { eScaleToFitScreen, eScaleDown, eFullSize };

class CImageViewUI_base
{
    friend CImageViewUI_base & UI();
    static boost::scoped_ptr<CImageViewUI_base> s_pUI;
protected:
    bool bGlobalDebuggingMode;
    double dClockStartTime; //For debugging--framerates etc.
    eDisplaySize displaySize = eFullSize;
public:
    CImageViewUI_base();
    virtual ~CImageViewUI_base() {}

    virtual void setDisplayMode(const eDisplayMode newDispMode) = 0;
    void setDisplaySize(const eDisplaySize newDispSize) { displaySize=newDispSize;};

    virtual char showOneImage(const cv::Mat & image, const std::string TITLE, const int nWaitMS=0) = 0;
    virtual char mergeAndShow(const cv::Mat & M1, const cv::Mat & M2, const cv::Mat & M3_optional, const int nWaitMS = 0) = 0;
    virtual char showText(const std::string text) = 0;

    //saves image in current output dir
    virtual void saveImage(const cv::Mat & image, const std::string & filename) = 0;

    /* Set a global debugging mode which will turn on more verbose output, slow the display, etc. Should only be set programatically when running code specifically to find an error */
    virtual void setDebugging() = 0;
    
    bool isDebugging() const 
    { 
        return bGlobalDebuggingMode; 
    }

    virtual bool isMT() const = 0; //Is the caller in a thread other than the 'main' thread? (in which case we don't display anything as OpenCV's display fn doesn't work well with threads)
    virtual bool makingVideo() const = 0;
    virtual void joinFrame(const cv::Mat & frame, const std::string & title, const int id) = 0;
    virtual std::string getLogPath() = 0;

    virtual void setNumBytesPrinted(const size_t nMaxNumBytesPrinted_in) = 0;
    
    double getElapsedTime() const;

    virtual bool crop3rdFrame() const {
        return makingVideo();    //Control video output 2.5* or 3* wide
    }

    double getCurrentTimeSecs() const;
    
    virtual const cv::Size & getMaxSize() const=0;
};

CImageViewUI_base & UI();


//Macro for convenience:
#define SHOW(im) { if(bVerbose) UI().showOneImage(im, #im); }
#define SHOW2(im, str) { if(bVerbose) UI().showOneImage(im, std::string(#im) + "-" + str); }
#define FASTSHOW(im) { if(bVerbose) UI().showOneImage(im, #im, 10); }
#define FASTSHOW2(im, str) { if(bVerbose) UI().showOneImage(im, std::string(#im) + "-" + str, 10); }

//Some functions should always have bVerbose true (turn them off by not calling them), and some places verbose is used to display some output before throwing some exception

//if bVerbose, creates a blank cv::Mat of appropriate size.
#define MAT(name) cv::Mat name; { if(bVerbose) name = cv::Mat(UI().getMaxSize(), CV_8UC3, cv::Scalar()); }
#define MAT2(name, frameId) cv::Mat name; { if(bVerbose) name = scene().getImageForOutput(frameId); }

void outputVideoTitles(const cv::Size & size);

#define int_TIME2(label, line, ...) const double dStartTimeSeconds_ ## line = UI().getElapsedTime();\
                                    __VA_ARGS__; \
                                    cout << "TIME_FUNCTION: " << label << " took " << (UI().getElapsedTime() - dStartTimeSeconds_ ## line) << "s" << endl;

#define int_TIME1(label, line, ...) int_TIME2(label, line, __VA_ARGS__)

#define TIME(label, ...) int_TIME1(label, __LINE__, __VA_ARGS__)


//Creates spreadsheet in log file dir. Adds headers if file doesn't exist yet
void logToTSV(const std::string label, const std::ostringstream & ss, const std::ostringstream & ss_headers, const bool bLogGlobal);

#endif /* PRUNINGUI_H_ */
