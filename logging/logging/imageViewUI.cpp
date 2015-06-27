/*
 * pruningUI.cpp
 *
 *  Created on: 1 Dec 2010
 *      Author: tom
 */
#define _CRT_SECURE_NO_WARNINGS
#include <logging/redirectCout.h>
#include <boost/filesystem.hpp>
#include <iomanip>

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/highgui/highgui_c.h>

#include <boost/date_time/posix_time/ptime.hpp>

#include "function_timer.h"
#include "PR.h"
#include <boost/scoped_ptr.hpp>
#include "imageViewUI_int.h"
#include <util/convert.h>

using namespace std;

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


//Watch order of static initialisation: top to bottom
boost::scoped_ptr<CImageViewUI_base> CImageViewUI_base::s_pUI ( new CImageViewUI() ); //construct 1st, destroy last
CEvalResults CEvalResults::s_evalResults;
CEvalFailure CEvalFailure::s_evalFailures;


CImageViewUI_base & UI()
{
    return *CImageViewUI_base::s_pUI;
}

std::string makeNiceFolderNameHACKTESTS()
{
    // Setup logging, windows, output directory, etc.

    std::string strFolderName = "logs/";

    const bool bRemoveLogsOnStartup = false;
    if(bRemoveLogsOnStartup) {
        if (boost::filesystem::exists(strFolderName)) {
            try {
                boost::filesystem::remove_all(strFolderName);
            } catch (...) {
                cout << "Error removing " << strFolderName << endl;
            }
        }
    }

    if (!boost::filesystem::exists(strFolderName)) {
        cout << "Creating directory " << strFolderName << endl;
        boost::filesystem::create_directory(strFolderName);
    }

    strFolderName.append(getTimeAndDate(false));
    strFolderName.append("/");

#ifdef _WIN32
    /**************  rdcm: parsing "/" and ":" to underlines "_" for windows *******************/
    char* cStrFolderName =  new char[strFolderName.size()+1];
    strcpy_s(cStrFolderName, strFolderName.size()+1, strFolderName.c_str());
    char* strFolderNameParts = strtok(cStrFolderName,"/:");
    strFolderName = "logs/";
    while(strFolderNameParts!=NULL) {
        strFolderName.append(strFolderNameParts);
        strFolderName.append("_");
        strFolderNameParts = strtok(NULL,"/:");
    }
    strFolderName.append("/");
    /************** END    rdcm: parsing "/" and ":" to underlines "_" for windows ****************/
#endif

    cout << "Creating directory " << strFolderName << endl;
    boost::filesystem::create_directory(strFolderName);

    return strFolderName;
}


double CImageViewUI_base::getCurrentTimeSecs() const
{
    //static boost::mutex mxGetTime;
    //boost::mutex::scoped_lock lock(mxGetTime);
    boost::posix_time::ptime currentTime = boost::posix_time::microsec_clock::local_time();

    double dCurrentTimeSecs = currentTime.time_of_day().total_milliseconds()/1000.0;

    while(dCurrentTimeSecs < dClockStartTime) {
        //We've passed midnight
        dCurrentTimeSecs += 24*60*60; //seconds per day; Keep it strictly increasing
    }

    return dCurrentTimeSecs;
}

double CImageViewUI_base::getElapsedTime() const
{
    const double dElapsedTime = getCurrentTimeSecs() - dClockStartTime;
    CHECK(dElapsedTime < 0, "Negative elapsed time. This could well wrap at midnight (but then the difference should be large, not a fraction of a second)");
    return dElapsedTime;
}

CImageViewUI_base::CImageViewUI_base() : bGlobalDebuggingMode(false), dClockStartTime(-1)
{
    dClockStartTime = getCurrentTimeSecs();//Must initialise dClockStartTime first
}

char waitKey(int nMS)
{
    if(nMS != CImageViewUI::SHORT_DISPLAY)
        CEventTimer::pauseLogging();

    return (char)cv::waitKey(nMS);
}

cv::Scalar getContrastingCol(const cv::Scalar & mean)
{
    if(mean[0 /*blue*/] < 220)
        return CV_RGB(255,255,255);
    else
        return CV_RGB(0,0,0);


    cv::Scalar col;
    for(int i=0; i<4; i++)
        col[i] = mean[i] < 128 ? 255 : 0;

    return col;
}

bool CImageViewUI::fastDisplay() const
{
    return displayMode == eFastDisplay || displayMode == eMakeFastVideo ||  displayMode == eMakeSplitVideo;
}

bool CImageViewUI::makingVideo() const
{
    return displayMode == eMakeVideo || displayMode == eMakeFastVideo ||  displayMode == eMakeSplitVideo;
}

void CImageViewUI::setNumBytesPrinted(const size_t nMaxNumBytesPrinted_in)
{
    pRedir->setNumBytesPrinted(nMaxNumBytesPrinted_in);
}

void CImageViewUI::setDebugging()
{
    if(makingVideo())
        throw CException("Not entering debug mode--debug mode off by macro"); //die right away (causes an abort)

    if(bGlobalDebuggingMode)
        return;

    bGlobalDebuggingMode = true;
    tryStartGedit();
    displayMode = eFullDisplay;
    skipList.clear();
    fastList.clear();
    cout << "\n\n############################## ERROR: entering debug mode #########################################\n" << endl;

    showText("ERROR: entering debug mode");
}

char CImageViewUI::showText(const std::string text)
{
    cv::Mat M;
    if(aImageHistory.size()>0) {
        M=aImageHistory.back();
    } else {
        M=cv::Mat(1000,1000, CV_8UC3, cv::Scalar());
    }

    cout << text << endl;
    cv::putText(M, text, cv::Point(20, 150), CV_FONT_HERSHEY_PLAIN, 2.5, CV_RGB(255,255,255), 3);

    return showOneImage(M, "_____", 0);
}
void CImageViewUI::setDisplayMode(const eDisplayMode newDispMode)
{
    const bool bVerbose = true;

    if(!makingVideo()) { //Don't change display mode when working on a video
        displayMode = newDispMode;
        COUT2("Changed display mode to ", displayMode);
    } else {
        COUT2("Didn't change display mode, ", displayMode);
    }
}
CImageViewUI::CImageViewUI() : pRedir(0), nSave(0), displayMode(eFullDisplay), mainThreadId(boost::this_thread::get_id()), aImageHistory(HISTORY_SIZE), videoSize(0,0)
{
    strFolderName = makeNiceFolderName();
    
    strLogfileName = strFolderName + "output.log";
    cerr << "Creating logfile " << strLogfileName << endl;

    strVidFrameDir = strFolderName + "/vidFrames/";

    pRedir = new redirectCout(strLogfileName.c_str(), false /*do not open gedit*/);

    //Just use 1 window
}

void CImageViewUI::setupWindow()
{
    CHECK(displayMode == eNoDisplay, "Shouldn't be displaying anything");
    if(aImageHistory.size()==0) {
        
        WINDOW_TITLE = "Output: q=quit, s=save, click=colour, 1,2,3=display speed, l=open logfile, d=open directory, p=skip this title, f/b=forward/back, i=flip, w=down, a=accel. this title, r=reset lists ";    
        if(IS_DEBUG)
            WINDOW_TITLE += "D";
        const std::string timeAndDate = getTimeAndDate(false);
        WINDOW_TITLE += timeAndDate.substr(timeAndDate.length()-13, 8);
        const bool bScaleToFit = displaySize == eScaleToFitScreen;
        const int nWindowPropFlags = bScaleToFit ? (int)cv::WND_PROP_FULLSCREEN : (int)cv::WINDOW_AUTOSIZE;
        cv::namedWindow(WINDOW_TITLE, nWindowPropFlags);
        cv::moveWindow(WINDOW_TITLE, 0, 0);
        if(bScaleToFit) 
            cv::resizeWindow(WINDOW_TITLE, 2000, 1000);
    }
}

CImageViewUI::~CImageViewUI()
{
    cv::destroyWindow(WINDOW_TITLE);
    delete pRedir;
}

void execCommand(std::string strCommand)
{
    cout << "Executing command: " << strCommand << endl;
    const int nResult = system(strCommand.c_str());
    if(IS_DEBUG) CHECK(-1 == nResult , "Execute command failed");
}

void tryOpen(std::string progname, std::string args)
{
    std::ostringstream ss;

    ss << "pgrep " << progname << " && " << progname << " " << args << std::endl;

    execCommand(ss.str());
}

void CImageViewUI::tryStartGedit() const
{
    if(!boost::filesystem::exists(strLogfileName))
        cout << "Log file " << strLogfileName << " no longer exists" << endl;

    tryOpen("geany", strLogfileName);
}
void CImageViewUI::tryStartNautilus() const
{
    if(!boost::filesystem::exists(strFolderName))
        cout << "Log file " << strFolderName << " no longer exists" << endl;

    //tryOpen("nautilus", strFolderName);
    std::ostringstream ss;
    ss << "nautilus \"" << strFolderName << "\"";
    execCommand(ss.str());
}

void CImageViewUI::saveImage(const cv::Mat & image, const std::string & filename)
{
    std::ostringstream ss;
    ss << strFolderName << "/" << filename;
    const string dirAndFilename = ss.str();
    cout << "Saving " << dirAndFilename << endl;
    cv::imwrite(dirAndFilename, image);
}

//Details on mouse handlers: http://www.cs.iit.edu/~agam/cs512/lect-notes/opencv-intro/index.html
void CImageViewUI::mouseHandlerMat(int event, int x, int y, int flags, void* pParam)
{
    cv::Mat & image = *((cv::Mat *) pParam);
    if(x<0||y<0|| x>=image.cols||y>=image.rows)
        return;

    switch (event) {
    case CV_EVENT_LBUTTONDOWN:
        cout << "(" << x << ", " << y << ") colour: ";

        if (image.depth() == CV_8U) {
            if (image.channels() == 1)
                cout << (int) image.at<uchar> (y, x);
            else if (image.channels() == 3) {
                //CRedGreenBlue rgb;
                cv::Vec3b rgb = image.at<cv::Vec3b > (y, x);
                cout << "RGB: " << (int) rgb[2] << ',' << (int) rgb[1] << ',' << (int) rgb[0] << ' ';

                //cout << colourLabel(image, cv::Point(x,y));
            } else {
                PP(image);
                THROW("number of channels in im not handled");
            }
        } else {
            if (image.channels() == 1)
                if(image.type() == CV_32FC1)
                    cout << image.at<float> (y, x);
                else if(image.type() == CV_64FC1)
                    cout << image.at<double> (y, x);
                else
                    cout << "unknown type" << endl;
            else {
                if(image.type() == CV_32FC3) {
                    //CRedGreenBlue rgb;
                    cv::Vec3f rgb = image.at<cv::Vec3f > (y, x);
                    cout << "RGB: " << rgb[2] << ',' << rgb[1] << ',' << rgb[0] << '\n';
                } else
                    cout << "unknown type" << endl;
            }
        }
        cout << endl;
        break;

    case CV_EVENT_RBUTTONDOWN:

        /*aCorners.push_back(C2dImagePointPx(x, y));
        if (aCorners.size() == 4) {
            cv::Mat Xpoint = cropRect(image, aCorners[0], aCorners[1], aCorners[2], aCorners[3], 50, 50);
            static int s_nId = 0;
            ostringstream ssFilename, ssFilenameSource;
            ssFilename << "crop" << s_nId << ".png";
            string filename = ssFilename.str();
            UI().saveImage(Xpoint, filename);

            for (int i = 0; i < 4; i++) {
                cv::line(image, aCorners[i], aCorners[(i + 1) % 4], CV_RGB(255, 0, 0), 1);
            }

            ssFilenameSource << "source" << s_nId << ".png";
            string filenameSource = ssFilenameSource.str();
            UI().saveImage(image, filenameSource);

            s_nId++;
            aCorners.clear();
        }*/
        break;
    }
}

void addTitleTextToIm(cv::Mat im, const std::string text, const int nRow)
{
    cv::Scalar mean = cv::mean(im.rowRange(im.rows/4,im.rows/4+std::min<int>(im.rows/4, 20)));

    const int nScale = clip<int>(im.rows/160, 3, 8);

    double hScale = 2.5;

    if(hScale<1)
        hScale=1;

    double dTextScale = 1;
    if(text.length() > 25)
        dTextScale -= (text.length() - 25)*0.02;
    if(dTextScale < 0.25)
        dTextScale = 0.25;

    cv::putText(im, text, cv::Point(3*nScale, 4*nScale+6*nScale*nRow), CV_FONT_HERSHEY_PLAIN, clip<double>(0.3*dTextScale*nScale, 0.8, 10), getContrastingCol(mean), clip<int>(cvRound(0.4*dTextScale*nScale), 1, 10));
    cout << text << endl;
}

/* convert for display */
void scaleAndConvert(const cv::Mat & src, cv::Mat & dest, const int nChannel = 0)
{
    double dMax = -1, dMin = -1;
    cv::minMaxLoc(src, &dMin, &dMax);
    cv::Mat imGrey;
    double alpha = 1, beta = 0;
    ostringstream ss;
    if ((dMax > 255 || dMax - dMin < 40 || dMin < 0) && dMin != dMax) {
        //Rescale image
        alpha = 255 / (dMax - dMin);
        beta = -alpha*dMin;

        ss << "Range: " << dMin << " to " << dMax << flush;
    }

    src.convertTo(dest, CV_8U, alpha, beta);
    if(ss.str().length() > 0)
        addTitleTextToIm(dest, ss.str(), nChannel+1);
}

char CImageViewUI::mergeAndShow(const cv::Mat & M1, const cv::Mat & M2, const cv::Mat & M3_optional, const int nWaitMS)
{
    const cv::Mat & M3 = (M3_optional.rows == 0) ? M2 : M3_optional;
    vector<cv::Mat> aMats(3);

    scaleAndConvert(M1, aMats[0],0);
    scaleAndConvert(M2, aMats[1],1);
    scaleAndConvert(M3, aMats[2],2);

    cv::Mat outIm;
    cv::merge(aMats, outIm);
    return showOneImage(outIm, "", nWaitMS);
}
bool CImageViewUI::isMT() const
{
    return boost::this_thread::get_id() != mainThreadId;
}

std::string trimTitle(const std::string & TITLE)
{
    const size_t nDashPos = TITLE.find('-');
    if(nDashPos == std::string::npos || nDashPos == 0 || TITLE.length() < 3 || nDashPos == TITLE.length()-1)
        return TITLE;

    std::string trimmed(TITLE.begin(), TITLE.begin()+(nDashPos-1));
    return trimmed;
}

char CImageViewUI::showOneImage(const cv::Mat & im_in, const string TITLE, const int nWaitMS)
{
    CHECK(!im_in.data || !im_in.rows, "showOneImage: image is zero-sized");

    if(displayMode == eNoDisplay)
        return 0;
        
    if(isMT()) {
        cout << "NOT displaying frame in call from thread " << boost::this_thread::get_id() << endl;
        return 0;
    }

    if(skipList.find(trimTitle(TITLE)) != skipList.end())
        return 0;

    int nWait = nWaitMS;
    if(fastList.find(trimTitle(TITLE)) != fastList.end())
        nWait = SHORT_DISPLAY;

    if (nWaitMS == 0) pp(TITLE, im_in);

    cv::Mat im;
    if(displaySize == eScaleDown && im_in.rows >= getMaxSize().height) //Can segfault here if construction of params() object is incomplete (some parameter exception thrown during contruction)
        cv::resize(im_in, im, cv::Size(), 0.5, 0.5, cv::INTER_NEAREST);
    else
        im = im_in;

    if(fastDisplay())
        nWait = SHORT_DISPLAY;
    else if (displayMode == eSlowDisplay || displayMode == eMakeVideo)
        nWait = (nWaitMS==0) ? 1000 : nWaitMS;

    setupWindow(); //Only once we're going to display something

    //aCorners.clear();
    if (im.channels() == 1) {
        cv::Mat imGrey;
        scaleAndConvert(im, imGrey);

        return showOneImage_int(imGrey, TITLE, nWait, im);
    } else if (im.channels() == 3 && im.depth() == CV_32F) {
        std::vector<cv::Mat> aChannels;
        cv::split(im, aChannels);
        double dTotalMin = HUGE, dTotalMax = -HUGE;
        for (int i = 0; i < 3; i++) {
            double dMin, dMax;
            cv::minMaxIdx(aChannels[i], &dMin, &dMax);
            if (dMin < dTotalMin)
                dTotalMin = dMin;
            if (dMax > dTotalMax)
                dTotalMax = dMax;
        }

        //Now scale to 0...255
        double alpha = 1, beta = 0;

        if (dTotalMax > 255 || dTotalMax - dTotalMin < 40 || dTotalMin < 0) {
            //Rescale image
            alpha = 255 / (dTotalMax - dTotalMin);
            beta = -alpha*dTotalMin;
            //ss << "Range: " << dMin << " to " << dMax << flush;
        }
        cv::Mat imBGR;
        im.convertTo(imBGR, CV_8U, alpha, beta);

        return showOneImage_int(imBGR, TITLE, nWait, im);
    } else {
        return showOneImage_int(im, TITLE, nWait, im);
    }
}

const char * CImageViewUI::getVidFramePath() const
{
    CHECK_P(strVidFrameDir.length() < 2, strVidFrameDir, "strVidFrameDir not init");
    boost::filesystem::create_directory(strVidFrameDir.c_str());
    CHECK_P(!boost::filesystem::exists(strVidFrameDir), strVidFrameDir, "strVidFrameDir does not exist");
    return strVidFrameDir.c_str();
}
//Im may have been scaled/etc. actualIm is the original, so we can click and get actual vals
char CImageViewUI::showOneImage_int(const cv::Mat & image, const string TITLE, const int nWaitMS, const cv::Mat & actualIm)
{
    const bool bVerbose = true;
    
    cv::Mat im = image.clone(); //as we are about to draw on it

    if(maxSize.height < im.rows)
        maxSize=im.size();

    addTitleTextToIm(im, TITLE, 0);

    aImageHistory.push_back(im);

    cv::imshow(WINDOW_TITLE, im);
    cv::setMouseCallback(WINDOW_TITLE, &CImageViewUI::mouseHandlerMat, (void *) &actualIm);

    int nBufferPosition = (int)aImageHistory.size()-1;

    char key = 0;
    for(;;) {
        enum { eForward, eBack, eWait, eDone } eState = eWait;

        do {
            key = (char)std::tolower(waitKey(nWaitMS));

            if(key>0 && key != ' ') COUT2("Key press: ", key);

            if ('q' == key) {
                COUT("Exiting");
                throw CException();
            } else if ('s' == key) {
                {
                    int nFrame = 0;
                    
                    std::ostringstream saveFN;
                    saveFN << strFolderName << "/" << "frame" << std::setw(8) << std::setfill('0') << nFrame << "_" << TITLE << nSave++ << ".png";
                    cv::imwrite(saveFN.str(), (nBufferPosition < (int(aImageHistory.size())-1)) ? aImageHistory[nBufferPosition] : im);
                }
            } else if ('1' == key) {
                setDisplayMode(eFullDisplay);
                COUT2("Changed display mode to ", displayMode);
                eState = eDone;
            } else if ('2' == key) {
                setDisplayMode(eSlowDisplay);
                COUT2("Changed display mode to ", displayMode);
                eState = eDone;
            } else if ('3' == key) {
                setDisplayMode(eFastDisplay);
                COUT2("Changed display mode to ", displayMode);
                eState = eDone;
            } else if ('p' == key) {
                skipList.insert(trimTitle(TITLE));
                COUT2("Added to skip list: ", trimTitle(TITLE));
                eState = eDone;
            } else if ('a' == key) {
                fastList.insert(trimTitle(TITLE));
                COUT2("Added to fastShow list: ", trimTitle(TITLE));
                eState = eDone;
            } else if ('r' == key) {
                fastList.clear();
                skipList.clear();
            } else if ('b' == key) {
                eState = eBack;
                nBufferPosition--;
                if(nBufferPosition<0)
                    nBufferPosition = 0;
            } else if ('f' == key) {
                eState = eForward;
                nBufferPosition++;
                if(nBufferPosition >= (int)aImageHistory.size() )
                    nBufferPosition = (int)aImageHistory.size()-1;
            } else if ('i' == key) {
                eState = eForward;
                COUT("Flipping image");
                cv::flip(aImageHistory[nBufferPosition], aImageHistory[nBufferPosition], 0);
            } else if (('w' == key) && aImageHistory[nBufferPosition].rows>1000 ) {
                eState = eForward;
                COUT("Moving image down");
                aImageHistory[nBufferPosition] = aImageHistory[nBufferPosition](cv::Range(aImageHistory[nBufferPosition].rows-1000, aImageHistory[nBufferPosition].rows), cv::Range::all());
            } else if ('l' == key) {
                COUT("Opening logfile if geany is open...");
                tryStartGedit();
            } else if ('d' == key) {
                COUT("Opening directory in Nautilus...");
                tryStartNautilus();
            } else if(key==-1 || ' ' == key || (key >= 'a' && key <= 'z') || (key >= '0' && key <= '9')) {
                if(key > 0 && key != ' ') COUT2("Returning key: ", key)
                    eState = eDone; //We want to return the key
            }

        } while (eState == eWait);

        if(eState == eDone)
            break;

        cout << "Showing image " << (nBufferPosition+1) << "/" << aImageHistory.size() << " from buffer" << endl;
        cv::imshow(WINDOW_TITLE, aImageHistory[nBufferPosition]);
        cv::setMouseCallback(WINDOW_TITLE, &CImageViewUI::mouseHandlerMat, (void *) &aImageHistory[nBufferPosition]);
    }

    if (makingVideo()) {
        static int s_nFrameID = 100000;
        const int nLongInterval = (displayMode == eMakeVideo) ? 60 : 1;
        const int nFrames = (nWaitMS >= LONG_DISPLAY) ? nLongInterval : 1;
        if(videoSize.area() == 0)
            videoSize = im.size();

        //breaks in Windows old opencv const bool bVerbose = true;
        //COUT(videoSize);
        //COUT(im.size());

        if(im.cols > videoSize.width)
            im.cols = videoSize.width; //crop larger (3x) frames

        if(im.size() == videoSize) {
            for (int nFrame = 0; nFrame < nFrames; nFrame++) {
                std::ostringstream fileName;
                fileName << getVidFramePath();
                if(displayMode == eMakeSplitVideo)
                    fileName << TITLE;
                fileName << setw(6) << setfill('0') << s_nFrameID << ".png";
                std::string strName = fileName.str();
                cv::imwrite(strName, im);

                if(pVideoWriter)
                    *pVideoWriter << im;

                s_nFrameID++;
            }
        }
    }

    return key;
}

void CFrameJoiner::addFrame(const cv::Mat & frame, const std::string & title, const int cam /*assumes 3 cameras*/)
{
    if(UI().isMT())
        return;

    if(joinedFrame.rows == 0)
    {
        joinedFrame = cv::Mat(frame.rows, cvRound(frame.cols*(UI().crop3rdFrame() ? 2.5 : 3)), CV_8UC3, cv::Scalar());
        nNumFramesAdded = 0;
    }
    
    nNumFramesAdded++;

    if(!UI().crop3rdFrame() || cam != 2) {
        cv::Mat range = joinedFrame(cv::Range::all(), cv::Range(cam * frame.cols,(cam+1)*frame.cols));
        frame.copyTo(range);
    } else {
        cv::Mat range = joinedFrame(cv::Range::all(), cv::Range(cam * frame.cols,cvRound((cam+0.5)*frame.cols)));
        cv::Range rangeOfFrame(0, frame.cols/2);
        cv::Mat subFrame = frame(cv::Range::all(), rangeOfFrame);
        subFrame.copyTo(range);
    }
    
    if(cam != 0) 
    {
        cv::Mat whiteStripe = joinedFrame(cv::Range::all(), cv::Range(cam * frame.cols, cvRound(cam*frame.cols+10)));
        whiteStripe.setTo(cv::Scalar(255,255,255));
    }    

    if(nNumFramesAdded == 3)
    {
        show(title);
        nNumFramesAdded = 0;
    }
}

void CFrameJoiner::show(const string & title) const
{
    UI().showOneImage(joinedFrame, title, 0);
}

inline /*duplicate fn*/ void drawText(cv::Mat & image, const cv::Point location, const std::string label, const cv::Scalar colour= CV_RGB(255,255,255), const int nSize=1, const bool bOutline=false)
{
    if(bOutline)
        cv::putText(image, label, location, CV_FONT_HERSHEY_PLAIN, 1.2*nSize, cv::Scalar(), nSize+6);
    cv::putText(image, label, location, CV_FONT_HERSHEY_PLAIN, 1.2*nSize, colour, nSize);
} 

void outputVideoTitles(const cv::Size & size) //match alphabetically with videoRenamerScript
{
    //ALWAYS_VERBOSE;
    //cv::Size size = sys().getCameraData(sys().eCameraAtStereoOrigin).getImSize();
    cv::Mat TITLE(size.height, cvRound(size.width*2.5), CV_8UC3, cv::Scalar(255,255,255));
    const double dSize = 4;

    drawText(TITLE, cv::Point (250,400), "Vision-based automated pruning", cv::Scalar(), cvRound(dSize+1));
    //drawText(TITLE, cv::Point (250,550), "3D reconstruction pipeline", cv::Scalar(), cvRound(dSize+1));
    drawText(TITLE, cv::Point (250,750), "Tom Botterill et al.", cv::Scalar(), cvRound(dSize));
    drawText(TITLE, cv::Point (250,900), "University of Canterbury, Christchurch, NZ", cv::Scalar(), cvRound(dSize));
    drawText(TITLE, cv::Point (250,1050), "May 2015", cv::Scalar(), cvRound(dSize));
    for(int i=0; i<2; i++)
        UI().showOneImage(TITLE, "A");



    TITLE.setTo(cv::Scalar(255,255,255));

    drawText(TITLE, cv::Point (250,400), "Vine pruning robot", cv::Scalar(), cvRound(dSize+1));
    for(int i=0; i<2; i++)
        UI().showOneImage(TITLE, "A");



    TITLE.setTo(cv::Scalar(255,255,255));
    drawText(TITLE, cv::Point (250,400), "3D vine reconstruction pipeline", cv::Scalar(), cvRound(dSize+1));
    drawText(TITLE, cv::Point (250,700), "May 2015", cv::Scalar(), cvRound(dSize));
    for(int i=0; i<25; i++)
        UI().showOneImage(TITLE, "A");
    UI().showOneImage(TITLE, "Z");

    TITLE.setTo(cv::Scalar(255,255,255));
    drawText(TITLE, cv::Point (250,400), "Match 2D vines between views", cv::Scalar(), cvRound(dSize+1));
    drawText(TITLE, cv::Point (250,700), "Green = matches selected by minimum", cv::Scalar(), cvRound(dSize));
    drawText(TITLE, cv::Point (250,850), "        expected loss corresponder", cv::Scalar(), cvRound(dSize));
    for(int i=0; i<25; i++)
        UI().showOneImage(TITLE, "B");

    TITLE.setTo(cv::Scalar(255,255,255));
    drawText(TITLE, cv::Point (250,400), "Detect 2D canes in each image", cv::Scalar(), cvRound(dSize+1));
    //drawText(TITLE, cv::Point (250,700), "", cv::Scalar(), cvRound(dSize));
    for(int i=0; i<25; i++)
        UI().showOneImage(TITLE, "C");

    TITLE.setTo(cv::Scalar(255,255,255));
    drawText(TITLE, cv::Point (250,400), "Incrementally reconstruct 3D vines", cv::Scalar(), cvRound(dSize+1));
    //drawText(TITLE, cv::Point (250,700), "", cv::Scalar(), cvRound(dSize));
    for(int i=0; i<25; i++)
        UI().showOneImage(TITLE, "R");
}

void logToOneTSV(const std::string & path, const std::ostringstream & ss, const std::ostringstream & ss_headers) 
{
    const bool bNew = !boost::filesystem::exists(path);
    std::ofstream out(path.c_str(), std::ios_base::app);
    if(bNew)
        out << ss_headers.str();
    out << ss.str();
    
    //cout << path << endl << ss_headers.str() << ss.str();
}

//Adds headers if file doesn't exist yet
void logToTSV(const string label, const std::ostringstream & ss, const std::ostringstream & ss_headers, const bool bLogGlobal)
{
    const std::string pathOne = UI().getLogPath() + "/" + label + ".tsv";

    logToOneTSV(pathOne, ss, ss_headers);
    
    if(bLogGlobal)
    {
        const std::string pathAll = "all" + label + ".tsv";
        logToOneTSV(pathAll, ss, ss_headers);
    }
}
