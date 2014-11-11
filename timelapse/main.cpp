#include <opencv2/opencv.hpp>
#include <opencv/highgui.h>
//#include <opencv2/highgui/highgui.hpp>
#include <boost/circular_buffer.hpp>

typedef boost::circular_buffer<cv::Mat> TImBuffer;

static const int HISTORY_SIZE=300;

int main(int argc, char **argv)
{
    cv::VideoCapture vidReader;
    //vidReader.open("/media/sdb1/timelapse/eth/ewap_dataset/seq_hotel/seq_hotel.avi"); //"/media/sdb1/goodVideos/P1040913-Working.MOV"
    ///media/sdb1/timelapse$ mencoder Waves.wmv -ovc raw -nosound -o waves.avi
    vidReader.open("/media/sdb1/timelapse/waves.avi");
    TImBuffer aFrames(HISTORY_SIZE);

    cv::Mat averageFrame8u, averageFrame32f;
    vidReader >> averageFrame8u; //For initialisation;

    averageFrame8u.convertTo(averageFrame32f, CV_32FC3);
    
/*cv::VideoWriter vidWriter("/media/sdb1/goodVideos/test.avi",0, 60, averageFrame.size());
    const int fourcc = CV_FOURCC('P', 'I', 'M', '1');//CV_FOURCC('M','J','P','G')CV_FOURCC_DEFAULT
    vidWriter.open("/media/sdb1/goodVideos/timelapse.MOV", fourcc, 60, averageFrame.size()); */
    cv::Mat timeLapsedFrame32f;
    
    enum { eVidFrame, eAverageFrame, eTimelapseFrame, eAveragePlusLatest } display = eVidFrame;
    
    const float dNewWeight = 0.03;
    
    
    for(int nIndex = 0; ;nIndex++)
    {
        cv::Mat vidFrame8u;
        vidReader >> vidFrame8u;
        if(vidFrame8u.rows == 0)
            break;
            
        aFrames.push_back(cv::Mat());
        vidFrame8u.convertTo(aFrames.back(), CV_32FC3);
        averageFrame32f = averageFrame32f*(1.0f-dNewWeight) + aFrames.back() * dNewWeight;
        
        
        if(display == eTimelapseFrame)
        {
            timeLapsedFrame32f = 0.4f*aFrames.back();
            for(float i=0; i<3.0f; i++)
            {
                int nIndex = (i*aFrames.size())/4;
                timeLapsedFrame32f += (0.1f*(i+1.0f)) * aFrames[nIndex];
            }
        }
            
        cv::Mat useFrame;
        switch(display)
        {
            case eAverageFrame:
                averageFrame32f.convertTo(useFrame, CV_8UC3);
                break;
            case eAveragePlusLatest:
                {
                    cv::Mat temp = (0.8*averageFrame32f + 0.2*aFrames.back());
                    temp.convertTo(useFrame, CV_8UC3);
                }
                break;
            case eVidFrame:
                useFrame = vidFrame8u;
                break;
            case eTimelapseFrame:
                timeLapsedFrame32f.convertTo(useFrame, CV_8UC3);;
                break;
        }

        std::ostringstream ss; ss << nIndex;
        
        if(nIndex%400==0)
        {
            std::string outfile = "/media/sdb1/timelapse/eth_av" + ss.str() + ".png";
            cv::imwrite(outfile, useFrame);
        }
        
        cv::Mat(useFrame, cv::Range(0,60), cv::Range(0,200)).setTo(cv::Scalar());
        cv::putText(useFrame, ss.str(), cv::Point(20,50), cv::FONT_HERSHEY_PLAIN, 4, CV_RGB(255,255,255));
        cv::imshow("Frame", useFrame);
        if(cv::waitKey(1) == 'q')
            break;
            
        //vidWriter << *pUseFrame;
    }
	
	return 0;
}
