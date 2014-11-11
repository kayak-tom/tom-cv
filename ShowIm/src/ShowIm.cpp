// The ShowIm program displays images. It is run by BoWSLAM to avoid the issues with 
// OpenCV's showim function not working in multithreaded applications.
//
// Sits in a poll loop, and displays bmp files saved by BoWSLAM.
//

#include <iostream>
#include <fstream>
#include "util/opencv.h"
#include "util/opencv_highgui.h"
#include <boost/filesystem/operations.hpp>

#ifndef USE_OLD_OPENCV
#include "opencv2/opencv.hpp"
#endif

using namespace std;

void getWinPos(int nWindows, int & nWinPosX, int & nWinPosY) {
    //Tile them 3x3
    const int nTile = 3;
    int dx = 1680 / nTile;
    int dy = 1050 / nTile;
    nWinPosX = dx * (nWindows % nTile);
    nWinPosY = dy * (nWindows / nTile);
}

int main(int argc, char* argv[]) {
    cout << "Image display process started" << endl; // prints !!!Hello World!!!

    ofstream myLog("showImlog"); //Kill any other display threads
    ofstream signalDone("done"); //Kill any other display threads
    signalDone.close();
    int count = 3;
    while (boost::filesystem::exists("done") && count--) {
        cvWaitKey(50);
        //cout << "waiting for im thread to die...\n";
    }
    remove("done");
    remove("quit");

    int nWindows = 0;
    char path[20];
    bool bQuit = false;
    while (!boost::filesystem::exists(boost::filesystem::path("done")) && !bQuit) {
        for (int nWindow = 0; nWindow <= nWindows; nWindow++) {
            sprintf(path, /*20,*/ "%d.bmp", (int) nWindow);

            if (boost::filesystem::exists(boost::filesystem::path(path))) {
                if (nWindows == nWindow) {
                    cvNamedWindow(path);
                    int nWinPosX, nWinPosY;
                    getWinPos(nWindows, nWinPosX, nWinPosY);
                    cvMoveWindow(path, nWinPosX, nWinPosY);

                    nWindows++;
                }
                IplImage * pIm = cvLoadImage(path);
                cvShowImage(path, pIm);

                if (argc > 1) cvWaitKey(0);

                remove(path);

                cvReleaseImage(&pIm);
            }
        }
        int nKey = cvWaitKey(50);
        if (nKey == 'q' || nKey == 'Q') {
            bQuit = true;
            ofstream signalQuit("quit"); //Kill main process
            signalQuit.close();
        } else if (nKey == 'p' || nKey == 'P')
            cvWaitKey(0);

        //Kill image loader once BoWSLAM has died:
        if (-1 == system("ps x | egrep 'BoWSLX?AM' > temp; if [ -s temp ]; then rm temp; else mv temp done; fi"))
            myLog << "Error killing image loader thread\n";
    }
    if (!bQuit) {
        char ch;
        ifstream donemsg("done");
        while (donemsg && donemsg.get(ch))
            myLog.put(ch);
    }
    remove("done");
    cvDestroyAllWindows();

    myLog.close();

    return 0;
}
