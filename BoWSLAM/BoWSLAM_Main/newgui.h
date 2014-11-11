/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once
#include <boost/filesystem.hpp>
#include "logging/prfxstream.h"
#include "FixEPS.h"

#ifndef USE_OLD_OPENCV
#include "opencv2/opencv.hpp"
#endif

#include "opencv/highgui.h"

class CRunGuiAp {

    enum eGuiMode {
        eNoGui, eSaveOnly, eSaveDisplay
    };
    const eGuiMode eDisplayMode;

    void tryKill() {
        ofstream signalDone("done");
        signalDone.close();
        int count = 3;
        while (boost::filesystem::exists("done") && count--) {
            usleep(50000);
            //cout << "waiting for im thread to die...\n";
        }
    };
public:

    enum eWindows {
        eWindowMap = 0, eWindowIm = 1, eWindowCorresp = 2, eWindowSimIm = 3, eAnglePlot = 4, eHeightPlot = 5, eNumWindows
    }; //check matches with aszSaveFrame below...
private:
    int anWinNums[eNumWindows], nextWinNum;
    static const int NONE = -1;

    void getPath(int winNum, char * pcBuff, int len) {
        sprintf_s(pcBuff, len, "%d.bmp", (int) winNum);
    }
    char * szSaveLoc;
public:

    bool showImages() const {
        return eDisplayMode != eNoGui;
    };

    CRunGuiAp(bool display = true) : eDisplayMode(display ?
#ifdef RUN100
    eSaveOnly
#else
    eSaveDisplay
#endif
    : eNoGui), szSaveLoc(0) {
        szSaveLoc = strdup(getTimeAndDate(false));

        char szMapDir[100];
        sprintf_s(szMapDir, 100, "maps/%s", szSaveLoc);
        boost::filesystem::create_directory("maps");
        boost::filesystem::create_directory(szMapDir);
        sprintf_s(szMapDir, 100, "maps/%s/aa_terminal", szSaveLoc);
        boost::filesystem::create_directory(szMapDir); //to open terminal in thunar

        if (eDisplayMode == eSaveDisplay) {
            tryKill();
            remove("done");
            remove("quit");
            const char * szShowImPath = "./ShowIm";
 
            const int NUM_SHOWIM_PATHS = 2;
            const char * aszPaths[] = {
                "../ShowIm/Release/ShowIm", //Eclipse
                "../ShowIm/dist/Release/GNU-Linux-x86/showim"
            }; //Netbeans

            if (!boost::filesystem::exists(szShowImPath)) {

                for (int i = 0; i < NUM_SHOWIM_PATHS; i++) {
                    if (boost::filesystem::exists(aszPaths[i])) {
                        cout << "Found ShowIm program at " << aszPaths[i] << endl;
                        boost::filesystem::copy_file(aszPaths[i], szShowImPath);
                    }
                }
            }

            cout << "Starting " << szShowImPath << "...";

            if (!boost::filesystem::exists(szShowImPath)) {
                cout << "ERROR: Image display program not found. Images will not be displayed" << endl;
                cout << "Failed to start image loader thread. Build ShowIm (release) and check it is in the correct place" << endl;
                cout << "Searching for: " << aszPaths[0] << ", " << aszPaths[1] << ", " << szShowImPath << endl;
            } else {

                std::string showImCommand(szShowImPath);
                showImCommand.append(" &");
                if (0 != system(showImCommand.c_str()))
                    cout << "Error spawning image display thread :" << showImCommand << endl;
                else {
                    cout << "Success\n";
                }
            }
        }
        for (int i = 0; i < eNumWindows; i++) {
            anWinNums[i] = NONE;

            char path[20];
            getPath((int) i, path, 20);
            remove(path);
        }

        nextWinNum = 0;
    };

    ~CRunGuiAp() {
        if (showImages()) tryKill();
        free(szSaveLoc);
    };

    void showIm(const IplImage * pIm, eWindows window, int nSaveFrame = 0) {
        if (!showImages()) return;

        //boost::mutex::scoped_lock scoped_lock(mlocationIds);
        //Wait until file gone
        if (anWinNums[window] == NONE)
            anWinNums[window] = nextWinNum++;
        int winNum = anWinNums[window];
        char path[20];
        getPath((int) winNum, path, 20);

        if (eDisplayMode == eSaveDisplay) {
            int timeout = 100;
            while (boost::filesystem::exists(boost::filesystem::path(path)) && timeout > 0) {
                usleep(50000);
                timeout--;
            }
        }
        if (!nSaveFrame) {
            if (eDisplayMode == eSaveDisplay)
                cvSaveImage(path, pIm);
        } else {
            char frameSavePath[100];
            static const char * aszSaveFrame[eNumWindows] = {"Map", "Im", "Corresp", "", "AnglePlot", "HeightPlot"};
            sprintf_s(frameSavePath, 100, "maps/%s/%s-%06d.png", szSaveLoc, aszSaveFrame[window], nSaveFrame);
            cvSaveImage(frameSavePath, pIm);
            if (eDisplayMode == eSaveDisplay) {
                if (!boost::filesystem::exists(boost::filesystem::path(path))) {
#ifdef __GNUC__
                    boost::filesystem::create_symlink(frameSavePath, path);
#else
                    boost::filesystem::copy_file(frameSavePath, path);
#endif
                } else
                    cout << "Warning: " << path << " already exists\n";
            }
        }
    };

    void saveEps(LibBoard::Board & board, int nSaveFrame) {
        if (anWinNums[eWindowMap] == NONE)
            anWinNums[eWindowMap] = nextWinNum++;
        int winNum = anWinNums[eWindowMap];
        char path[20];
        getPath((int) winNum, path, 20);

        char frameSavePath[100];
        static const char * aszSaveFrame[eNumWindows] = {"Map", "Im", "Corresp", "", "AnglePlot", "HeightPlot"};
        sprintf_s(frameSavePath, 100, "maps/%s/%s-%d.eps", szSaveLoc, aszSaveFrame[eWindowMap], nSaveFrame);
        LibBoard::Rect bb = board.boundingBox();
        //cout << "Bounding box: " << bb.left << ',' << bb.width << ',' << bb.top << ',' << bb.height << endl;
        const float MARGIN = max<float>(bb.width, bb.height) * 0.1;
        board.setPenColor(LibBoard::Color(255, 255, 255));
        board.drawRectangle(bb.left - MARGIN, bb.top + MARGIN, bb.width + 2 * MARGIN, bb.height + 2 * MARGIN);
        board.saveEPS(frameSavePath, 300, 300, 40);

        /*sprintf_s(frameSavePath, 100, "maps/%s/%s-%d-old.eps", szSaveLoc, aszSaveFrame[eWindowMap], nSaveFrame);
        board.saveEPS(frameSavePath, LibBoard::Board::BoundingBox, 10);
        if(boost::filesystem::exists(frameSavePath))
                CFixEPS::fixEpsBB(frameSavePath);
        else
                cout << "Error saving EPS drawing\n";*/
    };

    void getClusterFilename(int nClusterNum, int nLevel, char * frameSavePath) const {
        sprintf_s(frameSavePath, 1000, "maps/%s/FirstClusters%d-%d.png", szSaveLoc, nClusterNum, nLevel);
    }

    void getToroMapFilename(int nSaveFrame, int nToroConnectivity, char * frameSavePath) const {
        sprintf_s(frameSavePath, 100, "maps/%s/ToroMap%d-%d.graph", szSaveLoc, nToroConnectivity, nSaveFrame);
    }

    void getToroMapFilename2d(int nSaveFrame, int nToroConnectivity, char * frameSavePath) const {
        sprintf_s(frameSavePath, 100, "maps/%s/ToroMap2d%d-%d.graph", szSaveLoc, nToroConnectivity, nSaveFrame);
    }

    void getLogFilename(char * frameSavePath) const {
        sprintf_s(frameSavePath, 1000, "maps/%s/%s.log", szSaveLoc, szSaveLoc);
    }

    void getRDFilename(char * frameSavePath) const {
        sprintf_s(frameSavePath, 1000, "maps/%s/CorrectRD.png", szSaveLoc);
    }

    const char * outputDir() const {
        return szSaveLoc;
    }
};

