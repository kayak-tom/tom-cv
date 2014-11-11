#include "exception.h"
#include "calibration.h"
#include <fstream>
#include "convert.h"
#include "location.h"
#include <boost/filesystem/operations.hpp>
#include <string>
#include <stdio.h>

using namespace boost::filesystem;
using namespace std;

void CCamCalibMatrix::setId()
{
    adK[0] = adK[4] = adK[8] = 1;
    adK[1] = adK[2] = adK[3] = adK[5] = adK[6] = adK[7] = 0;

    nRDcoeffs=0;
    init_inv();
}
void CCamCalibMatrix::init_inv()
{
    K_inv.init(*this);
}

//Surely this is in boost or something... Windows doesn't like c:\path\..\path.cfg
void goUp1Level(const char * szPath, char * szPathUp)
{
    char const * pcEnd = szPath+strlen(szPath);
    while (pcEnd > szPath && (*pcEnd == '\\' || *pcEnd == '/'))
        pcEnd--;
    if(pcEnd == szPath)
    {
        //Length 0 or 1
        sprintf_s(szPathUp, 4, "../");
    }
    else
    {
        while (pcEnd > szPath && (*pcEnd != '\\' && *pcEnd != '/'))
            pcEnd--;

        int len = (int)(pcEnd-szPath);
        for(int i=0; i < len; i++)
            szPathUp[i] = szPath[i];

        szPathUp[len]=0;
    }
}

void CCamCalibMatrix::init(const char * szCalibFolderName, const int SCALE_DOWN, const int CROP)
{
#define CALIB_FORMAT "Camera calibration file format: 3x3 matrix (intrinsic calibration matrix), then number of radial distortion coefficients, then list of radial distortion coefficients. All seperated by whitespace. Example:\n"\
    "739.03483  0         403.31196\n"\
    "0          735.00067 257.79306\n"\
    "0          0         1\n\n"\
    "1\n"\
    "-0.06887\n"

    char szCalibFilename[1000],  szCalibFilename2[1000];
    std::ifstream calibFile;

    if(!exists(szCalibFolderName))
    { 
        cout << "Camera calibration filename location does not exist; looking for: " << szCalibFolderName << "" << endl;
        cout << "(" CALIB_FORMAT ")";
        THROW("Calibration file location not found");
    }
    else if(is_directory(szCalibFolderName))
    {
        sprintf_s(szCalibFilename, 1000, "%s/calib.txt", szCalibFolderName);
        boost::filesystem::path FN(szCalibFilename);
        calibFile.open(FN.string().c_str());

        //calibFile.open(szCalibFilename);

        if(!calibFile.is_open() || !calibFile.good())
        {
            //sprintf_s(szCalibFilename2, 1000, "%s/../calib.txt", szCalibFolderName);
            //boost::filesystem::path FN2(szCalibFilename2);
            //calibFile.open(FN2.native_file_string().c_str());

            char szFolderNameUp[1000];
            goUp1Level(szCalibFolderName, szFolderNameUp);
            sprintf_s(szCalibFilename2, 1000, "%s/calib.txt", szFolderNameUp);
            //calibFile.open(szCalibFilename2);
            boost::filesystem::path FN2(szCalibFilename2);
            calibFile.open(FN2.string().c_str());

            if(!calibFile.is_open() || !calibFile.good())
            {
                cout << "Cannot find camera calibration file 'calib.txt'; looking in " << szCalibFilename << " and " << szCalibFilename2 << endl;
                cout << "(" CALIB_FORMAT ")";

                sprintf_s(szCalibFilename, 1000, "%s", szCalibFolderName);
                calibFile.open(szCalibFilename);
                if(!calibFile.is_open() || !calibFile.good())
                {
                    cout << "Looking for calibration file " << szCalibFilename << endl;
                    THROW("Calibration filename supplied but file not found");
                }
                //THROW("Camera calibration file not found");
            }
            //else
              //  sprintf_s(szCalibFilename, 1000, "%s", szCalibFilename2);
        }
    }
    else //We have been given a path to a file--either a calibration file or a video file, with calib.txt in the same folder.
    {
        boost::filesystem::path file(szCalibFolderName);
        string path = szCalibFolderName;
        if(file.extension() != ".txt")
        {
            path=file.remove_filename().string();
            path.append("/calib.txt");
        }
        calibFile.open(path.c_str());

        if(!calibFile.is_open() || !calibFile.good())
        {
            cout << "Looking for calibration file " << path << endl;
            THROW("Filename supplied but calibration file not found in same location");
        }
    }

    if(!calibFile.good() || calibFile.eof())
    {
        THROW( "Error opening calibration file");
    }

    for(int nRow=0;nRow<3;nRow++)
        for(int nCol=0;nCol<3;nCol++)
        {
            calibFile >> adK[3*nRow + nCol];
            if(calibFile.eof())
                THROW("Error reading camera calibration file 'calib.txt'\n" CALIB_FORMAT);
        }

    calibFile >> nRDcoeffs;
    for(int i=0;i<nRDcoeffs; i++)
        calibFile >> adRD[i];

    calibFile.close();

    if(IS_DEBUG) CHECK(adK[0] <= 0 || adK[2] <= 0 || adK[4] <= 0 || adK[5] <= 0, "Elements 00, 02, 11, 12 of K must be positive");
    if(IS_DEBUG) CHECK(adK[6] != 0 || adK[7] != 0 || adK[8] != 1 , "currently bottom row of K must be 0 0 1");

    if(SCALE_DOWN > 1)
    {
        cout << "Adjusting calibration matrix for scaled-down images...\n";
        for(int i=0;i<6;i++)
            adK[i] /= SCALE_DOWN;
    }
    if(CROP > 0)
    {
        cout << "Cropping calibration matrix...\n";
        adK[2] -= CROP;
        adK[5] -= CROP;
    }

    init_inv();
}

void CCamCalibMatrix::testK(const int nWidth, const int nHeight) const
{
    int nEstCentreXdisp = nWidth/2 - (int)adK[0+2];
    int nEstCentreYdisp = nHeight/2- (int)adK[3+2];

    double dDisp = sqrt((double)(sqr(nEstCentreXdisp) + sqr(nEstCentreYdisp)));

    if(IS_DEBUG) CHECK( dDisp > nWidth / 5, "(Heuristic check) Cam centre is a bit far from calibration matrix's camera centre, is the right file in use?" )
}

void CCamCalibMatrix::CInvCamCalibMatrix::init(const CCamCalibMatrix & K) 
{
    double determinant =    +K.adK[0*3 + 0]*(K.adK[1*3 + 1]*K.adK[2*3 + 2]-K.adK[2*3 + 1]*K.adK[1*3 + 2])
                            -K.adK[0*3 + 1]*(K.adK[1*3 + 0]*K.adK[2*3 + 2]-K.adK[1*3 + 2]*K.adK[2*3 + 0])
                            +K.adK[0*3 + 2]*(K.adK[1*3 + 0]*K.adK[2*3 + 1]-K.adK[1*3 + 1]*K.adK[2*3 + 0]);

    double invdet = 1/determinant;
    adK_inv[0 + 0*3] =  (K.adK[1*3 + 1]*K.adK[2*3 + 2]-K.adK[2*3 + 1]*K.adK[1*3 + 2])*invdet;
    adK_inv[1 + 0*3] = -(K.adK[0*3 + 1]*K.adK[2*3 + 2]-K.adK[0*3 + 2]*K.adK[2*3 + 1])*invdet;
    adK_inv[2 + 0*3] =  (K.adK[0*3 + 1]*K.adK[1*3 + 2]-K.adK[0*3 + 2]*K.adK[1*3 + 1])*invdet;
    adK_inv[0 + 1*3] = -(K.adK[1*3 + 0]*K.adK[2*3 + 2]-K.adK[1*3 + 2]*K.adK[2*3 + 0])*invdet;
    adK_inv[1 + 1*3] =  (K.adK[0*3 + 0]*K.adK[2*3 + 2]-K.adK[0*3 + 2]*K.adK[2*3 + 0])*invdet;
    adK_inv[2 + 1*3] = -(K.adK[0*3 + 0]*K.adK[1*3 + 2]-K.adK[1*3 + 0]*K.adK[0*3 + 2])*invdet;
    adK_inv[0 + 2*3] =  (K.adK[1*3 + 0]*K.adK[2*3 + 1]-K.adK[2*3 + 0]*K.adK[1*3 + 1])*invdet;
    adK_inv[1 + 2*3] = -(K.adK[0*3 + 0]*K.adK[2*3 + 1]-K.adK[2*3 + 0]*K.adK[0*3 + 1])*invdet;
    adK_inv[2 + 2*3] =  (K.adK[0*3 + 0]*K.adK[1*3 + 1]-K.adK[1*3 + 0]*K.adK[0*3 + 1])*invdet;
}

inline double radCorrectionFactor(double dRadius_sq, const double * adRD, int nRDcoeffs)
{
    double dFactor = 1;
    if(nRDcoeffs>0)
    {
        double dRadiusPow = dRadius_sq;
        for(int i=1; ; i++)
        {
            dFactor += adRD[i-1]*dRadiusPow; //11-2-10 Positive numbers should now correct barrel distortion.
            if(i==nRDcoeffs) break;
            dRadiusPow *= dRadius_sq;
        }
    }
    return dFactor;
}

void CSimple2dPoint::correctRD(const CCamCalibMatrix & K)
{
    double dRadCorrectionFactor = radCorrectionFactor(sqr(x)+sqr(y), K.adRD, K.nRDcoeffs);
    x *= dRadCorrectionFactor;
    y *= dRadCorrectionFactor;
}

void CSimple2dPoint::calibrate(const CCamCalibMatrix & K, const bool CORRECT_RD)
{
    double xNew = K.K_inv.adK_inv[0]*x + K.K_inv.adK_inv[1]*y + K.K_inv.adK_inv[2];
    double yNew = K.K_inv.adK_inv[3]*x + K.K_inv.adK_inv[4]*y + K.K_inv.adK_inv[5];
    x=xNew; y=yNew;
    if(CORRECT_RD)
        correctRD(K);
}
void CSimple2dPoint::uncorrectRD(const CCamCalibMatrix & K)
{
    double dRadCorrectionFactor = radCorrectionFactor(sqr(x)+sqr(y), K.adRD, K.nRDcoeffs);
    x /= dRadCorrectionFactor;
    y /= dRadCorrectionFactor;
}

void CSimple2dPoint::uncalibrate(const CCamCalibMatrix & K)
{
    uncorrectRD(K);

    double xNew = K(0,0)*x + K(0,1)*y + K(0,2);
    double yNew = K(1,0)*x + K(1,1)*y + K(1,2);
    x=xNew; y=yNew;
}

CSimple2dPoint::CSimple2dPoint(const CLocation & l) : x(l.dx()), y(l.dy()) {};

void CBoWCorrespondences::calibrate(const CCamCalibMatrix & K, T2dPoints & aCalibratedPoints1, T2dPoints & aCalibratedPoints2, CInlierProbs & adArrLikelihood, CPointIdentifiers & pointIds, const bool CORRECT_RD) const
{
    int c = size();
    aCalibratedPoints1.clear(); aCalibratedPoints1.reserve(c);
    aCalibratedPoints2.clear(); aCalibratedPoints2.reserve(c);
    adArrLikelihood.clear(); adArrLikelihood.reserve(c);
    pointIds.clear(); pointIds.reserve(c);

    if(IS_DEBUG) CHECK((int)aCalibratedPoints1.size() != 0 || (int)aCalibratedPoints2.size() != 0 , "calibrateBoWCorr: Point vectors should be empty");

    CBoWCorrespondences::const_iterator ppEnd = end();
    int nPoint = 0;

    pointIds.clear();

    for (CBoWCorrespondences::const_iterator pCorresp = begin(); pCorresp < ppEnd; pCorresp++, nPoint++)
    {
        const CLocation loc1 = pCorresp->Location1();
        const CLocation loc2 = pCorresp->Location2();

        CSimple2dPoint p1(loc1);
        CSimple2dPoint p2(loc2);

        p1.calibrate(K, CORRECT_RD);
        p2.calibrate(K, CORRECT_RD);

        aCalibratedPoints1.push_back(p1);
        aCalibratedPoints2.push_back(p2);

        adArrLikelihood.push_back(pCorresp->priorProb());

        pointIds.push_back(CPointIds(loc1.id(), loc2.id()));
    }

    if(IS_DEBUG) CHECK(aCalibratedPoints1.size() != c || aCalibratedPoints2.size() != c || adArrLikelihood.size() != c || pointIds.size() != c , "calibrateBoWCorr: Point vectors should now be the same size");
}

