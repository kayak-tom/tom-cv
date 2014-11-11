#pragma once
#ifndef _CALIBRATION_H
#define _CALIBRATION_H

#include "Simple2dPoint.h"

class CCamCalibMatrix
{
    friend class CSimple2dPoint;
    friend class CCamera;
    friend class CImParams;
    friend std::ostream& operator<<(std::ostream& s, const CCamCalibMatrix & X);

protected:
    double adK[9];

    double adRD[2];
    int nRDcoeffs;

    class CInvCamCalibMatrix
    {
    protected:

    public:
        double adK_inv[9];
        void init(const CCamCalibMatrix & K);
    };

    CInvCamCalibMatrix K_inv;
public:
    void init(const char * szCalibFilename, const int SCALE_DOWN = 1, const int CROP = 0);
    void init_inv();
    
    void setId();

    double focalLength() const { return adK[0]; }

    void testK(const int nWidth, const int nHeight) const;

    bool canCorrectRD() const { return nRDcoeffs>0 && (nRDcoeffs>1 || adRD[0] != 0); }

    void setRDCoeffs(double * adNewRD, int nNewRDCoeffs) 
    {
        nRDcoeffs = nNewRDCoeffs;
        for(int i=0;i<nRDcoeffs; i++)
            adRD[i] = adNewRD[i];
    }
    double operator()(int r, int c) const
    {
        return adK[3*r + c];
    }
    
    double & operator()(int r, int c)
    {
        return adK[3*r + c];
    }
    
    double inv(int r, int c) const
    {
        return K_inv.adK_inv[3*r + c];
    }
};

#endif