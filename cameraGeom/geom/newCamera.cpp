#include "newCamera.h"
#include <Eigen/Dense> // only needed here for the determinant
#include <util/pp.h> 

//TODO: move in fns from threeDpoint.cpp which should be here (mostly in cameras.ipp now)

CNormalisedCameraMat::CNormalisedCameraMat(const TEigen3dRotation & R, const TEigen3dPoint & T, const eFrameForRT eFrame)
{
    if(eFrame == eRTFromCalibration)
    {
        P.block<3,3>(0,0) = R;
        P.col(3) = T;
    }
    else
    {
        P.block<3,3>(0,0) = R.transpose();
        P.col(3) = -R.transpose()*T;
        
        if(IS_DEBUG)
        {
            const double dError = (T - camPositionRelToFrame()).squaredNorm();
            
            if(dError > 0.000001)
            {
                REPEAT(3,
                ALWAYS_VERBOSE;
                COUT(R);
                COUT(T.transpose());
                COUT(dError);
                COUT(camPositionRelToFrame().transpose());
                COUT(P);
                COUT(M());
                
                //THROW("Coordinate frame error");
                cout << "Warning: Coordinate frame error (probably a rotation error somewhere -- TODO)" << endl);
            }
        }
    }

    const bool bVerbose = false;
    COUT(R);
    COUT(T.transpose());
    COUT(P);
    
    checkOk();
}

CNormalisedCameraMat::CNormalisedCameraMat(const C3dRotation & R, const TEigen3dPoint & T, const eFrameForRT eFrame)
{
    TEigen3dRotation Rmat;
    R.asMat(Rmat);
    *this = CNormalisedCameraMat(Rmat, T, eFrame);
    checkOk();
}

void CCameraMat_base::checkForNan() const
{
    CHECK_P(std::isnan(P.sum()), P, "Camera is nan");
    CHECK_P((P.norm() < 1  || P.norm() > (includesK() ? 1000000 : 1000)), P, "Camera is probably uninitialised");
}

bool CCameraMat_base::includesK() const
{
    return !within(M().determinant(), 1, 0.0001);
}

/*

template<>
void doCout(const char * label, const CCameraMat_base & t);
{
    cout << label << " = \n" << t << endl;
}

template<>
void doCout(const char * label, const C2dImagePoint & t)
{
    cout << label << " = " << t.transpose() << endl;
}

template<>
void doCout(const char * label, const C2dImagePointPx & t)
{
    cout << label << " = " << t.transpose() << endl;
}

template<>
void doCout(const char * label, const C3dWorldPoint & t)
{
    cout << label << " = " << t.transpose() << endl;
}
template<>
void doCout(const char * label, const TEigen3dRotation & t)
{
    cout << label << " = \n" << t << endl;
}*/
