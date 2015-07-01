#ifndef NEWCAMERA_H
#define    NEWCAMERA_H

/*
 * Move to Eigen matrices for camera stuff. Use double for float type, as projection may be a significant (and easily vectorisable) cost.
 *
 * Created on 31 October 2011, 3:33 PM
 */

#include <Eigen/Core>
#include <geom/geom.h>
#include <opencv2/core/core.hpp> //TODO: should be able to forward declare stuff and get rid of these opencv headers
#include <opencv2/core/types_c.h>
#include <opencv2/imgproc/imgproc.hpp>
#include <util/calibration.h>

#include <boost/optional.hpp>
using boost::optional;

#include "vectorPoints.h"
#include "pointWithThickness.h"

class CCameraPose;


/**
 * @class CCameraMat_base
 * @brief Base for a world-to-image camera or a world-to-pixel camera.
 */
class CCameraMat_base
{
    friend std::ostream& operator<<(std::ostream& s, const CCameraMat_base & P);
protected:
    TCamMat P; //Camera matrix. [May be calibrated or uncalibrated]

    /* Must project in front of camera */
    const optional<const C2dImagePointPx> project_int(const C3dWorldPoint & X) const HOT;

    void checkForNan() const;

    bool includesK() const;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    CCameraMat_base();
    /**
     * @brief Left 3x3 block of P (only a rotation if !includesK())
     * @return
     */
    const TEigen3dRotation M() const HOT;

    /* Right 3x1 block of P (not actually a translation if includesK()) */
    TTransformMat getTransformMat() const;

    bool isInit() const {
        return P.squaredNorm() != 0;
    }
    void checkInit() const {
        if(IS_DEBUG) CHECK(!isInit(), "Camera matrix not initialised. Absolute path to folder containing camera calibration MUST be included in image folder's cameraData.cfg file");
    }
    
    /* UNSAFE projection--no checks (assume checks happen after optimisation/whatever) */
    const C2dImagePointPx fastProject(const C3dWorldPoint & X) const HOT;

    const C3dWorldPoint translationToCam() const;

    /**
     * @brief Returns depth, if point is in front, or empty (w.r.t. individual cameras)
     * @return
     */
    const optional<const double> depth(const C3dWorldPoint & X) const;
    const double depth_checked(const C3dWorldPoint & p3d) const;

    /**
     * @brief Returns depth in the z direction (we should only create new structure at the correct distance in front of the camera, e.g. wires should fit within the doors)
     * @return
     */
    //const optional<const double> depth_z(const C3dWorldPoint & X) const;

    /**
     * @brief could be faster
     * @param p3d
     * @return
     */
    const bool testInFront(const C3dWorldPoint & X) const;

    const TCamMat & getPMat() const {
        return P;
    }
};

std::ostream& operator<<(std::ostream& s, const CCameraMat_base & P);

class CNormalisedCameraMat : public CCameraMat_base
{
public:
    inline const optional<const C2dImagePointPx> project(const C3dWorldPoint & X) const {
        return project_int(X);
    }
    const TEigen3dRotation rotationMat() const {
        checkOk();
        return M();
    }

    const C3dWorldPoint inverseTransform(const C3dWorldPoint & camPoint) const;
    const C3dWorldPoint transform(const C3dWorldPoint & camPoint) const;

    CNormalisedCameraMat();

    enum eFrameForRT { eRTFromCalibration /* make this R,T directly into a camera */, eRTFromGlobalPosition /*flip*/};
    /**
     * @brief Construct camera from its position and orientation.
     *
     * Camera is at T pointing in direction R
     * World point X is viewed in local camera frame to be at Xc.
     *
     * P X -> Xc
     * R^T (X - T) = Xc
     * or X = R Xc + T
     * (not entirely sure of this -- TODO: test)
     * REVERSED IN LOCALPOSITION
     * @param R
     * @param T
     * @return
     */
    CNormalisedCameraMat(const TEigen3dRotation & R, const TEigen3dPoint & T, const eFrameForRT eFrame);

    CNormalisedCameraMat(const C3dRotation & R, const TEigen3dPoint & T, const eFrameForRT eFrame);

    /* REVERSE of c'tor with eRTFromGlobalPosition */const C3dWorldPoint camPositionRelToFrame() const {
        return (-rotationMat().transpose()*translationToCam()).eval();
    }

    void checkOk() const {
        if(IS_DEBUG) 
        {
            checkInit();
            checkForNan();
            CHECK(includesK(), "M is not a rotation");
        }
    }
};

/* 3D coordinates to pixel coordinates */
class CPixelCameraMat : public CCameraMat_base
{
public:
    /**
     * @brief project to px; throws an exception if projection fails (behind or far-from camera)
     * @param X
     * @return
     */
    //C2dImagePointPx projectToPx_checked(const C3dWorldPoint & X) const HOT; 

    const optional<const C2dImagePointPx> projectToPx_int(const C3dWorldPoint & X) const HOT;

    CPixelCameraMat(const TEigen3dRotation & K, const CNormalisedCameraMat & P_in);

    CPixelCameraMat(const CPixelCameraMat & localToPixel, const TTransformMat & worldToLocal);

    /*Converts this camera to one which can project to enlarged images (for drawing structure which is out of view)*/
    void setPxOffset(const int nRows, const int nCols);

    CPixelCameraMat();

    void checkOk() const {
        if(IS_DEBUG) 
        {
            checkForNan();
            CHECK(!includesK(), "M is a pure rotation");
        }
    }

    /**
     * @brief Returns point with depth dDepth wrt this camera.
     * @param px
     * @param dDepth
     * @return
     */
    C3dWorldPoint pxToWorld_depth(const C2dImagePointPx & px, const double dDepth) const;
    
    /**
     * @brief Returns point with z-coordinate dz
     * @param px
     * @param dDepth
     * @return
     */
    C3dWorldPoint pxToWorld_z(const C2dImagePointPx & px, const double dDepth) const;

    /**
     * @brief converts a width measurement (metres perp to cam) at a particular location to pixels
     * @param dWidthInMetres
     * @param location
     * @return
     */
    const optional<const double> widthToPx(const double dWidthInMetres, const C3dWorldPoint & location) const;
    const double widthToPx_checked(const double dWidthInMetres, const C3dWorldPoint & location) const;

    const double focalLength() const;
    /**
     * @brief inverse of measurementToPx
     * @param location
     * @param dSizePx
     * @return
     */
    optional<const double> pxToWidth(const C3dWorldPoint & location, const double dSizePx) const;

private:
    optional<const double> pxToWidth_int(const C3dWorldPoint & location, const double dSizePx) const;

public:

};

inline CNormalisedCameraMat operator|(const C3dRotation & q, const C3dWorldPoint & T);

class CNewCamCalibMatrix
{
    Eigen::Matrix3d K, K_inv;
    Eigen::VectorXd aRDCoeffs;
public:
    CNewCamCalibMatrix();
    CNewCamCalibMatrix(const std::string strCalibFileName);
    CNewCamCalibMatrix(const cv::Size & imSize, const double dFocalLength);

    C2dImagePointPx toPxCoords(const C2dImagePoint & imPoint) const;
    C2dImagePoint toImCoords(const C2dImagePointPx & imPointPx) const;

    const double focalLength() const {
        if(IS_DEBUG) CHECK(K(0,0) != K(1,1), "Camera has non-square pixels so need to be more careful using focal length");
        return K(0,0);
    }

    const Eigen::Matrix3d & mat() const {
        return K;
    }
    const Eigen::Matrix3d & mat_inv() const {
        return K_inv;
    }

    //void setPxOffset(const int nRows, const int nCols);
};

/**
 * @class CWorldCamera
 * @author Tom Botterill
 * @date 09/08/12
 * @file newCamera.h
 * @brief Camera projecting points in world to points in pixel coordinates. TODO: RD correction (not necessary with the nice lenses on the grasshopper)
 */
class CWorldCamera : public CPixelCameraMat
{
    bool bInit;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    CWorldCamera(const TTransformMat & P_worldToLocal, const CPixelCameraMat & P_localToPx);

    void checkInit() const {
        if(IS_DEBUG) CHECK(!bInit, "World camera not yet initialised");
        if(IS_DEBUG) CHECK(P.squaredNorm() == 0, "World camera is all 0");
        checkOk();
    }

    CWorldCamera() : bInit(false) { }

    //const CCameraPose & getCamPose() const;

    void setCamMatForAnimation(const TCamMat & Pnew);
    //double depth_z(const C3dWorldPoint & X) const;
    
    //Y coordinate relative to this camera frame. Used for checking wires etc aren't below the ground or above the door entrance
   //double relativeY(const C3dWorldPoint & X) const;
    
    const optional<const C2dImagePointPx> projectToPx(const C3dWorldPoint & X) const;
    const optional<const C3dPolylineControlPoint::T2dPoint> projectToPx(const C3dPolylineControlPoint & X) const;
    const optional<const C3dPolylineControlPointWithThickness::T2dPoint> projectToPx(const C3dPolylineControlPointWithThickness & X) const;
    
    /**
     * @brief 
     * @param X
     * @param checkMode eNoChecks (only throws if not in front etc) or eCheckIO (prints an informative error message and throws)
     * @param structureType
     * @return 
     */
    template<class T3dPoint>
    const typename T3dPoint::T2dPoint projectToPx_checked(const T3dPoint & X) const
    {
        const optional<const typename T3dPoint::T2dPoint> px = projectToPx(X);
        CHECK_P(!px, X, "Probably point is behind camera on projectToPx_checked");
        return *px;
    }
};

class COptimisableSceneComponent;

typedef std::vector<CWorldCamera, Eigen::aligned_allocator<CWorldCamera> > TCams;

/*class CCameraData
{
    boost::scoped_ptr<CCam> pCamera;
    boost::scoped_ptr<CBackgroundModel> pBackground;
    cv::Size imSize;
    const int id;

    void initialiseSize(const std::string & pathToImageDir);

public:
    std::string getIdStr() const;
    const int getCamId() const {
        return id;
    }

    CCameraData(const CImagingSystemData * pImagingSystemBeingCreated, const int id, const int eGlobalCoordSysCam);

    //const cv::Size & getImSize() const { if(IS_DEBUG) CHECK(imSize.area() == 0, "Uninit size"); return imSize; }
    //void setImSize(const cv::Mat & sizeSource) const { imSize = cv::Size(sizeSource.rows, sizeSource.cols); }

    const CCam & getCalibration() const;
    const CBackgroundModel & getBackground() const {
        if(IS_DEBUG) CHECK(!pBackground, "Background not initialised yet");
        return *pBackground;
    }

    const cv::Size & getImSize() const {
        if(IS_DEBUG) CHECK(imSize.area()<=0, "imSize not initialised");
        return imSize;
    }

    //TODO: Integrate with background/mask
    const cv::Rect getActiveArea() const {
        const int nMargin = 50;
        cv::Rect activeArea(nMargin, nMargin, imSize.width-nMargin, imSize.height-nMargin);
        return activeArea;
    }
};*/

/*bool widthIsSensible_metres(const double dMetres);
bool CCheckThickness::isSensiblePx(const double dSizePx);
void CCheckThickness::checkSensiblePx(const double dSizePx);*/

#endif    /* NEWCAMERA_H */
