#ifndef CAMERA_IPP
#define CAMERA_IPP
#include "newCamera.h"
#include <util/pp.h>

#if defined(__GNUC__) || defined(INSTANTIATE_CAMERAS_IPP)

#ifndef IPP_INLINE
#  ifndef INSTANTIATE_CAMERAS_IPP
#    define IPP_INLINE inline 
#  else //Force an instantiation of everything, but allow multiple definitions where the ipp is included (Windows isn't happy with making everything inline).
#    define IPP_INLINE
#  endif
#endif

#include <Eigen/Dense>

//////////// CCameraMat_base ////////////////////////

IPP_INLINE CCameraMat_base::CCameraMat_base()
{
    P.setZero();
}

IPP_INLINE const TEigen3dRotation CCameraMat_base::M() const
{
    /*checkForNan(); a bit slow even in debug*/ return P.block<3,3>(0,0);
}

IPP_INLINE const bool CCameraMat_base::testInFront(const C3dWorldPoint & X) const
{
    const optional<const double> pDepth = depth(X);
    return pDepth && *pDepth > 0.01; //1cm at least
}

IPP_INLINE TTransformMat CCameraMat_base::getTransformMat() const
{
    TTransformMat T = TTransformMat::Identity();
    T.block<3,4>(0,0) = getPMat();
    return T;
}

IPP_INLINE const C3dWorldPoint CCameraMat_base::translationToCam() const
{
    checkInit();
    checkForNan();
    return P.block<3,1>(0,3);
}

IPP_INLINE const C2dImagePointPx CCameraMat_base::fastProject(const C3dWorldPoint & X) const
{
    const TEigen2dPointHomog twoDHomog = P*X.homog();
    return twoDHomog.head(2)/twoDHomog(2);
}

IPP_INLINE const optional<const C2dImagePointPx> CCameraMat_base::project_int(const C3dWorldPoint & X) const
{
    const bool bVerbose = false;
    
    if(IS_DEBUG) checkForNan();

    if(!testInFront(X))
    {
        COUT2("X not in front of this camera:", X);
        return optional<const C2dImagePointPx>();
    }
    
    //TEigen2dPointHomog twoDHomog = P*X.homog();
    //C2dImagePointPx twoD = twoDHomog.head(2)/twoDHomog(2);
    return fastProject(X);
}

///////// CNormalisedCameraMat ////////////////////////////

IPP_INLINE CNormalisedCameraMat::CNormalisedCameraMat()
{
    P.setIdentity();
}

IPP_INLINE CNormalisedCameraMat operator|(const C3dRotation & q, const C3dWorldPoint & T)
{
    TEigen3dRotation R;
    q.asMat(R);
    return CNormalisedCameraMat(R, T, CNormalisedCameraMat::eRTFromCalibration);
}

IPP_INLINE const C3dWorldPoint CNormalisedCameraMat::inverseTransform(const C3dWorldPoint & camPoint) const
{
    checkInit();
    const C3dWorldPoint localPoint = (rotationMat().transpose() * (camPoint - translationToCam())).eval();
    if(IS_DEBUG)
    {
        const C3dWorldPoint camPoint_again = transform(localPoint);
        if(IS_DEBUG) CHECK(!zero(1000*(camPoint_again-camPoint).squaredNorm()), "inverse CNormalisedCameraMat::transform failed");
    }    
    return localPoint;
}
IPP_INLINE const C3dWorldPoint CNormalisedCameraMat::transform(const C3dWorldPoint & camPoint) const
{
    checkInit();
    C3dWorldPoint transformedPoint = (P*camPoint.homog()).eval();
    return transformedPoint;        
}

/////////// CPixelCameraMat /////////////////////
/*IPP_INLINE C2dImagePointPx CPixelCameraMat::projectToPx_checked(const C3dWorldPoint & X) const
{
    const optional<const C2dImagePointPx> px = projectToPx(X);
    CHECKOPTIONAL(px);
    return *px;
}*/

IPP_INLINE const optional<const C2dImagePointPx> CPixelCameraMat::projectToPx_int(const C3dWorldPoint & X) const
{
    return project_int(X);
}

IPP_INLINE CPixelCameraMat::CPixelCameraMat()
{
    P.setZero();
}

IPP_INLINE C3dWorldPoint CPixelCameraMat::pxToWorld_z(const C2dImagePointPx & px, const double dZ) const
{
    checkOk();
    
    /* Solve P*X=px,(P*X).z() = dZ
     * 
     * Reconstruct X1, X2, at different depths, interpolate to get X
     * */    
    const TEigen3dPoint X1 = pxToWorld_depth(px, 0.5);
    const TEigen3dPoint X2 = pxToWorld_depth(px, 0.6);
     /*(X1 + lambda*(X2-X1)).z() = dZ
      * lambda = (dZ-X1.z())/((X2-X1).z())
      */
    const double lambda = (dZ-X1.z())/((X2-X1).z());
    
    const TEigen3dPoint X = X1+lambda*(X2-X1);
    if(IS_DEBUG) CHECK(!zero(X.z()-dZ), "pxToWorld_z failed");
    return X;
}

IPP_INLINE C3dWorldPoint CPixelCameraMat::pxToWorld_depth(const C2dImagePointPx & px, const double dDepth) const
{
    checkOk();
    
    //if(IS_DEBUG) CCheck3D::checkDepthIsWithinRange(dDepth, eSTAny);// CHECK(!depthIsReasonable(dDepth), "dDepth outside reasonable range");
    
    TEigen3dPoint X_cam = dDepth * px.homog();
    TEigen3dPoint X_shifted = X_cam - translationToCam();
    C3dWorldPoint X_world = (M().inverse()*X_shifted).eval();
    
    if(IS_DEBUG) CHECK(!within(*(depth(X_world)), dDepth, 0.0001), "Reconstruct point failed");
    
    if(IS_DEBUG)
    {
        C2dImagePointPx worldBackToPx = *projectToPx_int(X_world);
        if(IS_DEBUG) CHECK(!zero((worldBackToPx - px).squaredNorm()), "Project back into cam failed");
    }
    
    return X_world;
}

IPP_INLINE CPixelCameraMat::CPixelCameraMat(const TEigen3dRotation & K, const CNormalisedCameraMat & P_in)
{
    P = K * P_in.getPMat();
    checkOk();
}

IPP_INLINE CPixelCameraMat::CPixelCameraMat(const CPixelCameraMat & localToPixel, const TTransformMat & worldToLocal)
{
    P = localToPixel.getPMat()*worldToLocal;
    
    const bool bVerbose = false;
    if (bVerbose)
        cout << "Making camera from:" << endl;
    COUT(worldToLocal);
    //COUT(P_worldToLocal_ext);
    COUT(localToPixel);
    COUT(P);
    
    checkOk();
}


//inverse of measurementToPx
IPP_INLINE optional<const double> CPixelCameraMat::pxToWidth(const C3dWorldPoint & location, const double dSizePx) const {
    checkInit();
    
    const optional<const double> pdDiam = pxToWidth_int(location, dSizePx);
    
    if(!pdDiam) 
        return pdDiam;

    if(IS_DEBUG) CHECK(!within(dSizePx, widthToPx_checked(*pdDiam, location), 0.0001), "pxToWidth Inverse check failed" );
    
    if(*pdDiam <= 0)
        return optional<const double>();
    
    return pdDiam;
}

IPP_INLINE optional<const double> CPixelCameraMat::pxToWidth_int(const C3dWorldPoint & location, const double dSizePx) const {
    const optional<const double> pdDepth = depth(location);
    if(!pdDepth)
        THROW("pxToWidth_int requires a 3D point in front of the camera (todo: optional version?)");

    if(dSizePx<=0)
        return optional<const double>();

    const double dPxPerRad = focalLength();        
    return dSizePx * *pdDepth / dPxPerRad;
}

IPP_INLINE const optional<const double> CCameraMat_base::depth(const C3dWorldPoint & p3d) const
{
    checkInit();
    
    if(IS_DEBUG && M().determinant() <= 0)
    {
        ALWAYS_VERBOSE;
        COUT(*this);
        THROW("Expect a +ve determinant here--Probably using a disabled camera");
    }
    
    double dDepthNumerator = P.row(2) * p3d.homog();
    
    if(dDepthNumerator<=0)
        return optional<const double>();
    
    const double dDenom = M().row(2).norm();
    //COUT(dDenom);
    
    return optional<const double>(dDepthNumerator/dDenom);
}

IPP_INLINE const double CCameraMat_base::depth_checked(const C3dWorldPoint & p3d) const
{
    optional<const double> pDepth = depth(p3d);
    CHECKNOTNULL(pDepth);
    return *pDepth;
}

IPP_INLINE const double CPixelCameraMat::focalLength() const
{
    checkInit();

    const double dFocalLength = sqrt(M().determinant()); 
    
    const bool bVerbose = false;
    COUT(dFocalLength);
    
    if(IS_DEBUG) CHECK(dFocalLength < 100 || dFocalLength > 3000, "dFocalLength outside reasonable range");
    
    return dFocalLength;
}

IPP_INLINE const double CPixelCameraMat::widthToPx_checked(const double dWidthInMetres, const C3dWorldPoint & location) const {
    CHECK_P(dWidthInMetres<=0, dWidthInMetres, "Width 0 or negative");
    checkInit();
    optional<const double> pWidth = widthToPx(dWidthInMetres, location);
    CHECKNOTNULL(pWidth);
    return *pWidth;    
}

//The call back to pxToWidth_int checks all the widths are good again
IPP_INLINE const optional<const double> CPixelCameraMat::widthToPx(const double dWidthInMetres, const C3dWorldPoint & location) const {

    checkInit();
    
    CHECK(dWidthInMetres <= 0, "Negative width");
    
    const optional<const double> pdDepth = depth(location);
    if(!pdDepth || *pdDepth<=0)
        return optional<const double>();
        
    const double dDepth = *pdDepth;
    
    const double dWidthInCalibratedImage = dWidthInMetres / dDepth;
    const double dFocalLength = focalLength();
    double dWidthInPx = dFocalLength * dWidthInCalibratedImage;
    /*
    if(checkMode == eCheckIO)
    {
        CCheckThickness::checkSensiblePx(dWidthInPx, structureType);
    }
    else  if(checkMode == eValidateIO)
    {
        if(!CCheckThickness::isSensiblePx(dWidthInPx, structureType))
            return optional<const double>();
    }*/
    
    const bool bVerbose = false;
    if(bVerbose) cout << "Pixel width " << dWidthInPx << " computed from " << dWidthInMetres * 1000 << "mm at depth " << dDepth << "m" << endl;
    if(IS_DEBUG)
    {
        optional<const double> pdWidth = pxToWidth_int(location, dWidthInPx);
        CHECK(!pdWidth, "Could not reproject width--probably the input width or point are just OOB");
        CHECK(!within(dWidthInMetres, *pdWidth, 0.0001), "pxToWidth Inverse check failed" );
    }

    return optional<const double>(dWidthInPx);
}


///////////////// C3dWorldPoint ////////////////

IPP_INLINE C3dWorldPoint::C3dWorldPoint(const double x, const double y, const double z) : TEigen3dPoint(x,y,z) 
{
}

///////// CWorldCamera ///////////////////////


IPP_INLINE CWorldCamera::CWorldCamera(const TTransformMat & P_worldToLocal, const CPixelCameraMat & P_localToPx) 
 : CPixelCameraMat(P_localToPx, P_worldToLocal), bInit(true)
{
}

IPP_INLINE const optional<const C2dImagePointPx> CWorldCamera::projectToPx(const C3dWorldPoint & X) const
{
    return CPixelCameraMat::projectToPx_int(X);
}

IPP_INLINE const optional<const C3dPolylineControlPoint::T2dPoint> CWorldCamera::projectToPx(const C3dPolylineControlPoint & X) const
{
    const optional<const C2dImagePointPx> px = projectToPx(X.getPoint());
    
    if(px)
        return C3dPolylineControlPoint::T2dPoint(*px);
    else
        return optional<const C3dPolylineControlPoint::T2dPoint>();
}

IPP_INLINE const optional<const C3dPolylineControlPointWithThickness::T2dPoint> CWorldCamera::projectToPx(const C3dPolylineControlPointWithThickness & X) const
{
    const optional<const C2dImagePointPx> px = projectToPx(X.getPoint());
    
    if(px)
    {
        const optional<const double> pdThicknessPx = widthToPx(X.getWidth(), X.getPoint());
        if(pdThicknessPx)
            return C3dPolylineControlPointWithThickness::T2dPoint(*px, *pdThicknessPx);
    }

    return optional<const C3dPolylineControlPointWithThickness::T2dPoint>();
}


#endif // #if defined(__GNUC__) || defined(INSTANTIATE_CAMERAS_IPP)
#endif // CAMERA_IPP
