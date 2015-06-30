#include "vectorAngles.h"
//Now included at the top of every cpp file #include "stdafx.h"
#include <Eigen/Dense>
#include "lines.h"
#include <util/pp.h>


/*Reconstruct a point to a line:

Point x in view P, line x'+td in P'

PX=ax -> PX-ax = 0 -> P_1-3 X_1-3 - ax = -P_4 
P'X=bx'+cd' -> P'X-bx'+cd'=0
 * 
 * Or P'X.h=0 -> h . P'  X  =  0
 *              1x3 3x4 4x1 
 * 
 * 

6x6 matrix
(P_33  | -x  0  0 ) (X_1 X_2 X_3 a b c)^T = ( -P_col4 )
(P'_33 | 0  -x' -d)                         ( -P'_col4 )


Then where is P'X on the line (what is t)
*/
optional<C3dWorldPoint> reconstructPointToLine(const CWorldCamera & P_point, const C2dImagePointPx & x, const CWorldCamera & P_line, const C2dLine & line)
{
    const bool bVerbose = false;
    
    CHECK(P_point.getPMat() == P_line.getPMat(), "Reconstructing from a point and line in the same camera (no unique solution)");
    
    typedef Eigen::Matrix<double, 6, 6> TMatrix6d;
    typedef Eigen::Matrix<double, 6, 1> TVector6d;

    TMatrix6d K = TMatrix6d::Zero();
    K.block<3,3>(0,0) = P_point.M();
    K.block<3,3>(3,0) = P_line.M();
    K.block<3,1>(0,3) = -x.homog();
    K.block<3,1>(3,4) = -line.getPointOnLine().homog();
    K.block<2,1>(3,5) = -line.getDirection();
    //K.block<3,1>(3,5) = -C2dImagePoint(line.getDirection()).homog();
    
    COUT(K);
    
    TVector6d Pvec = TVector6d::Zero();
    Pvec.head<3>() = -P_point.translationToCam();
    Pvec.tail<3>() = -P_line.translationToCam();
    
    COUT(Pvec);
    
    const TVector6d X = K.inverse()*Pvec;
    
    COUT(X);
    
    const C3dWorldPoint soln = X.head<3>(); //homogeneous normalisation
    if(P_point.testInFront(soln) && P_line.testInFront(soln))
    {
        if(IS_DEBUG)
        {
            const C2dImagePointPx reprojPoint = P_point.projectToPx_checked(soln);
            if(!reprojPoint.isApprox(x, 0.001))
            {
                ALWAYS_VERBOSE;
                COUT(reprojPoint);
                COUT(x);
                COUT(soln);
                THROW("Reprojection to point failed");
            }
            CHECK(line.closestDistance(P_line.projectToPx_checked(soln)) > 0.0001, "Reprojection to line failed");
        }
        
        return soln;
    }
    else
        return optional<C3dWorldPoint>();
}

void shortenResidual(TEigen3dPoint & resid, const double dRadius)
{
    const double dResidLength = resid.norm();
    if(dResidLength>dRadius)
        resid *= (dResidLength-dRadius)/dResidLength;
    else
        resid.setZero();
}

C3dLine::C3dLine(const C3dBoundedLine & boundedLine)
{
    from2Points(boundedLine.getStartPoint(), boundedLine.getFinishPoint());
}

double C2dLine::signedAngle(const C2dLine & line) const
{
    return signedAngleBetweenUnitVectors(direction, line.direction);
}

double C2dLine::signedAngle(const C2dBoundedLine & line) const
{
    return signedAngle(line.unboundedLine());
}

TEigen2dPoint perpendicular(const TEigen2dPoint & dir)
{
    return TEigen2dPoint(-dir.y(),dir.x());
}


C3dBoundedLine makeRay(const C2dImagePointPx & p, const CWorldCamera & P)
{
    const C3dWorldPoint pointOnRay1 = P.pxToWorld_z(p, 0.3);
    const C3dWorldPoint pointOnRay2 = P.pxToWorld_z(p, 1.6); /* unbounded line todo check this depth range is ok*/

    return C3dBoundedLine(pointOnRay1, pointOnRay2);
}

/**
 * @brief Project ray from p out the the depth of the bounded line then find the closest point.
 * @param p
 * @param P
 * @param line
 * @return return closest point on \line rather than on the ray, so that the depth is sensible
 */
const C3dWorldPoint threeDPointNearLineAndMeasurement(const C2dImagePointPx & p, const CWorldCamera & P, const C3dBoundedLine & line)
{
    const bool bVerbose = false;

    /*if(!line.getStartPoint().isSensible(eSTAny) || !line.getFinishPoint().isSensible(eSTAny))
    {
        ALWAYS_VERBOSE;
        COUT(line);
        MAT2(componentAsLineOutsideSensibleRange, P.getFrameId());
        line.draw(componentAsLineOutsideSensibleRange, P, colour("3DCane"));
        SHOW(componentAsLineOutsideSensibleRange);
        THROW("Point on component-as-line outside sensible range");
    }*/
    
    const C3dBoundedLine ray = makeRay(p, P);
    const C3dWorldPoint closestPointToRay = line.closestPointToBoundedLine(ray);

    COUT(closestPointToRay);

    //Sometimes quite deep behind scene (could constrain it?)--should just reject closestPointToRay.checkSensible();

    return closestPointToRay;
}


/**
 * @brief Project ray from p out the the depth of the line then find the closest point.
 * @param p
 * @param P
 * @param line
 * @return return closest point on \line rather than on the ray, so that the depth is sensible
 */
const C3dWorldPoint threeDPointNearLineAndMeasurement(const C2dImagePointPx & p, const CWorldCamera & P, const C3dLine & line)
{
    const bool bVerbose = false;

    const C3dBoundedLine ray = makeRay(p, P);
    const C3dWorldPoint closestPointToRay = line.closestPointToBoundedLine(ray);

    COUT(closestPointToRay);
    return closestPointToRay;
}


/**
 * @brief Project ray from p out the the depth of the bounded line then find the shortest line between the closest point and the ray.
 * @param p
 * @param P
 * @param line
 * @return
 */
const C3dBoundedLine threeDLineBetweenLineAndMeasurement(const C2dImagePointPx & p, const CWorldCamera & P, const C3dBoundedLine & line)
{
    const C3dBoundedLine ray = makeRay(p, P);
    C3dWorldPoint closestPointToRay = line.closestPointToBoundedLine(ray);
    const C3dWorldPoint closestPointToLine = ray.closestPointToBoundedLine(line);

    if(closestPointToRay == closestPointToLine)
        closestPointToRay.x() += 0.00001; //Avoid exception here
        
    return C3dBoundedLine(closestPointToRay, closestPointToLine);
}

/**
 * @brief Ensure that new wire is broken into chunks (about) dSegmentLength long
 * @param line
 * @return
 */
void splitLine(const C3dBoundedLine & line, T3dPointVector & line_split, const double dSegmentLength)
{
    const bool bVerbose = true;

    const double dLength = line.length();
    int nNumSegments = std::max<int>(cvRound(dLength/dSegmentLength), 1);

    for(int nSegment = 0; nSegment <= nNumSegments; nSegment++) {
        double dRelPos = (double)nSegment/(double)nNumSegments;
        const C3dWorldPoint pointOnLine = (1-dRelPos)*line.getStartPoint() + dRelPos*line.getFinishPoint();
        line_split.push_back(pointOnLine);
        
        COUT(nSegment);
        COUT(pointOnLine);
    }
    
    CHECK(!line.getFinishPoint().isApprox(line_split.back()), "splitLine failed");
    CHECK(!line.getStartPoint().isApprox(line_split.front()), "splitLine failed");
}

bool onLine(const double t)
{
    return t>=0 && t <= 1;
}

bool withinLine(const double t) //not endpoints
{
    return t>0 && t < 1;
}


/* This isn't perfect (e.g. 2 non-intersecting lines in a T) but should be good enough in practice */
bool incompatLines(const C2dBoundedLine & line1, const C2dBoundedLine & line2, const eEndpoints nStartIdx)
{
    C2dLine perpLineAtEndpoint;
    perpLineAtEndpoint.fromPointAndDir(line2.getEndPoint(nStartIdx), line2.perpendicular());
    return (bool)line1.intersection(perpLineAtEndpoint);
}

//incompat if approx pll and one is
bool containedBetween(const C2dBoundedLine & line1, const C2dBoundedLine & line2)
{
    return (withinLine(line1.get_t(line2.getStartPoint())) && withinLine(line1.get_t(line2.getFinishPoint())));
}

bool incompatLines(const C2dBoundedLine & line1, const C2dBoundedLine & line2)
{
    const double dPxThresh = 1; //Allow a tiny bit of overlap at the end
    const optional<T2dBoundedLine_base> pLine1_trimmed = line1.trim(dPxThresh);
    const optional<T2dBoundedLine_base> pLine2_trimmed = line2.trim(dPxThresh);
    
    if(!pLine1_trimmed || !pLine2_trimmed)
        return false;
        
    const C2dBoundedLine line1_trimmed(*pLine1_trimmed);
    const C2dBoundedLine line2_trimmed(*pLine2_trimmed);
    const bool bIncompat = incompatLines(line1_trimmed, line2_trimmed, eStart) || incompatLines(line1_trimmed, line2_trimmed, eFinish)
                           || incompatLines(line2_trimmed, line1_trimmed, eStart) || incompatLines(line2_trimmed, line1_trimmed, eFinish)
                           || line1_trimmed.intersection(line2_trimmed)
                           || containedBetween(line1_trimmed, line2_trimmed)
                           || containedBetween(line2_trimmed, line1_trimmed);
    
    /*const bool bVerbose = false;
    if(bVerbose) {
        COUT(bIncompat);
        const cv::Scalar col = !bIncompat ? colour("Assigned") :  colour("Unassigned");
        MAT(green_compatible_red_incompat);
        line1.draw(green_compatible_red_incompat, 1, col);
        line2.draw(green_compatible_red_incompat, 1, col);
        SHOW(green_compatible_red_incompat);
    }*/

    return bIncompat;
}

double truncateCosAngle(double dCosAngle)
{
    const double dAbsCosAngle = fabs(dCosAngle);
    if(dAbsCosAngle>1) {
        CHECK(dAbsCosAngle>1.00001, "Error computing angle between vectors--probably not normalised");
        dCosAngle=sign(dCosAngle);
    }
    return dCosAngle;
}

double directedInvCosAngle(double dCosAngle, const bool bDirected)
{
    if(!bDirected)
        dCosAngle = fabs(dCosAngle);

    const double dAngle = acos(truncateCosAngle(dCosAngle));

    CHECKNAN(dAngle);

    return dAngle;
}

template<class TPoint>
double angleBetweenUnitVectors(const TPoint & dir1, const TPoint & dir2, const bool bDirected /*to return angles in range 0...Pi vs 0...Pi/2*/)
{
    if(IS_DEBUG) CHECK(!zero(dir1.squaredNorm()-1), "dir1 not normalised");
    if(IS_DEBUG) CHECK(!zero(dir2.squaredNorm()-1), "dir2 not normalised");

    const double dCosAngle = dir1.dot(dir2);
    return directedInvCosAngle(dCosAngle, bDirected);
}

double safe_acos(const double dCosAngle)
{
    return directedInvCosAngle(dCosAngle, true);
}

template<class TPoint>
double angleBetween3Points(const TPoint & p1, const TPoint & p2, const TPoint & p3)
{
    const bool bVerbose = false;
    COUT(p1);
    COUT(p2);
    COUT(p3);
    
    const TPoint p12 = p2 - p1;
    const TPoint p23 = p3 - p2;

    return angleBetweenVectors(p12, p23, true);
}

double signedAngleBetweenUnitVectors(const TEigen2dPoint & direction1, const TEigen2dPoint & direction2)
{
    if(IS_DEBUG) CHECK(!zero(direction1.squaredNorm()-1), "dir1 not normalised");
    if(IS_DEBUG) CHECK(!zero(direction2.squaredNorm()-1), "dir2 not normalised");

    const double dAngle = angleBetweenUnitVectors(direction1, direction2, true);
    if(direction1.x() * direction2.y() - direction1.y() * direction2.x() < 0)
        return -dAngle;
    else
        return dAngle;
}

double signedAngleBetweenVectors(const TEigen2dPoint & direction1, const TEigen2dPoint & direction2)
{
    return signedAngleBetweenUnitVectors(direction1.normalized(), direction2.normalized());
}

template<class TPoint>
double angleBetweenVectors(const TPoint & direction1, const TPoint & direction2, const bool bDirected /*to return angles in range 0...Pi vs 0...Pi/2*/ )
{
    const bool bVerbose = false;
    COUT(direction1);
    COUT(direction2);

    const double dDenom = sqrt(direction1.squaredNorm() * direction2.squaredNorm());
    COUT(dDenom);

    //if(IS_DEBUG) CHECK(dDenom == 0, "A direction vector has length 0--can't compute angle--this happens when polyline points are equal");
    if(dDenom == 0)
    {
        REPEAT(100, cout << "A direction vector has length 0--can't compute angle" << endl);
        return M_PI;
    }
    
    COUT(direction1.dot(direction2));
    const double dCosAngle = direction1.dot(direction2) / dDenom;
    
    COUT(dCosAngle);
    
    return directedInvCosAngle(dCosAngle, bDirected);
}

template double angleBetween3Points<C3dWorldPoint>(const C3dWorldPoint & p1, const C3dWorldPoint & p2, const C3dWorldPoint & p3);
template double angleBetween3Points<TEigen3dPoint>(const TEigen3dPoint & p1, const TEigen3dPoint & p2, const TEigen3dPoint & p3);
template double angleBetween3Points<TEigen2dPoint>(const TEigen2dPoint & p1, const TEigen2dPoint & p2, const TEigen2dPoint & p3);
template double angleBetween3Points<C2dImagePointPx>(const C2dImagePointPx & p1, const C2dImagePointPx & p2, const C2dImagePointPx & p3);

template double angleBetweenVectors<TEigen2dPoint>(const TEigen2dPoint & direction1, const TEigen2dPoint & direction2, const bool bDirected /*to return angles in range 0...Pi vs 0...Pi/2*/ );
template double angleBetweenVectors<TEigen3dPoint>(const TEigen3dPoint & direction1, const TEigen3dPoint & direction2, const bool bDirected /*to return angles in range 0...Pi vs 0...Pi/2*/ );
template double angleBetweenVectors<C3dWorldPoint>(const C3dWorldPoint & direction1, const C3dWorldPoint & direction2, const bool bDirected /*to return angles in range 0...Pi vs 0...Pi/2*/ );
template double angleBetweenVectors<C2dImagePointPx>(const C2dImagePointPx & direction1, const C2dImagePointPx & direction2, const bool bDirected /*to return angles in range 0...Pi vs 0...Pi/2*/ );

template double angleBetweenUnitVectors<TEigen2dPoint>(const TEigen2dPoint & dir1, const TEigen2dPoint & dir2, const bool bDirected /*to return angles in range 0...Pi vs 0...Pi/2*/);
template double angleBetweenUnitVectors<TEigen3dPoint>(const TEigen3dPoint & dir1, const TEigen3dPoint & dir2, const bool bDirected /*to return angles in range 0...Pi vs 0...Pi/2*/);
