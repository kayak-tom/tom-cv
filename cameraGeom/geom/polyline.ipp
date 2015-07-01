#ifndef POLYLINE_IPP
#define POLYLINE_IPP

#include "polyline.h"
#include "lines.ipp"
#include <boost/foreach.hpp>
#include <Eigen/Geometry>
#include "vectorAngles.h"
#include <util/pp.h>

inline double penalisedTanAngle(const double dDot, const double dCross)
{
    if(dDot <= 0)
        return sign(dCross) * CLMFunction::HUGE_RESIDUAL();

    return dCross/dDot;
}

template<class TVecType>
inline double crossProd(const TVecType & v1, const TVecType & v2);

template<>
inline double crossProd<C3dWorldPoint>(const C3dWorldPoint & v1, const C3dWorldPoint & v2)
{
    return v1.cross(v2).norm();
}

template<>
inline double crossProd<C2dImagePointPx>(const C2dImagePointPx & v1, const C2dImagePointPx & v2)
{
    return v1.y() * v2.x() - v1.x() * v2.y(); //Match signs with 3D cross prod
}

template<class TVecType>
double tanAngleBetweenVectors(const TVecType & seg1vec, const TVecType & seg2vec)
{
    const double dDot = seg1vec.dot(seg2vec);
    const double dCross = crossProd(seg1vec, seg2vec);

    return penalisedTanAngle(dDot, dCross);
}


template<class TControlPoint>
bool CPolyline_base<TControlPoint>::tooShort() const //a threshold used in a few places for getting rid of very short polylines--e.g. shorter than 1.5*thickness
{
    return (numPoints() < 2 || length() < minLength()); // short compared to its thickness, or single point
}

template<class TControlPoint>
typename CPolyline_base<TControlPoint>::TLineType CPolyline_base<TControlPoint>::approxAsLine() const
{
    return CPolyline_base<TControlPoint>::TLineType(getStartPoint(), getFinishPoint());
}

//Allow iteration over endpoints
template<class TControlPoint>
const typename CPolyline_base<TControlPoint>::TVecType & CPolyline_base<TControlPoint>::getEndPoint(const eEndpoints end) const
{
    if(end == eStart)
        return getStartPoint();
    else if(end == eFinish)
        return getFinishPoint();
    else
        THROW("Unhandled endpoint");
}

//Allow iteration over endpoints
template<class TControlPoint>
const typename CPolyline_base<TControlPoint>::TControlPoint & CPolyline_base<TControlPoint>::getEndControlPoint(const eEndpoints end) const
{
    if(end == eStart)
        return aControlPoints[0];
    else if(end == eFinish)
        return aControlPoints[numPoints()-1];
    else
        THROW("Unhandled endpoint");
}

template<class TControlPoint>
const typename CPolyline_base<TControlPoint>::TVecType CPolyline_base<TControlPoint>::midPoint() const
{
    if(numPoints()%2==1)
        return aControlPoints[numPoints()/2].getPoint();
    else
        return segment((numPoints()/2) - 1).midPoint();
}

template<class TControlPoint>
const typename CPolyline_base<TControlPoint>::TVecType CPolyline_base<TControlPoint>::endDirection(const eEndpoints end) const
{
    return ((end==eStart) ? -1.0 : 1.0) * endSegment(end).direction();
}

template<class TControlPoint>
double CPolyline_base<TControlPoint>::tanKinkAngleAtPoint(const int nPoint) const
{
    if(IS_DEBUG) CHECK(nPoint < 1 || nPoint > numPoints()-2, "nPoint OOB");
    const TVecType seg1vec = segment(nPoint-1).startToFinish();
    const TVecType seg2vec = segment(nPoint).startToFinish();

    return tanAngleBetweenVectors(seg1vec, seg2vec);
}

template<class TControlPoint>
double CPolyline_base<TControlPoint>::kinkAngleAtPoint(const int nPoint) const
{
    if(IS_DEBUG) CHECK_P(nPoint < 1 || nPoint > numPoints()-2, nPoint, "nPoint OOB");

    return angleBetween3Points<TVecType>(aControlPoints[nPoint-1].getPoint(), aControlPoints[nPoint].getPoint(), aControlPoints[nPoint+1].getPoint());
}

template<class TControlPoint>
void CPolyline_base<TControlPoint>::reverseDirection()
{
    for(int i=0; i<(int)numPoints()/2; i++) {
        std::swap(aControlPoints[i], aControlPoints[numPoints() - (i+1)]);
    }
}

template<class TControlPoint>
double C2dPolyline_base<TControlPoint>::getCurvature() const
{
    if(this->numSegments() <= 1)
        return 0;

    const double dAngleChange = this->segment(0).signedAngle(this->segment(this->numSegments()-1));
    const double dLength = this->length();
    return dAngleChange/dLength;
}


template<class TControlPoint>
const int CPolyline_base<TControlPoint>::closestSegmentIdx(const typename CPolyline_base<TControlPoint>::TVecType & p) const
{
    if(IS_DEBUG) CHECK(numPoints() < 2, "Too short polyline");

    double dMinDist = HUGE;
    int nClosestSegment = -1;

    for(int i=0; i<numSegments(); i++) {
        const double dDist = segment(i).closestDistance(p);
        if(dDist < dMinDist) {
            dMinDist = dDist;
            nClosestSegment = i;
        }
    }
    CHECKOOB(nClosestSegment, numSegments());
    return nClosestSegment;
}

template<class TControlPoint>
const typename CPolyline_base<TControlPoint>::TLineType CPolyline_base<TControlPoint>::closestSegment(const typename CPolyline_base<TControlPoint>::TVecType & p) const
{
    return segment(closestSegmentIdx(p));
}

template<class TControlPoint>
const typename CPolyline_base<TControlPoint>::TVecType CPolyline_base<TControlPoint>::closestPoint(const typename CPolyline_base<TControlPoint>::TVecType & p) const
{
    int nClosestSeg=-1;
    double dDistToPoly_sq = HUGE;
    return closestPoint_segIdx_distSq(p, nClosestSeg, dDistToPoly_sq).getPoint();
}

template<class TControlPoint>
const typename CPolyline_base<TControlPoint>::TVecType CPolyline_base<TControlPoint>::closestPointAndPosition(const typename CPolyline_base<TControlPoint>::TVecType & p, double & dPosition_index) const
{
    int nClosestSeg=-1;
    double dDistToPoly_sq = HUGE;
    return closestPoint_segIdx_distSq(p, nClosestSeg, dDistToPoly_sq, &dPosition_index).getPoint();
}

template<class TControlPoint>
const typename CPolyline_base<TControlPoint>::TVecType CPolyline_base<TControlPoint>::closestPointToBoundedLine(const typename TControlPoint::TLineType & boundedline) const
{
    TVecType closestPointOverall;
    double dDistToPoly = HUGE;
    for(int i=0; i<numSegments(); i++) {
        const double dDist = segment(i).closestDistanceToBoundedLine(boundedline);
        if(dDist < dDistToPoly) {
            dDistToPoly = dDist;
            closestPointOverall =  segment(i).closestPointToBoundedLine(boundedline);
        }
    }
    return closestPointOverall;
}

template<class TControlPoint>
const TControlPoint CPolyline_base<TControlPoint>::interpolateControlPoint(const typename CPolyline_base<TControlPoint>::TVecType & closestPointOverall, const int nClosestSeg, const double t_closest) const
{
    CHECKOOB(nClosestSeg, numPoints()-1);
    if(TControlPoint::HAS_THICKNESS) {

        //Use end widths as interpolated widths for points off the end (otherwise we'll get negative widths)
        if(t_closest <= 0)
            return TControlPoint(closestPointOverall, aControlPoints[nClosestSeg].getWidth());

        if(t_closest >= 1)
            return TControlPoint(closestPointOverall, aControlPoints[nClosestSeg+1].getWidth());

        const double dWidth = t_closest*aControlPoints[nClosestSeg+1].getWidth() + (1-t_closest)*aControlPoints[nClosestSeg].getWidth();
        //Don't check here as we sometime are manipulating invalid polylines, e.g. to fix them, CCheckThickness::checkSensible(dWidth, structureType, getMeasurementType(TVecType::RowsAtCompileTime));

        if(dWidth <= 0)
            REPEAT(100, cout << "Warning: interpolated control point width less than 0 " << TO_STRING(dWidth) << TO_STRING(t_closest) << TO_STRING(nClosestSeg) << endl << toString() << endl);

        /*if(dWidth < TControlPoint::minWidth())
        {
            ALWAYS_VERBOSE;
            COUT(dWidth);
            COUT(t_closest);
            COUT(nClosestSeg);
            COUT(this->toString());
            COUT(t_closest*aControlPoints[nClosestSeg].getWidth());
            COUT((1-t_closest)*aControlPoints[nClosestSeg+1].getWidth());
            THROW("Negative width interpolated");
        }*/
        return TControlPoint(closestPointOverall, dWidth);
    } else {
        return TControlPoint(closestPointOverall);
    }
}

template<class TControlPoint>
typename CPolyline_base<TControlPoint>::TVecType CPolyline_base<TControlPoint>::directionAtPoint(const int nPoint) const
{
    if(nPoint==0)
        return endDirection(eStart);
    
    if(nPoint==numPoints()-1)
        return endDirection(eFinish);
        
    return (aControlPoints[nPoint+1].getPoint()-aControlPoints[nPoint-1].getPoint()).normalized();
}

template<class TControlPoint>
const TControlPoint CPolyline_base<TControlPoint>::closestPoint_segIdx_distSq(const typename CPolyline_base<TControlPoint>::TVecType & p, int & nClosestSeg, double & dDistToPoly_sq, double * pdPosition_index) const
{
    //if(IS_DEBUG)
    CHECK_P(numPoints() < 2, numPoints(), "Polyline has 0 or 1 point--probably this polyline should never be used");

    TVecType closestPointOverall;
    double t_closest=-1;
    dDistToPoly_sq = HUGE;

    for(int i=0; i<numSegments(); i++) {
        double t=-1;
        const TVecType closestPoint = segment(i).closestPointAndt(p, t, false);
        const double dDist_sq = (closestPoint - p).squaredNorm();
        if(dDist_sq < dDistToPoly_sq) {
            dDistToPoly_sq = dDist_sq;
            closestPointOverall = closestPoint;
            nClosestSeg = i;
            t_closest = t;
        }
    }

    if(pdPosition_index)
        *pdPosition_index = (double)nClosestSeg + t_closest;

    CHECK(dDistToPoly_sq >= HUGE, "Failed to find a closest point (size 0?)");
    CHECK_P(nClosestSeg < 0 || nClosestSeg >= numSegments(), nClosestSeg, "Found segment idx OOB");
    return interpolateControlPoint(closestPointOverall, nClosestSeg, t_closest);
}

template<class TControlPoint>
const optional<const typename CPolyline_base<TControlPoint>::TVecType> CPolyline_base<TControlPoint>::closestPointStrictlyOnPoly(const typename CPolyline_base<TControlPoint>::TVecType & p) const
{
    const TVecType closestPt = closestPoint(p);

    const bool bOffStart = closestPt==getStartPoint();
    const bool bOffEnd = closestPt==getFinishPoint();

    if(!bOffStart && !bOffEnd)
        return closestPt;

    const double dDistanceToPoly = (p-closestPt).norm();
    const double dNoiseMargin = 1; // usually in pixels
    if(bOffStart) {
        const double dDistanceToStartUnbounded = segment(0).unboundedLine().closestDistance(p);
        if(dDistanceToStartUnbounded+dNoiseMargin > dDistanceToPoly)
            return closestPt;
    }

    if(bOffEnd) {
        const double dDistanceToEndUnbounded = segment(numSegments()-1).unboundedLine().closestDistance(p);
        if(dDistanceToEndUnbounded+dNoiseMargin > dDistanceToPoly)
            return closestPt;
    }

    return optional<const TVecType>();
}

template<class TControlPoint>
typename CPolyline_base<TControlPoint>::TVecType CPolyline_base<TControlPoint>::direction() const
{
    if(IS_DEBUG) CHECK(numPoints() < 2, "Too short to have a direction");

    const TVecType dir = (getFinishPoint() - getStartPoint());
    return dir.normalized();
}

template<class TControlPoint>
const typename CPolyline_base<TControlPoint>::TVecType CPolyline_base<TControlPoint>::pointAtDistance(const double dDistanceAlongPoly) const
{
    double dCumulativeLength=0;
    for(int i=0; i<numSegments(); i++) {
        const double dLength = segment(i).length();
        if(dCumulativeLength+dLength > dDistanceAlongPoly) {
            const double t = (dDistanceAlongPoly - dCumulativeLength)/dLength;
            return segment(i).getStartPoint() + t*segment(i).startToFinish();
        }
        dCumulativeLength+=dLength;
    }
    return getFinishPoint();
}

template<class TControlPoint>
const typename CPolyline_base<TControlPoint>::TLineType CPolyline_base<TControlPoint>::segmentAtDistance(const double dDistanceAlongPoly) const
{
    double dCumulativeLength=0;
    for(int i=0; i<numSegments()-1; i++) {
        const double dLength = segment(i).length();
        if(dCumulativeLength+dLength > dDistanceAlongPoly) {
            return segment(i);
        }
        dCumulativeLength+=dLength;
    }
    return segment(numSegments()-1);
}

template<class TControlPoint>
const typename CPolyline_base<TControlPoint>::TVecType CPolyline_base<TControlPoint>::pointAtPosition(const double dPosition) const
{
    if(dPosition<0)
        return getStartPoint();

    if(dPosition < (double)numSegments()) {
        const int nSegment = (int)floor(dPosition);
        const double t = dPosition-nSegment;
        return segment(nSegment).pointFromt(t);
    }

    return getFinishPoint();
}

template<class TControlPoint>
double CPolyline_base<TControlPoint>::thicknessAtPosition_index(const double dPosition) const
{
    if(dPosition<0)
        return getEndControlPoint(eStart).getWidth();

    if(dPosition < (double)numSegments()) {
        const int nSegment = (int)floor(dPosition);
        const double t = dPosition-nSegment;
        return (1-t)*aControlPoints[nSegment].getWidth() + t*aControlPoints[nSegment+1].getWidth();
    }

    return getEndControlPoint(eFinish).getWidth();
}

template<class TControlPoint>
const TControlPoint CPolyline_base<TControlPoint>::controlPointAtPosition(const double dPosition) const
{
    if(dPosition<0)
        return getEndControlPoint(eStart);

    if(dPosition < (double)numSegments()) {
        const int nSegment = (int)floor(dPosition);
        const double t = dPosition-nSegment;
        TVecType vec = segment(nSegment).pointFromt(t);
        if(!TControlPoint::HAS_THICKNESS)
            return TControlPoint(vec);

        const double dWidth = (1-t)*aControlPoints[nSegment].getWidth() + t*aControlPoints[nSegment+1].getWidth();
        return TControlPoint(vec, dWidth);
    }

    return getEndControlPoint(eFinish);
}

template<class TControlPoint>
const TControlPoint CPolyline_base<TControlPoint>::closestPointAndWidth(const typename CPolyline_base<TControlPoint>::TVecType & p, double * pdPosition_idx) const
{
    int nClosestSeg = -1;
    double dDistSq = -1;
    return closestPoint_segIdx_distSq(p, nClosestSeg, dDistSq, pdPosition_idx);
}

template<class TControlPoint>
double C2dPolyline_base<TControlPoint>::signedDistanceToPoint(const C2dImagePointPx & p) const
{
    const C2dBoundedLine closestSeg = this->closestSegment(p);
    const C2dImagePointPx closestPointOnLine = closestSeg.closestPoint(p);
    const TEigen2dPoint perp = perpendicular(closestSeg.startToFinish());
    const TEigen2dPoint vecToClosestPoint = closestPointOnLine-p;
    const double dDirection = perp.dot(vecToClosestPoint) > 0 ? 1 : -1;
    const double dDistToCentre = dDirection*vecToClosestPoint.norm();
    return dDistToCentre;
}

template<class TControlPoint>
double CPolyline_base<TControlPoint>::distanceToPoint(const typename CPolyline_base<TControlPoint>::TVecType & p, const bool bIncludeThickness) const
{
    TVecType closest = closestPoint(p);
    const double dDistToCentre = (closest-p).norm();
    if(!bIncludeThickness)
        return dDistToCentre;

    const double dThicknessHere = getWidth(closest); //inefficient
    const double dDistToEdge = dDistToCentre - 0.5*dThicknessHere;
    return (dDistToEdge>0) ? dDistToEdge : 0;
}

template<class TControlPoint>
double CPolyline_base<TControlPoint>::positionOfPoint_index(const typename CPolyline_base<TControlPoint>::TVecType & p) const
{
    const bool bVerbose = false;

    COUT2("Computing position on polyline by index...", p);
    
    CHECK(numPoints() == 0, "No points in polyline");
    
    if(numPoints() == 1)
        return 0;

    //First find the closest section
    double dClosestDist = HUGE;
    int nSeg = -1;
    for(int i=0; i<numSegments(); i++) {
        const double dDist = segment(i).closestDistance(p);
        if(dDist < dClosestDist) {
            dClosestDist = dDist;
            nSeg = i;
        }
    }
    COUT(nSeg);
    COUT(segment(nSeg));
    COUT(segment(nSeg).get_t(p));

    return nSeg + segment(nSeg).get_t(p);
}

template<class TControlPoint>
double CPolyline_base<TControlPoint>::indexToDistance(const double dIndex) const
{
    const int nSeg = (int)floor(clip<double>(dIndex, 0, numSegments()-1));

    double dDist = 0;
    for(int i=0; i<nSeg; i++)
        dDist += segment(i).length();

    const double dPosOnSegment = dIndex - (double)nSeg;
    const double dSegLength = segment(nSeg).length();

    return dDist + dPosOnSegment*dSegLength;
}

template<class TControlPoint>
CRange CPolyline_base<TControlPoint>::indexRangeToDistanceRange(const CRange & range_index) const
{
    return CRange(indexToDistance(range_index.getMin()), indexToDistance(range_index.getMax()));
}

template<class TControlPoint>
CRange CPolyline_base<TControlPoint>::getIndexRange() const
{
    return CRange(0, numSegments());
}

template<class TControlPoint>
CRange CPolyline_base<TControlPoint>::getDistanceRange() const
{
    return CRange(0, length());
}

template<class TControlPoint>
double CPolyline_base<TControlPoint>::positionOfPoint_distance(const typename CPolyline_base<TControlPoint>::TVecType & p) const
{
    const bool bVerbose = false;

    COUT("Computing position on polyline by length...");
    
    CHECK(numPoints() == 0, "No points in polyline");
    
    if(numPoints() == 1)
        return 0;

    //First find the closest section
    double dClosestDist = HUGE;
    double dCumulativeLength = 0, dCumulativeLengthAtClosestDistance = 0;
    int nSeg = -1;
    for(int i=0; i<numSegments(); i++) {
        const double dDist = segment(i).closestDistance(p);
        if(dDist < dClosestDist) {
            dClosestDist = dDist;
            nSeg = i;

            dCumulativeLengthAtClosestDistance=dCumulativeLength;
        }
        dCumulativeLength += segment(i).length();
    }

    return dCumulativeLengthAtClosestDistance + segment(nSeg).length()*segment(nSeg).get_t(p);
}

template<class TControlPoint>
double CPolyline_base<TControlPoint>::length() const
{
    double dLength = 0;
    for(int i=0; i<numSegments(); i++) {
        //dLength += segment(i).length(); This throws if endpoints are equal--might be that we just haven't fixed the polyline yet.
        dLength += (aControlPoints[i].getPoint() - aControlPoints[i+1].getPoint()).norm();//
    }
    return dLength;
}

/*template<class TControlPoint>
void CPolyline_base<TControlPoint>::extendPolyEnd(const double dExtensionPx, const bool bVerbose)
{
    const TVecType endDir = (aControlPoints[numPoints()-1].getPoint() - aControlPoints[numPoints()-2].getPoint()).normalized();
    const TVecType newVirtualEndPoint = (getFinishPoint() + / *segment(numPoints()-2).direction()* / endDir * dExtensionPx).eval();
    aControlPoints.push_back(TControlPoint(newVirtualEndPoint, aControlPoints[numPoints()-1].getWidth() ));
    COUT2("Extending 2D end to ", newVirtualEndPoint);
    if(IS_DEBUG) CHECK(!zero(kinkAngleAtPoint(numPoints()-2)), "Extend end failed");
}

template<class TControlPoint>
void CPolyline_base<TControlPoint>::extendPolyStart(const double dExtensionPx, const bool bVerbose)
{
    const TVecType endDir = (aControlPoints[1].getPoint() - aControlPoints[0].getPoint()).normalized();
    const TVecType newVirtualStartPoint = (getStartPoint() - / *segment(0).direction()* / endDir * dExtensionPx).eval();
    aControlPoints.insert(aControlPoints.begin(), TControlPoint(newVirtualStartPoint, aControlPoints[0].getWidth() ));
    COUT2("Extending 2D start to ", newVirtualStartPoint);
    if(IS_DEBUG) CHECK(!zero(kinkAngleAtPoint(1)), "Extend start failed");
}*/


template<class TControlPoint>
double CPolyline_base<TControlPoint>::getWidth(const typename CPolyline_base<TControlPoint>::TVecType & x) const
{
    if(IS_DEBUG) CHECK(!zero(distanceToPoint(x, false)), "getWidth is only for points *on* the polyline");

    if(numPoints()==1) //Relying on 'if(IS_DEBUG) CHECK' here
        return aControlPoints[0].getWidth();

    for(int i=0; i<numSegments(); i++) {
        const double dDist = segment(i).closestDistance(x);
        if(zero(dDist)) {
            //interpolate width from i and i+1
            double dRelPos = (aControlPoints[i].getPoint() - x).norm() - segment(i).length();

            return (1-dRelPos)*aControlPoints[i].getWidth() + dRelPos*aControlPoints[i].getWidth();
        }
    }
    THROW("x should lie on this polyline");
}

template<class TControlPoint>
typename CPolyline_base<TControlPoint>::TLineType CPolyline_base<TControlPoint>::segment(const int nSegStart) const
{
    CHECK_P(nSegStart >= numSegments() || nSegStart < 0, numSegments()+10000*nSegStart, "Segment start point OOB");
    return TLineType((*this)[nSegStart].getPoint(), (*this)[nSegStart+1].getPoint());
}

template<class TControlPoint>
typename CPolyline_base<TControlPoint>::TLineType CPolyline_base<TControlPoint>::endSegment(const eEndpoints end) const
{
	CHECK(bClosed, "Closed polylines don't really have 'ends'");
    if(IS_DEBUG) CHECK(numSegments()<=0, "Segment start point OOB");

    if(end==eFinish)
        return segment(numSegments()-1);
    else
        return segment(0);
}

template<class TControlPoint>
void CPolyline_base<TControlPoint>::toVector(Eigen::VectorXd & x) const // For LM optimisation
{
    if(IS_DEBUG) CHECK(x.rows() != dimension(), "Param vector dimension mismatch");

    int idx=0;
    BOOST_FOREACH(const TControlPoint & p, *this) {
        x.segment<TControlPoint::PARAMS_PER_POLY_CONTROL_POINT>(idx) = p.asVector();
        idx += TControlPoint::PARAMS_PER_POLY_CONTROL_POINT;
    }
}

template<class TControlPoint>
eSuccessStatus CPolyline_base<TControlPoint>::onePointFromVector(const Eigen::VectorXd & x, const int nIndex) // For LM optimisation
{
    TControlPoint & p = aControlPoints[nIndex];

    const Eigen::Matrix<double, TControlPoint::PARAMS_PER_POLY_CONTROL_POINT, 1> asVec = x.segment<TControlPoint::PARAMS_PER_POLY_CONTROL_POINT>(nIndex*TControlPoint::PARAMS_PER_POLY_CONTROL_POINT);

    p = TControlPoint(asVec, CConstructFromLMVector());

    return eSuccess;
}

template<class TControlPoint>
eSuccessStatus CPolyline_base<TControlPoint>::fromVector(const Eigen::VectorXd & x) // For LM optimisation
{
    const int nNumParamsPerControlPoint = TControlPoint::PARAMS_PER_POLY_CONTROL_POINT;
    const int nNumPoints = (int)x.size()/nNumParamsPerControlPoint;
    aControlPoints.resize(nNumPoints);

    for(int nPoint=0; nPoint < nNumPoints; nPoint++) {
        if(onePointFromVector(x, nPoint) == eFail)
            return eFail;
    }

    return eSuccess;
}

template<class TControlPoint>
void CPolyline_base<TControlPoint>::clearAndReserve(const int nReserveSize)
{
    aControlPoints.clear();
    aControlPoints.reserve(nReserveSize);
}

template<class TControlPoint>
void CPolyline_base<TControlPoint>::resize(const int nNewSize)
{
    aControlPoints.resize(nNewSize);
}

template<class TControlPoint>
void CPolyline_base<TControlPoint>::erase(typename CPolyline_base<TControlPoint>::iterator start, typename CPolyline_base<TControlPoint>::iterator finish )
{
    aControlPoints.erase(start, finish);
}
template<class TControlPoint>
void CPolyline_base<TControlPoint>::insert(typename CPolyline_base<TControlPoint>::iterator pos, typename CPolyline_base<TControlPoint>::const_iterator start, typename CPolyline_base<TControlPoint>::const_iterator finish )
{
    aControlPoints.insert(pos, start, finish);
}
template<class TControlPoint>
void CPolyline_base<TControlPoint>::push_front(const TControlPoint & controlPoint )
{
    aControlPoints.insert(aControlPoints.begin(), controlPoint);
}

template<class TControlPoint>
void CPolyline_base<TControlPoint>::push_end(const TControlPoint & p, const eEndpoints end)
{
    if(end==eStart)
        push_front(p);
    else
        push_back(p);
}

template<class TControlPoint>
void CPolyline_base<TControlPoint>::concatenate(const CPolyline_base<TControlPoint> & other)
{
    aControlPoints.insert(aControlPoints.end(), other.aControlPoints.begin(), other.aControlPoints.end() );
}


template<class TControlPoint>
void CPolyline_base<TControlPoint>::push_back(const TControlPoint & p)
{
    /* This happens normally when we interleave polylines: if(IS_DEBUG && aControlPoints.size() > 0)
        CHECK(p.getPoint() == aControlPoints.back().getPoint(), "Adding identical point to a polyline"); */

    aControlPoints.push_back(p);
}

template<class TControlPoint>
void CPolyline_base<TControlPoint>::pop_front()
{
    aControlPoints.erase(aControlPoints.begin());
}
template<class TControlPoint>
void CPolyline_base<TControlPoint>::pop_back()
{
    aControlPoints.pop_back();
}

template<class TControlPoint>
void CPolyline_base<TControlPoint>::pop_end(const eEndpoints end)
{
    if(end == eStart)
        pop_front();
    else
        pop_back();
}


template<class TControlPoint>
int CPolyline_base<TControlPoint>::closestControlPoint(const TVecType & p) const
{
    int nClosest=-1, n=0;
    double dClosest=HUGE;

    BOOST_FOREACH(const TControlPoint & px, *this) {
        double dDistSq = (px.getPoint() - p).squaredNorm();
        if(dDistSq < dClosest) {
            dClosest = dDistSq;
            nClosest = n;
        }

        n++;
    }
    return nClosest;
}


template<class TThickPoly, class TThinPoly>
void thickToThin(const TThickPoly & thick, TThinPoly & thin)
{
    if(thick.numPoints()==0)
        return;

    thin.clearAndReserve(thick.numPoints());
    BOOST_FOREACH(const typename TThickPoly::TControlPoint & p, thick) {
        thin.push_back(typename TThinPoly::TControlPoint(p.getPoint()));
    }
}

template<class TThinPoly, class TThickPoly>
void thinToThick(const TThinPoly & thin, TThickPoly & thick, const double dThickness)
{
    if(thin.numPoints()==0)
        return;

    thick.clearAndReserve(thin.numPoints());
    BOOST_FOREACH(const typename TThinPoly::TControlPoint & p, thin) {
        thick.push_back(typename TThickPoly::TControlPoint(p.getPoint(), dThickness));
    }
}

template<class TControlPoint>
double CPolyline_base<TControlPoint>::minLength() const //TODO: fix and make arguments (+ add maxLength?)
{
    if(TControlPoint::HAS_THICKNESS)
        return 1.95*averageThickness();
    /*if(TVecType::RowsAtCompileTime == 3)
        return std::min<double>(1.95*averageThickness(), 0.025);
    else
        return std::min<double>(1.95*averageThickness(), 4 );*/
    else
        return (TVecType::RowsAtCompileTime == 3) ? 0.02 : 5;
}

template<class TControlPoint>
double CPolyline_base<TControlPoint>::getDPEpsilon() const
{
    if(TControlPoint::HAS_THICKNESS) {
        double epsilon = std::max<double>(DPMinEpsilon(), 0.3*averageThickness());
        epsilon = std::min<double>(DPMaxEpsilon(), epsilon);
        return epsilon;
    } else {
        return DPMinEpsilon();
    }
}


template<class TLineType_in>
std::ostream & CPolylineControlPoint<TLineType_in>::pp(std::ostream & s) const
{
    s << "ControlPoint=" << getPoint().transpose();
    return s;
}

template<class TLineType_in>
std::ostream & CPolylineControlPointWithThickness<TLineType_in>::pp(std::ostream & s) const
{
    s << "ControlPoint=" << getPoint().transpose() << "-width=" << getWidth();
    return s;
}

template<class TLineType_in>
CPolylineControlPointWithThickness<TLineType_in>::CPolylineControlPointWithThickness(const TVecType & x, const double dWidth) : x(x), dWidth(dWidth)
{
    //checkOk(); The polylines are all over the place during optimisation so we don't want these strict checks
    checkBadNum();
}

template<class TLineType_in>
CPolylineControlPointWithThickness<TLineType_in>::CPolylineControlPointWithThickness(const TLMVector & lmVec, const CConstructFromLMVector &) : x(lmVec.template head<TVecType::RowsAtCompileTime>()), dWidth(lmVec.template tail<1>()(0))
{
    //checkOk();
    checkBadNum();
}

template<class TLineType_in>
double CPolylineControlPoint<TLineType_in>::distanceToPoint(const TVecType & vec) const
{
    const double dDist = (vec-getPoint()).norm();
    return dDist;
}

template<class TLineType_in>
double CPolylineControlPointWithThickness<TLineType_in>::distanceToPoint(const TVecType & vec) const
{
    const double dDist = (vec-getPoint()).norm();
    return std::max<double>(dDist-getRad(), 0);
}


#endif //POLYLINE_IPP
