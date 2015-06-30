#ifndef LINES_IPP
#define LINES_IPP

#include "lines.h"
#include "vectorAngles.h"

/* 
* This header implements the templatised functions for C2dLine etc.
* 
* #include it *only if* we want to be able to inline these functions.
*/

/////////// Unbounded line //////////////

template<class TVecType>
CUnboundedLine_base<TVecType>::CUnboundedLine_base()
{
    direction.setZero();
    pointOnLine.setZero();
}

template<class TVecType>
void CUnboundedLine_base<TVecType>::from2Points(const TVecType & p1, const TVecType & p2)
{
    if(IS_DEBUG) CHECK(p1 == p2, "Line end points are equal");
    direction = (p2 - p1).normalized();
    pointOnLine = p1;

    //direction.checkInit();
}

template<class TVecType>
void CUnboundedLine_base<TVecType>::fromPointAndDir(const TVecType & p, const typename TVecType::Base & dir)
{
    if(IS_DEBUG) CHECK(dir.squaredNorm() == 0, "Line has no direction");
    direction = dir;
    direction.normalize();
    pointOnLine = p;

    //direction.checkInit();
}

//Returns a point on the line with param lambda
template<class TVecType>
TVecType CUnboundedLine_base<TVecType>::point(const double & lambda) const
{
    return pointOnLine + direction * lambda;
}


/*
   pointOnLine + lambda * direction = other.pointOnLine + mu * other.direction
   Solve for lambda, mu
    *
    * Checked version -- will throw for pll lines

template<class TVecType>
TVecType CUnboundedLine_base<TVecType>::findIntersection(const TUnboundedLine & other) const
{
    return *findIntersection_robust(other);
} */

template<class TVecType>
double CUnboundedLine_base<TVecType>::undirectedAngle(const TUnboundedLine & other) const
{
    return angleBetweenUnitVectors(direction, other.direction, false);
}

template<class TVecType>
TVecType CUnboundedLine_base<TVecType>::closestPoint(const TVecType & p) const
{
    TVecType pointToPointOnLine = p - pointOnLine;
    double t = pointToPointOnLine.dot(direction); // AP.x*AB.x + AP.y*AB.y;
    return point(t);
}

template<class TVecType>
double CUnboundedLine_base<TVecType>::closestDistance(const TVecType & p) const
{
    return (p - closestPoint(p)).norm();
}

/*template<class TVecType>
template<class TBoundedLine>
TVecType CUnboundedLine_base<TVecType>::closestPointToBoundedLine_int(const TBoundedLine & boundedLine) const
{
    TVecType intersection = findIntersection(boundedLine.unboundedLine());

    //If intersection is on the bounded line, return intersection, else return whichever point on this line is closest to the bounded line endpoint
    TVecType closestPointOnBounded = boundedLine.closestPoint(intersection);

    return closestPoint(closestPointOnBounded);
}

template<class TVecType>
template<class TBoundedLine>
double CUnboundedLine_base<TVecType>::closestDistanceToBoundedLine_int(const TBoundedLine & boundedLine) const
{
    return boundedLine.closestDistance(findIntersection(boundedLine.unboundedLine()));
}*/

///////// C3dLine //////////////


///////// C2dLine //////////////


///////// C3dBoundedLine //////////////

//C3dBoundedLine::C3dBoundedLine() : T3dBoundedLine_base(C3dWorldPoint::Zeros(), C3dWorldPoint::Zeros());

///////// C2dBoundedLine //////////////


//////// Template base //////////////

template <class TVecType>
CBoundedLine_base<TVecType>::CBoundedLine_base(const TVecType & p1, const TVecType & p2)
{
    aEndPoints[eStart] = p1;
    aEndPoints[eFinish] = p2;

    //Skip check when running online
    CHECK_P(getStartPoint() == getFinishPoint() && getStartPoint() != TVecType::uninit(), getStartPoint(), "Line end points are equal (may happen in optimisations with very low probability)");
}

template <class TVecType>
optional< typename CBoundedLine_base<TVecType>::TBoundedLine> 
CBoundedLine_base<TVecType>::trim(const double dTrimThresh) const
{
    const double dLength = length();
    
    if(dLength <= 2*dTrimThresh)
        return optional< CBoundedLine_base<TVecType> > ();
    
    const typename TVecType::Base dir = direction();

    return CBoundedLine_base<TVecType>(getStartPoint() + dir*dTrimThresh, getFinishPoint() - dir*dTrimThresh);
}

template <class TVecType>
CBoundedLine_base<TVecType>::CBoundedLine_base()
{
    aEndPoints[eStart] = aEndPoints[eFinish] = TVecType::Zero();
}

template <class TVecType>
typename CBoundedLine_base<TVecType>::TUnboundedLine CBoundedLine_base<TVecType>::unboundedLine() const { 
    TUnboundedLine lineUnbounded; 
    lineUnbounded.from2Points(getStartPoint(), getFinishPoint());
    return lineUnbounded;
}

template <class TVecType>
double CBoundedLine_base<TVecType>::length() const
{
    return startToFinish().norm();
}

template <class TVecType>
const TVecType CBoundedLine_base<TVecType>::midPoint() const
{
    return 0.5*(getStartPoint() + getFinishPoint());
}

template <class TVecType>
double CBoundedLine_base<TVecType>::get_t(const TVecType & point) const
{
    const TVecType AP = point - getStartPoint();
    const TVecType AB = startToFinish();
    const double ab2 = AB.squaredNorm();
    const double ap_ab = AP.dot(AB); // AP.x*AB.x + AP.y*AB.y;
    const double t = ap_ab / ab2;
    return t;
}

template <class TVecType>
TVecType CBoundedLine_base<TVecType>::pointFromt(const double t) const 
{
    return getStartPoint() + startToFinish() * t;
}

template <class TVecType>
TVecType CBoundedLine_base<TVecType>::closestPointAndt(const TVecType & p, double & t, const bool bClamp_t) const 
{
    const double t_unclamped = get_t(p);
    const double t_clamped = segmentClamp(t_unclamped);
    
    t = bClamp_t ? t_clamped : t_unclamped;
    
    return pointFromt(t_clamped);
}

template <class TVecType>
TVecType CBoundedLine_base<TVecType>::closestPoint(const TVecType & point) const
{
    double t;
    return closestPointAndt(point, t, false);
}
/*template <class TVecType>
template<class TLineType>
TLineType CBoundedLine_base<TVecType>::reversed() const
{
    //return CBoundedLine_base<TVecType>::TBoundedLine(getFinishPoint(), getStartPoint());
    return TLineType(getFinishPoint(), getStartPoint());
}*/

template <class TVecType>
typename TVecType::Base CBoundedLine_base<TVecType>::direction() const
{
    return unboundedLine().getDirection();
}

template <class TVecType>
double CBoundedLine_base<TVecType>::closestDistance(const TVecType & point) const
{
    return (point - closestPoint(point)).norm();
}

template <class TVecType>
double CBoundedLine_base<TVecType>::undirectedAngle(const TBoundedLine & line) const
{
    return undirectedAngle(line.unboundedLine());
}

template <class TVecType>
double CBoundedLine_base<TVecType>::directedAngle(const TBoundedLine & line) const
{
    return directedAngle(line.unboundedLine());
}

template <class TVecType>
double CBoundedLine_base<TVecType>::undirectedAngle(const TUnboundedLine & line) const
{
    return unboundedLine().undirectedAngle(line);
}

template <class TVecType>
double CBoundedLine_base<TVecType>::directedAngle(const TUnboundedLine & line) const
{
    return angleBetweenUnitVectors(line.getDirection(), unboundedLine().getDirection(), true);
}

#endif //LINES_IPP
