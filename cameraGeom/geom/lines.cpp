//#include "stdafx.h"

#include "newCamera.h"
#include "lines.ipp"
#include <Eigen/Dense>
#include "vectorAngles.h"

///////// C3dLine /////////////

/**
 * @brief findIntersection for 2D lines. Can handle pll lines (return nothing)
 * @param other
 * @return
 */
optional<const C2dImagePointPx> C2dLine::findIntersection(const C2dLine & other) const
{
    Eigen::Matrix<double, 2, 2> directions;
    directions.block<2,1>(0,0) = direction;
    directions.block<2,1>(0,1) = -other.direction;

    if(fabs(directions.determinant()) < 1e-8)//Should usually be about 1
        return optional<const C2dImagePointPx>();

    const C2dImagePointPx diff = other.pointOnLine - pointOnLine;

    const TEigen2dPoint lambda_mu = directions.inverse() * diff;
    const C2dImagePointPx pointOnLine = point(lambda_mu(0));
    if(IS_DEBUG) CHECK(!zero((other.point(lambda_mu(1)) - pointOnLine).squaredNorm()), "findIntersection failed (parallel lines?)");

    return optional<const C2dImagePointPx>(pointOnLine);
}

//joinPoint is the point on the sphere closest to this line.
void C3dBoundedLine::closestPointOnSphere(const C3dPolylineControlPointWithThickness & closestHeadPart, C3dWorldPoint & joinPoint, double & dDistance) const
{
    const C3dWorldPoint closestOnLine = closestPoint(closestHeadPart.getPoint());
    const TEigen3dPoint vectorToCentre = closestHeadPart.getPoint() - closestOnLine;
    const double dDistCentreToLine = vectorToCentre.norm();
    
    const double dRad = closestHeadPart.getRad();
    
    //2 cases: 1: line intersects sphere
    if(dDistCentreToLine < dRad)
    {
        //The line might poke a little way into the sphere. Change to unbounded to simplify computations:
        const C3dWorldPoint closestOnUnboundedLine = unboundedLine().closestPoint(closestHeadPart.getPoint());
        const double dDistCentreToUnboundedLine = (closestHeadPart.getPoint() - closestOnUnboundedLine).norm();
        
        //The closest on the sphere, the closest point on the line, and the intersection point form a right angle triangle.
        //Select the point on the line where the hypotenuse has length dRad
        // dRad^2 = dDistCentreToLine^2 + offset^2
        const double dOffset = sqrt(sqr(dRad) - sqr(dDistCentreToUnboundedLine));
        joinPoint = closestOnUnboundedLine - dOffset*direction();
        dDistance = 0; 
        
    } else {   
        //2: line doesn't intersect sphere
        joinPoint = closestHeadPart.getPoint() - dRad*vectorToCentre.normalized();
        dDistance = dDistCentreToLine - dRad;
    }
    
    const double dDistJoinPointToSphere = (joinPoint - closestHeadPart.getPoint()).norm();
    CHECK_P(!zero(dDistJoinPointToSphere-dRad), dDistJoinPointToSphere, "Error computing closestPointOnSphere");
}


optional<const C2dLine> C3dLine::projectUnbounded(const CWorldCamera & P) const
{
    const C3dWorldPoint point2 = 0.01*direction+pointOnLine;
    const C3dBoundedLine bounded(pointOnLine,point2);
    return bounded.projectUnbounded(P);
}

C3dWorldPoint C3dLine::closestPointToLine(const C3dLine & other) const // from http://objectmix.com/graphics/133793-coordinates-closest-points-pair-skew-lines.html
{
    const C3dWorldPoint P21 = other.pointOnLine - pointOnLine;
    const C3dWorldPoint M = other.direction.cross(direction);
    const double m2 = M.squaredNorm();
    if(m2==0) {
        cout << "Parallel lines" << endl;
        return pointOnLine;
    }

    C3dWorldPoint R = P21.cross(M) / m2;

    double t1 = R.dot(other.direction);
    C3dWorldPoint Q1 = point(t1);

    /*if (IS_DEBUG) {
        double t2 = R.dot(direction);
        C3dWorldPoint Q2 = other.point(t2);

        double linelineDist = (Q2 - Q1).norm();
        if(IS_DEBUG) CHECK(!zero(linelineDist - closestDistance(other)), "Closest point calculation failed to find 2 points which are closestDist apart");
    }*/

    return Q1;
}

double C3dLine::closestDistanceToLine(const C3dLine & other) const
{
    return other.closestDistance(closestPointToLine(other));
}

C3dWorldPoint C3dLine::closestPointToBoundedLine(const C3dBoundedLine & boundedLine) const
{
    //return closestPointToBoundedLine_int<C3dBoundedLine>(boundedLine);
    const C3dWorldPoint closestUnbounded = closestPointToLine(boundedLine.unboundedLine());
    
    C3dWorldPoint closestPointOnBounded = boundedLine.closestPoint(closestUnbounded);

    return closestPoint(closestPointOnBounded);
}

double C3dLine::closestDistanceToBoundedLine(const C3dBoundedLine & boundedLine) const
{
    return boundedLine.closestDistance(closestPointToBoundedLine(boundedLine));
}

///////// C2dLine /////////////

void C2dLine::draw(cv::Mat & M, const cv::Scalar colour, const int nThickness) const
{
    CHECK_MAT_INIT(M);
    C2dImagePointPx p1 = point(-1000);
    C2dImagePointPx p2 = point(1000);
    cv::line(M, p1, p2, colour, nThickness);
}

C2dImagePointPx C2dLine::closestPointToBoundedLine(const C2dBoundedLine & boundedLine) const
{
    //return closestPointToBoundedLine_int<C2dBoundedLine>(boundedLine);
    
    optional<const C2dImagePointPx> intersection = findIntersection(boundedLine.unboundedLine());
    if(intersection)
    {
        //If intersection is on the bounded line, return intersection, else return whichever point on this line is closest to the bounded line endpoint
        C2dImagePointPx closestPointOnBounded = boundedLine.closestPoint(*intersection);

        return closestPoint(closestPointOnBounded);
    }
    else
    {
        return getPointOnLine();
    }    
}
double C2dLine::closestDistanceToBoundedLine(const C2dBoundedLine & boundedLine) const
{
    return boundedLine.closestDistance(closestPointToBoundedLine(boundedLine));
    //return closestDistanceToBoundedLine_int<C2dBoundedLine>(boundedLine);
}

///////// C3dBoundedLine //////////////

C3dWorldPoint C3dBoundedLine::closestPointToBoundedLine(const C3dBoundedLine & boundedLine) const
{
    return closestPoint(C3dLine(unboundedLine()).closestPointToBoundedLine(boundedLine));
}
C3dWorldPoint C3dBoundedLine::closestPointToLine(const C3dLine & otherLine) const
{
    const C3dWorldPoint closestOnOther = otherLine.closestPointToBoundedLine(*this);
    return closestPoint(closestOnOther);
}

double C3dBoundedLine::closestDistanceToBoundedLine(const C3dBoundedLine & boundedLine) const
{
    return boundedLine.closestDistance(closestPointToBoundedLine(boundedLine));
}

double C3dBoundedLine::closestDistanceToLine(const C3dLine & unboundedLine) const
{
    return unboundedLine.closestDistance(closestPointToLine(unboundedLine));
}

void C3dBoundedLine::draw(cv::Mat & M, const CWorldCamera & P, const cv::Scalar & col, const int nThickness) const
{
    CHECK_MAT_INIT(M);
    const optional<const C2dImagePointPx> p0 = P.projectToPx(getStartPoint());
    const optional<const C2dImagePointPx> p1 = P.projectToPx(getFinishPoint());

    if(p0 && p1)
        cv::line(M, *p0, *p1, col, nThickness);
}

optional<const C2dBoundedLine> C3dBoundedLine::projectBounded(const CWorldCamera & P) const
{
    const optional<const C2dImagePointPx> p1 = P.projectToPx(getStartPoint());
    const optional<const C2dImagePointPx> p2 = P.projectToPx(getFinishPoint());
    if(!p1 || !p2)
        return optional<const C2dBoundedLine>();

    return optional<const C2dBoundedLine>(C2dBoundedLine(*p1,*p2));
}

const C2dBoundedLine C3dBoundedLine::fastProjectBounded(const CWorldCamera & P) const
{
    return C2dBoundedLine(P.fastProject(getStartPoint()), P.fastProject(getFinishPoint()));
}

optional<const C2dLine> C3dBoundedLine::projectUnbounded(const CWorldCamera & P) const
{
    optional<const C2dBoundedLine> bounded = projectBounded(P);
    if(bounded)
        return optional<const C2dLine>(bounded->unboundedLine());
    else
        return optional<const C2dLine>();
}

bool C3dBoundedLine::planeIntersect(const C3dWorldPoint & n, const double d, C3dWorldPoint & intersection) const
{
    double t = (d-n.dot(getStartPoint()))/n.dot(getFinishPoint()-getStartPoint());

    if(t>1 || t<0)
        return false;

    intersection = (getStartPoint() + t*(getFinishPoint()-getStartPoint())).eval();

    //cout << intersection.transpose() << " intersects plane at t=" << t << endl;

    if(IS_DEBUG) CHECK(!zero(intersection.dot(n) - d), "plane intersect computation failed");

    return true;
}


C3dBoundedLine::C3dBoundedLine(const C3dWorldPoint & p1, const C3dWorldPoint & p2) : T3dBoundedLine_base(p1, p2)
{
    if(IS_DEBUG)
    {
        p1.checkInit();
        p2.checkInit();
        CHECK(getStartPoint() == getFinishPoint() && getStartPoint() != C3dWorldPoint::uninit(), "Line end points are equal");
    }
}

C3dLine C3dBoundedLine::unboundedLine() const { 
    C3dLine lineUnbounded; 
    lineUnbounded.from2Points(getStartPoint(), getFinishPoint());
    return lineUnbounded;
}




///////// C2dBoundedLine //////////////

C2dImagePointPx C2dBoundedLine::closestPointToBoundedLine(const C2dBoundedLine & boundedLine) const
{
    return closestPoint(C2dLine(unboundedLine()).closestPointToBoundedLine(boundedLine));
}
C2dImagePointPx C2dBoundedLine::closestPointToLine(const C2dLine & otherLine) const
{
    const C2dImagePointPx intersect = *(unboundedLine().findIntersection(otherLine));
    return closestPoint(intersect);
}

double C2dBoundedLine::closestDistanceToBoundedLine(const C2dBoundedLine & boundedLine) const
{
    return boundedLine.closestDistance(closestPointToBoundedLine(boundedLine));
}

double C2dBoundedLine::closestDistanceToLine(const C2dLine & unboundedLine) const
{
    return unboundedLine.closestDistance(closestPointToLine(unboundedLine));
}

C2dBoundedLine::eLineSide C2dBoundedLine::side(const C2dImagePointPx & point) const
{
    C2dImagePointPx bottom, top;
    if(getStartPoint().y() > getFinishPoint().y()) {
        bottom=getStartPoint();
        top=getFinishPoint();
    } else {
        bottom=getFinishPoint();
        top=getStartPoint();
    }

    C2dImagePointPx pll = bottom-top;
    C2dImagePointPx perp(-pll.y(), pll.x());

    return (perp.dot(point-bottom) < 0) ? eLeftOfLine : eRightOfLine;
}

TEigen2dPoint C2dBoundedLine::perpendicular() const
{
    return ::perpendicular(direction());
}

double C2dBoundedLine::signedAngle(const C2dLine & line) const
{
    return unboundedLine().signedAngle(line);
}

double C2dBoundedLine::signedAngle(const C2dBoundedLine & line) const
{
    return signedAngle(line.unboundedLine());
}

C2dLine C2dBoundedLine::unboundedLine() const { 
    C2dLine lineUnbounded; 
    lineUnbounded.from2Points(getStartPoint(), getFinishPoint());
    return lineUnbounded;
}

optional<const C2dImagePointPx> C2dBoundedLine::intersection(const C2dBoundedLine & line) const
{
    optional<const C2dImagePointPx> pointOfIntersection = unboundedLine().findIntersection(line.unboundedLine());
    if(!pointOfIntersection)
        return pointOfIntersection;

    if(zero(closestDistance(*pointOfIntersection)) && zero(line.closestDistance(*pointOfIntersection)))
        return pointOfIntersection;
    else
        return optional<const C2dImagePointPx>();
}

optional<const C2dImagePointPx> C2dBoundedLine::intersection(const C2dLine & line) const
{
    optional<const C2dImagePointPx> pointOfIntersection = unboundedLine().findIntersection(line);

    if(pointOfIntersection && zero(closestDistance(*pointOfIntersection))) //TODO: What about pll lines that do intersect? Prob should just check it this ever happens.
        return optional<const C2dImagePointPx>(pointOfIntersection);
    else
        return optional<const C2dImagePointPx>();
}

/*bool C2dBoundedLine::intersectsAtMidpoint(const C2dBoundedLine & otherLine, const double dTrimThresh) const
{
    const optional<T2dBoundedLine_base> trim1 = trim(dTrimThresh);
    const optional<T2dBoundedLine_base> trim2 = otherLine.trim(dTrimThresh);
    
    if(trim1 && trim2)
        return (bool)C2dBoundedLine(*trim1).intersection(C2dBoundedLine(*trim2));    
    else
        return false;
}*/

void C2dBoundedLine::draw(cv::Mat & M, const int nThickness, const cv::Scalar colour) const
{
    //const bool bVerbose = false;
    CHECK_MAT_INIT(M);
    const C2dImagePointPx p1 = getStartPoint();
    const C2dImagePointPx p2 = getFinishPoint();
    /*if(!inIm(p1, M, -2000) || !inIm(p1, M, -2000))
    {
        COUT2("Warning: not drawing 2D line " , *this);
        COUT(colour);
        return;
    }*/
    cv::line(M, p1, p2, colour, nThickness);
}


std::ostream& operator<<(std::ostream& s, const C3dLine & X) {
    s << "Point: " << X.getPointOnLine().transpose() << ", dir: " << X.getDirection().transpose() << " ";
    return s;
}
std::ostream& operator<<(std::ostream& s, const C2dLine & X) {
    s << "Point: " << X.getPointOnLine().transpose() << ", dir: " << X.getDirection().transpose() << " ";
    return s;
}
std::ostream& operator<<(std::ostream& s, const C3dBoundedLine & X) {
    s << "Point 1: " << X.getStartPoint().transpose() << ", Point 2: " << X.getFinishPoint().transpose() << " ";
    return s;
}
std::ostream& operator<<(std::ostream& s, const C2dBoundedLine & X) {
    s << "Point 1: " << X.getStartPoint().transpose() << ", Point 2: " << X.getFinishPoint().transpose() << " ";
    return s;
}


template class CBoundedLine_base<C2dImagePointPx>;// = T2dBoundedLine_base;
template class CBoundedLine_base<C3dWorldPoint>; // = T3dBoundedLine_base;
template class CUnboundedLine_base<C2dImagePointPx>;// = T2dUnboundedLine_base;
template class CUnboundedLine_base<C3dWorldPoint>;  // = T3dUnboundedLine_base;
