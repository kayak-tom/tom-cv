/* Template classes for line segments
 * 
 * Implementations of these functions are in lines.ipp (you should normally #include lines.h and not lines.ipp)
 * 
 * Directions ARE important.
 * 
 * Each bounded line has 2 *endpoints* labelled *start* and *finish*. These match *direction*, *startToFinish* and match all the equivalent concepts for polylines.
 * 
 * Each bounded line is parameterised by t. t=0 is the start point, t=1 is the endpoint.
 * 
 * */

#ifndef LINES_H
#define LINES_H

 
class C2dBoundedLine;
class C2dLine;
class C3dBoundedLine;
class C3dLine;
class CWorldCamera;
class C3dPolylineControlPointWithThickness;

enum eSuccessStatus { eFail, eSuccess };

enum eEndpoints { eStart=0, eFinish=1, NUM_ENDPOINTS=2, UNINIT_ENDPOINT=-1 }; //Better way to allow enumeration over endpoints:  for(eEndpoints end=eStart; end < NUM_ENDPOINTS; end++)
inline eEndpoints otherEnd(const eEndpoints end) { return (end==eStart) ? eFinish : eStart; }
inline void operator++(eEndpoints & e, int)
{
    e=(eEndpoints)(((int)e)+1);
}

template <class TVecType>
class CUnboundedLine_base {
protected:
    typename TVecType::Base direction; //normalised line direction
    TVecType pointOnLine;

    typedef CUnboundedLine_base<TVecType> TUnboundedLine;
public:
    CUnboundedLine_base();
    
    //Specifically the line will be in the direction from p1 to p2
    void from2Points(const TVecType & p1, const TVecType & p2);

    void fromPointAndDir(const TVecType & p, const typename TVecType::Base & dir);

    //Returns a point on the line with param lambda
    TVecType point(const double & lambda) const;

    double undirectedAngle(const TUnboundedLine & other) const;
    
    //TVecType closestPoint(const TUnboundedLine & line) const;
    TVecType closestPoint(const TVecType & p) const;
    double closestDistance(const TVecType & p) const;
    //double closestDistance(const TUnboundedLine & line) const;
    
    //optional<const TVecType> findIntersection_robust(const TUnboundedLine & other) const; //TODO: rename me
    
    const typename TVecType::Base & getDirection() const { return direction; }
    const TVecType & getPointOnLine() const { return pointOnLine; }
    
/*protected:
    //Prevent instantiation problems by keeping these protected, then instantiating from derived class 
    template<class TBoundedLine>
    TVecType closestPointToBoundedLine_int(const TBoundedLine & boundedLine) const;
    
    template<class TBoundedLine>
    double closestDistanceToBoundedLine_int(const TBoundedLine & boundedLine) const;*/
};

typedef CUnboundedLine_base<C2dImagePointPx> T2dLine;

class C2dLine : public T2dLine {
public:
    C2dLine() {}
    explicit C2dLine(const T2dLine & base) { (T2dLine & )(*this) = base; }

    double signedAngle(const C2dLine & line) const;
    double signedAngle(const C2dBoundedLine & line) const;

    void draw(cv::Mat & M, const cv::Scalar colour, const int nThickness = 1) const;
    
    /*
       pointOnLine + lambda * direction = other.pointOnLine + mu * other.direction 
       Solve for lambda, mu
     */
    optional<const C2dImagePointPx> findIntersection(const C2dLine & other) const;
    
    C2dImagePointPx closestPointToBoundedLine(const C2dBoundedLine & boundedLine) const;
    double closestDistanceToBoundedLine(const C2dBoundedLine & boundedLine) const;
};

typedef CUnboundedLine_base<C3dWorldPoint> T3dLine;

class C3dLine : public T3dLine {
public:
    C3dLine() {}
    C3dLine(const C3dBoundedLine & boundedLine);
    explicit C3dLine(const T3dLine & base) { (T3dLine &)(*this) = base; }

    optional<const C2dLine> projectUnbounded(const CWorldCamera & P) const;
    /**
     * @brief closest point on this line to the other line. Pll not handled properly at the moment
     * @param line
     * @return 
     */
    C3dWorldPoint closestPointToLine(const C3dLine & line) const;
    double closestDistanceToLine(const C3dLine & line) const;
    
    C3dWorldPoint closestPointToBoundedLine(const C3dBoundedLine & boundedLine) const;
    double closestDistanceToBoundedLine(const C3dBoundedLine & boundedLine) const;
};


/////////////////////////////////////////////////////////// Bounded /////////////////////////////////////////

template <class TVecType_in>
class CBoundedLine_base
{
public:
    typedef TVecType_in TVecType;
protected:
    TVecType aEndPoints[NUM_ENDPOINTS];
    
    typedef CUnboundedLine_base<TVecType> TUnboundedLine;
    typedef CBoundedLine_base<TVecType> TBoundedLine;
public:
    CBoundedLine_base(const TVecType & p1, const TVecType & p2);
    CBoundedLine_base();

    /**
     * @brief Return a line trimmed by dTrimThresh off each end
     * @param dTrimThresh
     * @return 
     */
    optional< TBoundedLine> trim(const double dTrimThresh) const;
    
    double length() const;
    //double directionRads() const;

    TVecType closestPoint(const TVecType & point) const HOT;
    TVecType closestPointAndt(const TVecType & p, double & t, const bool bClamp_t) const HOT;
    double closestDistance(const TVecType & point) const HOT;

    //double closestDistanceRENAMEME(const C2dBoundedLine & otherLine) const HOT;

    double undirectedAngle(const TUnboundedLine & line) const;
    double undirectedAngle(const TBoundedLine & line) const;
    double directedAngle(const TUnboundedLine & line) const;
    double directedAngle(const TBoundedLine & line) const;

    //Unbounded line will be in the direction from p1 to p2
    TUnboundedLine unboundedLine() const;
    
    /**
     * @brief t from the closest-point-on-clipped-line equation, with t=0 being startPoint and t=1 being the endpoint
     * @param point
     * @return 
     */
    double get_t(const TVecType & point) const;
    //Turn t from get_t into a point on the line 
    TVecType pointFromt(const double t) const;
    
    const TVecType midPoint() const;
    const TVecType & getStartPoint() const { return aEndPoints[eStart]; }
    const TVecType & getFinishPoint() const { return aEndPoints[eFinish]; }
    const TVecType & getEndPoint(const eEndpoints whichEnd) const { return aEndPoints[whichEnd]; }
    
    /**
     * @brief Normalised perpendicular vector. Used for signed distances so direction shouldn't be changed.
     * @return 
     */
    //TVecType perpendicular() const;
    /**
     * @brief Normalised direction vector in start -> end direction.
     * @return 
     */
    typename TVecType::Base direction() const;
    /**
     * @brief vector in start -> end direction, NOT normalised
     * Parameterised by t
     * @return 
     */
    TVecType startToFinish() const { return getFinishPoint() - getStartPoint(); }
    
    template<class TLineType>
    TLineType reversed() const
    {
        return TLineType(getFinishPoint(), getStartPoint());
    }
};

typedef CBoundedLine_base<C2dImagePointPx> T2dBoundedLine_base;
    
class C2dBoundedLine : public T2dBoundedLine_base
{
public:
    C2dBoundedLine(const C2dImagePointPx & p1, const C2dImagePointPx & p2) : T2dBoundedLine_base(p1, p2) {}
    C2dBoundedLine() {}
    explicit C2dBoundedLine(const T2dBoundedLine_base & boundedLine) { (T2dBoundedLine_base &)(*this) = boundedLine; }
    void operator=(const T2dBoundedLine_base & boundedLine) { (T2dBoundedLine_base &)(*this) = boundedLine; }

    typedef T2dBoundedLine_base::TVecType TVecType;

    enum eLineSide {eLeftOfLine, eRightOfLine};
    eLineSide side(const C2dImagePointPx & point) const;
    void draw(cv::Mat & M, const int nThickness, const cv::Scalar colour) const;    
    
    double signedAngle(const C2dLine & line) const;
    double signedAngle(const C2dBoundedLine & line) const;
    
    //Normalised perpendicualr vector
    TEigen2dPoint perpendicular() const;    
    C2dLine unboundedLine() const;

    optional<const C2dImagePointPx> intersection(const C2dLine & line) const;
    optional<const C2dImagePointPx> intersection(const C2dBoundedLine & line) const;
    //bool intersectsAtMidpoint(const C2dBoundedLine & line, const double dTrimThresh) const;

    C2dImagePointPx closestPointToBoundedLine(const C2dBoundedLine & boundedLine) const;
    C2dImagePointPx closestPointToLine(const C2dLine & unboundedLine) const;
    double closestDistanceToBoundedLine(const C2dBoundedLine & boundedLine) const;
    double closestDistanceToLine(const C2dLine & unboundedLine) const;
    
    typedef C2dLine TUnboundedLine;
};

typedef CBoundedLine_base<C3dWorldPoint> T3dBoundedLine_base;

class C3dBoundedLine : public T3dBoundedLine_base {
public:
    C3dBoundedLine(const C3dWorldPoint & p1, const C3dWorldPoint & p2);
    C3dBoundedLine() {}
    explicit C3dBoundedLine(const T3dBoundedLine_base & boundedLine) { (T3dBoundedLine_base &)(*this) = boundedLine; }

    typedef T3dBoundedLine_base::TVecType TVecType;
    
    /**
     * @brief Find intersection x on line where n.x=d
     * @param n
     * @param d
     * @param intersection
     * @return true on success
     */
    bool planeIntersect(const C3dWorldPoint & n, const double d, C3dWorldPoint & intersection) const;
    
    void draw(cv::Mat & M, const CWorldCamera & P, const cv::Scalar & col, const int nThickness = 1) const;
    
    optional<const C2dBoundedLine> projectBounded(const CWorldCamera & P) const HOT;
    optional<const C2dLine> projectUnbounded(const CWorldCamera & P) const;
    const C2dBoundedLine fastProjectBounded(const CWorldCamera & P) const HOT;
    
    C3dLine unboundedLine() const;
    
    C3dWorldPoint closestPointToBoundedLine(const C3dBoundedLine & boundedLine) const;
    C3dWorldPoint closestPointToLine(const C3dLine & unboundedLine) const;
    double closestDistanceToBoundedLine(const C3dBoundedLine & boundedLine) const;
    double closestDistanceToLine(const C3dLine & unboundedLine) const;
//joinPoint is the point on the sphere closest to this line.
    void closestPointOnSphere(const C3dPolylineControlPointWithThickness & closestHeadPart, C3dWorldPoint & joinPoint, double & dDistance) const;

    typedef C3dLine TUnboundedLine;
};


std::ostream& operator<<(std::ostream& s, const C3dLine & X);
std::ostream& operator<<(std::ostream& s, const C2dLine & X);
std::ostream& operator<<(std::ostream& s, const C3dBoundedLine & X);
std::ostream& operator<<(std::ostream& s, const C2dBoundedLine & X);

//For clamping t to 0...1 (which finds the closest point on a bounded rather than unbounded line)
inline double segmentClamp(const double t)
{
    if (t<0)
        return 0;
    if (t>1)
        return 1;
    return t;
}


#endif
