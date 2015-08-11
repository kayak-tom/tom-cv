#ifndef POLYLINE_H
#define POLYLINE_H

#include <util/exception.h>
#include <util/convert.h>

/**
 * @class CRange
 * @brief Represents a range of allowed values for a parameter.
 */
class CRange
{
    static const int UNINIT = -1000000;
    double dMin, dMax;
public:
    CRange() : dMin(UNINIT), dMax(UNINIT) {}

    CRange(const double dMin, const double dMax) : dMin(dMin), dMax(dMax) {
       CHECK_P(dMin > dMax, this, "Bad range");//Skip check when running online
    }
    
    static CRange makeRangeFromUnordered(const double dMin, const double dMax)
    {
        return (dMin < dMax) ? CRange(dMin, dMax) : CRange(dMax, dMin);
    }
    
    bool isInit() const {
        return dMin != (double)UNINIT && dMax != (double)UNINIT;
    }

    bool within_absolute(const double x, const double dMargin /* dMargin > 0 increases the isSensible accepted volume */) const {
        return x>=(dMin-dMargin) && x<=(dMax+dMargin);
    }

    bool within_relative(const double x, const double dRelMargin /* dMargin > 0 increases the isSensible accepted volume */) const {
        CHECK(dRelMargin < 0, "dRelMargin must be +ve");
        return x>=(1.0-dRelMargin)*dMin && x<=(1.0+dRelMargin)*dMax;
    }

    double clip(const double x, const double dMargin = 0) const {
        return ::clip<double>(x, dMin-dMargin, dMax+dMargin);
    }

    double getMin() const {
        return dMin;
    }
    double getMiddle() const {
        return 0.5*(dMin+dMax);
    }
    double getMax() const {
        return dMax;
    }
    
    double getRange() const { return getMax() - getMin(); }
    
    static CRange intersect(const CRange & R1, const CRange & R2);
    static CRange makeUnion(const CRange & R1, const CRange & R2);
};

std::ostream & operator<<(std::ostream & s, const CRange & range);

/* Thickness-specific functionality is implemented for all polylines, but will throw exceptions for C2dPolyline and C3dPolyline
 *
 *
 *
 *
 *                                                              CPolyLine<TControlPoint>
 *                                                           /                            \
 *                                                          /                              \
 *                                                         /                                \
 *                                                        /                                  \
 *                                     C2dPolyline_base                                         C3dPolyline_base
 *                           [2D, with or without thickness]                                     [3D, with or without thickness]
 *                         /                          \                                            /                          \
 *                        /                           \                                           /                           \
 *           C2dPolyline                         C2dPolylineWithThickness                 C3dPolyline                         C3dPolylineWithThickness
 *          [no thickness]                     [with thickness]                     [no thickness]                           [with thickness]
 *
 * 
 * Directions ARE important.
 * 
 * Each polyline is a wrapper around a vector of control points aControlPoints
 * 
 * Each control point has a point. polylines 'with thickness' have a thickness as well.
 * 
 * Each polyline has 2 *endpoints* labelled *start* (aControlPoints.begin) and *finish* (aControlPoints.back) These match *direction*, *startToFinish*, the direction of 'segment', and match all the equivalent concepts for bounded lines.
 * 
 */
#include "newCamera.h"

class CShortSection;

static const int DEFAULT_THICKNESS = -1,  //Constant for drawing thick polylines with their correct thickness, thin polylines with thickness 1
                DRAW_OUTLINE=-2;

template<class TControlPoint_in>
class CPolyline_base
{
public:
    typedef TControlPoint_in TControlPoint;
    typedef typename TControlPoint::TVecType TVecType;
    typedef typename TControlPoint::TLineType TLineType;
    typedef CPolyline_base<TControlPoint> TPolyline;
protected:
    typedef std::vector<TControlPoint, Eigen::aligned_allocator<TControlPoint> > TPolylineVector;
public:
    typedef typename TPolylineVector::const_iterator const_iterator;
    typedef typename TPolylineVector::iterator iterator;
protected:
    TPolylineVector aControlPoints;
    bool bClosed;
public:
    CPolyline_base(const bool bClosed=false) : bClosed(bClosed) {}
    
    /* size==numPoints()
     * std::size_t size() const { 
		return aControlPoints.size();
	}*/

    const TControlPoint & front() const { return aControlPoints.front(); }
    TControlPoint & front() { return aControlPoints.front(); }
    const TControlPoint & back() const { return aControlPoints.back(); }
    TControlPoint & back() { return aControlPoints.back(); }
    
    int numPoints() const {
        return (int)aControlPoints.size();
    }

    const_iterator begin() const { return aControlPoints.begin(); }
    const_iterator end() const { return aControlPoints.end(); }
    iterator begin() { return aControlPoints.begin(); } //Needed to make BOOST_FOREACH work properly?
    iterator end() { return aControlPoints.end(); }
    void push_back(const TControlPoint & p);
    void clearAndReserve(const int nReserveSize);
    void resize(const int nNewSize);
    void pop_front();
    void pop_back();
    void pop_end(const eEndpoints end);
    void erase(iterator start, iterator finish );
    void insert(iterator pos, const_iterator start, const_iterator finish );
    void push_front(const TControlPoint & controlPoint );

    void push_end(const TControlPoint & p, const eEndpoints end) ;

    void concatenate(const TPolyline & other);
    void removeDuplicates(); //remove exact duplicates in-place
    
    TControlPoint & operator[] (const int nIdx) //'Ref' is important as g2o optimiser refers back to polyline points.
    { 
        if(IS_DEBUG && !bClosed) CHECKOOB(nIdx, aControlPoints.size());
        return aControlPoints[nIdx % numPoints()]; 
    }    

    const TControlPoint & operator[] (const int nIdx) const
    { 
        if(!bClosed) /*if(IS_DEBUG) restore for a while*/ CHECKOOB(nIdx, aControlPoints.size());
        return aControlPoints[nIdx% numPoints()];
    }
    
    //For conversion to/from vectors in LM optimisation
    int dimension() const { return TControlPoint::PARAMS_PER_POLY_CONTROL_POINT * numPoints(); }
    
    //Convert an LM parameter index to a control point index
    int indexToCPIndex(const int nLMVecIndex) const { return (nLMVecIndex>=0) ? (nLMVecIndex/TControlPoint::PARAMS_PER_POLY_CONTROL_POINT) : -1; }
    
    TLineType approxAsLine() const;
    
    /**
     * @brief Distance from centre to p. May be incorrect for very curved polylines.
     * @param p
     * @return
     */
    double distanceToPoint(const TVecType & p, const bool bIncludeThickness) const;
    
    const TVecType closestPointOnSurface(const TVecType & p) const;
    const TVecType closestPoint(const TVecType & p) const;
    const TVecType closestPointAndPosition(const TVecType & p, double & dPosition_index) const;
    const TVecType closestPointToBoundedLine(const TLineType & boundedline) const;
    const double closestCollisionDistanceToPoly(const TPolyline & otherPoly) const;
    /**
     * @brief Closest point; also return idx of the closest segment, and the distance (squared) to the closest point.
     * @param p
     * @return
     */
    const TControlPoint closestPoint_segIdx_distSq(const TVecType & p, int & nClosestSeg, double & dDistToPoly_sq, double * pdPosition_index = 0) const;
    
    //For returning control points: interpolate width at t if needed
    const TControlPoint interpolateControlPoint(const typename CPolyline_base<TControlPoint>::TVecType & point, const int nClosestSeg, const double t_closest) const;
    /**
     * @brief If the closest point is off the end then return nothing.
     * @param p
     * @return
     */
    const optional<const TVecType> closestPointStrictlyOnPoly(const TVecType & p) const;
    const TLineType closestSegment(const TVecType & p) const;
    double tanKinkAngleAtPoint(const int nPoint) const;
    double kinkAngleAtPoint(const int nPoint) const;
    TVecType getMaxKinkPoint() const;
    
    TVecType directionAtPoint(const int nPoint) const;
    
    double maxDistanceFromPoly(const TPolyline & other) const;

    double minLength() const;
    
    double getDPEpsilon() const;
    
    /**
     * @brief Closest midpoint and interpolated width (a control point with no width if polyline doesn't have it)
     * @param p
     * @return
     */
    const TControlPoint closestPointAndWidth(const TVecType & p, double * pdPosition_idx=0) const;

    //dPosition is strictly in the range 0...length(). Return point at that position.
    const TVecType pointAtDistance(const double dDistanceAlongPoly) const;
    
    //dPosition is strictly in the range 0...length(). Return point at that position.
    const TLineType segmentAtDistance(const double dDistanceAlongPoly) const;
    
    //dPositionAlongPoly is strictly in the range 0...numPoints() (even off the end), and integer values match corresponding control point positions.
    const TVecType pointAtPosition(const double dPositionAlongPoly) const;
    const TControlPoint controlPointAtPosition(const double dPosition) const;
    
    const int closestSegmentIdx(const TVecType & p) const;
    
    double thicknessAtPosition_index(const double dPosition) const;
    
private:
    int closestControlPoint(const TVecType & p) const;
public:
    CPolyline_base(const typename TPolylineVector::const_iterator & begin, const typename TPolylineVector::const_iterator & end) : aControlPoints(begin, end), bClosed(false) {}

    /**
     * @brief Position along polyline of p. Integer positions correspond to control points, and the fractional part corresponds to the relative position along the closest segment. NOT clipped to 0...numPoints
     * @param p
     * @return Position along polyline of p.
     */
    double positionOfPoint_index(const TVecType & p) const;
    
    
    /**
     * @brief Polyline is parameterised 0...length. Points off the polyline are *not* clipped, so can take any value.
     * @param p
     * @return Length (metres) along polyline of p
     */
    double positionOfPoint_distance(const TVecType & p) const;

    double indexToDistance(const double dIndex) const;
    CRange indexRangeToDistanceRange(const CRange & range_index) const;
    
    CRange getIndexRange() const;
    CRange getDistanceRange() const;
    
    const TVecType & getStartPoint() const {
        CHECK(numPoints() == 0, "Empty polyline");
        CHECK(bClosed, "Closed polylines don't really have 'ends'");
        return aControlPoints[0].getPoint();
    }
    const TVecType & getFinishPoint() const {
        CHECK(numPoints() == 0, "Empty polyline");
        CHECK(bClosed, "Closed polylines don't really have 'ends'");
        return aControlPoints[numPoints() - 1].getPoint();
    }
    //Allow iteration over endpoints
    const TVecType & getEndPoint(const eEndpoints end) const;

    //The index of the control point at end
    int getEndIdx(const eEndpoints end) const {
        CHECK(bClosed, "Closed polylines don't really have 'ends'");
        return (end == eStart) ? 0 : (numPoints() - 1);
    }
    
    //Allow iteration over endpoints
    const TControlPoint & getEndControlPoint(const eEndpoints end) const;

    const TVecType midPoint() const;
    
    const TVecType endDirection(const eEndpoints end) const;

    std::string toString() const;
    
    TVecType direction() const;
    void reverseDirection();

    const int numSegments() const {
        return bClosed ? numPoints() : (numPoints() - 1);
    }
    TLineType segment(const int nSegStart) const;
    TLineType endSegment(const eEndpoints end) const;

    double getWidth(const TVecType & x) const;

    double length() const;

    /**
     * @brief Test if the point p is a midpoint, i.e. not a tip
     * @param p
     * @return
     */
    bool isMidPoint(const TVecType & p, const double dREThresh) const;

    void addLocalisationErrors(const double dGaussianNoiseSD);

    /**
     * @brief remove super-short sections, colinear points
     * @return
     */
    optional<TPolyline> truncate(const double dMinLength, const double dMaxLength, const double dEpsilon, const bool bVerbose) const;
    
    optional<TPolyline> truncate_defaults(const bool bVerbose) const
    {
        return truncate(minLength(), HUGE, getDPEpsilon(), bVerbose);
    }

    static TPolyline douglasPeucker(const TPolyline & polyLine, const double dMaxLength, const double dEpsilon, const bool bVerbose);
    TPolyline removeShortSections(const double dMinLength, const bool bVerbose) const;
    TPolyline splitLongSections(const double dMaxLength, const bool bVerbose) const;
    
    void splitAtPoint(const TVecType & p, TPolyline * aPolylines) const;
    void subPolyline(const CRange & indexRange, TPolyline & subPoly) const;

    //Find a polyline close to this one with maxKinkAngle below dNewKinkAngle;
    void dropKinkAngle(const double dNewKinkAngle);

    //Move up a pyramid level: points *= 2; extra points inserted.
    void doubleUp();
    
    void transform(const Eigen::Matrix<double, TControlPoint_in::TVecType::RowsAtCompileTime + 1, TControlPoint_in::TVecType::RowsAtCompileTime + 1> & T);
private:
    /**
     * @brief Douglas-Peucker minimum epsilon
     * @return TODO fix truncation etc.
     */
    static double DPMinEpsilon() {
        return TVecType::RowsAtCompileTime == 2 ? 1 : 0.01;
    }
    static double DPMaxEpsilon() {
        return TVecType::RowsAtCompileTime == 2 ? 5 : 0.05;
    }
    double maxKinkAngle_int(int & nMaxKinkAnglePos) const;
    
public:
    
    /**
     * @brief Average thickness in pixels
     * @return
     */
    double averageThickness() const;

    bool tooShort() const; //a threshold used in a few places for getting rid of very short polylines--e.g. shorter than 1.5*thickness
    
    
    /**
     * @brief Check polyline isn't too kinked.
     */
    void checkNotTooKinked(const double dMaxKinkAngle = M_PI) const;

    void splitSegment(const int nSegmentToSplit);
    
    double maxKinkAngle() const;
    
    
    /**
     * @brief bad assignments cause the polyline to zigzag
     * @return Maximum difference between the polyline direction and the average direction of the two neighbouring segments. If smoothly curved it should be 0. A bit sensitive to segment lengths...
     */
    double maxDislocationAngle() const;
    double dislocationAngleAtSegment(const int nPoint) const;

    void toVector(Eigen::VectorXd & x) const; // For LM optimisation
    eSuccessStatus fromVector(const Eigen::VectorXd & x); // For LM optimisation
    eSuccessStatus onePointFromVector(const Eigen::VectorXd & x, const int nPointIndex); // For LM optimisation
    
    
    /**
     * @brief Used for checking how much polylines overlap (for checking test reconstructions match the original), checking recon. polylines aren't the same as existing polylines, etc.
     * 
     * Function will tend to overestimate the overlap a little, especially for kinked/zig-zag polylines
     * 
     * @param otherPoly
     * @return Proportion of otherPoly's area which overlaps with this poly. 1 if this contains the other poly (i.e. its control points), 0 if they don't intersect
     */
    double contains(const TPolyline & otherPoly, const bool bSlowAndAccurate=false) const;
    
    double area() const;
    double volume() const;
};


//If max > 2*min then we can always find a path through the polyline, *provided we subdivide the polyline sufficiently*
class CMinMaxApproxSettings
{
    double dMin, dMax;
public:
    CMinMaxApproxSettings(const double dMin, const double dMax) : dMin(dMin), dMax(dMax) 
    {
        CHECK(dMin < 0, "dMin too short");
        CHECK(dMax <= 0, "dMax too short");
        CHECK(dMax < 2*dMin, "dMax too short"); 
    }
    
    double getMin() const { return dMin; }
    double getMax() const { return dMax; }

    bool tooShort(const double dDistFromPrevCP) const { return dDistFromPrevCP < dMin; }
    bool tooLong(const double dDistFromPrevCP) const { return dDistFromPrevCP > dMax; }
};

class CPolyApproxSettings
{
    friend std::ostream & operator<<(std::ostream & s, const CPolyApproxSettings & settings);
protected:
    CMinMaxApproxSettings minMaxApprox;
    double dTargetSegLength, //Penalty 0 for segments this long.
                 dEpsilon; //A penalty per CP in proportion to the Dougles Peucker error (max distance). Normalised by target segment length for area error
public:
    CPolyApproxSettings(const double dMin, const double dMax_in, const double dEpsilon=0, const double dTargetSegLength_in=-1 /* default = max */);
    
    const double penalty(const double dDistFromPrevCP) const;

    const double getEpsilon() const { return dEpsilon; }
    const double getTargetLength() const { return dTargetSegLength; }
    const CMinMaxApproxSettings & getMinMax() const { return minMaxApprox; }
};

std::ostream & operator<<(std::ostream & s, const CPolyApproxSettings & settings);

class C2dPolyApproxSettings : public CPolyApproxSettings
{
    const double dBudScale; //buds smaller than this are removed first.
    const bool bMaxDistanceCost; //Douglas-Peucker-like max distance vs. area cost
public:
    C2dPolyApproxSettings(const double dMin, const double dMax_in, const double dEpsilon=0, const double dTargetSegLength_in=-1, const double dBudScale=-1, const bool bMaxDistanceCost_in=true)
        : CPolyApproxSettings(dMin, dMax_in, dEpsilon, dTargetSegLength_in), dBudScale(dBudScale), bMaxDistanceCost(bMaxDistanceCost_in)
    {
        
    }
    
    double getBudScale() const { return dBudScale; }
    bool useMaxDistanceCost() const { return bMaxDistanceCost; }
    const double penalty(const double dDistFromPrevCP) const 
    {
        const double dPenalty = CPolyApproxSettings::penalty(dDistFromPrevCP);
        if(useMaxDistanceCost())
            return dPenalty;
        else
            return dTargetSegLength*dPenalty;
    }
};

/* Implementation of functions only relevent for 2d polylines
 * TControlPoint = CPolylineControlPoint<2d> or CPolylineControlPointWithThickness<2d>
 * */
template<class TControlPoint>
class C2dPolyline_base : public CPolyline_base<TControlPoint>   //=this::TPolyline
{
    void drawDetails(cv::Mat & M) const;
    void drawLine(cv::Mat & M, const cv::Scalar colour, const int nThickness) const;
    void drawOutline(cv::Mat & M, const cv::Scalar col, const int nOutlineThickness) const;
public:
    typedef C2dPolyline_base<TControlPoint> T2dPolyline;
    typedef CPolyline_base<TControlPoint> TPolyline_base;
    typedef std::vector<C2dPolyline_base<TControlPoint> > TCandidateApproximations;
    typedef C2dPolyApproxSettings TPolyApproxSettings;

    C2dPolyline_base(const bool bClosed=false) : TPolyline_base(bClosed) {}
    C2dPolyline_base(const TPolyline_base & base) : TPolyline_base(base) {}

    /**
     * @brief Finds intersection between epiline and polygon, starting from section startIndex
     * Multiple intersections are resolved by finding the first index equal to or after startIndex. This only works if the 2D polylines are traversed in the same direction, and are not too kinked.
     * 
     * BETTER: use pointForRecon
     * 
     * @param epiline
     * @return
     */
    optional<const C2dImagePointPx> pointOfIntersection(const C2dBoundedLine & epiline) const;

    optional<const C2dImagePointPx> pointForRecon(const C2dBoundedLine & epiline, const double dMinAngleForWellConditionedRecon) const;
    
    /**
     * @brief *Signed* distance from centre to p. May be incorrect for very curved polylines.
     * @param p
     * @return
     */
    double signedDistanceToPoint(const C2dImagePointPx & p) const;
    
    void sortByPolyline(const T2dPolyline & polyToSortWRT);
    
    /**
     * @brief Used as a metric for deciding whether to join cane segments
     *
     * Units: radians per pixel (signed)
     *
     * @return
     */
    double getCurvature() const;

    void dropKinkAngle_signed(const double dNewKinkAngle, const bool bFromBelow);

    /**
     * @brief Draw 2D polyline
     * @param M
     * @param colour
     * @param bDrawDetails Add control point labels
     * @param nThickness DEFAULT_THICKNESS draws with control point thickness (1 if control point doesn't have thickness), otherwise specify a thickness. DRAW_OUTLINE draws the outline
     */ 
    void draw(cv::Mat & M, const cv::Scalar colour, const bool bDrawDetails, const int nThickness=DEFAULT_THICKNESS) const;    
    //template compatible draw-and-show method
    //void show(std::string label) const;

    eSuccessStatus approximate(C2dPolyline_base & polyline_out, const C2dPolyApproxSettings & polyApproxSettings, const bool bVerbose) const;
    
    //New truncation
    enum ePolyTruncationJoinType { eNoJoin, eReplaceWithIntersection, eReplaceWithOneEdge };
    ePolyTruncationJoinType isBud(const int nBudSegmentEnd, int & nBudEnd, typename TPolyline_base::TVecType & interpControlPoint, const double dMinSegmentLength, const double dBudScale, const bool bVerbose) const;

    double maxKinkAngle_signed(int& nMaxKinkAnglePos, const bool bFromBelow) const;
    double signedKinkAngleAtPoint(const int nPoint) const;

protected:
    //For each chain which is in total shorter than dMinLength, either deletes a point or interpolates a join as appropriate
    T2dPolyline removeShortSections_new(const double dMinLength, const bool bVerbose) const;
    T2dPolyline decideHowToTruncateShortSection(const CShortSection & shortSection, const double dMinLength, const bool bVerbose) const;
    double measureApproximationQuality(const CShortSection & shortSection, T2dPolyline & candidateApproximation, const double dMinLength, const bool bVerbose) const;
    
public:    
    //Assumes this polyline is not self-intersecting
    double polygonalArea_nonSelfIntersecting() const;
};

typedef C2dPolyline_base<C2dPolylineControlPoint > T2dPolyline;
class C2dPolyline : public T2dPolyline
{
public:
    C2dPolyline(const bool bClosed=false) : T2dPolyline(bClosed){}
    C2dPolyline(const T2dPolyline & base) : T2dPolyline(base) {} //{ T2dPolyline::operator=(base); }
    C2dPolyline(const T2dPolyline::TPolyline_base & base) : T2dPolyline(base) {} //{ T2dPolyline::TPolyline_base::operator=(base); }
    
    //C2dPolyline approximate(const double dMinLength, const double dMaxLength, const double dEpsilon, const double dBudScale, const eApproxMethod approxMethod, const bool bVerbose) const;
};

class C3dPolylineWithThickness;

typedef C2dPolyline_base<C2dPolylineControlPointWithThickness > T2dPolylineWithThickness;

class C2dPolylineWithThickness : public T2dPolylineWithThickness
{
public:
    C2dPolylineWithThickness(const bool bClosed=false) : T2dPolyline(bClosed){}
    
    C2dPolylineWithThickness(const T2dPolylineWithThickness & base) : T2dPolylineWithThickness(base) {} //{ T2dPolylineWithThickness::operator=(base); }
    //C2dPolylineWithThickness(const T2dPolylineWithThickness::TPolyline_base & base)  { T2dPolylineWithThickness::TPolyline_base::operator=(base); }
    //void operator=(const T2dPolylineWithThickness & base) { T2dPolylineWithThickness::operator=(base); }

    optional<const C3dPolylineWithThickness> polyToWorld(const CWorldCamera & P, const double dDepthPrior) const;
};

/* Declare functions only relevent for 3d polylines
 * TControlPoint = CPolylineControlPoint<3d> or CPolylineControlPointWithThickness<3d>
 * */
template<class TControlPoint>
class C3dPolyline_base : public CPolyline_base<TControlPoint>   //=this::TPolyline
{
public:
    typedef CPolyline_base<TControlPoint> TPolyline; 
    typedef CPolyApproxSettings TPolyApproxSettings;
protected:
    void drawDetails(cv::Mat & M, const CWorldCamera & P) const;
    
    template<class TProjectToPolyType>
    void drawLine(cv::Mat & M, const CWorldCamera & P, const cv::Scalar colour, const int nThickness) const;
    
    eSuccessStatus projectOne(const CWorldCamera & P, C2dPolyline & thinPoly, const int nControlPointIdx, const bool bVerboseOnFailure) const;
    eSuccessStatus projectOne(const CWorldCamera & P, C2dPolylineWithThickness & thinPoly, const int nControlPointIdx, const bool bVerboseOnFailure) const;

    void projectOne_fast(const CWorldCamera & P, C2dPolyline & thinPoly, const int nControlPointIdx) const;
    //void projectOne_fast(const CWorldCamera & P, C2dPolylineWithThickness & thinPoly, const int nControlPointIdx) const;
public:
    C3dPolyline_base(const bool bClosed=false) : TPolyline(bClosed){}
    C3dPolyline_base(const TPolyline & base) : TPolyline(base) {} //{ TPolyline::operator=(base); }

    eSuccessStatus project(const CWorldCamera & P, C2dPolyline & thinPoly, const bool bVerboseOnFailure=false, const int nControlPointIdx = -1 /*Project all points*/) const;
    eSuccessStatus project(const CWorldCamera & P, C2dPolylineWithThickness & poly, const bool bVerboseOnFailure=false, const int nControlPointIdx = -1 /*Project all points*/) const;

    void project_incomplete(const CWorldCamera & P, C2dPolyline & thinPoly) const;
    void project_incomplete(const CWorldCamera & P, C2dPolylineWithThickness & poly) const;
    
    bool isVisible(const CWorldCamera & P) const;
};

typedef C3dPolyline_base<C3dPolylineControlPoint> T3dPolyline;
class C3dPolyline : public T3dPolyline
{
public:
    C3dPolyline(const bool bClosed=false) : T3dPolyline(bClosed){}
    C3dPolyline(const T3dPolyline & base) : T3dPolyline(base) {} // { T3dPolyline::operator=(base); }
    //C3dPolyline(const T3dPolyline::TPolyline & base)  { T3dPolyline::TPolyline::operator=(base); }    
    /**
     * @brief Project and draw 3D polyline, should be compatible with C3dPolylineWithThickness::draw for template substitution.
     * @param M
     * @param P
     * @param colour
     * @param bDrawDetails Add control point labels
     * @param nThickness DEFAULT_THICKNESS draws with thickness 1, otherwise specify a thickness. 
     */ 
    void draw(cv::Mat & M, const CWorldCamera & P, const cv::Scalar colour, const bool bDrawDetails, const int nThickness=DEFAULT_THICKNESS) const;

    //template compatible draw-and-show method
    //void show(std::string label) const;

    eSuccessStatus approximate(C3dPolyline & polyline_out, const CPolyApproxSettings & polyApproxSettings, const bool bVerbose) const;
    
    eSuccessStatus project_mask(const CWorldCamera & P, C2dPolyline & thinPoly, std::vector<int> & aCorresponding3dPolyIndices) const;
};


typedef C3dPolyline_base<C3dPolylineControlPointWithThickness> T3dPolylineWithThickness;

class C3dPolylineWithThickness : public T3dPolylineWithThickness
{
public:
    C3dPolylineWithThickness(const bool bClosed=false) : T3dPolylineWithThickness(bClosed) {}
    C3dPolylineWithThickness(const T3dPolylineWithThickness & base) : T3dPolylineWithThickness(base) {} //{ T3dPolylineWithThickness::operator=(base); }
    //C3dPolylineWithThickness(const T3dPolylineWithThickness::TPolyline & base)  { T3dPolylineWithThickness::TPolyline::operator=(base); }
    
    /**
     * @brief Project and draw 3D polyline, should be compatible with C3dPolyline::draw for template substitution.
     * @param M
     * @param P
     * @param colour
     * @param bDrawDetails Add control point labels
     * @param nThickness DEFAULT_THICKNESS draws with projected thickness, otherwise specify a thickness. DRAW_OUTLINE draws the outline
     */ 
    void draw(cv::Mat & M, const CWorldCamera & P, const cv::Scalar colour, const bool bDrawDetails, const int nThickness=DEFAULT_THICKNESS) const;

    //template compatible draw-and-show method
    //void show(std::string label) const;

    eSuccessStatus approximate(C3dPolylineWithThickness & polyline_out, const CPolyApproxSettings & polyApproxSettings, const bool bVerbose) const;
    
    TVecType closestPointToLineOnSurface(const C3dBoundedLine & line) const;
    
    eSuccessStatus projectOneThickness(const CWorldCamera & P, C2dPolylineWithThickness & thinPoly, const int nControlPointIndex) const;
    eSuccessStatus projectThickness(const CWorldCamera & P, C2dPolylineWithThickness & thinPoly) const;
};

class CReconSVMFeature;

template<class TThickPoly, class TThinPoly>
void thickToThin(const TThickPoly & thick, TThinPoly & thin);

template<class TThinPoly, class TThickPoly>
void thinToThick(const TThinPoly & thin, TThickPoly & thick, const double dThickness);

/**
 * @brief Retu
 * @param poly1
 * @param poly2
  * @return 0 (no overlap) to 1 (equal polylines). Uses 'overlap' methods
*/
template<class TThickPoly>
double intersection(const TThickPoly & poly1, const TThickPoly & poly2, const bool bSlowAndAccurate);


std::ostream & operator<<(std::ostream& s, const C3dPolylineControlPointWithThickness & X);
std::ostream & operator<<(std::ostream& s, const C2dPolylineControlPointWithThickness & X);
std::ostream & operator<<(std::ostream& s, const C3dPolylineControlPoint & X);
std::ostream & operator<<(std::ostream& s, const C2dPolylineControlPoint & X);

std::ostream & operator<<(std::ostream& s, const C2dPolyline & X);
std::ostream & operator<<(std::ostream& s, const C2dPolylineWithThickness & X);
std::ostream & operator<<(std::ostream& s, const C3dPolyline & X);
std::ostream & operator<<(std::ostream& s, const C3dPolylineWithThickness & X);

template<class TVecType>
double tanAngleBetweenVectors(const TVecType & seg1vec, const TVecType & seg2vec);

    
//template-compatible method for debugging -- no frameId = latest frame (3D) or blank (2D)
/*template<typename TPolyType>
void show2dPolyline(const TPolyType & polyline, std::string label="", CFrameId frameId = CFrameId());

template<typename TPolyType>
void show3dPolyline(const TPolyType & polyline, std::string label="", CFrameId frameId = CFrameId());*/

//Used for measuring reconstruction and approximation quality.
template<typename TPolyline>
double polyMaxSeparation(const TPolyline & poly1, const TPolyline & poly2);

//Used for detecting possible assignments, etc.
template<typename TPolyline>
double polyMinSeparation(const TPolyline & poly1, const TPolyline & poly2);

#endif // POLYLINE_H
