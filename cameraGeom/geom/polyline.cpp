#include "polyline.h"
//#include "threeDPoint.h"
#include <util/random.h>
//#include "reconstruct3d.h"
#include "imageUtil.h"
#include <geom/levMarNumerical.h>
#include <Eigen/Geometry>
#include "lines.h"
#include <math.h>

#include "polyline.ipp"
#include "cameras.ipp"

#include "vector_hash.h"
#include "fixedPolylineOptimisation.h"
#include <boost/scoped_array.hpp>

enum ePolyApproxMethod {
    ePolyApproxBB /* Branch+bound (exponential time but can do things like minimise max kink angle or fix number of segments) */,
    ePolyApproxDijkstra,
    ePolyApproxDP /*Douglas-Peucker*/
};

template <class TPolyline, class TErrorFunction>
void extendAndApproximate(const TPolyline& polyline_in,
                          TPolyline& polyline_out,
                          const typename TPolyline::TPolyApproxSettings& approxSettings,
                          const ePolyApproxMethod approxMethod,
                          const bool bVerbose);

template <class TPolyline>
void removeBuds(const TPolyline& polyline,
                TPolyline& truncatedPolyline,
                const double dMinLength,
                const double dBudScale,
                const bool bVerbose)
{
    typedef typename TPolyline::TControlPoint TControlPoint;

    CHECK(dBudScale <= 0, "Bad bud scale");

    truncatedPolyline.clearAndReserve(polyline.numPoints());

    for(int i = 0; i < polyline.numPoints(); i++) {
        int nBudEnd = -1;
        typename TPolyline::TVecType interpControlPoint;
        const typename TPolyline::ePolyTruncationJoinType bestJoinType =
            polyline.isBud(i, nBudEnd, interpControlPoint, dMinLength, dBudScale, bVerbose);
        switch(bestJoinType) {
        case TPolyline::eNoJoin:
            truncatedPolyline.push_back(polyline[i]);
            break;
        case TPolyline::eReplaceWithIntersection: {
            TControlPoint cp;
            if(TControlPoint::HAS_THICKNESS) {
                const double dInterpThickness = 0.5 * (polyline[nBudEnd].getWidth() + polyline[i].getWidth());
                cp = TControlPoint(interpControlPoint, dInterpThickness);
            } else {
                cp = TControlPoint(interpControlPoint);
            }
            truncatedPolyline.push_back(cp);
            i = nBudEnd;
        } break;
        case TPolyline::eReplaceWithOneEdge:
            truncatedPolyline.push_back(polyline[i]);
            i = nBudEnd - 1;
            break;
        }
    }

    //const int nNumChanges = polyline.numPoints() - truncatedPolyline.numPoints();

/*    if(bVerbose && nNumChanges > 0) {
        polyline.show("withBuds");
        truncatedPolyline.show("noBuds");
    }*/

}

template <class T2dPolyline>
double areaBetweenLineAndPolyline(const T2dPolyline& polylineToApproximate,
                                  const int nFrom,
                                  const int nTo,
                                  const C2dBoundedLine& line)
{
    const bool bVerbose=false;
    typedef std::map<double, C2dImagePointPx> TPointPositionMap;
    TPointPositionMap aPositions;

    for(int nSegmentCandidate = nFrom; nSegmentCandidate < nTo; nSegmentCandidate++) {
        C2dBoundedLine seg = polylineToApproximate.segment(nSegmentCandidate);
        const double eps = 0.0000001;
        optional<const C2dImagePointPx> pIntersect = line.intersection(seg); // TODO check this returns nothing for pll?
        if(pIntersect && fabs(seg.direction().dot(line.direction()) - 1) > eps) // not pll
        {
            const double dPos = polylineToApproximate.positionOfPoint_index(*pIntersect);
            if(fabs(dPos - (double)nFrom) > eps && fabs(dPos - (double)nTo) > eps) // It's not the start or the end
                                                                                   // point
                aPositions[dPos] = *pIntersect;
        }
    }

    T2dPolyline subArea;

    subArea.clearAndReserve((nTo - nFrom) + 3);

    if(aPositions.size() == 0) {
        // Normal case: no sUI().isDebugging()elf-intersection

        // Copy the current polyline into subArea backwards, from a point beyond and a point earlier than this section
        for(int i = nFrom; i <= nTo; i++)
            subArea.push_back(polylineToApproximate[i]);

        subArea.push_back(line.getFinishPoint());
        subArea.push_back(line.getStartPoint());

        return subArea.polygonalArea_nonSelfIntersecting();
    }

    aPositions[nTo] = line.getFinishPoint();
    double dTotalArea = 0;
    int nStartOnThis = nFrom;
    C2dImagePointPx lastIntersectionPoint = line.getStartPoint();

    BOOST_FOREACH(const TPointPositionMap::value_type& intersectionEnd, aPositions) {
        subArea.clearAndReserve(0);

        int nEndOnThis /* index of point on this before the crossing */ = (int)std::floor(intersectionEnd.first);

        for(int i = nStartOnThis; i <= nEndOnThis; i++)
            subArea.push_back(polylineToApproximate[i]);

        subArea.push_back(intersectionEnd.second);
        subArea.push_back(lastIntersectionPoint);
        lastIntersectionPoint = intersectionEnd.second;

        dTotalArea += subArea.polygonalArea_nonSelfIntersecting();

        nStartOnThis = (int)std::ceil(intersectionEnd.first); // = nEndOnThis + 1;
    }

    COUT2("Sum of multiple areas = ", dTotalArea);
    return dTotalArea;
}

template <class TPolyline> class C2dAreaError
{
public:
    double
    operator()(const TPolyline& polylineToApproximate, const int nFrom, const int nTo, const C2dBoundedLine& line) const
    {
        return areaBetweenLineAndPolyline(polylineToApproximate, nFrom, nTo, line);
    }
};

CPolyApproxSettings::CPolyApproxSettings(const double dMin,
                                         const double dMax_in,
                                         const double dEpsilon,
                                         const double dTargetSegLength_in /*=-1 default = max */)
    : minMaxApprox(dMin, dMax_in)
    , dTargetSegLength((dTargetSegLength_in > 0) ? dTargetSegLength_in : dMax_in)
    , dEpsilon(dEpsilon)
{
    CHECK(dTargetSegLength <= dMin, "Target length too short");
}

const double CPolyApproxSettings::penalty(const double dDistFromPrevCP) const
{
    const double dRelPenalty =
        fabs((dDistFromPrevCP - dTargetSegLength) /
             (dTargetSegLength - minMaxApprox.getMin())); // 0 for target distance, 1 for min length
    return dRelPenalty; // * sqr(dEpsilon);
}

std::ostream& operator<<(std::ostream& s, const CPolyApproxSettings& settings)
{
    s << "PolyApproxSettings: " << TO_STRING(settings.dEpsilon) << TO_STRING(settings.dTargetSegLength)
      << TO_STRING(settings.getMinMax().getMin()) << TO_STRING(settings.getMinMax().getMax());
    return s;
}

/**
 * @class CMaxDistanceError
 * @brief maximum distance between any point on the polyline section and the line (not the same as Douglas-Peucker
 * error, which is a cost = infinity or 0 threshold). Always occurs at a polyline vertex
 */
template <class TPolyLine> class CMaxDistanceError
{
public:
    double operator()(const TPolyLine& polylineToApproximate,
                      const int nFrom,
                      const int nTo,
                      const typename TPolyLine::TLineType& line) const
    {
        double dMaxDist = 0;
        for(int i = nFrom; i <= nTo; i++) {
            const typename TPolyLine::TControlPoint& cp = polylineToApproximate[i];
            const double dDist = line.closestDistance(cp.getPoint());
            if(dDist > dMaxDist)
                dMaxDist = dDist;
        }
        return dMaxDist;
    }
};

/**
 * @class CMaxDistanceError
 * @brief exp(maximum distance^2) between any point on the polyline section and the line. Always occurs at a polyline
vertex
 *
template<class TPolyLine>
class CExpMaxDistanceError : public CMaxDistanceError<TPolyLine>
{
public:
    double operator()(const TPolyLine & polylineToApproximate, const int nFrom, const int nTo, const typename
TPolyLine::TLineType & line) const
    {
        const double dMaxDist = CMaxDistanceError<TPolyLine>(polylineToApproximate, nFrom, nTo, line);
        return exp(dMaxDist^2);
    }
};
 */
template <class TControlPoint> void CPolyline_base<TControlPoint>::removeDuplicates()
{
    int nCPOut = 0;
    for(int nCP = 0; nCP < numPoints(); nCP++) {
        if(nCP == 0 || !aControlPoints[nCP].getPoint().isApprox(aControlPoints[nCPOut - 1].getPoint()))
            aControlPoints[nCPOut++] = aControlPoints[nCP];
    }
    aControlPoints.resize(nCPOut);
}

template <class TPolyline> class CNoBudRemover
{
public:
    CNoBudRemover()
    {
    }

    void operator()(const TPolyline& polyline,
                    TPolyline& polyline_out,
                    const CPolyApproxSettings&,
                    const bool bVerbose) const
    {
        // Nothing to do but remove exact duplicates
        polyline_out = polyline;
        polyline_out.removeDuplicates();
    }
};

// typename C2dPolyline_base<TControlPoint>::T2dPolyline C2dPolyline_base<TControlPoint>::removeBuds(const double
// dMinLength, const double dBudScale, const bool bVerbose) const
template <class TPolyline> class CBudRemover
{
public:
    CBudRemover()
    {
    }

    void operator()(const TPolyline& polyline,
                    TPolyline& polyline_out,
                    const C2dPolyApproxSettings& polyApproxSettings,
                    const bool bVerbose) const
    {
        // First remove exact duplicates
        CNoBudRemover<TPolyline> removeDuplicates;
        TPolyline polylineTemp;
        removeDuplicates(polyline, polylineTemp, polyApproxSettings, bVerbose);
        removeBuds<TPolyline>(polylineTemp,
                              polyline_out,
                              polyApproxSettings.getMinMax().getMin(),
                              polyApproxSettings.getBudScale(),
                              bVerbose);
    }
};

template <class TControlPoint>
double
CPolyline_base<TControlPoint>::maxDistanceFromPoly(const typename CPolyline_base<TControlPoint>::TPolyline& other) const
{
    double dMaxDist = 0;
    BOOST_FOREACH(const TControlPoint& cp, aControlPoints) {
        const double dDist = other.distanceToPoint(cp.getPoint(), false);

        if(dDist > dMaxDist)
            dMaxDist = dDist;
    }

    return dMaxDist;
}

template <class TPolyline, class TBudRemover, class TErrorFn>
eSuccessStatus approximate_int(const TPolyline& polyline_in,
                               TPolyline& polyline_out,
                               const typename TPolyline::TPolyApproxSettings& polyApproxSettings,
                               const bool bVerbose)
{
    //const bool bShowPolys = false;

    COUT2("Before approximation: ", polyline_in.numPoints());
    COUT2("Before approximation: ", polyline_in.length());
    //if(bShowPolys)
      //  polyline_in.show("Original");

    TPolyline polyline_temp;

    const TBudRemover budRemover;
    budRemover(polyline_in, polyline_temp, polyApproxSettings, bVerbose);
    if(polyline_temp.numPoints() < 0)
        return eFail;

    COUT2("After bud removal: ", polyline_temp.toString());
    //if(bShowPolys)
      //  polyline_temp.show("NoBuds");

    polyline_out.clearAndReserve(polyline_in.numPoints());

    // cout << polyline_temp.toString() << endl;

    extendAndApproximate<TPolyline, TErrorFn>(
        polyline_temp, polyline_out, polyApproxSettings, ePolyApproxDijkstra, bVerbose);

    COUT2("After extendAndApproximate: ", polyline_out.numPoints());
    COUT2("After extendAndApproximate: ", polyline_out.length());
    //if(bShowPolys)
      //  polyline_out.show("afterApproximate");

    return (polyline_out.numPoints() > 0) ? eSuccess : eFail;
}

template <class TPolyline, class TErrorFn>
eSuccessStatus approximate_int2(const TPolyline& polyline_in,
                                TPolyline& polyline_out,
                                const typename TPolyline::TPolyApproxSettings& polyApproxSettings,
                                const bool bVerbose)
{
    if(polyApproxSettings.getBudScale() > 0)
        return approximate_int<TPolyline, CBudRemover<TPolyline>, TErrorFn>(
            polyline_in, polyline_out, polyApproxSettings, bVerbose);
    else
        return approximate_int<TPolyline, CNoBudRemover<TPolyline>, TErrorFn>(
            polyline_in, polyline_out, polyApproxSettings, bVerbose);
}

template <class TControlPoint>
eSuccessStatus C2dPolyline_base<TControlPoint>::approximate(C2dPolyline_base<TControlPoint>& polyline_out,
                                                            const C2dPolyApproxSettings& polyApproxSettings,
                                                            const bool bVerbose) const
{
    if(polyApproxSettings.useMaxDistanceCost())
        return approximate_int2<C2dPolyline_base<TControlPoint>, CMaxDistanceError<C2dPolyline_base<TControlPoint> > >(
            *this, polyline_out, polyApproxSettings, bVerbose);
    else
        return approximate_int2<C2dPolyline_base<TControlPoint>, C2dAreaError<C2dPolyline_base<TControlPoint> > >(
            *this, polyline_out, polyApproxSettings, bVerbose);
}

eSuccessStatus C3dPolyline::approximate(C3dPolyline& polyline_out,
                                        const CPolyApproxSettings& polyApproxSettings,
                                        const bool bVerbose) const
{
    return approximate_int<C3dPolyline, CNoBudRemover<C3dPolyline>, CMaxDistanceError<C3dPolyline> >(
        *this, polyline_out, polyApproxSettings, bVerbose);
}

eSuccessStatus C3dPolylineWithThickness::approximate(C3dPolylineWithThickness& polyline_out,
                                                     const CPolyApproxSettings& polyApproxSettings,
                                                     const bool bVerbose) const
{
    return approximate_int<C3dPolylineWithThickness,
                           CNoBudRemover<C3dPolylineWithThickness>,
                           CMaxDistanceError<C3dPolylineWithThickness> >(
        *this, polyline_out, polyApproxSettings, bVerbose);
}

template <typename TPolyline> double polyMaxSeparation(const TPolyline& poly1, const TPolyline& poly2)
{
    return std::max<double>(poly1.maxDistanceFromPoly(poly2), poly2.maxDistanceFromPoly(poly1));
}

template double polyMaxSeparation<C3dPolyline>(const C3dPolyline& poly1, const C3dPolyline& poly2);

template <typename TPolyline> double polyMinSeparation(const TPolyline& poly1, const TPolyline& poly2)
{
    return std::min<double>(poly1.maxDistanceFromPoly(poly2), poly2.maxDistanceFromPoly(poly1));
}

template double polyMinSeparation(const C2dPolyline& poly1, const C2dPolyline& poly2);

template double polyMinSeparation(const C2dPolylineWithThickness& poly1, const C2dPolylineWithThickness& poly2);

//

bool approxArea(const double dArea1, const double dArea2)
{
    return fabs(dArea1 - dArea2) < 1e-6;
}

void testAreaComputation()
{
    const bool bVerbose = false;

    C2dPolyline test, linePllToTest;
    test.push_back(C2dImagePointPx(100, 100)); // 0
    test.push_back(C2dImagePointPx(150, 100)); // 1
    test.push_back(C2dImagePointPx(200, 150)); // 2 a triangle area 2500
    test.push_back(C2dImagePointPx(250, 100)); // 3
    test.push_back(C2dImagePointPx(300, 100)); // 4
    test.push_back(C2dImagePointPx(300, 200)); // 5 a square, area 10000
    test.push_back(C2dImagePointPx(400, 200)); // 6
    test.push_back(C2dImagePointPx(400, 100)); // 7
    test.push_back(C2dImagePointPx(500, 100)); // 8

    BOOST_FOREACH(const C2dPolylineControlPoint& cp, test) {
        linePllToTest.push_back(C2dImagePointPx(cp.getPoint().x(), test[0].getPoint().y()));
    }

    const double adAreas[] = {
        0, 0, 1250, 2500, /*4*/ 2500, 2500, /*6*/ 12500, 12500, 12500
    }; // Cumulative area up-to index i
    const double adOffsetAreas[] = {
        0, 2500, 3750, 5000, /*4*/ 7500, 7500, /*6*/ 12500, 12500, 17500
    }; // Area up-to index i

    for(int i = 1; i < test.numPoints(); i++) {
        const C2dBoundedLine line(linePllToTest[0].getPoint(), linePllToTest[i].getPoint());
        CHECK(!approxArea(areaBetweenLineAndPolyline(test, 0, i, line), adAreas[i]), "Area computation test failed");

        if(i > 1) {
            const C2dBoundedLine line(linePllToTest[1].getPoint(), linePllToTest[i].getPoint());
            CHECK(!approxArea(areaBetweenLineAndPolyline(test, 1, i, line), adAreas[i]),
                  "Area computation test failed");
        }

        // And an offset line
        TEigen2dPoint offset(0, -50);
        const C2dBoundedLine line_offset(linePllToTest[0].getPoint() + offset, linePllToTest[i].getPoint() + offset);
        const double dArea = adAreas[i] + line_offset.length() * 50;
        CHECK(!approxArea(areaBetweenLineAndPolyline(test, 0, i, line_offset), dArea),
              "Offset area computation test failed");

        // And an offset line 50 px away, so it intersects the square
        TEigen2dPoint offset2(0, 50);
        const C2dBoundedLine line_offset2(linePllToTest[0].getPoint() + offset2, linePllToTest[i].getPoint() + offset2);
        CHECK(!approxArea(areaBetweenLineAndPolyline(test, 0, i, line_offset2), adOffsetAreas[i]),
              "Offset area computation test failed");

        // And an offset line 50 px away, and sloping slightly back so it intersects the top of the triangle and the
        // square, and isn't pll anywhere
        TEigen2dPoint offset3(-0.1, -1);
        const C2dBoundedLine line_offset3(linePllToTest[0].getPoint() + offset2 + offset3,
                                          linePllToTest[i].getPoint() + offset2 + i * offset3);
        CHECK(!within(areaBetweenLineAndPolyline(test, 0, i, line_offset3), adOffsetAreas[i], 0.1),
              "Offset area computation test failed");
    }

    // for(ePolyApproxMethod method = ePolyApproxBB /* Branch+bound */; method <= ePolyApproxDijkstra; ((int&)method)++
    // )
    {
        C2dPolyline polyApprox;

        const double eps = 1e-6;
        C2dPolyApproxSettings approxSettings0(1, 1000);
        CHECK(test.approximate(polyApprox, approxSettings0, bVerbose) == eFail, "Error approximating polyline");
        CHECK(polyApprox.numPoints() < test.numPoints(), "This approximation shouldn't have lost any points");
        CHECK(polyMaxSeparation<C2dPolyline>(polyApprox, test) > eps, "Poly approiximation failed");

        C2dPolyApproxSettings approxSettings1(40, 110);
        CHECK(test.approximate(polyApprox, approxSettings1, bVerbose) == eFail, "Error approximating polyline");
        CHECK(polyApprox.numPoints() < test.numPoints(), "This approximation shouldn't have lost any points");
        CHECK(polyMaxSeparation<C2dPolyline>(polyApprox, test) > eps, "Poly approiximation failed");

        C2dPolyApproxSettings approxSettings2(1, 90);
        CHECK(test.approximate(polyApprox, approxSettings2, bVerbose) == eFail, "Error approximating polyline");
        CHECK(polyApprox.numPoints() < test.numPoints() + 3, "This approximation should have gained at least 3 points");
        CHECK(polyMaxSeparation<C2dPolyline>(polyApprox, test) > eps, "Poly approximation failed");

        C2dPolyApproxSettings approxSettings3(110, 220);
        CHECK(test.approximate(polyApprox, approxSettings3, bVerbose) == eFail, "Error approximating polyline");
        CHECK(polyApprox.getStartPoint() != test.getStartPoint(), "Poly approx start point failed");
        CHECK(polyApprox.getFinishPoint() != test.getFinishPoint(), "Poly approx start point failed");

        C2dPolyApproxSettings approxSettings4(220, 440);
        CHECK(test.approximate(polyApprox, approxSettings4, bVerbose) == eFail, "Error approximating polyline");
        CHECK(test.approxAsLine().length() != polyApprox.approxAsLine().length(), "Error approximating polyline");

        C2dPolyApproxSettings approxSettings5(1, 2);
        CHECK(test.approximate(polyApprox, approxSettings5, bVerbose) == eFail, "Error approximating polyline");
        const double dError = polyMaxSeparation<C2dPolyline>(polyApprox, test);
        CHECK_P(dError > 1 /* eps */, dError, "Poly approximation with many points failed");
        if(dError > eps)
            COUT2("Warning: error greater than 0. What is going on with kink angles here?", dError);
    }
}

//////////// Dijkstra DAG optimal approximation ///////////////////

/*Setup edges. For each vertex we can iterate over in-edges (from earlier vertices at distances between min and max) and
*out-edges (to later vertices at distances between min and max)
*
* Each vertex has a lowestCostSource index (initially not set) and a lowestCost variable (initially infinity)
* For each vertex in turn: update shortest path to all out-vertices in future (also on 'extra' polyline)
*
*/

//////////// B+B optimal approximation ///////////////////
// typedef C2dPolyline_base<TControlPoint>::T2dPolyline T2dPolyline;

/**
 * @brief
 *
 * Must first split the start and end sections enough that there exists a point between min and max distance from the
 *end TODO
 *
 * @param polylineToApproximate //the original polyline, plus some extra control points. Choose a subset of these XOR
 *extraPoints as the approximation
 * @param extraPoints //If non-zero, these points are alternatives to the corresponding points in polylineToApproximate
 *(i.e. midpoint OR point interpolated from 2 neighbours, endpoint OR interpolated endpoint)
 * @param approximation
 * @param dCumulativeError
 * @param nIdxIntoPolylineWithExtraPoints
 * @param nIndexToSecondLast We've computed the error up to here (-1 for empty polylines, 0 for 1 point polylines)
 * @param bestApproximation
 * @param dLeastCumulativeError
 * @param dMinLength
 * @param dMaxLength
 */
// template<class TControlPoint,>
template <class TPolyline, class TErrorFunction>
void approxPolyBB(const TPolyline& polylineToApproximate,
                  const TPolyline& extraPoints,
                  TPolyline& approximation,
                  double dCumulativeError,
                  const int nIdxIntoPolyline,
                  const int nIndexToSecondLast,
                  TPolyline& bestApproximation,
                  double& dLeastCumulativeError,
                  const double dMinLength,
                  const double dMaxLength)
{
    // Update the error: we just added the back point to \approximation at position \nIdxIntoPolyline
    if(nIdxIntoPolyline > 0) {
        // Compute the new added area. Might be duplicate control points. Segments might cross

        CHECK(nIndexToSecondLast >= nIdxIntoPolyline || nIndexToSecondLast < 0 || nIdxIntoPolyline < 0,
              "Index mismatch");

        // const double dExtraError = areaBetweenLineAndPolyline(polylineToApproximate, nIndexToSecondLast,
        // nIdxIntoPolyline, approximation.segment(approximation.numSegments()-1));
        TErrorFunction errorFn;
        const double dExtraError = errorFn(polylineToApproximate,
                                           nIndexToSecondLast,
                                           nIdxIntoPolyline,
                                           approximation.segment(approximation.numSegments() - 1));
        dCumulativeError += dExtraError;
    }

    // If we're at the end then consider how good an approximation it is:
    if(nIdxIntoPolyline == polylineToApproximate.numPoints() - 1) {
        if(dCumulativeError < dLeastCumulativeError) {
            dLeastCumulativeError = dCumulativeError;
            bestApproximation = approximation;
        }
    } else {
        for(int i = nIdxIntoPolyline + 1; i < polylineToApproximate.numPoints(); i++) {
            if(approximation.numPoints() > 0) {
                const double dLengthToThisCP =
                    (approximation.getFinishPoint() - polylineToApproximate[i].getPoint()).norm();
                if(dLengthToThisCP > dMaxLength)
                    break;

                if(dLengthToThisCP < dMinLength)
                    continue;
            } else {
                if(i > 0)
                    break;
            }

            // Add either polylineToApproximate[i] or extraPoints[i] as the next point
            approximation.push_back(polylineToApproximate[i]);
            approxPolyBB<TPolyline, TErrorFunction>(polylineToApproximate,
                                                    extraPoints,
                                                    approximation,
                                                    dCumulativeError,
                                                    i,
                                                    nIdxIntoPolyline,
                                                    bestApproximation,
                                                    dLeastCumulativeError,
                                                    dMinLength,
                                                    dMaxLength);
            approximation.pop_back();

            if(extraPoints[i].getPoint() != TPolyline::TControlPoint::uninit()) {
                approximation.push_back(polylineToApproximate[i]);
                approxPolyBB<TPolyline, TErrorFunction>(polylineToApproximate,
                                                        extraPoints,
                                                        approximation,
                                                        dCumulativeError,
                                                        i,
                                                        nIdxIntoPolyline,
                                                        bestApproximation,
                                                        dLeastCumulativeError,
                                                        dMinLength,
                                                        dMaxLength);
                approximation.pop_back();
            }
        }
    }
}

// State of a vertex
class CDijkstraState
{
    double dClosest;
    int nIndexOfPredecessor, nCandidateSource;

public:
    CDijkstraState()
        : dClosest(HUGE)
        , nIndexOfPredecessor(-1)
        , nCandidateSource(-1)
    {
    }

    void update(const double dDistance, const int nIndex, const int nCandidateSource_in)
    {
        if(dDistance < dClosest) {
            dClosest = dDistance;
            nIndexOfPredecessor = nIndex;
            nCandidateSource = nCandidateSource_in;
        }
    }

    bool isReachable() const
    {
        return dClosest < HUGE;
    }

    double getLeastError() const
    {
        return dClosest;
    }

    void getPrevIndices(int& nIndex, int& nCandidateSource_out)
    {
        CHECK(!isReachable(), "Unreachable vertex selected when reconstructing shortest path");
        nIndex = nIndexOfPredecessor;
        nCandidateSource_out = nCandidateSource;
    }
};

class CDijkstraStatePair
{
    CDijkstraState aPolyVertexStates[2]; // original pouyline, and extra points
public:
    CDijkstraStatePair()
    {
    }

    CDijkstraState& getState(const int nCandidateSource)
    {
        return aPolyVertexStates[nCandidateSource];
    }
};

/* Dijkstra's algorithm on a Directed Acyclic Graph:
 * http://en.wikipedia.org/wiki/Dijkstra%27s_algorithm#Running_time (see 'note')
 * */
template <class TPolyline, class TErrorFn>
void approxPolyDijkstra(const TPolyline& thisPolylineWithExtraControlPoints,
                        const TPolyline& extraPoints,
                        TPolyline& approximation,
                        const typename TPolyline::TPolyApproxSettings& polyApproxSettings)
{
    const bool bVerbose = false;
    // make all this more verbose. Afterwards un-comment the other tests
    std::vector<CDijkstraStatePair> aStates(thisPolylineWithExtraControlPoints.numPoints());
    const TPolyline* apCandidateSources[2] = { &thisPolylineWithExtraControlPoints, &extraPoints };

    for(int nCandidateSource = 0; nCandidateSource < 2; nCandidateSource++) {
        if((*apCandidateSources[nCandidateSource])[0].getPoint() != TPolyline::TControlPoint::uninit())
            aStates[0].getState(nCandidateSource).update(0, -1, -1);
    }

    for(int i = 0; i < thisPolylineWithExtraControlPoints.numPoints(); i++) {
        for(int nCandidateSource = 0; nCandidateSource < 2; nCandidateSource++) {
            if(!aStates[i].getState(nCandidateSource).isReachable())
                continue;

            const double dShortestPathToThisVertex = aStates[i].getState(nCandidateSource).getLeastError();

            typename TPolyline::TControlPoint lastCPonPoly = (*apCandidateSources[nCandidateSource])[i];
            if(lastCPonPoly.getPoint() == TPolyline::TControlPoint::uninit())
                continue;

            if(bVerbose)
                cout << "Completed path up to CP " << i << "(" << nCandidateSource << ")"
                     << " length=" << dShortestPathToThisVertex << endl;

            int nMaxIdx = -1;
            for(int nCandidateDest = 0; nCandidateDest < 2; nCandidateDest++) {
                const int nUpperLimit = (nMaxIdx < 0) ? thisPolylineWithExtraControlPoints.numPoints() : nMaxIdx;
                for(int nFuture = i + 1; nFuture < nUpperLimit; nFuture++) {
                    if(bVerbose)
                        cout << "Future CP " << nFuture << "(" << nCandidateDest << ")"
                             << ": ";

                    typename TPolyline::TControlPoint nextCPonPoly = (*apCandidateSources[nCandidateDest])[nFuture];

                    if(nextCPonPoly.getPoint() == TPolyline::TControlPoint::uninit()) {
                        if(bVerbose)
                            cout << "uninit" << endl;
                        continue;
                    }

                    const double dDistToNextCP = (lastCPonPoly.getPoint() - nextCPonPoly.getPoint()).norm();

                    if(dDistToNextCP == 0)
                        cout << "Warning: duplicate point in polyline: " << lastCPonPoly << " in "
                             << apCandidateSources[nCandidateSource]->toString() << endl;

                    if(polyApproxSettings.getMinMax().tooLong(dDistToNextCP)) {
                        if(bVerbose)
                            cout << "too distant" << endl;
                        nMaxIdx = nFuture + 1;
                        break;
                    }

                    if(polyApproxSettings.getMinMax().tooShort(dDistToNextCP)) {
                        if(bVerbose)
                            cout << dDistToNextCP << " too close" << endl;
                        continue;
                    }

                    typename TPolyline::TLineType candidateSeg(
                        lastCPonPoly.getPoint(), nextCPonPoly.getPoint()); // This fails for self-intersecting polylines
                    // const double dDistToNextCP = candidateSeg.length();

                    TErrorFn errorFn;
                    const double dExtraError =
                        sqr(errorFn(thisPolylineWithExtraControlPoints, i, nFuture, candidateSeg) /
                            ((polyApproxSettings.getEpsilon() > 0) ? polyApproxSettings.getEpsilon() : 1));

                    const double dTotalError =
                        dShortestPathToThisVertex + dExtraError + polyApproxSettings.penalty(dDistToNextCP);
                    COUT(dDistToNextCP);
                    COUT(polyApproxSettings.penalty(dDistToNextCP));

                    aStates[nFuture].getState(nCandidateDest).update(dTotalError, i, nCandidateSource);

                    if(bVerbose)
                        cout << "new shortest path length = "
                             << aStates[nFuture].getState(nCandidateDest).getLeastError() << endl;
                }
            }
        }
    }
    double dBestEndpointError = HUGE;
    const int nLastPoint = thisPolylineWithExtraControlPoints.numPoints() - 1;
    int nBestEndpoint = -1;
    for(int nCandidateEndpointSource = 0; nCandidateEndpointSource < 2; nCandidateEndpointSource++) {
        if(aStates[nLastPoint].getState(nCandidateEndpointSource).isReachable()) {
            if(aStates[nLastPoint].getState(nCandidateEndpointSource).getLeastError() < dBestEndpointError) {
                dBestEndpointError = aStates[nLastPoint].getState(nCandidateEndpointSource).getLeastError();
                nBestEndpoint = nCandidateEndpointSource;
            }
        }
    }

    if(nBestEndpoint == -1 || dBestEndpointError >= HUGE) {
        ALWAYS_VERBOSE;
        COUT(thisPolylineWithExtraControlPoints.toString());
        COUT(polyApproxSettings);
        
        std::cerr << "No shortest path found: Probably opposite-join has joined polylines which give a kinked polyline when concatenated. TODO restore exception (?). "<< thisPolylineWithExtraControlPoints.toString() << endl << polyApproxSettings << endl;

        //THROW("No shortest path found");
        
        //Yuk
        const auto temp = thisPolylineWithExtraControlPoints.removeShortSections(polyApproxSettings.getMinMax().getMin(), bVerbose);
        const auto temp2 = temp.splitLongSections(polyApproxSettings.getMinMax().getMax(), bVerbose);
        
        for(const auto & cp : temp2)
        {
            approximation.push_back(cp);
        }
        
        std::cerr << "After approximation: " << approximation.toString() << endl;
        return;
    }

    int nIndex = nLastPoint, nCandidateSource = nBestEndpoint;
    for(;;) {
        approximation.push_back((*apCandidateSources[nCandidateSource])[nIndex]);

        aStates[nIndex].getState(nCandidateSource).getPrevIndices(nIndex, nCandidateSource);

        if(nIndex < 0)
            break;

        CHECK(nCandidateSource < 0, "Error reconstructing polyline shortest path");
    }

    approximation.reverseDirection();
}

template <class TPolyline, class TErrorFn>
void extendAndApproximate(const TPolyline& polyline_in,
                          TPolyline& polyline_out,
                          const typename TPolyline::TPolyApproxSettings& approxSettings,
                          const ePolyApproxMethod approxMethod,
                          const bool bVerbose)
{
    //const bool bShowPolys = false;

    if(polyline_in.approxAsLine().length() < approxSettings.getMinMax().getMin())
        return;

    TPolyline thisPolylineWithExtraControlPoints; // Add midpoints as well
    TPolyline interpolatedPoints, extraPoints; // Indices correspond to thisPolylineWithExtraControlPoints. Add all the
                                               // possible extra points to this. We'll choose a subset of these for the
                                               // approximation

    /*interpolatedPoints.clearAndReserve(polyline_in.numPoints());

    //Setup interpolated points
    for(int i=0; i<polyline_in.numPoints(); i++ )
    {
        const C2dBoundedLine segi=polyline_in.segment(std::min<int>(i, polyline_in.numSegments()-1);
        const typename TPolyline::TVecType diri=segi.direction();
        const typename TPolyline::TVecType extraPoint = typename TPolyLine::TControlPoint::uninit();
        if(i==0)
        {
            if(polyline_in.numPoints()>2)//we can interpolate an alternative start point
            {
                const C2dLine seg1=polyline_in.segment(1).unboundedLine();
                const double dDot = seg1.direction().dot(diri);
                if(dDot != 1 && dDot > 0)
                     extraPoint=seg1.closestPoint(polyline_in.getStartPoint());
            }

        }
        else if (i==polyline_in.numPoints()-1)
        {
            if(polyline_in.numPoints()>2)//we can interpolate an alternative end point
            {
                const C2dLine seg_1=polyline_in.segment(polyline_in.numSegments()-2).unboundedLine();
                const double dDot = seg_1.direction().dot(diri);
                if(dDot != 1 && dDot > 0)
                     extraPoint=seg_1.closestPoint(polyline_in.getFinishPoint());
            }
        }
        else
        {

        }
        interpolatedPoints.push_back(extraPoint);
    }*/

    for(int i = 0; i < polyline_in.numPoints(); i++) {
        thisPolylineWithExtraControlPoints.push_back(polyline_in[i]);

        if(i < polyline_in.numSegments()) {
            const typename TPolyline::TControlPoint::TLineType thisSeg = polyline_in.segment(i);
            const double dThisSegLength = thisSeg.length();

            // We want 2 splits between the min and max length from each CP:
            double dTargetLength = 0.5 * (approxSettings.getMinMax().getMax() - approxSettings.getMinMax().getMin());

            dTargetLength = std::min<double>(dTargetLength, 0.5 * approxSettings.getTargetLength());

            int nNumSplits = std::max<int>(
                2,
                (int)std::ceil(dThisSegLength /
                               dTargetLength)); // number of extra points. TODO: make sure there's enough points here.

            // These checks should rarely ever be hit now...

            // If i>0 and the start of this segment is within dMinLength of the startpoint then increase the number of
            // splits until the first split is less than dMaxLength of the start (otherwise BB will get stuck)
            // Repectively the end

            // This can still fail if min is too close to max
            // CHECK(dMinLength > 0.5*dMaxLength, "Min too close to max--BB may fail to step past a short end segment");
            if(i >= 1) {
                const double dDistToThisSegStart = (thisSeg.getStartPoint() - polyline_in.getStartPoint()).norm();
                if(dDistToThisSegStart < approxSettings.getMinMax().getMin()) {
                    for(;;) {
                        const double t = 1.0 / nNumSplits;
                        const double dDistToThisSegFirstSplitpoint =
                            (thisSeg.pointFromt(t) - polyline_in.getStartPoint()).norm();
                        if(dDistToThisSegFirstSplitpoint < approxSettings.getMinMax().getMax())
                            break;
                        nNumSplits++;
                    }
                }
            }
            if(i < polyline_in.numSegments() - 1) {
                std::vector<CDijkstraStatePair> aStates(thisPolylineWithExtraControlPoints.numPoints());

                const double dDistToThisSegEnd = (thisSeg.getFinishPoint() - polyline_in.getFinishPoint()).norm();
                if(dDistToThisSegEnd < approxSettings.getMinMax().getMin()) {
                    std::vector<CDijkstraStatePair> aStates(thisPolylineWithExtraControlPoints.numPoints());

                    for(;;) {
                        const double t = 1.0 - 1.0 / nNumSplits;
                        const double dDistToThisSegLastSplitpoint =
                            (thisSeg.pointFromt(t) - polyline_in.getFinishPoint()).norm();
                        if(dDistToThisSegLastSplitpoint < approxSettings.getMinMax().getMax())
                            break;

                        nNumSplits++;
                    }
                }
            }
            for(int nSplitPoint = 1; nSplitPoint < nNumSplits; nSplitPoint++) {
                const double t = (double)nSplitPoint / (double)nNumSplits;
                const typename TPolyline::TVecType point = thisSeg.pointFromt(t);
                typename TPolyline::TControlPoint cp;
                if(TPolyline::TControlPoint::HAS_THICKNESS) {
                    const double dStartThickness = polyline_in[i].getWidth();
                    const double dEndThickness = polyline_in[i + 1].getWidth();
                    cp = typename TPolyline::TControlPoint(point, t * dEndThickness + (1 - t) * dStartThickness);
                } else {
                    cp = typename TPolyline::TControlPoint(point);
                }
                thisPolylineWithExtraControlPoints.push_back(cp);
            }
        }
    }

    // extraPoints add interp start point TODO
    for(int i = 0; i < thisPolylineWithExtraControlPoints.numPoints(); i++) {
        typename TPolyline::TControlPoint uninit = TPolyline::TControlPoint::uninitCP();
        extraPoints.push_back(uninit);
    }
    // extraPoints add interp end point TODO

    /*if(bShowPolys) {
        polyline_in.show("polyBeforeApprox");
        thisPolylineWithExtraControlPoints.show("polyExtraPoints");
        extraPoints.show("alternativeControlPoints");
    }*/

    CHECK(!zero(thisPolylineWithExtraControlPoints.maxKinkAngle() - polyline_in.maxKinkAngle()),
          "Adding points has changed poly max kink angle");
    REPEAT(100,
           CHECK(!zero(polyMaxSeparation<TPolyline>(thisPolylineWithExtraControlPoints, polyline_in)),
                 "Adding points has changed poly")); // This can trip thickness constraints e.g. before a bad 3d
                                                     // component is removed from the optimisation

    // TPolyline bestApproximation;

    if(approxMethod == ePolyApproxBB) // Doesn't actually bound at the moment, so its slow.
    {
        const double dCumulativeError = 0;
        double dLeastCumulativeError = HUGE;

        TPolyline approximation;
        approxPolyBB<TPolyline, TErrorFn>(thisPolylineWithExtraControlPoints,
                                          extraPoints,
                                          approximation,
                                          dCumulativeError,
                                          -1,
                                          -1,
                                          polyline_out,
                                          dLeastCumulativeError,
                                          approxSettings.getMinMax().getMin(),
                                          approxSettings.getMinMax().getMax());

        CHECK(dLeastCumulativeError >= HUGE, "Failed to find *any* approximation--this should be impossible");
        CHECKBADNUM(dLeastCumulativeError);
    } else if(approxMethod == ePolyApproxDijkstra) {
        approxPolyDijkstra<TPolyline, TErrorFn>(
            thisPolylineWithExtraControlPoints, extraPoints, polyline_out, approxSettings);
    } else
        THROW("polyline approx method not handled");

    /*if(bShowPolys) {
        polyline_out.show("bestApproximatn");
    }*/
}

///////////////////////////////////////////////////////////////////////////

optional<const C3dPolylineWithThickness> C2dPolylineWithThickness::polyToWorld(const CWorldCamera& P,
                                                                               const double dDepthPrior) const
{
    C3dPolylineWithThickness polyProjOutOtWorld;
    polyProjOutOtWorld.clearAndReserve(numPoints());
    BOOST_FOREACH(const C2dPolylineControlPointWithThickness& cp, *this) {
        const C3dWorldPoint p3d = P.pxToWorld_z(cp.getPoint(), dDepthPrior);
        optional<const double> dWidthMetres = P.pxToWidth(p3d, cp.getWidth());
        if(!dWidthMetres)
            return optional<const C3dPolylineWithThickness>();
        polyProjOutOtWorld.push_back(C3dPolylineControlPointWithThickness(p3d, *dWidthMetres));
    }
    return polyProjOutOtWorld;
}

template <class TControlPoint>
void CPolyline_base<TControlPoint>::transform(const Eigen::Matrix<double,
                                                                  TControlPoint::TVecType::RowsAtCompileTime + 1,
                                                                  TControlPoint::TVecType::RowsAtCompileTime + 1>& T)
{
    BOOST_FOREACH(TControlPoint& X, *this) {
        Eigen::Matrix<double, TControlPoint::TVecType::RowsAtCompileTime + 1, 1> X_trans_homog =
            T * X.getPoint().homog();
        X.getPoint() = (X_trans_homog.head(TControlPoint::TVecType::RowsAtCompileTime));
    }
}

template <class TControlPoint>
void C3dPolyline_base<TControlPoint>::project_incomplete(const CWorldCamera& P, C2dPolyline& thinPoly) const
{
    thinPoly.clearAndReserve(this->numPoints());
    BOOST_FOREACH(const TControlPoint& X, *this) {
        optional<const C2dImagePointPx> px = P.projectToPx(X.getPoint());

        if(px)
            thinPoly.push_back(C2dPolylineControlPoint(*px));
    }
}

template <class TControlPoint>
void C3dPolyline_base<TControlPoint>::project_incomplete(const CWorldCamera& P, C2dPolylineWithThickness& poly) const
{
    CHECK(TControlPoint::HAS_THICKNESS == false,
          "Don't use this projection fn (project to thick polyline) for thin 3D polylines");

    poly.clearAndReserve(this->numPoints());
    BOOST_FOREACH(const TControlPoint& X, this->aControlPoints) {
        optional<const C2dImagePointPx> px =
            P.projectToPx(X.getPoint()); // Restricts to nearby. Todo make a projectPointAndWidth method.

        if(px) {
            const optional<const double> pdWidth =
                P.widthToPx(X.getWidth(), X.getPoint());

            if(pdWidth) {
                poly.push_back(C2dPolylineControlPointWithThickness(*px, *pdWidth));
            }
        }
    }
}

template <class TControlPoint> bool C3dPolyline_base<TControlPoint>::isVisible(const CWorldCamera& P) const
{
    BOOST_FOREACH(const TControlPoint& X, this->aControlPoints) {
        optional<const C2dImagePointPx> px = P.projectToPx(X.getPoint());
        if(px && px->x() > -25)
            return true;
    }
    return false;
}

template <class TControlPoint>
void C3dPolyline_base<TControlPoint>::projectOne_fast(const CWorldCamera& P,
                                                      C2dPolyline& thinPoly,
                                                      const int nControlPointIdx) const
{
    const TControlPoint& X = this->aControlPoints[nControlPointIdx];
    thinPoly[nControlPointIdx] = C2dPolylineControlPoint(P.fastProject(X.getPoint()));
}

// nControlPointIdx is roughly the same as LM -- the idx of the only point that has changed when computing derivatives.
template <class TControlPoint>
eSuccessStatus C3dPolyline_base<TControlPoint>::projectOne(const CWorldCamera& P,
                                                           C2dPolyline& thinPoly,
                                                           const int nControlPointIdx,
                                                           const bool bVerboseOnFailure) const
{
    const TControlPoint& X = this->aControlPoints[nControlPointIdx];
    optional<const C2dImagePointPx> px = P.projectToPx(X.getPoint());

    if(!px) {
        if(bVerboseOnFailure) {
            REPEAT(5000, cout << "polyline.project failed because X=" << X << endl);
        }
        return eFail;
    }

    thinPoly[nControlPointIdx] = C2dPolylineControlPoint(*px);

    return eSuccess;
}

// nControlPointIdx is roughly the same as LM -- the idx of the only point that has changed when computing derivatives.
template <class TControlPoint>
eSuccessStatus C3dPolyline_base<TControlPoint>::project(const CWorldCamera& P,
                                                        C2dPolyline& polyLine2d,
                                                        const bool bVerboseOnFailure,
                                                        const int nControlPointIdx) const
{
    if(nControlPointIdx >= 0) {
        projectOne_fast(P, polyLine2d, nControlPointIdx);
        return eSuccess;
    }

    polyLine2d.resize(this->numPoints());

    for(int i = 0; i < this->numPoints(); i++) {
        if(projectOne(P, polyLine2d, i, bVerboseOnFailure) == eFail)
            return eFail;
    }

    return eSuccess;
}

// nControlPointIdx is roughly the same as LM -- the idx of the only point that has changed when computing derivatives.
template <class TControlPoint>
eSuccessStatus C3dPolyline_base<TControlPoint>::projectOne(const CWorldCamera& P,
                                                           C2dPolylineWithThickness& polyLine2d,
                                                           const int nControlPointIdx,
                                                           const bool bVerboseOnFailure) const
{
    const TControlPoint& X = this->aControlPoints[nControlPointIdx];

    optional<const C2dImagePointPx> px =
        P.projectToPx(X.getPoint()); // Restricts to nearby.

    if(!px) {
        if(bVerboseOnFailure) {
            REPEAT(5000, cout << "polyline.project failed because X=" << X << endl);
        }
        return eFail;
    }

    const optional<const double> pdWidth = P.widthToPx(X.getWidth(), X.getPoint());

    if(!pdWidth) {
        if(bVerboseOnFailure) {
            if(pdWidth)
                cout << "Width projects to " << *pdWidth << ", ";
            cout << "Failed to project polyline to thin 2D polyline because of width projecting point X=" << X << endl;
        }

        return eFail;
    }

    polyLine2d[nControlPointIdx] = C2dPolylineControlPointWithThickness(*px, *pdWidth);

    return eSuccess;
}

template <class TControlPoint>
eSuccessStatus C3dPolyline_base<TControlPoint>::project(const CWorldCamera& P,
                                                        C2dPolylineWithThickness& polyLine2d,
                                                        const bool bVerboseOnFailure,
                                                        const int nControlPointIdx) const
{
    CHECK(TControlPoint::HAS_THICKNESS == false,
          "Don't use this projection fn (project to thick polyline) for thin 3D polylines");

    if(nControlPointIdx >= 0)
        return projectOne(P, polyLine2d, nControlPointIdx, bVerboseOnFailure);

    polyLine2d.resize(this->numPoints());

    for(int i = 0; i < this->numPoints(); i++) {
        if(projectOne(P, polyLine2d, i, bVerboseOnFailure) == eFail)
            return eFail;
    }

    return eSuccess;
}

eSuccessStatus C3dPolyline::project_mask(const CWorldCamera& P,
                                         C2dPolyline& polyLine2d,
                                         std::vector<int>& aCorresponding3dPolyIndices) const
{
    aCorresponding3dPolyIndices.clear();
    aCorresponding3dPolyIndices.reserve(numPoints());
    polyLine2d.clearAndReserve(numPoints());
    for(int i = 0; i < numPoints(); i++) {
        const TControlPoint& X = aControlPoints[i];
        optional<const C2dImagePointPx> px =
            P.projectToPx(X.getPoint()); // Restricts to nearby.

        if(px) {
            polyLine2d.push_back(C2dPolylineControlPoint(*px));
            aCorresponding3dPolyIndices.push_back(i);
        }
    }
    return (polyLine2d.numPoints() > 1) ? eSuccess : eFail;
}

/*Translated from http://en.wikipedia.org/wiki/Ramer%E2%80%93Douglas%E2%80%93Peucker_algorithm
 * TODO: DP is not optimal--Im sure there's a better approach.
 *
 */
template <class TControlPoint>
typename CPolyline_base<TControlPoint>::TPolyline
CPolyline_base<TControlPoint>::douglasPeucker(const TPolyline& polyLine,
                                              const double dMaxLength,
                                              const double dEpsilon,
                                              const bool bVerbose)
{
    if(polyLine.numPoints() == 1)
        return polyLine;

    // Find the point with the maximum distance (furthest from the line between the endpoints)
    double dmax = -1;
    int index = 0;
    typename TControlPoint::TLineType boundedLine(
        polyLine.getStartPoint(),
        polyLine.getFinishPoint()); // TODO: possibly this should be unbounded but probably not.
    for(int i = 1; i < (int)polyLine.numPoints() - 1; i++) {
        const double d = boundedLine.closestDistance(polyLine[i].getPoint());
        if(d > dmax) {
            index = i;
            dmax = d;
        }
    }

    // If max distance is greater than epsilon, recursively simplify
    if(dmax >= dEpsilon) {
        // Recursive call
        const TPolyline lineFirstPart(polyLine.begin(), polyLine.begin() + index);
        TPolyline firstPartResults = douglasPeucker(lineFirstPart, dMaxLength, dEpsilon, bVerbose);
        const TPolyline lineSecondPart(polyLine.begin() + index, polyLine.end());
        const TPolyline secondPartResults = douglasPeucker(lineSecondPart, dMaxLength, dEpsilon, bVerbose);

        firstPartResults.aControlPoints.insert(firstPartResults.aControlPoints.end(),
                                               secondPartResults.aControlPoints.begin(),
                                               secondPartResults.aControlPoints.end());
        return firstPartResults;
    }

    TPolyline startAndEnd;
    startAndEnd.push_back(polyLine[0]);

    if((polyLine[polyLine.numPoints() - 1].getPoint() - polyLine[0].getPoint()).squaredNorm() > sqr(dMaxLength)) {
        // SPECIAL CASE: we don't want to delete all the midpoints here because that will leave us with too long a
        // segment
        // We could either prefer to keep the midpoint with the largest deviation, or the midpoint which is most central
        // We really want to keep a minimal number of points
        // Sometimes the section might already be too long

        for(int i = 1; i < polyLine.numPoints() - 1; i++) {
            if((startAndEnd.getFinishPoint() - polyLine[i + 1].getPoint()).squaredNorm() >
               sqr(dMaxLength)) // The next will be too far away
                startAndEnd.push_back(polyLine[i]); // could end up with 2 close together still
        }
    }

    startAndEnd.push_back(polyLine[polyLine.numPoints() - 1]);

    return startAndEnd;
}

/**
 * @brief This is not always correct either
 * @param
 * @return
 */
template <class TControlPoint>
typename CPolyline_base<TControlPoint>::TPolyline
CPolyline_base<TControlPoint>::removeShortSections(const double dMinSectionLength, const bool bVerbose) const
{
    COUT(dMinSectionLength);

    TPolyline poly_noShortSections;

    if(numPoints() < 2 || (getStartPoint() - getFinishPoint()).norm() < dMinSectionLength)
        return poly_noShortSections;

    poly_noShortSections.push_back(getEndControlPoint(eStart));

    typename TPolylineVector::value_type endPoint = getEndControlPoint(eFinish);
    for(int i = 1; i < numPoints() - 1; i++) {
        const typename TPolylineVector::value_type newBackPoint =
            poly_noShortSections[poly_noShortSections.numPoints() - 1];
        if((aControlPoints[i].getPoint() - newBackPoint.getPoint()).norm() > dMinSectionLength &&
           (aControlPoints[i].getPoint() - endPoint.getPoint()).norm() > dMinSectionLength) {
            poly_noShortSections.push_back(aControlPoints[i]);
            if(IS_DEBUG)
                CHECK(poly_noShortSections.segment((int)poly_noShortSections.numPoints() - 2).length() <
                          dMinSectionLength,
                      "Segment length is too short after truncation.");
        }
    }
    poly_noShortSections.push_back(endPoint);
    if(IS_DEBUG)
        CHECK(poly_noShortSections.segment((int)poly_noShortSections.numPoints() - 2).length() < dMinSectionLength,
              "Last segment length is too short after truncation.");

    return poly_noShortSections;
}

/*double widthToMidpointThreshold(const double dWidth)
{
    return dWidth / * twice actual distance to polyline edge * / +
params().Correspondence.CANE_CORRESPONDENCE.TIP_TO_MIDPOINT_THRESHOLD;
}*/

/*template<class TControlPoint>
bool isCloseToControlPoint(const typename TControlPoint::TVecType & p, const TControlPoint & controlPoint, const double
dREThresh)
{
    return (controlPoint.getPoint() - p).norm() < dREThresh;// widthToMidpointThreshold(controlPoint.getWidth());
}*/

// If a tip reconstructed w.r.t. a polyline is a midpoint on that polyline then the reconstruction is clearly incorrect.
template <class TControlPoint>
bool CPolyline_base<TControlPoint>::isMidPoint(const TVecType& p, const double dREThresh) const
{
    const bool bVerbose = false;

    const double dPos = positionOfPoint_distance(p);

    const double dLength = length();
    if(dLength <= 2 * dREThresh)
        return false;

    const bool bResult = (dPos > dREThresh) && (dPos < (dLength - dREThresh));

    COUT2("isMidPoint=", bResult);
    COUT(dREThresh);
    COUT(dPos);

    return bResult;
}

template <class TControlPoint>
typename CPolyline_base<TControlPoint>::TPolyline
CPolyline_base<TControlPoint>::splitLongSections(const double dMaxSectionLength, const bool bVerbose) const
{
    // const bool bVerbose = false;

    COUT(dMaxSectionLength);
    COUT(numSegments());

    TPolyline poly_noLongSections;

    if(numPoints() < 2)
        return *this;

    poly_noLongSections.push_back(aControlPoints[0]);
    for(int nSection = 0; nSection < numSegments(); nSection++) {
        const double dLength = segment(nSection).length();
        if(dLength > dMaxSectionLength) {
            const int nNumSegmentsReq = (int)ceil(dLength / dMaxSectionLength);
            const double dNewSegLength = dLength / nNumSegmentsReq;
            CHECK(dNewSegLength > dMaxSectionLength, "dMaxSectionLength failed");
            for(int i = 1; i < nNumSegmentsReq; i++) {
                const double dProp = (double)i / (double)nNumSegmentsReq;
                const typename TControlPoint::TVecType point =
                    (1 - dProp) * aControlPoints[nSection].getPoint() + dProp * aControlPoints[nSection + 1].getPoint();
                poly_noLongSections.push_back(interpolateControlPoint(point, nSection, dProp));
            }
        }
        poly_noLongSections.push_back(aControlPoints[nSection + 1]);
    }

    COUT(poly_noLongSections.numSegments());

    COUT(this->toString());
    COUT(poly_noLongSections.toString());

    CHECK(poly_noLongSections.numSegments() < numSegments(), "Too many splits");

    return poly_noLongSections;
}

template <class TControlPoint> double CPolyline_base<TControlPoint>::averageThickness() const
{
    if(IS_DEBUG)
        CHECK(numPoints() == 0, "0-length polyline");
    double dAvThickness = 0;
    BOOST_FOREACH(const TControlPoint& cp, *this) {
        dAvThickness += cp.getWidth();
    }

    dAvThickness /= numPoints();
    return dAvThickness;
}

template <class TControlPoint> double CPolyline_base<TControlPoint>::maxKinkAngle_int(int& nMaxKinkAnglePos) const
{
    double dMaxKinkAngle = 0;
    for(int i = 1; i < numPoints() - 1; i++) {
        const double dKinkAngle = kinkAngleAtPoint(i);
        if(dKinkAngle > dMaxKinkAngle) {
            dMaxKinkAngle = dKinkAngle;
            nMaxKinkAnglePos = i;
        }
    }
    return dMaxKinkAngle;
}

template <class TControlPoint> double CPolyline_base<TControlPoint>::maxKinkAngle() const
{
    int nMaxKinkAnglePos = -1;
    return maxKinkAngle_int(nMaxKinkAnglePos);
}

template <class TControlPoint>
typename CPolyline_base<TControlPoint>::TVecType CPolyline_base<TControlPoint>::getMaxKinkPoint() const
{
    TVecType maxKinkPoint = getStartPoint();
    double dMaxKinkAngle = 0;
    for(int i = 1; i < numPoints() - 1; i++) {
        const double dKinkAngle = kinkAngleAtPoint(i);
        if(dKinkAngle > dMaxKinkAngle) {
            dMaxKinkAngle = dKinkAngle;
            maxKinkPoint = (*this)[i].getPoint();
        }
    }
    return maxKinkPoint;
}

// template<class TControlPoint>
// double CPolyline_base<TControlPoint>::Curvature() const
//{
//	double curv = 0;
//	for(int i=1; i<numPoints()-1; i++)
//	{
//		double dist = (TPolylineVector[i+1] - TPolylineVector[i]).norm() + (TPolylineVector[i] -
//TPolylineVector[i-1]).norm()
//        curv += 2*(180.0 - kinkAngleAtPoint(i))/dist;
//	}
//	return curv;
//
//}

template <class TControlPoint>
void CPolyline_base<TControlPoint>::checkNotTooKinked(const double dMaxKinkAngleThresh) const
{
    if(maxKinkAngle() > dMaxKinkAngleThresh) {
        ALWAYS_VERBOSE;
        COUT(maxKinkAngle());
        COUT2("Threshold: ", dMaxKinkAngleThresh);

        COUT(toString());

        THROW("Polyline getting too kinked")
    }
}

template <class TControlPoint> void CPolyline_base<TControlPoint>::splitSegment(const int nSegmentToSplit)
{
    const bool bVerbose = true;

    COUT2("Before", maxKinkAngle());
    COUT2("Before", toString());
    COUT(nSegmentToSplit);

    checkNotTooKinked(M_PI / 2);

    if(TControlPoint::HAS_THICKNESS) {
        const double dNewWidth = 0.5 * ((*this)[nSegmentToSplit].getWidth() + (*this)[nSegmentToSplit + 1].getWidth());
        TControlPoint newMidpoint(segment(nSegmentToSplit).midPoint(), dNewWidth);
        aControlPoints.insert(aControlPoints.begin() + nSegmentToSplit + 1, newMidpoint);
    } else {
        const TControlPoint newMidpoint(segment(nSegmentToSplit).midPoint());
        aControlPoints.insert(aControlPoints.begin() + nSegmentToSplit + 1, newMidpoint);
    }

    COUT2("After", maxKinkAngle());
    COUT2("After", toString());

    checkNotTooKinked(M_PI / 2);
}

/* Return what fraction of otherPoly is contained within this
 * */
template <class TControlPoint>
double CPolyline_base<TControlPoint>::contains(const TPolyline& otherPoly, const bool bSlowAndAccurate) const
{
    const bool bVerbose = false;

    if(bSlowAndAccurate) {
        const optional<TPolyline> otherPoly_moreSections =
            otherPoly.truncate(0, otherPoly.averageThickness() * 3, 0, bVerbose);
        if(!otherPoly_moreSections)
            return contains(otherPoly, false); // Truncation failed (probably too short), just use the non-truncated one

        return contains(*otherPoly_moreSections, false);
    }

    Eigen::ArrayXd adLengths = Eigen::ArrayXd::Zero(otherPoly.numSegments());

    Eigen::ArrayXd adSeparations = Eigen::ArrayXd::Zero(otherPoly.numPoints());
    Eigen::ArrayXd adR1 = Eigen::ArrayXd::Zero(otherPoly.numPoints());
    Eigen::ArrayXd adR2 = Eigen::ArrayXd::Zero(otherPoly.numPoints());

    for(int i = 0; i < otherPoly.numPoints(); i++)
    // BOOST_FOREACH(const TControlPoint & cp, otherPoly)
    {
        const TControlPoint& cp = otherPoly[i];
        const TControlPoint closestHere = closestPointAndWidth(cp.getPoint());

        const double dSeparation = (closestHere.getPoint() - cp.getPoint()).norm(); //'d' on mathworld

        adSeparations(i) = dSeparation;
        adR1(i) = 0.5 * cp.getWidth();
        adR2(i) = 0.5 * closestHere.getWidth();

        /*const double dOverlap = min3(cp.getWidth(), closestHere.getWidth(),
        std::max<double>(0.5*(closestHere.getWidth() + cp.getWidth()) - dSeparation, 0));
        adOverlaps(i) = dOverlap;*/

        if(i < otherPoly.numSegments())
            adLengths(i) = otherPoly.segment(i).length();
    }
    COUT(adSeparations);
    COUT(adR1);
    COUT(adR2);

    const Eigen::ArrayXd adOverlaps = Eigen::ArrayXd::Zero(otherPoly.numPoints()).max((adR1 + adR2) - adSeparations);

    // Average overlap along each segment, overestimate because we linearly interpolate down to 0 overlap, when 0 means
    // 'no overlap'
    const Eigen::ArrayXd adOverlaps_min = adOverlaps.min(2 * adR1).min(2 * adR2);
    const Eigen::ArrayXd adAverageOverlap_overestimate =
        0.5 * (adOverlaps_min.head(otherPoly.numSegments()) + adOverlaps_min.tail(otherPoly.numSegments()));

    if(TControlPoint::TVecType::RowsAtCompileTime == 2) {
        const double dOverlap2D = (adAverageOverlap_overestimate * adLengths).sum();
        const double dArea = otherPoly.area();
        COUT(adAverageOverlap_overestimate);
        COUT(adLengths);
        COUT(dOverlap2D);
        COUT(dArea);
        return dOverlap2D / dArea;
    } else {
        // TODO use http://mathworld.wolfram.com/Circle-CircleIntersection.html
        // For now just assume the overlapping volume is a cylinder/has a circular cross-section (sometimes it is--when
        // one encloses the other)
        const double dOverlap3DOverestimate = M_PI * ((0.5 * adAverageOverlap_overestimate).square() * adLengths).sum();
        const double dVolume = otherPoly.volume();

        return dOverlap3DOverestimate / dVolume;
    }
}

template <class TControlPoint> double CPolyline_base<TControlPoint>::area() const
{
    CHECK(TControlPoint::TVecType::RowsAtCompileTime != 2,
          "No area for a 3D polyline (surface area isn't very useful)");

    double dArea = 0;
    for(int i = 0; i < numSegments(); i++) {
        const double dAvThickness = 0.5 * (aControlPoints[i].getWidth() + aControlPoints[i + 1].getWidth());
        dArea += dAvThickness * segment(i).length();
    }
    return dArea;
}

template <class TControlPoint> double CPolyline_base<TControlPoint>::volume() const
{
    // http://en.wikipedia.org/wiki/Frustum#Volume "volume of a circular cone frustum"
    CHECK(TControlPoint::TVecType::RowsAtCompileTime != 3, "No volume for a 2D polyline");

    double dVolume = 0;
    for(int i = 0; i < numSegments(); i++) {
        const double R1 = 0.5 * aControlPoints[i].getWidth();
        const double R2 = 0.5 * aControlPoints[i + 1].getWidth();
        const double dHeronianBase = sqr(R1) + R1 * R2 + sqr(R2);
        dVolume += (M_PI / 3) * dHeronianBase * segment(i).length();
    }
    return dVolume;
}

template <class TThickPoly>
double intersection(const TThickPoly& poly1, const TThickPoly& poly2, const bool bSlowAndAccurate)
{
    return std::min<double>(poly1.contains(poly2, bSlowAndAccurate), poly2.contains(poly1, bSlowAndAccurate));
}
template <class TControlPoint> double CPolyline_base<TControlPoint>::maxDislocationAngle() const
{
    double dMaxDislocation = 0;
    for(int i = 1; i < numSegments() - 1; i++) {
        const double dDislocation = dislocationAngleAtSegment(i);
        if(dDislocation > dMaxDislocation)
            dMaxDislocation = dDislocation;
    }
    return dMaxDislocation;
}

template <class TControlPoint> double CPolyline_base<TControlPoint>::dislocationAngleAtSegment(const int nPoint) const
{
    CHECK(nPoint <= 0 || nPoint >= numSegments() - 1, "Dislocation angle only defined for middle segments");
    const TVecType averageDir = (segment(nPoint - 1).direction() + segment(nPoint + 1).direction());
    const TVecType segmentDir = segment(nPoint).startToFinish();
    return angleBetweenVectors(averageDir, segmentDir, true);
}

template <class TControlPoint>
optional<typename CPolyline_base<TControlPoint>::TPolyline>
CPolyline_base<TControlPoint>::truncate(const double dMinLength,
                                        const double dMaxLength,
                                        const double dEpsilon,
                                        const bool bVerbose) const
{
    //    newTruncate(dMinLength, dMaxLength, dEpsilon, 10, bVerbose); //Just testing
    //#pragma message("TB:  remove newTruncate testing code from here")

    COUT2("Size before", numPoints());
    COUT2("Length before", length());

    TPolyline poly_noShortSections = removeShortSections(dMinLength, bVerbose);

    COUT2("Size after shorten", poly_noShortSections.numPoints());
    COUT2("Length after shorten", poly_noShortSections.length());

    if(poly_noShortSections.length() < dMinLength) {
        COUT("Too short");
        return optional<TPolyline>();
    }

    TPolyline poly_reducedDP = douglasPeucker(poly_noShortSections, dMaxLength, dEpsilon, bVerbose);

    COUT2("Size after", poly_reducedDP.numPoints());
    COUT2("Length after", poly_reducedDP.length());

    if(poly_reducedDP.length() < dMinLength) {
        COUT("Too short");
        return optional<TPolyline>();
    }

    const TPolyline poly_noLongSections = poly_reducedDP.splitLongSections(dMaxLength, bVerbose);
    return poly_noLongSections;
}

/* Find 2 segments separated by a distance of at most dBudScale with a short chain of edges between.
 *
 * We're going to replace the bud by joining at the intersection of the 2 ends
 *
 * if the short chain goes *outside* the angle of the intersection (on *either side* for shallow angles) then replace by
 *joining at the intersection of the 2 ends
 *
 * if the short chain goes *inside* the the line between the 2 sections then
 *    if (separation > dMinLength)
 *       replace with an edge between the 2 ends
 *    else
 *       replace by joining at the intersection of the 2 ends
 *
 * If ambiguous, choose the longest join
 *
 * If length approx 0 delete points
 * **/
template <class TControlPoint>
typename C2dPolyline_base<TControlPoint>::ePolyTruncationJoinType C2dPolyline_base<TControlPoint>::isBud(
    const int nBudSegmentEnd,
    int& nBudEnd,
    typename C2dPolyline_base<TControlPoint>::TPolyline_base::TVecType& interpControlPoint,
    const double dMinSegmentLength,
    const double dBudScale,
    const bool bVerbose) const
{
    nBudEnd = -1;
    if(nBudSegmentEnd == 0 || nBudSegmentEnd + 3 >= this->numPoints())
        return eNoJoin;

    const C2dBoundedLine startLine = this->segment(nBudSegmentEnd - 1);
    const C2dLine startLineUnbounded = startLine.unboundedLine();

    ePolyTruncationJoinType bestJoinType = eNoJoin;

    const double dAngleFromStart =
        signedAngleBetweenUnitVectors(startLineUnbounded.getDirection(), this->segment(nBudSegmentEnd).direction());

    for(int nBudEndSeg = nBudSegmentEnd + 2; nBudEndSeg < this->numSegments(); nBudEndSeg++) {
        const C2dBoundedLine endLine = this->segment(nBudEndSeg);

        C2dBoundedLine lineBetweenStartAndEnd(startLine.getFinishPoint(), endLine.getStartPoint());
        const double dSeparation = lineBetweenStartAndEnd.length();
        if(dSeparation > dBudScale)
            continue; // might be lower later... We'll remove meanders with this bud truncation.

        const C2dLine endLineUnbounded = endLine.unboundedLine();

        optional<const C2dImagePointPx> pIntersection = startLineUnbounded.findIntersection(endLineUnbounded);

        const C2dImagePointPx intersection = pIntersection ? *pIntersection : lineBetweenStartAndEnd.midPoint();
        const double t_start = startLine.get_t(intersection);
        const double t_end = endLine.get_t(intersection);

        // const bool bIntersectionIsInMiddle = (t_start < 1 || t_end > 0);
        const bool bIntersectionIsOutsideRegion = (t_start <= 0 || t_end >= 1);
        if(bIntersectionIsOutsideRegion)
            continue;

        const double dAngleStartToFinish =
            signedAngleBetweenUnitVectors(startLineUnbounded.getDirection(), endLineUnbounded.getDirection());
        const double dAngleToEnd =
            signedAngleBetweenUnitVectors(this->segment(nBudEndSeg - 1).direction(), endLineUnbounded.getDirection());

        // Greedy: overwrite any decisions we made earlier and join as many as possible
        // If the signs are different, eReplaceWithOneEdge (even if the replacement segment is very short)
        if(sign(dAngleFromStart) != sign(dAngleToEnd)) {
            if(!bIntersectionIsOutsideRegion) {
                bestJoinType = eReplaceWithOneEdge;
                nBudEnd = nBudEndSeg;
            }
        } else { // Same sign
            const double nSign = sign(dAngleToEnd);
            if(sign(dAngleStartToFinish) == nSign && dSeparation > dMinSegmentLength) {
                bestJoinType = eReplaceWithOneEdge;
                nBudEnd = nBudEndSeg;
            } else if(!bIntersectionIsOutsideRegion) {
                bestJoinType = eReplaceWithIntersection;
                nBudEnd = nBudEndSeg;
                interpControlPoint = intersection;
            }
        }
    }
    return bestJoinType;
}

/**
 * @brief
 * @param dMinLength
 * @param dMaxLength
 * @param dEpsilon
 * @param dBudScale Remove buds smaller than this. A bud is a chain where the endpoints are closer than this, and where
replacing the endpoints with one intersection point will TODO define.
 * @param bVerbose
 * @return
 *
template<class TControlPoint>
typename C2dPolyline_base<TControlPoint>::T2dPolyline C2dPolyline_base<TControlPoint>::newTruncate(const double
dMinLength, const double dMaxLength, const double dEpsilon, const double dBudScale, const bool bVerbose) const
{
    //Remove *very* short sections
    const C2dPolyline_base<TControlPoint>::TPolyline_base polyNoZeroSegments =
(this->removeShortSections(0.1*dMinLength, bVerbose));
    COUT(polyNoZeroSegments.numPoints());

    const C2dPolyline_base<TControlPoint> polyNoBuds =
((C2dPolyline_base<TControlPoint>&)polyNoZeroSegments).removeBuds(dMinLength, dBudScale, bVerbose);
    COUT(polyNoBuds.numPoints());

    const C2dPolyline_base<TControlPoint> polyNoShortSections = polyNoBuds.removeShortSections_new(dMinLength,
bVerbose);
    COUT(polyNoShortSections.numPoints());

    return polyNoShortSections;
}*/

void test2dCaneEdgeTruncation()
{
    const bool bVerbose = false;

    C2dPolyline polyTest, polyOut;

    polyTest.push_back(C2dImagePointPx(10, 10));
    polyTest.push_back(C2dImagePointPx(10, 20));

    C2dPolyApproxSettings approxSettings0(6, 12, 0, -1, 8);
    CHECK(polyTest.approximate(polyOut, approxSettings0, bVerbose) != eSuccess, "Truncation failed");
    CHECK(polyOut.numPoints() != 2, "Truncated a cane which shouldn't have been truncated");
    C2dPolyApproxSettings approxSettings1(12, 24, 0);
    CHECK(polyTest.approximate(polyOut, approxSettings1, bVerbose) != eFail,
          "Should fail to truncate a cane shorter than min length");
    CHECK(polyOut.numPoints() != 0, "Returned a test2dCcane shorter than min length");

    polyTest.push_back(C2dImagePointPx(12, 20));

    polyTest.push_back(C2dImagePointPx(10, 22));

    polyTest.push_back(C2dImagePointPx(10, 30));

    CHECK(polyTest.approximate(polyOut, approxSettings0, bVerbose) != eSuccess, "Truncation failed");
    CHECK(polyOut.numPoints() != 3, "Didn't remove a bud");

    polyTest.clearAndReserve(10);
    polyTest.push_back(C2dImagePointPx(10, 10));
    polyTest.push_back(C2dImagePointPx(10, 30));
    polyTest.push_back(C2dImagePointPx(10, 30));
    polyTest.push_back(C2dImagePointPx(10, 50));
    C2dPolyApproxSettings approxSettings2(6, 60, 0, 20 /* target length */, 8);
    CHECK(polyTest.approximate(polyOut, approxSettings2, bVerbose) != eSuccess, "Truncation failed");
    CHECK(polyOut.numPoints() != 3, "Didn't remove a duplicate point");

    polyTest.clearAndReserve(10);
    polyTest.push_back(C2dImagePointPx(10, 10));
    polyTest.push_back(C2dImagePointPx(10, 30));
    polyTest.push_back(C2dImagePointPx(10, 30.01));
    polyTest.push_back(C2dImagePointPx(10, 50));
    CHECK(polyTest.approximate(polyOut, approxSettings2, bVerbose) != eSuccess, "Truncation failed");
    CHECK(polyOut.numPoints() != 3, "Didn't remove a near-duplicate point");

    polyTest.clearAndReserve(10);
    polyTest.push_back(C2dImagePointPx(10, 10));
    polyTest.push_back(C2dImagePointPx(10, 30));
    polyTest.push_back(C2dImagePointPx(10, 50));
    C2dPolyApproxSettings approxSettings3(6, 14, 0, -1, 8);
    CHECK(polyTest.approximate(polyOut, approxSettings3, bVerbose) != eSuccess, "Truncation failed");
    CHECK(polyOut.numPoints() != 5, "Error splitting a polyline into shorter sections");
}

void test2dCaneEdgeTruncation2()
{
    const bool bVerbose = false;

    C2dPolyline polyTest, polyOut;

    C2dImagePointPx centre(500, 500);

    for(int nPoints = 0; nPoints < 360; nPoints++) {
        const double dAngle = nPoints * (M_PI / 180);

        C2dPolylineControlPoint cp(centre + 100 * TEigen2dPoint(sin(dAngle), cos(dAngle)));
        polyTest.push_back(cp);
    } // Circumferance 628, section lengths 1.8

    C2dPolyApproxSettings approxSettings(1, 20, 1, 2 /* target length */);
    CHECK(polyTest.approximate(polyOut, approxSettings, bVerbose) != eSuccess, "Truncation failed");
    CHECK(polyOut.numPoints() != 360, "Approximation error");
    CHECK(polyMaxSeparation(polyTest, polyOut) > 1e-8, "Approx failed");

    C2dPolyApproxSettings approxSettings2(1, 20, 10, 3.7 /* target length */);
    CHECK(polyTest.approximate(polyOut, approxSettings2, bVerbose) != eSuccess, "Truncation failed");
    CHECK_P(polyOut.numPoints() != 180,
            polyOut.numPoints(),
            "Approximation error (actually anywhere from 178-182 would be ok)");

    C2dPolyApproxSettings approxSettings3(1, 20, 0.00001, 4 /* target length */);
    CHECK(polyTest.approximate(polyOut, approxSettings3, bVerbose) != eSuccess, "Truncation failed");
    CHECK_P(polyOut.numPoints() != 360, polyOut.numPoints(), "Approximation error");
    CHECK(polyMaxSeparation(polyTest, polyOut) > 1e-8, "Approx failed");

    for(double dEpsilon = 0; dEpsilon < 5; dEpsilon++) {
        C2dPolyApproxSettings approxSettings4(1, 20, dEpsilon, 4 * dEpsilon + 5 /* target length */);
        CHECK(polyTest.approximate(polyOut, approxSettings4, bVerbose) != eSuccess, "Truncation failed");

        const double dMaxSeparation = polyMaxSeparation(polyTest, polyOut);
        COUT(polyOut.numSegments());
        COUT(polyOut.length() / polyOut.numSegments());
        COUT(dMaxSeparation);
        CHECK_P(dMaxSeparation > dEpsilon + 1, dMaxSeparation, "Approx failed (a bit of a trade-off here)");
    }
}

// A section which will be truncated somehow
class CShortSection
{
    int nStart, nEnd;
    double dTotalLength;

public:
    CShortSection(const int nStart, const int nEnd, const double dTotalLength)
        : nStart(nStart)
        , nEnd(nEnd)
        , dTotalLength(dTotalLength)
    {
    }

    void shiftIndices(const int nNumPointsDeleted)
    {
        nStart -= nNumPointsDeleted;
        nEnd -= nNumPointsDeleted;
    }

    double getTotalLength() const
    {
        return dTotalLength;
    }
    int getStart() const
    {
        return nStart;
    }
    int getFinishIdx() const
    {
        return nEnd;
    }
};

// Area between approximated polyline and original polyline
// Should be ordered the same as original, and neither should be self-intersecting with themselves
template <class TControlPoint>
double C2dPolyline_base<TControlPoint>::measureApproximationQuality(
    const CShortSection& shortSection,
    typename C2dPolyline_base<TControlPoint>::T2dPolyline& candidateApproximation,
    const double dMinLength,
    const bool bVerbose) const
{
    // Split into a list of separate areas at intersection points. Sum areas

    typedef std::map<double, C2dImagePointPx> TPointPositionMap;
    TPointPositionMap aPositions;

    for(int nSegmentCandidate = 0; nSegmentCandidate < candidateApproximation.numSegments(); nSegmentCandidate++) {
        C2dBoundedLine seg = candidateApproximation.segment(nSegmentCandidate);

        for(int nOtherSideOnOriginal = shortSection.getStart(); nOtherSideOnOriginal < shortSection.getFinishIdx();
            nOtherSideOnOriginal++) {
            optional<const C2dImagePointPx> pIntersect = this->segment(nOtherSideOnOriginal).intersection(seg);
            if(pIntersect) {
                aPositions[candidateApproximation.positionOfPoint_index(*pIntersect)] = *pIntersect;
            }
        }
    }

    typename C2dPolyline_base<TControlPoint>::T2dPolyline subArea;

    subArea.clearAndReserve(candidateApproximation.numPoints() +
                            (shortSection.getFinishIdx() - shortSection.getStart()) + 3);

    int nStartOnThis = (shortSection.getStart() > 0) ? (shortSection.getStart() - 1) : shortSection.getStart();
    int nEnd = (shortSection.getFinishIdx() < this->numPoints() - 1) ? (shortSection.getFinishIdx() + 1) :
                                                                       shortSection.getFinishIdx();

    if(aPositions.size() == 0) {
        // Normal case: no self-intersection

        // Copy the current polyline into subArea backwards, from a point beyond and a point earlier than this section
        for(int i = nEnd; i >= nStartOnThis; i--)
            subArea.push_back(this->aControlPoints[i]);

        // And copy the candidate approximation in forwards
        subArea.insert(subArea.end(), candidateApproximation.begin(), candidateApproximation.end());

        return subArea.polygonalArea_nonSelfIntersecting();
    }

    aPositions[nEnd] = this->aControlPoints[nEnd].getPoint();
    double dTotalArea = 0;
    int nStartOnApprox = 0;

    BOOST_FOREACH(const TPointPositionMap::value_type& intersectionEnd, aPositions) {
        subArea.clearAndReserve(0);

        int nEndOnThis /* index of point on this before the crossing */ = (int)std::floor(intersectionEnd.first);
        int nEndOnApprox /* index of point on this before the crossing */ =
            (int)std::floor(candidateApproximation.positionOfPoint_index(intersectionEnd.second));

        for(int i = nEndOnThis; i >= nStartOnThis; i--)
            subArea.push_back(this->aControlPoints[i]);

        for(int i = nStartOnApprox; i <= nEndOnApprox; i++)
            subArea.push_back(candidateApproximation[i]);

        dTotalArea += subArea.polygonalArea_nonSelfIntersecting();

        nStartOnThis = nEndOnThis + 1;
        nStartOnApprox = nEndOnApprox + 1;
    }

    return dTotalArea;
}

// http://en.wikipedia.org/wiki/Polygonal_area#Area_and_centroid
template <class TControlPoint> double C2dPolyline_base<TControlPoint>::polygonalArea_nonSelfIntersecting() const
{
    double dTotalArea = 0;
    for(int i = 0; i < this->numPoints(); i++) {
        const C2dImagePointPx& x1 = this->aControlPoints[i].getPoint();
        const C2dImagePointPx& x2 = this->aControlPoints[(i + 1) % this->numPoints()].getPoint();

        dTotalArea += x1.x() * x2.y() - x1.y() * x2.x();
    }

    const double dArea = fabs(0.5 * dTotalArea);

    return dArea;
}

/*template<class TControlPoint>
void C2dPolyline_base<TControlPoint>::addAllCombinations(const CShortSection & shortSection, const double dMinLength,
typename C2dPolyline_base<TControlPoint>::TCandidateApproximations & aApproximations) const
{
    int nLength = 1+(shortSection.getFinishIdx()-shortSection.getStart());
    if(bIncludeEnd)

}*/

/*
 * 1) One or more short sections which together are too short -> delete 1 of the points OR replace with intersection
 *point
 *
 * 2) Several short sections which together are a good length (1 to 2 x minlength) -> approximate with 1 line OR delete
 *intermediate points OR replace with intersection point
 *
 * 3) Several short sections which together are long enough for multiple other sections -> approximate with multiple
 *lines somehow
 *
 * 4) The start and end section could be either moved, or have a point deleted, as well (probably do this first).
 *
 *
 * Approach: find a subsection where edges are all too short, with endpoints which are not too short.
 *
 * For each subsection: shorten
 *
 *
 * **/
template <class TControlPoint>
typename C2dPolyline_base<TControlPoint>::T2dPolyline
C2dPolyline_base<TControlPoint>::decideHowToTruncateShortSection(const CShortSection& shortSection,
                                                                 const double dMinLength,
                                                                 const bool bVerbose) const
{
    C2dPolyline_base<TControlPoint> truncatedPolyline;

    const bool bIsStart = shortSection.getStart() == 0;
    const bool bIsEnd = shortSection.getFinishIdx() == this->numPoints() - 1;

    typedef std::vector<C2dPolyline_base<TControlPoint> > TCandidateApproximations;
    TCandidateApproximations aApproximations;

    if(shortSection.getTotalLength() < dMinLength) // case 1: (never stretches from the start to the end)
    {
        // Either keep 1 point or interpolate a new point
        C2dImagePointPx interpSection;
        C2dBoundedLine prevSection, nextSection;
        C2dLine prevSectionUnbounded, nextSectionUnbounded;
        if(!bIsStart) {
            prevSection = this->segment(shortSection.getStart() - 1);
            prevSectionUnbounded = prevSection.unboundedLine();
        }
        if(!bIsEnd) {
            nextSection = this->segment(shortSection.getFinishIdx());
            nextSectionUnbounded = nextSection.unboundedLine();
        }

        const C2dBoundedLine lineBetweenStartAndEnd(this->aControlPoints[shortSection.getStart()].getPoint(),
                                                    this->aControlPoints[shortSection.getFinishIdx()].getPoint());

        if(!bIsStart && !bIsEnd) {
            optional<const C2dImagePointPx> pIntersection = prevSectionUnbounded.findIntersection(nextSectionUnbounded);

            const C2dImagePointPx intersection = pIntersection ? *pIntersection : lineBetweenStartAndEnd.midPoint();
            const double t_start = prevSection.get_t(intersection);
            const double t_end = nextSection.get_t(intersection);

            // const bool bIntersectionIsInMiddle = (t_start < 1 || t_end > 0);
            const bool bIntersectionIsOutsideRegion = (t_start <= 0 || t_end >= 1);
            if(!bIntersectionIsOutsideRegion) {
                aApproximations.resize(1);
                aApproximations.back().push_back(TControlPoint(intersection));
            }

            // Also add each point as a candidate
            aApproximations.resize(1 + shortSection.getFinishIdx() - shortSection.getStart());
            for(int i = shortSection.getStart(); i <= shortSection.getFinishIdx(); i++) {
                aApproximations[i - shortSection.getStart()].push_back(this->aControlPoints[i]);
            }
        } else // At the start or the end... Either the start|end point or an interpolated start|end point.
        {
            aApproximations.resize(2);

            if(bIsStart) {
                aApproximations[0].push_back(this->aControlPoints[0]);

                TControlPoint cpStartInterp(nextSectionUnbounded.closestPoint(this->aControlPoints[0].getPoint()),
                                            this->aControlPoints[0].getWidth());
                aApproximations[1].push_back(cpStartInterp);
            } else if(bIsEnd) {
                aApproximations[0].push_back(this->aControlPoints[this->numPoints() - 1]);

                TControlPoint cpStartInterp(
                    nextSectionUnbounded.closestPoint(this->aControlPoints[this->numPoints() - 1].getPoint()),
                    this->aControlPoints[this->numPoints() - 1].getWidth());
                aApproximations[1].push_back(cpStartInterp);
            } else {
                THROW("Start and End should be impossible when total length less than min");
            }
        }

    } else // a longer-than-min-length section
    {
        aApproximations.resize(2);
        aApproximations[0].push_back(this->aControlPoints[shortSection.getStart()]);
        aApproximations[1].push_back(this->aControlPoints[shortSection.getFinishIdx()]);

        // Try all combinations of edges (TODO)
        // addAllCombinations(shortSection, dMinLength, aApproximations);
    }

    double dMinError = HUGE;
    C2dPolyline_base<TControlPoint> bestExtension;
    BOOST_FOREACH(const C2dPolyline_base<TControlPoint>& candidateExt, aApproximations) {
        const double dError = candidateExt.polygonalArea_nonSelfIntersecting();
        if(dError < dMinError) {
            bestExtension = candidateExt;
            dMinError = dError;
        }
    }

    return bestExtension;
}

template <class TControlPoint>
typename C2dPolyline_base<TControlPoint>::T2dPolyline
C2dPolyline_base<TControlPoint>::removeShortSections_new(const double dMinLength, const bool bVerbose) const
{
    C2dPolyline_base<TControlPoint> truncatedPolyline;
    truncatedPolyline.clearAndReserve(this->numPoints());

    // Eigen::ArrayXd lengths(this->numPoints());
    std::vector<CShortSection> aShortSectionsToTruncate;

    bool bInShortSection = false;
    int nSectionStart = -1;
    double dSectionCumulativeLength = 0;

    for(int i = 0; i < this->numSegments(); i++) {
        const double dThisLength = this->segment(i).length();
        // lengths(i) = dThisLength;

        if(dThisLength < dMinLength) {
            if(!bInShortSection) {
                bInShortSection = true;
                nSectionStart = i;
                dSectionCumulativeLength = 0;
            }

            dSectionCumulativeLength += dThisLength;
        } else {
            if(bInShortSection) {
                bInShortSection = false;
                aShortSectionsToTruncate.push_back(CShortSection(nSectionStart, i - 1, dSectionCumulativeLength));
            }
        }
    }
    if(bInShortSection) {
        aShortSectionsToTruncate.push_back(CShortSection(nSectionStart, this->numSegments(), dSectionCumulativeLength));
    }

    // int nNumPointsDeleted = 0;

    std::vector<CShortSection>::iterator pShortSections = aShortSectionsToTruncate.begin();

    const bool bNothingToTruncate = (pShortSections == aShortSectionsToTruncate.end());

    for(int i = 0; i < this->numPoints(); i++) // i is index into the *source*
    {
        const int nNumPointsDeleted = i - truncatedPolyline.numPoints();

        if(bNothingToTruncate || i < pShortSections->getStart()) {
            truncatedPolyline.push_back(this->aControlPoints[i]);
        } else {
            i = pShortSections->getFinishIdx(); // skip past this short section
            pShortSections->shiftIndices(nNumPointsDeleted);
            typename C2dPolyline_base<TControlPoint>::T2dPolyline newShortSectionPolyline =
                decideHowToTruncateShortSection(*pShortSections, dMinLength, bVerbose);
            truncatedPolyline.insert(
                truncatedPolyline.end(), newShortSectionPolyline.begin(), newShortSectionPolyline.end());
        }
    }



    return truncatedPolyline;
}

//#pragma message("TB: todo new truncation ")
/*
 * Operations: Delete, Merge, move to optimum [split?]
 *
 * Cost: integrated distance + too long penalty + too short penalty
 *
 *
 * * We have to split the longer sections
 *
 * Take original + add as many points as possible? [and remove duplicates]. Each is a 'split point'. Each section has
length > dMinLength.
 *
 * 'Active' polyline has same number of points, and a binary field with each 'on' or 'off'. After we change it, we
interpolate all the 'off' points
 *
 * Compute error by summing error in each section: http://mathworld.wolfram.com/Quadrilateral.html
 * [if p and q don't cross, find intersection and compute 2 triangle areas]
 *
 *
 *
 *
 *
 * Add all possible points, iterate:
 * * Delete lowest cost
 * * Optimise neighbours
 *
 *
 * **First find buds and remove**
 *
//truncatedPolyline
for(int i=0;i<this->numPoints(); i++)
{

}

return truncatedPolyline;
}
*/

/*void C3dCanePrecursor::setNodeStatus(const CVirtualNodeStatus & newStatus, const eEndpoints end)
{
    CVirtualNodeStatus & status = (end==eStart) ? startNodeStatus : endNodeStatus;
    CHECK(status.nodeStatusAtCaneEnd() != eVineUnknown, "Resetting or changing existing status");
    status = newStatus;
}


//return true if statuses turn out to be incompatible
bool C3dCanePrecursor::setNodeStatusOrReject(const C2dPolylineControlPoint & polyPoint, const CWorldCamera & P, const
eNodeStatus newStatus)
{
    const bool bVerbose=true;
    const C2dImagePointPx & p = polyPoint.getPoint();

    if(IS_DEBUG) CHECK((int)newStatus < 0 || newStatus >= NUM_NODE_STATUSES, "Bad input status");
    if(IS_DEBUG) CHECK((int)startNodeStatus < 0 || startNodeStatus >= NUM_NODE_STATUSES, "Bad start node status");
    if(IS_DEBUG) CHECK((int)endNodeStatus < 0 || endNodeStatus >= NUM_NODE_STATUSES, "Bad end node status");

    const double dDistToStart = (p - P.projectToPx_checked(polyline.getStartPoint(), eCheckIO, eSTVine)).norm();
    const double dDistToEnd = (p - P.projectToPx_checked(polyline.getFinishPoint(), eCheckIO, eSTVine)).norm();

    enum { eNodeIsStart, eNodeIsEnd, eNodeIsNeither } eWhichNode = eNodeIsNeither;

    double dThreshInPx = 2*polyPoint.getWidth(); //There may be 2 tips which reproject to slightly different places, so
allow some tolerence here

    if(dDistToStart < dDistToEnd && (dDistToStart < dThreshInPx))
        eWhichNode = eNodeIsStart;
    else if(dDistToEnd < dThreshInPx)
        eWhichNode = eNodeIsEnd;

    if(newStatus != eVineImageEdge && eWhichNode == eNodeIsNeither) {
        COUT("Rejecting correspondence with tip projecting to middle, or off the edge (should usually project to an
endpoint if the other end is off the edge, because other segment should have been extended appropriately)");
        COUT(dDistToStart);
        COUT(dDistToEnd);
        return true;
    }

    //if(IS_DEBUG) CHECK(newStatus == eVineTip && eWhichNode == eNodeIsNeither, "2D Tip should project exactly to an end
node"); //not a problem, EXCEPT we should have extended the 3D recon out to the tip--'neither' just means that the tip
is off the edge

    COUT(getVineNodeName(newStatus));

    if(eWhichNode == eNodeIsNeither) {
        COUT("Not setting node status--endpoint doesn't project to node");
        COUT(dDistToStart);
        COUT(dDistToEnd);
        return false;
    }

    eNodeStatus & statusToSet = (eWhichNode == eNodeIsStart) ? startNodeStatus : endNodeStatus;

    COUT(getVineNodeName(statusToSet));

    const eNodeStatus otherStatus = (eWhichNode == eNodeIsStart) ? endNodeStatus : startNodeStatus;
    CHECK(otherStatus == eVineTip && newStatus == eVineTip, "Both ends of cane are set as tips");

    if(statusToSet == eVineUnknown || statusToSet == newStatus || statusToSet == eVineImageEdge) {
        statusToSet = newStatus;
        COUT2("Changed status to ", getVineNodeName(statusToSet));
        return false;
    }

    //Don't want to change tip (or similar) to edge
    if(newStatus == eVineImageEdge) { //e.g. edge and tip correspondence at end
        COUT("Not changing previously set status to edge status");
        return false;
    }

    COUT(newStatus);
    COUT(statusToSet);
    THROW("Bad status change");
}*/

// For drawing polyline segments
template <class TControlPoint> double getLineWidth(const TControlPoint& controlPoint, const int nThickness)
{
    return clip<double>(
        controlPoint.getWidth_default(nThickness), 1.0, 255.0 /* stop opencv complaining about thick lines>255 */);
}
template <class TControlPoint> int getCircleRad(const TControlPoint& controlPoint, const int nThickness)
{
    return clip<int>(cvRound(0.5 * controlPoint.getWidth_default(nThickness)), 1, 1000);
}

template <class TControlPoint>
void drawSegment(cv::Mat& M,
                 const cv::Scalar colour,
                 const int nThickness,
                 const TControlPoint& start,
                 const TControlPoint& end)
{
    CHECK(nThickness == DRAW_OUTLINE, "shouldn't be calling this fn for outlines");
    
    //CHECK(start.getPoint().x() < 0 && end.getPoint().x() > 2000 || end.getPoint().x() < 0 && start.getPoint().x() > 2000, "Cane projection error");

    const double dStartWidth = getLineWidth(start, nThickness);
    const double dEndWidth = getLineWidth(end, nThickness);

    if(start.getPoint() != end.getPoint()) {
        C2dBoundedLine seg(start.getPoint(), end.getPoint());

        //SHOW(M);
        if(fabs(dStartWidth - dEndWidth) < 3 && dStartWidth < 25)
            seg.draw(M, cvRound(0.5 * (dStartWidth + dEndWidth)), colour);
        else {
            const int nNumPoints = 4;
            cv::Point aPoints[nNumPoints];

            const TEigen2dPoint perp = seg.perpendicular();

            aPoints[0] = C2dImagePointPx(start.getPoint() + 0.5 * dStartWidth * perp);
            aPoints[1] = C2dImagePointPx(end.getPoint() + 0.5 * dEndWidth * perp);
            aPoints[2] = C2dImagePointPx(end.getPoint() - 0.5 * dEndWidth * perp);
            aPoints[3] = C2dImagePointPx(start.getPoint() - 0.5 * dStartWidth * perp);

            cv::fillConvexPoly(M, aPoints, nNumPoints, colour);
        }
    }
    //SHOW(M);
}

// template-compatible method for debugging -- no frameId = latest frame (3D)
// template<class TControlPoint>
// void C2dPolyline_base<TControlPoint>::show(std::string label, CFrameId frameId) const
/*template <typename TPolyType> void show2dPolyline(const TPolyType& polyline, std::string label, CFrameId frameId)
{
    ALWAYS_VERBOSE;
    if(frameId == CFrameId())
        frameId = CFrameId(scene().getLatestProcessedFrameTime(), sys().eCameraAtStereoOrigin);

    MAT2(debugPoly, frameId);
    polyline.draw(debugPoly, colour("2DCane"), true);
    SHOW2(debugPoly, label);
}

// template-compatible method for debugging -- no frameId = latest frame (3D)
// template<class TControlPoint>
// void C3dPolyline_base<TControlPoint>::show(std::string label, CFrameId frameId) const
template <typename TPolyType> void show3dPolyline(const TPolyType& polyline, std::string label, CFrameId frameId)
{
    ALWAYS_VERBOSE;
    if(frameId == CFrameId())
        frameId = CFrameId(scene().getLatestInitPositionTime(), sys().eCameraAtStereoOrigin);

    MAT2(debugPoly3D, frameId);
    polyline.draw(debugPoly3D, scene().getCamera(frameId), colour("3DCane"), true);
    SHOW2(debugPoly3D, label);
}*/

template <class TControlPoint>
void C2dPolyline_base<TControlPoint>::drawLine(cv::Mat& M, const cv::Scalar colour, const int nThickness) const
{
    if(nThickness == DRAW_OUTLINE)
        drawOutline(M, colour, 1);
    else {
        for(int i = 0; i < this->numSegments(); i++) {
            drawSegment(M, colour, nThickness, (*this)[i], (*this)[i + 1]);
            if(i > 0 && i < this->numSegments())
                cv::circle(M,
                           this->aControlPoints[i].getPoint(),
                           getCircleRad(this->aControlPoints[i], nThickness),
                           0.5*colour,
                           -1);
        }
    }

    const bool bMasking = (M.channels() == 1);
    if(this->numPoints() > 0 && !bMasking)
        cv::circle(M, this->getStartPoint(), getCircleRad(this->aControlPoints[0], nThickness), 0.5 * colour);
}

template <class TControlPoint>
void C2dPolyline_base<TControlPoint>::drawOutline(cv::Mat& M, const cv::Scalar col, const int nOutlineThickness) const
{
    CHECK_MAT_INIT(M);
    const int nNumPoints = 2 * this->numPoints();
    boost::scoped_array<cv::Point> aPoints(new cv::Point[nNumPoints]);
    for(int i = 0; i < this->numPoints(); i++) {

        const int nSegBefore = std::max<int>(i - 1, 0);
        const int nSegAfter = std::min<int>(i, this->numSegments() - 1);
        C2dImagePointPx perpAtJoin =
            this->segment(nSegBefore).perpendicular() + this->segment(nSegAfter).perpendicular();
        const double dNorm = perpAtJoin.norm();
        if(dNorm > 0)
            perpAtJoin /= dNorm;
        else {
            perpAtJoin = this->segment(nSegBefore).perpendicular();
            if(perpAtJoin.norm() == 0) {
                cout << "Warning: polyline has 2 0-length segments\nthis=" << this << endl;
                return;
            }
        }

        perpAtJoin *= 0.5 * this->aControlPoints[i].getWidth();

        const C2dImagePointPx p1 = (this->aControlPoints[i].getPoint() + perpAtJoin).eval();
        const C2dImagePointPx p2 = (this->aControlPoints[i].getPoint() - perpAtJoin).eval();
        aPoints[i] = p1;
        aPoints[nNumPoints - (i + 1)] = p2;
    }

    const cv::Point* pCorners = aPoints.get();

    cv::polylines(M, &pCorners, &nNumPoints, 1, true, col, nOutlineThickness);
}

template <class TControlPoint>
void C2dPolyline_base<TControlPoint>::draw(cv::Mat& M,
                                           const cv::Scalar colour,
                                           const bool bDrawDetails,
                                           const int nThickness) const
{
    drawLine(M, colour, nThickness);
    if(bDrawDetails)
        drawDetails(M);
}


template <class TControlPoint>
void C3dPolyline_base<TControlPoint>::drawDetails(cv::Mat& M, const CWorldCamera& P) const
{
    cout << "Drawing 3D polyline with details:" << endl;
    for(int i = 0; i < this->numPoints(); i++) {

        optional<const C2dImagePointPx> pp = P.projectToPx(this->aControlPoints[i].getPoint());

        std::ostringstream ss;
        ss << "i=" << i << " x=(" << this->aControlPoints[i].getPoint().transpose() << ")";

        if(TControlPoint::HAS_THICKNESS)
            ss << " thickness=" << this->aControlPoints[i].getWidth();

        if(i > 0 && i < (this->numPoints() - 1)) {
            const double dLength0 = this->segment(i - 1).length();

            if(this->segment(i).length() > 0 && dLength0 > 0)
                ss << " kink angle=" << this->kinkAngleAtPoint(i);
        }
        if(i < this->numSegments())
            ss << " length=" << this->segment(i).length();

        if(pp) {
            const C2dImagePointPx& p2d = *pp;
            drawText(M, p2d, ss.str());
        }

        cout << ss.str() << endl;
    }
}

template <class TControlPoint> void C2dPolyline_base<TControlPoint>::drawDetails(cv::Mat& M) const
{
    CHECK_MAT_INIT(M);

    cout << "Drawing 2D polyline with details:" << endl;
    for(int i = 0; i < this->numPoints(); i++) {

        const C2dImagePointPx& p2d = this->aControlPoints[i].getPoint();
        std::ostringstream ss;
        ss << "i=" << i << " x=(" << p2d.transpose() << ")";

        if(TControlPoint::HAS_THICKNESS)
            ss << " thickness=" << this->aControlPoints[i].getWidth();

        if(i > 0 && i < (this->numPoints() - 1))
            ss << " kink angle=" << this->kinkAngleAtPoint(i);

        drawText(M, p2d, ss.str());

        cout << ss.str() << endl;
    }
}

template <class TControlPoint>
template <class TProjectToPolyType>
void C3dPolyline_base<TControlPoint>::drawLine(cv::Mat& M,
                                               const CWorldCamera& P,
                                               const cv::Scalar colour,
                                               const int nThickness) const
{
    CHECK_MAT_INIT(M);

    const bool bVerbose = true;

    TProjectToPolyType p2;
    if(project(P, p2) == eFail) {
        REPEAT(10, COUT2("WARNING (3D Polyline draw): projection failed drawing polyline", this->toString()));
        return;
    }
    
    //try{
    p2.draw(M, colour, false, nThickness);
        
    /*} catch(CException & ex)
    {
        COUT("Rendering/projection error");
        COUT(p2.toString());
        COUT(this->toString());
        COUT(P);
        COUT(colour);
    }
    */
}

void C3dPolylineWithThickness::draw(cv::Mat& M,
                                    const CWorldCamera& P,
                                    const cv::Scalar colour,
                                    const bool bDrawDetails,
                                    const int nThickness) const
{
    drawLine<C2dPolylineWithThickness>(M, P, colour, nThickness);
    if(bDrawDetails)
        drawDetails(M, P);
}

void C3dPolyline::draw(cv::Mat& M,
                       const CWorldCamera& P,
                       const cv::Scalar colour,
                       const bool bDrawDetails,
                       const int nThickness) const
{
    drawLine<C2dPolyline>(M, P, colour, nThickness);
    if(bDrawDetails)
        drawDetails(M, P);
}

template <class TControlPoint> std::string CPolyline_base<TControlPoint>::toString() const
{
    std::ostringstream ss;
    ss << "Polyline size " << numPoints() << " length " << length() << " max kink angle=" << maxKinkAngle() << endl;
    for(int i = 0; i < numPoints(); i++) {
        ss << i << ": " << aControlPoints[i];

        if(i > 0 && i < (this->numPoints() - 1))
            ss << " kink angle=" << this->kinkAngleAtPoint(i);

        if(i < this->numSegments())
            ss << " length="
               << (aControlPoints[i + 1].getPoint() - aControlPoints[i].getPoint())
                      .norm(); // this->segment(i).length();

        ss << endl;
    }
    return ss.str();
}

std::ostream& operator<<(std::ostream& s, const C3dPolylineControlPointWithThickness& X)
{
    return X.pp(s);
}

std::ostream& operator<<(std::ostream& s, const C2dPolylineControlPointWithThickness& X)
{
    return X.pp(s);
}

std::ostream& operator<<(std::ostream& s, const C3dPolylineControlPoint& X)
{
    return X.pp(s);
}

std::ostream& operator<<(std::ostream& s, const C2dPolylineControlPoint& X)
{
    return X.pp(s);
}

std::ostream& operator<<(std::ostream& s, const C2dPolyline& X)
{
    return s << X.toString();
}

std::ostream& operator<<(std::ostream& s, const C2dPolylineWithThickness& X)
{
    return s << X.toString();
}

std::ostream& operator<<(std::ostream& s, const C3dPolyline& X)
{
    return s << X.toString();
}

std::ostream& operator<<(std::ostream& s, const C3dPolylineWithThickness& X)
{
    return s << X.toString();
}

void testPolylineClosestPoint()
{
    const bool bVerbose = false;

    C2dPolylineWithThickness poly;

    for(int i = 0; i < 5; i++) {
        poly.push_back(
            C2dPolylineControlPointWithThickness((10 * C2dImagePointPx(i, sqr(i - 2))).eval(), 2 /* width */));
    }

    for(int i = 0; i < (int)poly.numPoints(); i++) {
        const C2dImagePointPx p = poly[i].getPoint();
        COUT(poly.positionOfPoint_index(p));
        CHECK(poly.positionOfPoint_index(p) != i, "Position failed");

        CHECK(!zero((poly.closestPoint(p) - p).norm()), "Closest point failed");
        CHECK(!zero((poly.closestPointAndWidth(p).getPoint() - p).norm()), "Closest point failed");
        CHECK(!zero((*(poly.closestPointStrictlyOnPoly(p)) - p).norm()), "Closest point failed");

        const C2dImagePointPx p_below = (p + C2dImagePointPx(0, -1)).eval();

        const double dDist1 = (poly.closestPoint(p_below) - p_below).norm();
        CHECK(dDist1 <= 0 || dDist1 > 1, "Dist to poly failed");
        const double dDist2 = (*poly.closestPointStrictlyOnPoly(p_below) - p_below).norm();
        CHECK(dDist2 <= 0 || dDist2 > 1, "Dist to poly failed");
        CHECK(dDist1 != dDist2, "Dist to poly failed");

        const double dDist3 = (poly.closestPointAndWidth(p_below).getPoint() - p_below).norm();
        CHECK(dDist1 != dDist3, "Dist to poly failed");
        CHECK(poly.closestPointAndWidth(p_below).getWidth() != 2, "compute width at point failed");
    }

    for(int i = 0; i < (int)poly.numPoints() - 1; i++) {
        const C2dImagePointPx pointAt25 = (0.75 * poly[i].getPoint() + 0.25 * poly[i + 1].getPoint()).eval();
        COUT(poly.positionOfPoint_index(pointAt25));
        CHECK(!zero(poly.positionOfPoint_index(pointAt25) - (i + 0.25)), "Position failed");

        CHECK(!zero((poly.closestPoint(pointAt25) - pointAt25).norm()), "Closest point failed");
        CHECK(!zero((*(poly.closestPointStrictlyOnPoly(pointAt25)) - pointAt25).norm()), "Closest point failed");
    }

    C2dImagePointPx offFrontBack[2];
    offFrontBack[0] = 2 * poly[0].getPoint() - poly[1].getPoint();
    offFrontBack[1] = 2 * poly[poly.numPoints() - 1].getPoint() - poly[poly.numPoints() - 2].getPoint();

    // const double adDists[2] = { poly.segment(0).length()

    for(int i = 0; i < 2; i++) {
        CHECK(poly.closestPoint(offFrontBack[i]) != poly[i * (poly.numPoints() - 1)].getPoint(),
              "Closest point failed");
        CHECK(poly.closestPointStrictlyOnPoly(offFrontBack[i]), "Closest point strictly on poly should return nothing");
    }

    COUT("Polyline closest point test complete");
}

void testPolylineTruncation(const double dThickness)
{
    const bool bVerbose = false;
    COUT(dThickness);

    C3dPolylineWithThickness poly;
    const double root3inv = pow(1.0 / 3.0, 0.5), root2inv = sqrt(0.5);
    TEigen3dPoint ones_start(0, 0, 1.0);
    TEigen3dPoint ones_dir(root3inv, root3inv, root3inv), onesPerp(root2inv, -root2inv, 0);
    // Should only truncate to nothing if total length is about the same as the width
    poly.push_back(C3dPolylineControlPointWithThickness(ones_start, dThickness));
    poly.push_back(C3dPolylineControlPointWithThickness(ones_start + ones_dir * dThickness, dThickness));

    CHECK(poly.truncate_defaults(bVerbose), "Truncation should return nothing for a poly as short as it is wide");

    poly.push_back(C3dPolylineControlPointWithThickness((ones_start + ones_dir * (2 * dThickness)).eval(), dThickness));
    CHECK(!poly.truncate_defaults(bVerbose),
          "Truncation should return something for a poly twice as long as it is wide");
    CHECK_P(poly.truncate_defaults(bVerbose)->numPoints() != 2,
            poly.truncate_defaults(bVerbose)->toString(),
            "Truncation should return something for a poly twice as short as it is wide");

    poly.push_back(
        C3dPolylineControlPointWithThickness(ones_start + ones_dir * 0.1 + onesPerp * dThickness * 0.25, dThickness));
    poly.push_back(C3dPolylineControlPointWithThickness(ones_start + ones_dir * 0.2, dThickness));

    CHECK(poly.truncate_defaults(bVerbose)->numPoints() != 2, "DP should truncate everything for deviation this small");

    poly.push_back(
        C3dPolylineControlPointWithThickness(ones_start + ones_dir * 0.3 + onesPerp * dThickness * 2, dThickness));
    poly.push_back(C3dPolylineControlPointWithThickness(ones_start + ones_dir * 0.4, dThickness));

    CHECK(poly.truncate_defaults(bVerbose)->numPoints() != 4, "DP should not truncate kink for deviation this large");
}

void testPolylineIntersection(const bool bSlowAndAccurate)
{
    C3dPolylineWithThickness poly;
    C3dPolylineWithThickness poly_short;
    C3dPolylineWithThickness poly_thin;
    C3dPolylineWithThickness poly_disjoint;
    for(int i = 0; i < 4; i++) {
        const C3dPolylineControlPointWithThickness cp(C3dWorldPoint(0.1 * i, 0.01 * sqr(i), 1), 0.02);
        poly.push_back(cp);

        if(i == 1 || i == 2)
            poly_short.push_back(cp);

        const C3dPolylineControlPointWithThickness cp_thin(cp.getPoint(), 0.01);
        poly_thin.push_back(cp_thin);

        const C3dPolylineControlPointWithThickness cp_disjoint(C3dWorldPoint(0.1 * i, 0.01 * sqr(i), 1.5), 0.01);
        poly_disjoint.push_back(cp_disjoint);
    }

    CHECK(!zero(poly.length() * M_PI * sqr(0.01) - poly.volume()), "Poly volume check");

    CHECK(!zero(poly_disjoint.volume() - poly.volume()), "Poly volume check");

    CHECK(!zero(4 * poly_thin.volume() - poly.volume()), "Poly volume check");

    CHECK(!zero(intersection(poly, poly, bSlowAndAccurate) - 1), "Poly intersection with itself should be 1");
    CHECK(!zero(intersection(poly_short, poly_short, bSlowAndAccurate) - 1),
          "Poly intersection with itself should be 1");
    CHECK(!zero(intersection(poly_thin, poly_thin, bSlowAndAccurate) - 1), "Poly intersection with itself should be 1");

    CHECK(!zero(intersection(poly_disjoint, poly, bSlowAndAccurate)),
          "Poly intersection with poly_disjoint should be 0");
    CHECK(!zero(intersection(poly_disjoint, poly_short, bSlowAndAccurate)),
          "Poly intersection with poly_disjoint should be 0");
    CHECK(!zero(intersection(poly_disjoint, poly_thin, bSlowAndAccurate)),
          "Poly intersection with poly_disjoint should be 0");

    CHECK(!zero(intersection(poly, poly_thin, bSlowAndAccurate) - 0.25),
          "Poly intersection with poly_thin should be 0.25");
    CHECK(!zero(poly_thin.contains(poly, bSlowAndAccurate) - 0.25), "Poly intersection with poly_thin should be 0.25");
    CHECK(!zero(poly.contains(poly_thin, bSlowAndAccurate) - 1), "Poly should contain poly_thin");

    CHECK(!zero(poly.contains(poly_short, bSlowAndAccurate) - 1), "Poly should contain poly_short");
    const double dVolofpolyInPoly_short = poly_short.contains(poly, bSlowAndAccurate);
    CHECK_P(dVolofpolyInPoly_short < 0.1 || dVolofpolyInPoly_short > 0.5,
            dVolofpolyInPoly_short,
            "Poly_short should overlap poly");

    CHECK(!zero(poly_thin.contains(poly_short, bSlowAndAccurate) - 0.25), "Poly should contain 25% of poly_short");
    const double dVolOfPoly_thinInPoly_short = poly_short.contains(poly_thin, bSlowAndAccurate);
    CHECK_P(!within(dVolOfPoly_thinInPoly_short, dVolofpolyInPoly_short, 0.15),
            dVolOfPoly_thinInPoly_short,
            "Poly_thin should overlap poly");
}

void testFastPolylineDist()
{

    
    C2dPolylineWithThickness poly;

    CRandom::fast_srand(1);

    for(int i = 0; i < 8; i++)
        poly.push_back(C2dPolylineControlPointWithThickness(C2dImagePointPx(100 + 20 * i, 500 + sqr(i)), 5));

    CPolylinePrecomputed_base<C2dPolylineWithThickness> precomputedPoly(poly);

    for(int i = 0; i < 10000; i++) {
        const bool bRecompute = (i % 10 == 5);
        const int nToRecompute = i % poly.numPoints();

        if(bRecompute) {
            poly[nToRecompute].getPoint().x() += CRandom::Uniform(-0.5, 0.5);
            precomputedPoly.update(nToRecompute);

        }

        C2dImagePointPx x(CRandom::Uniform(0, 200), CRandom::Uniform(400, 600));
        C2dImagePointPx x_fast = precomputedPoly.closestPoint(x);
        C2dImagePointPx x_old = poly.closestPoint(x);

        CHECK(!x_fast.isApprox(x_old), "Point-polyline distance computation failed");
    }
}

void testIsMidPoint()
{
    C2dPolylineWithThickness poly;
    for(int i = 0; i < 8; i++)
        poly.push_back(C2dPolylineControlPointWithThickness(C2dImagePointPx(100 + 20 * i, 500 + sqr(i)), 5));

    CHECK(poly.isMidPoint(poly[0].getPoint(), 20), "testIsMidPoint");
    CHECK(!poly.isMidPoint(poly[1].getPoint(), 20), "testIsMidPoint");
    CHECK(!poly.isMidPoint(poly[6].getPoint(), 20), "testIsMidPoint");
    CHECK(poly.isMidPoint(poly[7].getPoint(), 20), "testIsMidPoint");

    CHECK(poly.isMidPoint(poly[0].getPoint(), 0.1), "testIsMidPoint");
    CHECK(poly.isMidPoint(poly[1].getPoint(), 30), "testIsMidPoint");
    CHECK(poly.isMidPoint(poly[6].getPoint(), 40), "testIsMidPoint");

    CHECK(poly.isMidPoint(poly[4].getPoint(), 0.5 * poly.length()), "testIsMidPoint");

    CHECK(poly.isMidPoint(C2dImagePointPx(80, 500), 20), "testIsMidPoint");
    CHECK(poly.isMidPoint(C2dImagePointPx(400, 500), 20), "testIsMidPoint");
}

void testPolyline()
{
    testPolylineClosestPoint();

    for(double dThickness = 0.005; dThickness < 0.03; dThickness += 0.005)
        testPolylineTruncation(dThickness);

    testPolylineIntersection(false);
    testPolylineIntersection(true);

    testAreaComputation();

    test2dCaneEdgeTruncation();
    test2dCaneEdgeTruncation2();

    testFastPolylineDist();

    testIsMidPoint();
}

template <class TControlPoint> void CPolyline_base<TControlPoint>::addLocalisationErrors(const double dGaussianNoiseSD)
{
    BOOST_FOREACH(TControlPoint& p, aControlPoints) {
        p.addLocalisationErrors(dGaussianNoiseSD);
    }
}

template <class TVecType> void addLocalisationErrors_int(const double dGaussianNoiseSD, TVecType& x)
{
    TVecType vec;
    for(int i = 0; i < vec.size(); i++)
        vec(i) = CRandom::Normal(0, dGaussianNoiseSD);
    x += vec;
}

template <class TLineType_in>
void CPolylineControlPointWithThickness<TLineType_in>::addLocalisationErrors(const double dGaussianNoiseSD)
{
    addLocalisationErrors_int<TVecType>(dGaussianNoiseSD, x);
}

template <class TLineType_in>
void CPolylineControlPoint<TLineType_in>::addLocalisationErrors(const double dGaussianNoiseSD)
{
    addLocalisationErrors_int<TVecType>(dGaussianNoiseSD, x);
}

/**
 * @class CSort2dPolyBy2dPoly
 * @brief Defines an order on points by the position of their closest point on a 2D polyline
 */

template <class T2dPolyline>
class CSort2dPolyBy2dPoly : public std::binary_function<const typename T2dPolyline::TControlPoint&,
                                                        const typename T2dPolyline::TControlPoint&,
                                                        bool>
{
    typedef typename T2dPolyline::TControlPoint TControlPoint;
    const T2dPolyline& poly2d;

public:
    CSort2dPolyBy2dPoly(const T2dPolyline& poly2d)
        : poly2d(poly2d)
    {
    }

    double position(const TControlPoint& a) const
    {
        const double dPos = poly2d.positionOfPoint_index(a.getPoint());

        return dPos;
    }

    bool operator()(const TControlPoint& a, const TControlPoint& b) const
    {
        if(position(a) < position(b))
            return true;
        else if(position(a) > position(b) || a.getPoint() == b.getPoint())
            return false;

        // We (very rarely) get 2 endpoints and one point perpendicular to the endpoint (so it has the same
        // positionOfPoint_index). Following the sort, the repeated points are not adjacent. Solution: Add a hash to the
        // sort to ensure that repeated points are always adjacent when sorted.
        // Points are not equal but 'position' is the same
        return vectorHash(a.getPoint()) < vectorHash(b.getPoint());
    }

private:
    void drawOne(cv::Mat& test2dPolySorting, const T2dPolyline& cane1) const
    {
        BOOST_FOREACH(const TControlPoint& p, cane1) {
            drawNumber(test2dPolySorting, C2dImagePointPx(p.getPoint()), position(p));
        }
    }

public:
    void test(const T2dPolyline& caneToSort) const
    {
        BOOST_FOREACH(const TControlPoint& p, caneToSort) {
            cout << p << " pos=" << position(p) << endl;
        }
    }
};

template <class TControlPoint>
void C2dPolyline_base<TControlPoint>::sortByPolyline(
    const typename C2dPolyline_base<TControlPoint>::T2dPolyline& polyToSortWRT)
{
    CSort2dPolyBy2dPoly<typename C2dPolyline_base<TControlPoint>::T2dPolyline> sort2d(polyToSortWRT);

    // sort2d.test(*this);
    // sort2d.test(polyToSortWRT);

    std::sort(this->aControlPoints.begin(), this->aControlPoints.end(), sort2d);

    // sort2d.test(*this);
}

template <class TControlPoint>
void CPolyline_base<TControlPoint>::splitAtPoint(const typename CPolyline_base<TControlPoint>::TVecType& p,
                                                 CPolyline_base<TControlPoint>* aPolylines) const
{
    TControlPoint cpAtJoin = closestPointAndWidth(p);

    const double dPos = positionOfPoint_index(cpAtJoin.getPoint());
    const double eps = 0.01; // 1% of a segment length

    aPolylines[eFinish].push_back(cpAtJoin);

    for(int i = 0; i < numPoints(); i++) {
        double diPos = (double)i;
        if(diPos < dPos - eps)
            aPolylines[eStart].push_back(this->aControlPoints[i]);
        else if(diPos > dPos + eps)
            aPolylines[eFinish].push_back(this->aControlPoints[i]);
    }

    aPolylines[eStart].push_back(cpAtJoin);
    const double dL1 = aPolylines[eStart].length();
    const double dL2 = aPolylines[eFinish].length();
    if(dL1 > 0 && dL2 > 0) // otherwise because of the eps they don't exactly sum to 1
    {
        if(fabs(dL1 + dL2 - length()) > 2.1 * eps * segment((int)std::floor(dPos)).length()) {
            const bool bVerbose = true;COUT(aPolylines[eStart].toString());COUT(aPolylines[eFinish].toString());COUT(this->toString());COUT(cpAtJoin);COUT(dPos);
            COUT(dL1);
            COUT(dL2);
            COUT(length());
            COUT(eps);
            THROW("Polyline parts don't sum to the same length");
        }
    }
}

template <class TControlPoint>
void CPolyline_base<TControlPoint>::subPolyline(const CRange& indexRange, TPolyline& subPoly) const
{
    subPoly.clearAndReserve(numPoints());

    const bool bVerbose = false;

    COUT2("subPolyline", this->toString());
    COUT(indexRange);

    if(!indexRange.isInit())
        return;

    if(indexRange.getMin() > 0 && fabs(indexRange.getMin() - (double)cvRound(indexRange.getMin())) > 0.01) {
        subPoly.push_back(controlPointAtPosition(indexRange.getMin()));
    }

    for(int i = 0; i < numPoints(); i++) {
        if(indexRange.within_absolute((double)i, 0))
            subPoly.push_back(aControlPoints[i]);
    }

    if(indexRange.getMax() < ((double)numPoints() - 1.0) &&
       fabs(indexRange.getMax() - (double)cvRound(indexRange.getMax())) > 0.01) {
        subPoly.push_back(controlPointAtPosition(indexRange.getMax()));
    }

    COUT(subPoly.toString());
}
// Find a polyline close to this one with maxKinkAngle below dNewKinkAngle;
template <class TControlPoint> void CPolyline_base<TControlPoint>::doubleUp()
{
    double dLength = length();

    const int numPointsNew = numPoints() * 2 - 1;

    TPolylineVector aNewPoints(numPointsNew);

    for (int nPoint = 0; nPoint<numPoints(); nPoint++)
    {
        aNewPoints[nPoint * 2] =  aControlPoints[nPoint].scale(2);
    }

    for (int nPoint = 1; nPoint<numPointsNew - 1; nPoint += 2)
    {
        if (TControlPoint::HAS_THICKNESS)
        {
            aNewPoints[nPoint] = TControlPoint(0.5*(aNewPoints[nPoint - 1].getPoint() + aNewPoints[nPoint + 1].getPoint()), 0.5*(aNewPoints[nPoint - 1].getWidth() + aNewPoints[nPoint + 1].getWidth()));
        }
        else
        {
            aNewPoints[nPoint] = TControlPoint(0.5*(aNewPoints[nPoint - 1].getPoint() + aNewPoints[nPoint + 1].getPoint()));
        }
    }

    aControlPoints = aNewPoints;

    CHECK_P(!within(dLength * 2, length(), 0.00001), length(), "Length upscale check failed");
}



// Find a polyline close to this one with maxKinkAngle below dNewKinkAngle;
template <class TControlPoint> void CPolyline_base<TControlPoint>::dropKinkAngle(const double dNewKinkAngle)
{
    bool bVerbose = false;

    for(;;) {
        int nNewMKIndex = -1;
        CHECK(dNewKinkAngle <= 0, "Need dNewKinkAngle>0");
        for(int nIter = 0; nIter < 200; nIter++) {
            int nMaxKinkAnglePos = -1;
            const double dMaxKinkAngle = maxKinkAngle_int(nMaxKinkAnglePos);

            if(dNewKinkAngle > dMaxKinkAngle)
                return;

            TControlPoint& cp = aControlPoints[nMaxKinkAnglePos];

            const TLineType lineAround(aControlPoints[nMaxKinkAnglePos - 1].getPoint(),
                                       aControlPoints[nMaxKinkAnglePos + 1].getPoint());
            TVecType zeroAnglePoint = lineAround.closestPoint(cp.getPoint());
            // zeroAnglePoint can be far too close to one end for near-right angle corners. Average it with the midpoint
            // to keep it sensible
            zeroAnglePoint = 0.5 * (zeroAnglePoint + lineAround.midPoint()).eval(/*eval needed*/);

            const double dShift = /* The 0.975 prevents steps being so small that it never converges */ 0.975 *
                                  dNewKinkAngle / dMaxKinkAngle;
            const TVecType newPoint = (1 - dShift) * zeroAnglePoint + dShift * cp.getPoint();

            cp.getPoint() = newPoint;

            const double dNewMKA = maxKinkAngle_int(nNewMKIndex);
            CHECK(dNewMKA > dMaxKinkAngle && nNewMKIndex == nMaxKinkAnglePos,
                  "Max kink angle increased (without being transferred to a different control point)");

            // if(bVerbose) show3dPolyline<C3dPolyline>(*this, ::toString(dNewMKA) + '-' + ::toString(nNewMKIndex));

            COUT2("After an iteration dropping kink angle", toString());
        }

        erase(begin() + nNewMKIndex,
              begin() + nNewMKIndex + 1); // Should always converge to a straight line if necessary (not ideal though)
        // if(bVerbose) THROW("Failed to truncate polyline");

        // bVerbose = true;
    }
    // cout << "Polyline which can't be truncated: " << toString() << endl;
    // THROW("Failed to truncate polyline");
}

template <class TControlPoint>
const double CPolyline_base<TControlPoint>::closestCollisionDistanceToPoly(const TPolyline& otherPoly) const
{
    // This is not perfect because it ignores thicknesses until later, but should be close enough for all nice polylines

    double dClosestSeparation = HUGE;
    TVecType overallClosestHere = TVecType::uninit(), overallClosestOther = TVecType::uninit();
    for(int nSeg = 0; nSeg < numSegments(); nSeg++) {
        const TLineType seg = segment(nSeg);
        const TVecType closestOnOther = otherPoly.closestPointToBoundedLine(seg);

        const TVecType closestHere = seg.closestPoint(closestOnOther);

        const double dSeparation = (closestHere - closestOnOther).norm();

        if(dSeparation < dClosestSeparation) {
            dClosestSeparation = dSeparation;
            overallClosestHere = closestHere;
            overallClosestOther = closestOnOther;
        }
    }

    if(overallClosestHere == TVecType::uninit())
        return HUGE;

    const TControlPoint cp1 = closestPointAndWidth(overallClosestHere);
    const TControlPoint cp2 = otherPoly.closestPointAndWidth(overallClosestOther);
    return distanceBetween(cp1, cp2);
}

/*  Possibilitiess: multiple intersections => fail
                    One intersection (or one closest point) =>
                      if ill-conditioned => fail
                      else success
*/
template <class TControlPoint>
optional<const C2dImagePointPx>
C2dPolyline_base<TControlPoint>::pointForRecon(const C2dBoundedLine& epiline,
                                               const double dMinAngleForWellConditionedRecon) const
{
    const bool bVerbose = false;

    CHECK(this->numSegments() < 1, "Not enough points in cane");

    double dClosestDist = HUGE;
    C2dImagePointPx closestPointOverall;
    TEigen2dPoint closestSegDirection;

    for(int i = 0; i < this->numSegments(); i++) {

        const C2dImagePointPx closestPoint = epiline.closestPointToBoundedLine(this->segment(i));
        const double dDist = this->segment(i).closestDistance(closestPoint);

        if(dDist <= dClosestDist) {
            // If we have already found an intersection (and it was elsewhere):
            if(dClosestDist == 0 && !closestPoint.isApprox(closestPointOverall)) {
                COUT("Warning: multiple distinct intersection points between epiline and 2D cane"); // This is common
                                                                                                    // because we try
                                                                                                    // all pairs T-E for
                                                                                                    // recon
                return optional<const C2dImagePointPx>();
            }
            closestPointOverall = closestPoint;
            dClosestDist = dDist;
            closestSegDirection = this->segment(i).startToFinish();
        }
    }

    const double dAngle = angleBetweenVectors<TEigen2dPoint>(closestSegDirection, epiline.startToFinish(), false);
    if(dAngle < dMinAngleForWellConditionedRecon)
        return optional<const C2dImagePointPx>();

    return closestPointOverall;
}

/* Returns fail for *either* duplicate intersection, or no intersection */
template <class TControlPoint>
optional<const C2dImagePointPx>
C2dPolyline_base<TControlPoint>::pointOfIntersection(const C2dBoundedLine& epiline) const
{
    optional<C2dImagePointPx> closestIntersection;
    for(int i = 0; i < this->numSegments(); i++) {
        optional<const C2dImagePointPx> intersection = this->segment(i).intersection(epiline);
        if(intersection) {
            if(closestIntersection && !zero(0.5 * (*intersection - *closestIntersection).norm())) {
                cout << "Warning: multiple distinct intersection points between epiline and 2D cane (very unlikely "
                        "unless cane is very kinked)" << endl;
                return optional<const C2dImagePointPx>();
            }

            closestIntersection = *intersection;
        }
    }
    if(closestIntersection)
        return optional<const C2dImagePointPx>(*closestIntersection);
    else
        return optional<const C2dImagePointPx>();
}

void C2dPolylineControlPoint::draw(cv::Mat& image, const cv::Scalar& col, const int nRad, const bool bDrawX) const
{
    CHECK_MAT_INIT(image);

    if(bDrawX)
        drawX(image, getPoint(), col, nRad);
    else
        cv::circle(image, getPoint(), nRad, col, -1);
}

void C2dPolylineControlPointWithThickness::draw(cv::Mat& image, const cv::Scalar& col, const bool bDrawX) const
{
    const int nRad = getCircleRad(*this, -1);
    C2dPolylineControlPoint(getPoint()).draw(image, col, nRad, bDrawX);
}

void C3dPolylineControlPoint::draw(cv::Mat& image,
                                   const CWorldCamera& P,
                                   const cv::Scalar& col,
                                   const double dRad,
                                   const bool bDrawX) const
{
    C3dPolylineControlPointWithThickness(getPoint(), dRad).draw(image, P, col, bDrawX);
}

void C3dPolylineControlPointWithThickness::draw(cv::Mat& image,
                                                const CWorldCamera& P,
                                                const cv::Scalar& col,
                                                const bool bDrawX) const
{
    optional<const C2dPolylineControlPointWithThickness> pcp = P.projectToPx(*this);
    if(pcp)
        pcp->draw(image, col, bDrawX);
}

double C3dPolylineControlPointWithThickness::volume() const
{
    return (4.0 / 3.0) * M_PI * cube(getRad());
}

C3dWorldPoint C3dPolylineControlPointWithThickness::closestPointOnSurface(const C3dWorldPoint& p) const
{
    TEigen3dPoint vecToHeadCentre = p - getPoint();

    if(vecToHeadCentre.norm() == 0) {
        cout << "Warning: closest point to surface called on a centrepoint (undefined)" << endl;
        //THROW("testing");
        vecToHeadCentre.z() = 1;
    }

    const TEigen3dPoint closestPointOnSphere = getPoint() + vecToHeadCentre.normalized() * getRad();
    return closestPointOnSphere;
}

// The line is directed, so the point will be the first intersection of the line into the polyline
C3dWorldPoint
C3dPolylineWithThickness::closestPointToLineOnSurface(const C3dBoundedLine& line /*, double & dDistance*/) const
{
    // Same as closestPointOnSphere: 2 cases.
    const C3dWorldPoint closestOnPoly = closestPointToBoundedLine(line);
    const C3dPolylineControlPointWithThickness closestCP = closestPointAndWidth(closestOnPoly);
    const C3dWorldPoint closestOnLine = line.closestPoint(closestOnPoly);
    const C3dWorldPoint vecToPolyCentre = closestOnPoly - closestOnLine;

    const double dDistToPolyCentre = vecToPolyCentre.norm();

    if(dDistToPolyCentre > closestCP.getRad()) {
        const C3dWorldPoint joinPoint = closestOnPoly - vecToPolyCentre.normalized() * closestCP.getRad();

        CHECK(!zero((joinPoint - closestOnPoly).norm() - closestCP.getRad()),
              "Error computing point where line is closest to cylinder");

        return joinPoint;
    } else {
        // Line intersects thick poly cylinder.

        // The closest on the cylinder, the closest point on the line, and the intersection point form a right angle
        // triangle.
        // Select the point on the line where the hypotenuse has length dRad
        // dRad^2 = dDistCentreToLine^2 + offset^2
        const double dOffset = sqrt(sqr(closestCP.getRad()) - sqr(dDistToPolyCentre));
        const C3dWorldPoint joinPoint = closestOnLine - dOffset * line.direction();
        // dDistance = 0;
        CHECK(!zero(line.unboundedLine().closestDistance(joinPoint)), "Error computing point where line intersects "
                                                                      "sphere (unbounded because line could be fully "
                                                                      "contained within cylinder)");
        return joinPoint;
    }
}

template <class TControlPoint>
const typename CPolyline_base<TControlPoint>::TVecType
CPolyline_base<TControlPoint>::closestPointOnSurface(const typename CPolyline_base<TControlPoint>::TVecType& p) const
{
    const TControlPoint closest = closestPointAndWidth(p);
    TVecType vecToSurface = (p - closest.getPoint()); //.normalized()
    if((p - closest.getPoint()).squaredNorm() == 0) {
        cout << "Warning: closestPointOnSurface called on centre point" << endl;
                //THROW("testing");

        vecToSurface = closestSegment(p).direction(); // now make perpendicular
        std::swap(vecToSurface.y(), vecToSurface.x());
        vecToSurface.y() *= -1;
    }
    const TVecType closestPointOnTrunkSurface = closest.getPoint() + vecToSurface.normalized() * closest.getRad();
    return closestPointOnTrunkSurface;
}

eSuccessStatus C3dPolylineWithThickness::projectThickness(const CWorldCamera& P,
                                                          C2dPolylineWithThickness& thinPoly) const
{
    for(int i = 0; i < numPoints(); i++) {
        if(projectOneThickness(P, thinPoly, i) == eFail)
            return eFail;
    }
    return eSuccess;
}

eSuccessStatus C3dPolylineWithThickness::projectOneThickness(const CWorldCamera& P,
                                                             C2dPolylineWithThickness& poly2d,
                                                             const int nControlPointIndex) const
{
    const C3dPolylineControlPointWithThickness& cp = aControlPoints[nControlPointIndex];
    const optional<const double> pdThickness2d = P.widthToPx(cp.getWidth(), cp.getPoint());
    if(!pdThickness2d)
        return eFail;

    C2dPolylineControlPointWithThickness& cp_2d = poly2d[nControlPointIndex];

    cp_2d = C2dPolylineControlPointWithThickness(cp_2d.getPoint(), *pdThickness2d);

    return eSuccess;
}

std::ostream & operator<<(std::ostream & s, const CRange & range)
{
    s << "Range [" << range.getMin() << ", " << range.getMax() << "]";
    return s;
}

template class CPolylineControlPointWithThickness<C2dBoundedLine>;
template class CPolylineControlPointWithThickness<C3dBoundedLine>;

template class CPolyline_base<C2dPolylineControlPoint>;
template class CPolyline_base<C2dPolylineControlPointWithThickness>;
template class CPolyline_base<C3dPolylineControlPoint>;
template class CPolyline_base<C3dPolylineControlPointWithThickness>;

template class C2dPolyline_base<C2dPolylineControlPoint>;
template class C2dPolyline_base<C2dPolylineControlPointWithThickness>;
template class C3dPolyline_base<C3dPolylineControlPoint>;
template class C3dPolyline_base<C3dPolylineControlPointWithThickness>;

template void thickToThin<C3dPolylineWithThickness, C3dPolyline>(C3dPolylineWithThickness const&, C3dPolyline&);
template void thickToThin<C2dPolylineWithThickness, C2dPolyline>(C2dPolylineWithThickness const&, C2dPolyline&);

template void thinToThick<C3dPolyline, C3dPolylineWithThickness>(C3dPolyline const&,
                                                                 C3dPolylineWithThickness&,
                                                                 const double dThickness);
template void thinToThick<C2dPolyline, C2dPolylineWithThickness>(C2dPolyline const&,
                                                                 C2dPolylineWithThickness&,
                                                                 const double dThickness);

template double intersection<C3dPolylineWithThickness>(const C3dPolylineWithThickness& poly1,
                                                       const C3dPolylineWithThickness& poly2,
                                                       const bool bSlowAndAccurate);
template double intersection<C2dPolylineWithThickness>(const C2dPolylineWithThickness& poly1,
                                                       const C2dPolylineWithThickness& poly2,
                                                       const bool bSlowAndAccurate);

template double tanAngleBetweenVectors<C3dWorldPoint>(const C3dWorldPoint& seg1vec, const C3dWorldPoint& seg2vec);
template double tanAngleBetweenVectors<C2dImagePointPx>(const C2dImagePointPx& seg1vec, const C2dImagePointPx& seg2vec);

