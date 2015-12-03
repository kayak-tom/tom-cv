/*
 * A wrapper around a polyline that is used repeatedly for computing closestPoint/closestDistance/length/direction. Memoises all of these
 */

#ifndef FIXEDPOLYLINE_H
#define FIXEDPOLYLINE_H

#include "polyline.ipp"
#include "vectorAngles.h"

/**
 * @class CPolylinePrecomputed_base
 * @brief Precompute all the segment directions and lengths. These routines are optimised versions of very similar routines in polyline.*
 */
template<class TPolyline>
class CPolylinePrecomputed_base
{
    typedef std::vector<typename TPolyline::TVecType, Eigen::aligned_allocator<typename TPolyline::TVecType> > TPointVec;
    typedef std::vector<typename TPolyline::TVecType, Eigen::aligned_allocator<typename TPolyline::TVecType::Base> > TDirectionVec;
    
    const TPolyline * pPolyline;
    
    TPointVec aPoints;
    TDirectionVec aDirections;
    std::vector<double> adLengths;
    
    double dTotalLength;
public:
    CPolylinePrecomputed_base(const TPolyline & polyline) : pPolyline(&polyline), dTotalLength(0)
    {
        if(polyline.numPoints()>1)
            update(-1);
    }
    
    CPolylinePrecomputed_base() : pPolyline(0), dTotalLength(0)
    {
    }
    
    const TPolyline & getPolyline() const { return *pPolyline; }
    
    inline const TEigen2dPoint endDirection(const eEndpoints end) const { return (end==eStart) ? (-aDirections[0]).eval() : aDirections.back(); }
    
    void update1(const int nIndexToUpdate)
    {
        const bool bVerbose = false;

        if(IS_DEBUG) CHECKNOTNULL(pPolyline);
        const TPolyline & polyline = *pPolyline;
        
        aPoints[nIndexToUpdate] = polyline[nIndexToUpdate].getPoint();
        typename TPolyline::TLineType line = polyline.segment(nIndexToUpdate);
        adLengths[nIndexToUpdate] = line.length();
        aDirections[nIndexToUpdate] = line.startToFinish()/adLengths[nIndexToUpdate];
        
        COUT2("update1 ", nIndexToUpdate);
    }
    
    void update(const int nIndexToUpdate = -1 /*recompute all points*/)
    {
        const bool bVerbose = false;
        COUT2("update ", nIndexToUpdate);
        
        if(IS_DEBUG) CHECKNOTNULL(pPolyline);
        const TPolyline & polyline = *pPolyline;

        if(nIndexToUpdate == -1)
        {
            aPoints.resize(polyline.numSegments());
            aDirections.resize(polyline.numSegments());
            adLengths.resize(polyline.numSegments());
            dTotalLength = 0;
            for(int i=0; i<polyline.numSegments(); i++)
            {
                update1(i);
                dTotalLength += adLengths[i];
            } 
        }
        else
        {
            const int nStart = (nIndexToUpdate > 0) ? (nIndexToUpdate-1) : 0;
            const int nEnd = (nIndexToUpdate<polyline.numSegments()) ? nIndexToUpdate : (polyline.numSegments()-1);
            for(int i=nStart; i<=nEnd; i++)
            {
                dTotalLength -= adLengths[i];
                update1(i);
                dTotalLength += adLengths[i];
            } 
        }
    }
    
    typename TPolyline::TVecType closestPoint_distanceSq_position(const typename TPolyline::TVecType & p, double & dClosestDistanceSq, double & dPosition_index, int & nClosestSeg) const
    {
        const bool bVerbose = false;

        typename TPolyline::TVecType closestPoint = TPolyline::TVecType::uninit();
        dClosestDistanceSq = HUGE;
        nClosestSeg = -1;
        double dClosest_t = -1;
        for(int nSegIndex=0; nSegIndex < (int)aDirections.size(); nSegIndex++)
        {
            //const typename TPolyline::TVecType::Base diff = p-aPoints[nSegIndex];
            double t=-1;// = segmentClamp(diff.dot(aDirections[nSegIndex]));
            typename TPolyline::TVecType closestOnThisSegment = closestPointOnSegment(p, nSegIndex, t);
            const double dDistanceSq = (closestOnThisSegment-p).squaredNorm();
            if(dDistanceSq < dClosestDistanceSq)
            {
                dClosestDistanceSq = dDistanceSq;
                nClosestSeg = nSegIndex;
                closestPoint = closestOnThisSegment;
                dClosest_t = t;
            }
        }
        
        dPosition_index = (double)nClosestSeg+dClosest_t;
        if(bVerbose)
        {
            const double dErr = (closestPoint - pPolyline->closestPoint(p)).norm();
            if(dErr > 0.00001)
            {
                COUT(dErr);
                COUT(closestPoint);
                COUT(pPolyline->closestPoint(p));
                COUT(dClosest_t);
                COUT(nClosestSeg);
                
                COUT(*pPolyline);
                for(int nSegIndex=0; nSegIndex < (int)aDirections.size(); nSegIndex++)
                    cout << TO_STRING(aPoints[nSegIndex]) << TO_STRING(aDirections[nSegIndex]) << endl;
                    
                COUT(dTotalLength);
                COUT(pPolyline->length());
                
                THROW("Closest point computation failed--polyline is out of date");
            }
        }
        
        return closestPoint;
    }
    
    typename TPolyline::TVecType closestPointOnSegment(const typename TPolyline::TVecType & p, const int nSegIndex, double & t) const
    {
        const typename TPolyline::TVecType::Base diff = p-aPoints[nSegIndex];
        
        const double dDistAlongSeg = clip<double>(diff.dot(aDirections[nSegIndex]), 0, adLengths[nSegIndex]); // segmentClamp(diff.dot(aDirections[nSegIndex])/adLengths[nSegIndex]);
        
        t = dDistAlongSeg/adLengths[nSegIndex];
        
        return aPoints[nSegIndex] + dDistAlongSeg * aDirections[nSegIndex];
    }
    
    typename TPolyline::TVecType closestPoint(const typename TPolyline::TVecType & p) const
    {
        double dClosestDistanceSq=-1, dPosition_index=-1;
        int nClosestSeg = -1;
        return closestPoint_distanceSq_position(p, dClosestDistanceSq, dPosition_index, nClosestSeg);
    }
	double signedDistanceToPoint(const typename TPolyline::TVecType & p) const
	{
		double dClosestDistanceSq = -1, dPosition_index = -1;
		int nClosestSeg = -1;
		typename TPolyline::TVecType closestPoint = closestPoint_distanceSq_position(p, dClosestDistanceSq, dPosition_index, nClosestSeg);

		const double dDist = sqrt(dClosestDistanceSq);

		if ((closestPoint-p).dot(perpendicular(nClosestSeg)) < 0)
			return -dDist;
		else
			return dDist;
	}

    void unboundedSegment(const int nSeg, typename TPolyline::TLineType::TUnboundedLine & line) const
    {
        line.fromPointAndDir(aPoints[nSeg], aDirections[nSeg]);
    }
    
    TEigen2dPoint perpendicular(const int nSeg) const
    {
        return ::perpendicular(aDirections[nSeg]); //compiler may complain fo 3D polylines...
    }
    
    const optional<const typename TPolyline::TVecType> closestPointStrictlyOnPoly(const typename TPolyline::TVecType & p) const
    {
        double dClosestDistanceSq=-1, dPosition_index=-1 /*this is clipped I think...*/;
        int nClosestSeg = -1;
        const typename TPolyline::TVecType closestPt = closestPoint_distanceSq_position(p, dClosestDistanceSq, dPosition_index, nClosestSeg);
        
        const int nLastSegIdx = ((int)aDirections.size() - 1);

        const bool bOffStart = dPosition_index <= 0;// closestPt.isApprox(pPolyline->getStartPoint());
        const bool bOffEnd = dPosition_index >= (double)aDirections.size();// closestPt.isApprox(pPolyline->getFinishPoint());

        if(!bOffStart && !bOffEnd)
            return closestPt;
            
        typename TPolyline::TLineType::TUnboundedLine line;

        const double dNoiseMargin = 1; // usually in pixels
        if(bOffStart) {
            unboundedSegment(0, line);
        } else if(bOffEnd) {
            unboundedSegment(nLastSegIdx, line);
        } 
        
        const double dDistanceToEndUnbounded = line.closestDistance(p);
        if(sqr(dDistanceToEndUnbounded+dNoiseMargin) > dClosestDistanceSq)
            return closestPt;
        else
            return optional<const typename TPolyline::TVecType>();
    }

};

/**
 * @class CClosestPoint
 * @brief A class for computing closest points which remembers which segment was closest.
 *
 * Usage CClosestPoint computeClosest(poly);
 * p_closest = computeClosest.closestPoint(p);
 * 
 * Currently only used in cane_reconstructer -- we compute closest points w.r.t. detected 2D canes, which are constant.
 */
template<class TPolyType>
class CClosestPoint
{
    typedef typename TPolyType::TVecType TVecType;
    typedef CPolylinePrecomputed_base<TPolyType> TPolyPrecomp;
    
    static const int SEG_NOT_SET=-1;

    TPolyPrecomp * pPoly;
    int nClosestSeg;
    double dDistToPoly_sq;
    
    //CPolylinePrecomputed_base<TPolyType> * pPolylinePrecomputed;

    TVecType recomputeClosestMemoise(const TVecType & p) HOT;
public:
    CClosestPoint(TPolyPrecomp & poly) : pPoly(&poly), nClosestSeg(SEG_NOT_SET), dDistToPoly_sq(HUGE) {

    }
    CClosestPoint() : pPoly(0), nClosestSeg(SEG_NOT_SET), dDistToPoly_sq(HUGE) {
        //Needed for c'tors of containers of CClosestPoint, but should never actually be used
    }

    TVecType closestPoint_simple(const TVecType & p, const bool bStepIsDelta) HOT;
    
    //Small speedup without affecting results, although there are lots of cases where the approximation is poor
    //typename TPolyType::TControlPoint closestPoint_approx(const TVecType & p, const bool bVerbose) HOT;
};

template<class TPolyType>
typename CClosestPoint<TPolyType>::TVecType CClosestPoint<TPolyType>::recomputeClosestMemoise(const typename TPolyType::TVecType & p)
{
    double dUnused=-1;
//  return closestPoint_distanceSq_position(p, dClosestDistanceSq, dPosition_index, nClosestSeg);
    return pPoly->closestPoint_distanceSq_position(p, dDistToPoly_sq, dUnused, nClosestSeg);
}

template<class TPolyType>
typename CClosestPoint<TPolyType>::TVecType CClosestPoint<TPolyType>::closestPoint_simple(const typename TPolyType::TVecType & p, const bool bStepIsDelta)
{
    const bool bVerbose = false;

    if(IS_DEBUG) CHECKNOTNULL(pPoly);

    if(!bStepIsDelta || nClosestSeg == SEG_NOT_SET) {
        return recomputeClosestMemoise(p);
    }

    double dt_unused = -1;
    
    COUT2("Return closestPointOnSegment", nClosestSeg);
    return pPoly->closestPointOnSegment(p, nClosestSeg, dt_unused);
}

#endif //FIXEDPOLYLINE_H