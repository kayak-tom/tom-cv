/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

/*
 * inlierCounter.h
 *
 *  Created on: 8/10/2009
 *      Author: tom
 */

#ifndef INLIERCOUNTER_H_
#define INLIERCOUNTER_H_

#include "util/exception.h"
#include "util/convert.h"
#include "util/Simple2dPoint.h"
#include "models.h"
#include "cRansacTerminator.h"
#include <Eigen/Core>

class CInlierCounter {
protected:
    const double dThreshold_sq; // threshhold used may be adapted when doing topdown refinement

    const T2dPoints & ap;
    const T2dPoints & ap_prime;

public:

    CInlierCounter(const double dThreshold, const T2dPoints & ap, const T2dPoints & ap_prime) : dThreshold_sq(sqr(dThreshold)), ap(ap), ap_prime(ap_prime) {
    }

    //void so can be boost::bound. 
    virtual void countInliers(const CModel & model, CRansacTerminator * pTerminator, const int nThread, const int nBGC, CMask & mask, double * adResiduals, int & nInlierCount, double dThresholdScale) const = 0;

    double threshold_sq() const {
        return dThreshold_sq;
    };
};

class CEssentialMatInlierCounter : public CInlierCounter {
    void timeInlierCounter() const;
    bool isInlier_old(const C3x3MatModel & f, const CSimple2dPoint & p1, const CSimple2dPoint & p2, double * pdResidual, const double dThreshold_sq_use) const HOT HARD_INLINE;
    bool isInlier_vectorise(const Eigen::Matrix3f & f, const CSimple2dPoint & p1, const CSimple2dPoint & p2, double * pdResidual, const float dThreshold_sq_use) const HOT HARD_INLINE;
    bool isInlier_vectoriseDouble(const Eigen::Matrix<double, 4, 2 > & f, const Eigen::Vector4d & fcol3, const CSimple2dPoint & p1, const CSimple2dPoint & p2, double * pdResidual, const double dThreshold_sq_use) const HOT HARD_INLINE;

public:
    template<typename TFloat>
    static bool isInlier(const TFloat * f, const CSimple2dPoint & p1, const CSimple2dPoint & p2, double * pdResidual, const TFloat dThreshold_sq_use) HOT HARD_INLINE;
    
    //Actual L2 Sampson's error (H+Z)
    template<typename TFloat>
    static bool isInlier_SE(const TFloat * f, const CSimple2dPoint & p1, const CSimple2dPoint & p2, double * pdResidual, const TFloat dThreshold_sq_use) HOT HARD_INLINE;

    CEssentialMatInlierCounter(const double dThreshold, const T2dPoints & ap, const T2dPoints & app) : CInlierCounter(dThreshold, ap, app) {
    }

    //void so can be boost::bound. NOT VIRTUAL also for bind
    virtual void countInliers(const CModel & model, CRansacTerminator * pTerminator, const int nThread, const int nBGC, CMask & mask, double * adResiduals, int & nInlierCount, double dThresholdScale) const HOT;
};

class CHomographyInlierCounter : public CInlierCounter {
protected:
    inline bool isInlier(const C3x3MatModel & H, const CSimple2dPoint & p1, const CSimple2dPoint & p2, double * pdResidual, const double dThreshold_sq_use) const HOT HARD_INLINE;

public:

    CHomographyInlierCounter(const double dThreshold, const T2dPoints & p1, const T2dPoints & p2) : CInlierCounter(dThreshold, p1, p2) {
    }

    virtual void countInliers(const CModel & model, CRansacTerminator * pTerminator, const int nThread, const int nBGC, CMask & mask, double * adResiduals, int & nInlierCount, double dThresholdScale) const HOT;
};
#endif /* INLIERCOUNTER_H_ */
