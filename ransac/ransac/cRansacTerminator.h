/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

/*
 * CRansacTerminator.h
 *
 *  Created on: 7/05/2009
 *      Author: tom
 */

#ifndef CRANSACTERMINATOR_H_
#define CRANSACTERMINATOR_H_

#include "util/exception.h"
#include "math.h"

class CRansacTerminator {
public:
	CRansacTerminator();
	virtual ~CRansacTerminator();
	virtual bool observeStatus(bool bStatus) { return false; }
	virtual bool observeStatus(int nThread, bool bStatus) { return false; }
	virtual void reset() { }
	virtual int bestGoodCount() { return -1; } //For debug output
	virtual bool isThreadsafe() const { return true; }
};

class CWaldSACTerminator : public CRansacTerminator {
    const int nCount, nSampleSize, nModelsPerSample;
    const double dPointsPer1ModelTime, dCountInv;

protected:
    double A, delta, epsilon;
    int nInliers;
    double dLambda_num, dLambda_denom;
    int nBestInlierCount, n2ndBestInlierCount, nStabilityCheckInterval, nTotalDataPoints, nTotalInliersToBadModels, nDataPoints;

    virtual void updateWaldSACThresh();
    bool observeStatus(bool bStatus);
public:
    CWaldSACTerminator(int nCount, int nSampleSize, int nModelsPerSample, double dPointsPer1ModelTime, double minPropInliers);

    virtual bool observeStatus(int nThread, bool bStatus) {
        if (nThread > 0)
            THROW("Not threadsafe");
        return observeStatus(bStatus);
    }
    virtual void reset();

    virtual int bestGoodCount() {
        return nBestInlierCount;
    } //For debug output

    virtual bool isThreadsafe() const {
        return false;
    }
};

class CLogWaldSACTerminator : public CWaldSACTerminator {
    double logLambda_on_A, logdelta_on_epsilon, log1takedelta_on_1takeepsilon;
    bool observeStatus(bool bStatus);

protected:
    virtual void updateWaldSACThresh();
    void updateLogs();
    void updateLogA();
public:

    CLogWaldSACTerminator(int nCount, int nSampleSize, int nModelsPerSample, double dPointsPer1ModelTime, double minPropInliers) : CWaldSACTerminator(nCount, nSampleSize, nModelsPerSample, dPointsPer1ModelTime, minPropInliers) {
        updateLogs();
    }

    virtual bool observeStatus(int nThread, bool bStatus) {
        if (nThread > 0)
            THROW("Not threadsafe");
        return observeStatus(bStatus);
    }

    void reset();
};

class CTerminatorLogger : public CRansacTerminator {
    CRansacTerminator * pTerminator;
    bool bTerminated;
    int nNumPointsAtTermination, nInliersAtTermination, nFalseTerminations, nRealBGC, nIterations, nTotalEvaluationsOnBadHyp, nMaxEvaluationsOnBadHyp, nCount, nTotalInliers;
public:

    CTerminatorLogger(CRansacTerminator * pTerminator, int nCount) : pTerminator(pTerminator), bTerminated(false), nNumPointsAtTermination(0), nInliersAtTermination(0), nFalseTerminations(0), nRealBGC(0), nIterations(-1),
    nTotalEvaluationsOnBadHyp(0), nMaxEvaluationsOnBadHyp(0), nCount(nCount), nTotalInliers(0) {
    }

    virtual ~CTerminatorLogger();

    virtual bool observeStatus(int nThread, bool bStatus);
    virtual void reset();

    virtual int bestGoodCount() {
        return pTerminator->bestGoodCount();
    }

    virtual bool isThreadsafe() const {
        return false;
    }

};


class CSimpleTerminator : public CRansacTerminator
{
	const int nTotalCount;
	const int nBelowParCutoff;
	int nBestInlierCount;

	static const int MAX_THREADS = 16; //Can be increased
	int anInliers[MAX_THREADS], anCount[MAX_THREADS];

	//void updateWaldSACThresh() ;
public:
	CSimpleTerminator(int nTotalCount) : nTotalCount(nTotalCount), nBelowParCutoff(20), nBestInlierCount(0) { resetCount(); }
	//CSimpleTerminator() : nTotalCount(0), nBelowParCutoff(20), nInliers(0), nBestInlierCount(0) { resetCount(); }
	//virtual bool observeStatus(bool bStatus);
	virtual bool observeStatus(int nThread, bool bStatus);
	virtual void reset();
	void resetCount();
	virtual int bestGoodCount() { return nBestInlierCount; } //For debug output
	virtual bool isThreadsafe() const { return true; }
};

class CBBTerminator : public CRansacTerminator
{
	double k, constant_expr;
protected:
	const int N;
	int I;
	int nBestInlierCount;
	double p;
	int n;

public:
	CBBTerminator(int nTotalCount, double dMinPropInliers, double dPropFail);

	virtual bool observeStatus(int nThread, bool bStatus);
	virtual void reset();
	virtual int bestGoodCount() { return nBestInlierCount; } //For debug output
	virtual bool isThreadsafe() const { return false; }
};

class CBBLinearTerminator : public CBBTerminator
{
	const double T;
	double linearThresh;
public:
	CBBLinearTerminator(int nTotalCount, double dMinPropInliers, double dPropFail) : CBBTerminator(nTotalCount, dMinPropInliers, dPropFail), T(sqrt(-0.5*log(dPropFail)))
	{
			

	}

	virtual bool observeStatus(int nThread, bool bStatus);
	virtual void reset();
};

#endif /* CRANSACTERMINATOR_H_ */
