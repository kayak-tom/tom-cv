/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * CRansacTerminator.cpp
 *
 *  Created on: 7/05/2009
 *      Author: tom
 */

#include "cRansacTerminator.h"
#include "util/exception.h"
#include "util/convert.h"
#include <iostream>
#include <boost/math/distributions/normal.hpp>
#include <boost/filesystem.hpp>
#include <fstream>

CRansacTerminator::CRansacTerminator() {
}

CRansacTerminator::~CRansacTerminator() {
}
        

CWaldSACTerminator::CWaldSACTerminator(int nCount, int nSampleSize, int nModelsPerSample, double dPointsPer1ModelTime, double minPropInliers)
: nCount(nCount), nSampleSize(nSampleSize), nModelsPerSample(nModelsPerSample), dPointsPer1ModelTime(dPointsPer1ModelTime), dCountInv(1.0 / nCount), A(HUGE), delta(0.5), epsilon(delta + 0.0001), //just to keep it stable (and stop any termination) before initialised
nInliers(0), nBestInlierCount((int) (nCount * (double) minPropInliers)), n2ndBestInlierCount(0), nStabilityCheckInterval(20),
nTotalDataPoints(0), nTotalInliersToBadModels(0), nDataPoints(0) {
    epsilon = (nBestInlierCount + 1) * dCountInv;
    delta = epsilon * 0.5;
    updateWaldSACThresh();
};

void CLogWaldSACTerminator::updateWaldSACThresh() {
    CWaldSACTerminator::updateWaldSACThresh();
    updateLogs();
}

void CLogWaldSACTerminator::updateLogs()
{
    log1takedelta_on_1takeepsilon = log(1 - delta) - log(1 - epsilon);
    logdelta_on_epsilon = log(delta) - log(epsilon);
    updateLogA();
}
void CLogWaldSACTerminator::updateLogA()
{
    logLambda_on_A = -log(A);
}

void CWaldSACTerminator::updateWaldSACThresh() {
    if (epsilon > 0.9999) {
        A = 1;
        return;
    }

    /*if (delta == 0) {
        A = HUGE;
        return;
    }*/

    double ms = nModelsPerSample;
    double tm = dPointsPer1ModelTime;
    double Pg = pow(epsilon, nSampleSize);
    double C = (1 - delta) * log((1 - delta) / (1 - epsilon)) + delta * log(
            delta / epsilon);
    double K1 = tm / Pg;
    double K2 = ms / (Pg * C);
    double K = K1 / K2 + 1;
    A = K;
    CHECKNAN(A);
    for (int ii = 0; ii < 40; ii++) { //Iterate to convergence
        double eps = K + log(A) - A;
        A += eps;
        //cout << "A=" << A << ",";
        if (fabs(eps) < 0.001) break;
    }
    //cout << endl;
}

bool CWaldSACTerminator::observeStatus(bool bStatus) {
    nDataPoints++;

    if (bStatus) {
        nInliers++;
        dLambda_num *= delta; //Prob inlier/outlier | Hbad
        dLambda_denom *= epsilon; //I think epsilon must be > delta
    } else {
        dLambda_num *= (1 - delta); //Prob inlier/outlier | Hbad
        dLambda_denom *= (1 - epsilon); //I think epsilon must be > delta
    }

    if (EXPECT(dLambda_num > A * dLambda_denom, 0)) {
        return true;
    }

    //Keep it stable
    if (EXPECT(nStabilityCheckInterval <= 0, 0)) {
        nStabilityCheckInterval = 35;
        const double tiny = 1e-10, big = 1e+10;
        if (dLambda_num > 0 && dLambda_denom > 0) {
            while ((dLambda_num < tiny) && (dLambda_denom < tiny)) {
                dLambda_num *= big;
                dLambda_denom *= big;
            }
            while ((dLambda_num > big) && (dLambda_denom > big)) {
                dLambda_num *= tiny;
                dLambda_denom *= tiny;
            }
        }
    }
    nStabilityCheckInterval--;
    return false;
}

bool CLogWaldSACTerminator::observeStatus(bool bStatus) {
    nDataPoints++;

    if (bStatus) {
        nInliers++;
        logLambda_on_A += logdelta_on_epsilon; //Prob inlier/outlier | Hbad
        //dLambda_denom += logepsilon ; //I think epsilon must be > delta
    } else {
        logLambda_on_A += log1takedelta_on_1takeepsilon; //Prob inlier/outlier | Hbad
        //dLambda_denom += log1takeepsilon; //I think epsilon must be > delta
    }

    if (EXPECT(logLambda_on_A > 0, 0)) {
        //cout << "Terminated after " << nDataPoints << ", bgc =" << bestGoodCount() << endl;
        return true;
    }

    return false;
}
void CLogWaldSACTerminator::reset() {
    CWaldSACTerminator::reset();
    updateLogA();
}

void CWaldSACTerminator::reset() {

    bool bUpdate = false;

    if (nInliers > nBestInlierCount) {
        nBestInlierCount = nInliers;
        epsilon = nBestInlierCount * dCountInv;
        bUpdate = true;
    } else if (nInliers > 0 && nInliers < nBestInlierCount) {
        //update delta
        nTotalDataPoints += nDataPoints;
        nTotalInliersToBadModels += nInliers;

        double dNewDelta = nTotalInliersToBadModels / (double) nTotalDataPoints;
        if (!within(delta, dNewDelta, 0.05)) {
            delta = dNewDelta;
            bUpdate = true;
        }
        if(IS_DEBUG) CHECK(delta>epsilon, "Delta calc failed");
    }

    if (bUpdate)
        updateWaldSACThresh();

    dLambda_num = 1, dLambda_denom = 1;

    if(IS_DEBUG) CHECK(delta > epsilon && delta > 0.000001 && epsilon < 1, "CWaldSACTerminator::reset: delta >= epsilon");

    nInliers = 0;
    nDataPoints = 0;
}

CBBTerminator::CBBTerminator(int nTotalCount, double dMinPropInliers, double dPropFail) : N(nTotalCount), I(0), nBestInlierCount((int)(N*dMinPropInliers)) 
{
	boost::math::normal ND(-M_PI*(2.0/3.0), M_PI/5);
	k = boost::math::quantile(ND, dPropFail);
	//try -2.22792, sqrt(0.36916)
}
bool CBBTerminator::observeStatus(int nThread, bool bStatus)
{
	n++;
	if(bStatus)
	{
		I++;
		return false;
	}
	else
	{
		double np = n*p;
		double LHS = I - np;
		double RHS = n*(N-n)*constant_expr;// constant_expr = k*sqrt(N * p*(1-p))/sqr(N);

		return LHS<RHS;
	}

}
	
void CBBTerminator::reset()
{
	if(I>nBestInlierCount)
		nBestInlierCount = I;

	p = (nBestInlierCount+1)/(double)N;
	n=0;
	I=0;
	constant_expr = k * sqrt(N * p*(1-p))/sqr(N);
}

bool CBBLinearTerminator::observeStatus(int nThread, bool bStatus)
{
	n++;
	if(bStatus)
	{
		I++;
		return false;
	}
	else
	{
		double np = n*p;
		double LHS = I - np;

		return LHS < linearThresh;
	}
}
	
void CBBLinearTerminator::reset()
{
	CBBTerminator::reset();
	linearThresh = T/sqrt(N*p*(1-p));
}

void CSimpleTerminator::reset()
{
	for(int i=0; i< MAX_THREADS; i++)
		if(anInliers[i] > nBestInlierCount) nBestInlierCount = anInliers[i];

	resetCount();
}

void CSimpleTerminator::resetCount()
{
	for(int i=0; i< MAX_THREADS; i++)
		anCount[i] = anInliers[i] = 0;
}

bool CSimpleTerminator::observeStatus(int nThread, bool bStatus)
{
	anCount[nThread]++;

	if(bStatus)
	{
		anInliers[nThread]++;
	}
	else
	{
		//once nStartChecking below par stop.
		//nInliers + nStartChecking < nBestInlierCount * nCount/nTotalCount
		if((anInliers[nThread] + nBelowParCutoff)*nTotalCount < nBestInlierCount * anCount[nThread])
		{
			//std::cout << "Terminating after " << anInliers[nThread] << " inliers from " << anCount[nThread] << " of " << nTotalCount << std::endl;
			return true;
		}
		else if(nTotalCount - anCount[nThread] < nBestInlierCount - anInliers[nThread])
			return true; //Not possible to be best now
	}
	return false;
}

CTerminatorLogger::~CTerminatorLogger() {
    reset(); /*log last time*/

    bool bNew = false;
    if (!boost::filesystem::exists("terminationRes.tsv"))
        bNew = true;
    std::ofstream results("terminationRes.tsv", std::ios_base::app);

    if (bNew)
        results << "nIterations" << '\t' << "dPropEvalsMadeOnBadHyps" << '\t' << "bBestMissed\tTotalSpeedup" << endl;

    double dPropEvalsMadeOnBadHyps = nTotalEvaluationsOnBadHyp / (double)nMaxEvaluationsOnBadHyp;

    int bBestMissed = (bestGoodCount() < nRealBGC) ? 1 : 0;

    int nGoodHyps = nIterations - (int)(nMaxEvaluationsOnBadHyp / (double) nCount);
    double nTotalSpeedup = (nTotalEvaluationsOnBadHyp + nGoodHyps * nCount) / ((double) nIterations * nCount);

    results << nIterations << '\t' << dPropEvalsMadeOnBadHyps << '\t' << bBestMissed << '\t' << nTotalSpeedup << endl;
    results.close();
}

bool CTerminatorLogger::observeStatus(int nThread, bool bStatus) {
    if (!bTerminated) {
        nNumPointsAtTermination++; //Count number of iterations before termination.
        if (bStatus)
            nInliersAtTermination++;
        bTerminated = pTerminator->observeStatus(nThread, bStatus);
    }

    if(bStatus)
        nTotalInliers++;

    return false;
}

void CTerminatorLogger::reset() {

    if (nNumPointsAtTermination > 0) //Not the first time
    {
        //Log data
        if (bTerminated) {
            bool bFalseTermination = (nTotalInliers > bestGoodCount());
            if (bFalseTermination)
                nFalseTerminations++;
        }

        //If not the best (whether or not it terminated):
        if(nInliersAtTermination < bestGoodCount())
        {
            nTotalEvaluationsOnBadHyp += nNumPointsAtTermination;
            nMaxEvaluationsOnBadHyp += nCount;
        }

        if (nTotalInliers > nRealBGC)
            nRealBGC = nTotalInliers;
    }

    pTerminator->reset();
    bTerminated = false;
    nNumPointsAtTermination = 0;
    nInliersAtTermination = 0;
    nTotalInliers = 0;
    nIterations++;
}