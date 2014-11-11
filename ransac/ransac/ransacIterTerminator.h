/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

/*
 * ransacIterTerminator.h
 *
 *  Created on: 7/10/2009
 *      Author: tom
 */

#ifndef RANSACITERTERMINATOR_H_
#define RANSACITERTERMINATOR_H_

#include "util/Simple2dPoint.h"
#include "util/exception.h"
#include "util/convert.h"
#include <algorithm>

//Decide when to stop iterating (this one stops after maxiters)
class CIterTerminator
{
protected:
	int nMaxIters;
	int nBGC, nIters;

	virtual int updateMaxIters() { return nMaxIters; };

public:
	virtual bool isThreadsafe() const { return false; } //False by default Doesn't actually matter.

	CIterTerminator(int nMaxIters) : nMaxIters(nMaxIters), nBGC(0), nIters(0) {}

	//True if best
	virtual bool updateBGC(int nNewGC)
	{
		if(nBGC < nNewGC)
		{
			nBGC = nNewGC;

			//Try reducing iters by whatever alg.
			nMaxIters = std::min<int>(nMaxIters, updateMaxIters());
			return true;
		}
		else
			return false;
	}
	void updateNumIters(int nNewIters) const {  } //NOT TS but shouldn't actually matter

	bool terminate(int nNewIters) { nIters += nNewIters; return nIters>=nMaxIters; }

	int BGC() const { return nBGC; }
	int numIters() const { return nIters; }
	int maxIters() const { return nMaxIters; }
};

//Stop when inlier set gets this big
class CPropIterTerminator : public CIterTerminator
{
protected:
	const int nNumInliersToTerminateOn;

	virtual int updateMaxIters() { if(nBGC >= nNumInliersToTerminateOn) return 0; else return nMaxIters; };

public:
	virtual bool isThreadsafe() const { return false; } //False by default Doesn't actually matter.

	CPropIterTerminator(int nMaxIters, int nNumInliersToTerminateOn) : CIterTerminator(nMaxIters), nNumInliersToTerminateOn(nNumInliersToTerminateOn) {}

	//True if best
};


/*class CTSIterTerminator
{
	boost
protected:
	virtual int updateMaxIters() { return nMaxIters; };

public:
	virtual bool isThreadsafe() const { return true; } //False by default

	CIterTerminator(int nMaxIters) : nMaxIters(nMaxIters), nBGC(0), nIters(0) {}
};*/

//Classic RANSAC update to reduce number of iters
class CRansacIterTerminator : public CIterTerminator
{
	const int nCount, nChoose;
	const double dProbFindNInliers;
public:
	CRansacIterTerminator(int nMaxIters, int nCount, int nChoose, double pProbSuccess)
	: CIterTerminator(nMaxIters), nCount(nCount), nChoose(nChoose), dProbFindNInliers(pProbSuccess)
	{
		//can actually be true as nCount is a sum of probabilities if(IS_DEBUG) CHECK(nCount < nChoose, "Insufficient points to choose nChoose of them!");
		//Updating num iters a priori is a bit dubious so don't--if we don't find a hypothesis it indicates our inlier proportion estimate was too high.
	}

	virtual int updateMaxIters()
	{
		if(IS_DEBUG) CHECK(nCount < nChoose, "Insufficient points to choose nChoose of them!");

		double dPropInliersFound = nBGC/(double)nCount;

		if(dPropInliersFound > 0.99) return 0; //dProbInlier >= 1 (greater if we're using a realistic true inlier estimate)

		double dProbChooseNInliers = pow(dPropInliersFound, nChoose);
		if((dProbChooseNInliers<0.001)) return MAX_INT;

		int nIters = std::max<int>(0, (int)(log(1-dProbFindNInliers)/log(1-dProbChooseNInliers)));

		//cout << nCount << " points, ";
		//cout << nChoose << " needed for model, ";
		//cout << dProbInlier << " of them are inliers, ";
		//cout << dProbFindNInliers << " is the prob of finding a good subset (old method) " << endl;

		return nIters;
	}
};

//Use prior probabilities to estimate a realistic inlier count--use this to terminate early, as above
class CPPIterTerminator : public CRansacIterTerminator
{
	static int getRealisticInlierCount(const CInlierProbs & adArrLikelihood)
	{
		double nRealisticInlierCount = 0;
		for (int i = 0; i < (int)adArrLikelihood.size(); i++)
			nRealisticInlierCount += adArrLikelihood[i];

		return (int)((nRealisticInlierCount+3) * 2.0); //Todo: 2.0 is a scale to make this compatible with what we observe.
	}

public:
	CPPIterTerminator(int nMaxIters, const CInlierProbs & adPriorProbs, int nChoose, double pProbSuccess) : CRansacIterTerminator(nMaxIters,  getRealisticInlierCount(adPriorProbs), nChoose, pProbSuccess) {};
};

#endif /* RANSACITERTERMINATOR_H_ */
