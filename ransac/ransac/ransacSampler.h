/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

/*
 * ransacSampler.h
 *
 *  Created on: 12/02/2009
 *      Author: tom
 */

#ifndef RANSACSAMPLER_H_
#define RANSACSAMPLER_H_

#include "util/Simple2dPoint.h"
#include "util/exception.h"
#include "util/dynArray.h"
#include <boost/thread/mutex.hpp>
#include "util/random.h"
#include <boost/pool/singleton_pool.hpp>
#include <boost/pool/pool_alloc.hpp>
#include "util/set2.h"
#include "util/smallHashTable.h"
#include "util/optimisation_attributes.h"

class CSampler
{
protected:
	const int n, N;

public:
	//virtual bool isThreadsafe() const { return false; } No longer need to be TS--locking done by RANSAC
	virtual bool choose(TSubSet & anSample) HOT = 0;
    CSampler(int nSampleSize, int nNumPoints) : n(nSampleSize), N(nNumPoints)
	{
		if(IS_DEBUG) CHECK(n>N, "CSampler: Sample size >= count");
		//if(IS_DEBUG) CHECK((int)aPoints0.size() == 0 || aPoints1.size() == 0, "CSampler: No points");
		//if(IS_DEBUG) CHECK(aPoints0.size() != aPoints1.size(), "CSampler: Different numbers of points");
	};
	virtual ~CSampler() {};

	int hypSetSize() const { return n; }
	int numPoints() const { return N; }
};

//Normal RANSAC
class CRANSACSampler : public CSampler
{
public:
	virtual bool choose(TSubSet & anSample) HOT;
	CRANSACSampler(int nSampleSize, int nNumPoints) : CSampler(nSampleSize, nNumPoints) {  };
};

//RANSAC where incompatable points are not chosen together
class CDisjointRANSACSampler: public CSampler
{
protected:
	const CPointIdentifiers & pointIds;

public:
	virtual bool choose(TSubSet & anSample) HOT HOT;
	CDisjointRANSACSampler(int nSampleSize, const CPointIdentifiers & pointIds_in) : CSampler(nSampleSize, pointIds_in.size()), pointIds(pointIds_in) {  };
};

class CGuidedMLESACSampler: public CDisjointRANSACSampler
{
	CInlierProbs adPriorProb_doNotUseDirectly;
protected:
	CInlierProbs & adPriorProb;
	const bool bDisjoint;
public:
	virtual bool choose(TSubSet & anSample) HOT;
	CGuidedMLESACSampler(int nSampleSize, const CPointIdentifiers & pointIds, const int * anArrLikelihood, const double dPriorProb, bool bDisjoint);
	CGuidedMLESACSampler(int nSampleSize, const CPointIdentifiers & pointIds, CInlierProbs & adPriorProb, bool bDisjoint) : CDisjointRANSACSampler(nSampleSize, pointIds), adPriorProb_doNotUseDirectly(0), adPriorProb(adPriorProb), bDisjoint(bDisjoint) {  };
	~CGuidedMLESACSampler() {};
};

//Don't use this--use CBaySACSampler
class CNaiveBayesSampler: public CGuidedMLESACSampler
{
	CDynArray<int> anOrderedIndices;
	void init();
public:
	virtual bool choose(TSubSet & anSample);
	CNaiveBayesSampler(int nSampleSize, const CPointIdentifiers & pointIds, const int * anArrLikelihood,const double dPriorProb, bool bDisjoint) : CGuidedMLESACSampler(nSampleSize, pointIds, anArrLikelihood, dPriorProb, bDisjoint), anOrderedIndices(pointIds.size()) { init(); };
    CNaiveBayesSampler(int nSampleSize, const CPointIdentifiers & pointIds, CInlierProbs & adPriorProb, bool bDisjoint): CGuidedMLESACSampler(nSampleSize, pointIds, adPriorProb, bDisjoint), anOrderedIndices(pointIds.size()) { init(); };
};

//Sample using BaySAC; efficient when prior probabilities come from a continuous set, otherwise use CBaySACSampler_discretePP.
//Does not support disjoint N-M correspondences any more (use CBaySACSampler_discretePP)
// **I haven't maintained this for a while--CBaySACSampler_discretePP is probably faster and more reliable
class CBaySACSampler: public CGuidedMLESACSampler
{
	typedef std::pair<int, double *> idxPair;
	struct probPred
	{
		bool operator()(const idxPair & p1, const idxPair & p2) const
		{
			return *(p1.second) > *(p2.second);
		}
	};
	typedef std::multiset<idxPair, probPred> TSortedProbs;
	TSortedProbs sortedProbs;
	void init() HOT;
public:
	virtual bool choose(TSubSet & anSample) HOT;
	CBaySACSampler(int nSampleSize, const CPointIdentifiers & pointIds, const int * anArrLikelihood, const double dPriorProb/*, bool bDisjoint*/) : CGuidedMLESACSampler(nSampleSize, pointIds, anArrLikelihood, dPriorProb, false) { init(); };
	CBaySACSampler(int nSampleSize, const CPointIdentifiers & pointIds, CInlierProbs & adPriorProb /*, bool bDisjoint*/): CGuidedMLESACSampler(nSampleSize, pointIds, adPriorProb, false) { init(); };
};

//Sample using BaySAC; efficient when prior probabilities come from a discrete set (e.g. 20 have p=0.6, 30 have p=0.3, 150 have p=0.15)
//Supports N-M correspondences
class CBaySACSampler_discretePP: public CGuidedMLESACSampler
{
    typedef CFixedArray<int, NN_MAX> TFastIntVec;
    typedef CFixedArray<int, 2*(NN_MAX-1)> TLongerFastIntVec; // a N-M match can be incompatible with (N-1) + (M-1) others

	typedef set2<int, std::less<int>/*, CIndividualPool_NoFree_Allocator<int, 32 >*/ > TIntSet;
	typedef map2<double, TIntSet, std::greater<double> /*, std::allocator<std::pair< double, TIntSet> > / *CIndividualPool_NoFree_Allocator<std::pair<double, TIntSet>, 512 >*/ > TSortedProbs;
	TSortedProbs sortedProbVectors;

	std::vector<TLongerFastIntVec> aaIncompatable;

	void init();

	typedef CSmallHashTable<int, 40, 1, intHash<40> > TExcludeSet;
	int chooseFromSet(int numAlreadyChosen, TSubSet & anSample, const TIntSet & options, TExcludeSet & excludeSet, int & nNumExcluded);
public:
	virtual bool choose(TSubSet & anSample) HOT;
	CBaySACSampler_discretePP(int nSampleSize, const CPointIdentifiers & pointIds, const int * anArrLikelihood, const double dPriorProb, bool bDisjoint) : CGuidedMLESACSampler(nSampleSize, pointIds, anArrLikelihood, dPriorProb, bDisjoint), aaIncompatable((int)pointIds.size()) { init(); }
	CBaySACSampler_discretePP(int nSampleSize, const CPointIdentifiers & pointIds, CInlierProbs & adPriorProb, bool bDisjoint): CGuidedMLESACSampler(nSampleSize, pointIds, adPriorProb, bDisjoint), aaIncompatable(pointIds.size()) { init(); }
};

class CPROSACSampler: public CGuidedMLESACSampler
{
	std::vector<int> anOrderedIndices;
	int T, t, g;
	double dT;
	bool bSystematic;
	std::vector<const int *> aSamples;
	void init();
public:
	virtual bool choose(TSubSet & anSample) HOT;
	CPROSACSampler(int nSampleSize, const CPointIdentifiers & pointIds, const int * anArrLikelihood,const double dPriorProb, bool bDisjoint, bool bSystematic) : CGuidedMLESACSampler(nSampleSize, pointIds, anArrLikelihood, dPriorProb, bDisjoint), anOrderedIndices(pointIds.size()), T(1), t(0), g(n), dT(1), bSystematic(bSystematic) { init(); };
	CPROSACSampler(int nSampleSize, const CPointIdentifiers & pointIds, CInlierProbs & adPriorProb, bool bDisjoint, bool bSystematic): CGuidedMLESACSampler(nSampleSize, pointIds, adPriorProb, bDisjoint), anOrderedIndices(pointIds.size()), T(1), t(0), g(n), dT(1), bSystematic(bSystematic) { init(); };
	~CPROSACSampler()
	{
		for(std::vector<const int *>::const_iterator ppSample = aSamples.begin(); ppSample != aSamples.end(); ppSample++)
			delete [] *ppSample;
	}
};

template<class bitType>
class CStatus//Implements the bitfield or whatever
{
protected:
	const int N;
	bitType * abBits;
public:
	CStatus(const int count) : N(count) { abBits = new bitType[N]; };
	~CStatus() { delete [] abBits; };
	inline void set(int n, bitType val) { abBits[n] = val; };
	inline bitType get(int n) { return abBits[n]; };
};

template<>
class CStatus<bool> //Implements the bitfield or whatever
{
protected:
	const int N;
	unsigned * abBits;
public:
	CStatus(const int count) : N(count) { abBits = new unsigned[N/32 + 1]; };
	~CStatus() { delete [] abBits; };
	void set(int n, bool val) {

		int idx = n/32;
		int bit = n % 32;
		unsigned bitMask = 1 << bit;
		if(val)
			abBits[idx] |= bitMask;
		else
			abBits[idx] ^= bitMask;
	};
	bool get(int n)
	{
		int idx = n/32;
		int bit = n % 32;
		unsigned bitMask = 1 << bit;
		return (abBits[idx] & bitMask) != 0;
	};
};

//SimSAC sampling, supporting N--M correspondences
class CSimSACSampler: public CGuidedMLESACSampler
{
	typedef int TStatus; //int and bool appear fastest...
	typedef std::vector<int *> THistory;
	typedef std::vector<int> THistoryIndices;

	THistory history;
	THistoryIndices anHistIndices;

	int T;
	const bool bExact;

	static const TStatus SS_UNSET = 2;
    static const TStatus SS_TRUE = 1;
    static const TStatus SS_FALSE = 0;

	void chooseML_exact(TSubSet & anSampleHist);
	void chooseML_approx(TSubSet & anSampleHist);
	void chooseML_approx2(TSubSet & anSampleHist);
	bool chooseML(TSubSet & anSample);

	CStatus<TStatus> simStatus;
	CDynArray<int> anAllIndices;
	const int MAX_ITERS;
public:
	void init();
	virtual bool choose(TSubSet & anSample) HOT;
	CSimSACSampler(int nSampleSize, const CPointIdentifiers & pointIds, const int * anArrLikelihood, const double dPriorProb, bool bDisjoint, bool bExact, const int T, const int MAX_ITERS = 100) : CGuidedMLESACSampler(nSampleSize, pointIds, anArrLikelihood, dPriorProb, bDisjoint), T(T), bExact(bExact), simStatus(pointIds.size()), anAllIndices(pointIds.size()), MAX_ITERS(MAX_ITERS) { init(); };
	CSimSACSampler(int nSampleSize, const CPointIdentifiers & pointIds, CInlierProbs & adPriorProb, bool bDisjoint, bool bExact, const int T, const int MAX_ITERS = 100): CGuidedMLESACSampler(nSampleSize, pointIds, adPriorProb, bDisjoint), T(T), bExact(bExact), simStatus(pointIds.size()), anAllIndices(pointIds.size()), MAX_ITERS(MAX_ITERS) { init(); };
	~CSimSACSampler() {
		for(THistory::iterator paHist = history.begin(); paHist<history.end(); paHist++)
			delete [] *paHist;
	};
};

#endif /* RANSACSAMPLER_H_ */
