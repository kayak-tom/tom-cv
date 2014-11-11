/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

/*
 * refiner.h
 *
 *  Created on: 7/10/2009
 *      Author: tom
 */

#ifndef REFINER_H_
#define REFINER_H_

#include "util/Simple2dPoint.h"
#include "util/exception.h"

class CModel;

class CModelRefiner //Model==e.g. fit an essential matrix to >=8 points. May also refine mask.
{
protected:
	const int nMinPoints;
public:
	CModelRefiner(int nMinPoints) : nMinPoints(nMinPoints)
	{
		if(IS_DEBUG) CHECK(nMinPoints<0, "CModelRefiner::fitModel: Invalid sizes");
	}

	virtual bool fitModel(CMask & mask, CModel & pModel, bool bVerbose) = 0;

	int minNumPoints() const { return nMinPoints; } //Cant refine with less than this many points
};

// Do not refine model (0 topdown iterations; for testing)
class CNoModelRefiner : public CModelRefiner
{
public:
    CNoModelRefiner(int nMinPoints) : CModelRefiner(nMinPoints) {}
    
    virtual bool fitModel(CMask & mask, CModel & pModel, bool bVerbose) { return true; }
};

// Any model working on image coordinates
class CImCorrModelRefiner : public CModelRefiner
{
protected:
	const T2dPoints & p1;
	const T2dPoints & p2;
	virtual bool fitModel_int(CMask & mask, CModel & pModel, bool bVerbose) = 0;
public:
	CImCorrModelRefiner(int nMinPoints, const T2dPoints & p1, const T2dPoints & p2) : CModelRefiner(nMinPoints), p1(p1), p2(p2)
	{
		if(IS_DEBUG) CHECK(nMinPoints<=0, "CImCorrModelRefiner::fitModel: Invalid sizes");
	}

	virtual bool fitModel(CMask & mask, CModel & pModel, bool bVerbose)
	{
		const int nSize = p1.size();
		if(IS_DEBUG) CHECK(nSize < nMinPoints, "CImCorrModelRefiner::fitModel: Insufficient points");
		if(IS_DEBUG) CHECK(nSize != (int)p2.size() || nSize != (int)mask.size(), "CImCorrModelRefiner::fitModel: Size mismatch");
		return fitModel_int(mask, pModel, bVerbose);
	}
};

#endif /* REFINER_H_ */
