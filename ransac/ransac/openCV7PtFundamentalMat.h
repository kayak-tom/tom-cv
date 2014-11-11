/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

/*
 * openCV8PtEssentialMatModified.h
 *
 *  Created on: 8/10/2009
 *      Author: tom
 */

#ifndef OPENCV7PT_H_
#define OPENCV7PT_H_

//#include "refiner.h"
#include "hypothesiser.h"

/*class COpenCVHomography: public CImCorrModelRefiner
{
	virtual bool fitModel_int(CMask & mask, CModel & pModel, bool bVerbose);
public:
	COpenCVHomography(const T2dPoints & p1, const T2dPoints & p2) : CImCorrModelRefiner(4, p1, p2) {}
	virtual ~COpenCVHomography() {}
};*/

class COpenCV7PtFundamentalMat : public CImCorrModelHypothesiser //Model==essential matrix from >=8 points
{
	//Return number of hypotheses
	virtual int getModels_int(const TSubSet & anHypSet, CModels & pModels);
public:
	virtual bool isThreadsafe() const { return true; } // as has no state

	COpenCV7PtFundamentalMat(const T2dPoints & p1, const T2dPoints & p2) : CImCorrModelHypothesiser(7, p1, p2) {}

        virtual int modelsPerIteration() const { return 2; /*average*/ }
        virtual int timePerIteration() const { return 700; /* times time to evaluate a hypothesis */ }
};

#endif /* OPENCV7PT_H_ */
