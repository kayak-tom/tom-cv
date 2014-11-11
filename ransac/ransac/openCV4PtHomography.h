/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

/*
 * openCV8PtEssentialMatModified.h
 *
 *  Created on: 8/10/2009
 *      Author: tom
 */

#ifndef OPENCV4PT_HOMOGRAPHY_H_
#define OPENCV4PT_HOMOGRAPHY_H_

#include "refiner.h"
#include "hypothesiser.h"

class COpenCVHomography: public CImCorrModelRefiner
{
	virtual bool fitModel_int(CMask & mask, CModel & pModel, bool bVerbose);
public:
	COpenCVHomography(const T2dPoints & p1, const T2dPoints & p2) : CImCorrModelRefiner(4, p1, p2) {}
	virtual ~COpenCVHomography() {}
};

class COpenCV4ptHomography : public CImCorrModelHypothesiser //Model==essential matrix from >=8 points
{
	//Return number of hypotheses
	virtual int getModels_int(const TSubSet & anHypSet, CModels & pModels);
public:
	virtual bool isThreadsafe() const { return true; } // as has no state

	COpenCV4ptHomography(const T2dPoints & p1, const T2dPoints & p2) : CImCorrModelHypothesiser(4, p1, p2) {}

        virtual int modelsPerIteration() const { return 1; }
        virtual int timePerIteration() const { return 100; /* times time to evaluate a hypothesis */ }
};

#endif /* OPENCV4PT_HOMOGRAPHY_H_ */
