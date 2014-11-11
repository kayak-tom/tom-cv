/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/* Classes derived from the CModelHypothesiser class compute RANSAC hypotheses from minimal
 * point sets, i.e. all or most of the Essential matrices compatible with sets of 5 or 7 points.
 * 
 * Choose a hypothesis algorithm by setting RANSAC.HypothesisAlgorithm, or use the derived classes 
 * individually. E.g.:
 *

 const T2dPoints points0, points1; //2 arrays (CDynArray's) of N 2D image points.

 CEssentialMatLMbasis levMarBasisHypothesisGenerator(points0, points1, 5)); //Make a 5-point hypothesis-generator (derived from CModelHypothesiser)

 
 TSubSet anHypSet; // An array of n (e.g. n=5) integers, which are the indices of the 
                   // points in points0, points1 forming the hypothesis set.
 ...

 CModels * pModels = levMarBasisHypothesisGenerator.getModels(anHypSet);
 // A set of the models (essential matrices) is returned (see models.h)
 

 *
 *
 *
 *
 */

#pragma once

#ifndef HYPOTHESISER_H_
#define HYPOTHESISER_H_

#include "models.h"
#include "util/Simple2dPoint.h"

// Generates a set of models compatible with a subset of datapoints. Base class for Essential matrix estimators,
// homography estimators, etc.
class CModelHypothesiser 
{
protected:
	const int nPoints; //Minimum number of points from which a hypothesis can be generated

public:
	virtual bool isThreadsafe() const { return false; } //False by default

	CModelHypothesiser(int nPoints) : nPoints(nPoints)
	{
		if(IS_DEBUG) CHECK(nPoints<0, "CModelHypothesiser::fitModel: Invalid num points");
	}

	//Return a set of models compatible with the hypothesis set anHypSet (a set of indices for data points).
	virtual const CModels * getModels(const TSubSet & anHypSet) = 0;

	int numPoints() const { return nPoints; }

    virtual int modelsPerIteration() const = 0; //These are used to compute the threshholds used in WaldSAC
    virtual int timePerIteration() const = 0; // In terms of the cost per test, i.e. 16ns for Sampson's distance 
};

class CImCorrModelHypothesiser : public CModelHypothesiser
{
protected:
	const T2dPoints & p1;
	const T2dPoints & p2;

	//Return number of hypotheses
	virtual int getModels_int(const TSubSet & anHypSet, CModels & pModels) = 0;
public:
	CImCorrModelHypothesiser(int nPoints, const T2dPoints & p1, const T2dPoints & p2) : CModelHypothesiser(nPoints), p1(p1), p2(p2)
	{
		if(IS_DEBUG) CHECK(nPoints<0, "CModelHypothesiser::fitModel: Invalid num points");
	}

	//Return models
	virtual const CModels * getModels(const TSubSet & anHypSet);
};

#endif /* HYPOTHESISER_H_ */
