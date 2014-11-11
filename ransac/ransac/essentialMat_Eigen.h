/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */
/* This code calls calibrated_fivepoint_helper */

// 5-point algorithm and notes in essentialMat_Eigen.cpp file

#pragma once

/*
 * essentialMat_Eigen.h
 *
 *  Created on: 8/04/2009
 *      Author: tom
 */

#ifndef ESSENTIALMAT_EIGEN_H_
#define ESSENTIALMAT_EIGEN_H_

#include "hypothesiser.h"

class C5ptEssentialMat : public CImCorrModelHypothesiser 
{
protected:
	const double dUprightThresh; // If > 0 only upright matrices (within this thresh) returned.
	//Return number of hypotheses
	virtual int getModels_int(const TSubSet & anHypSet, CModels & pModels);
public:
	virtual bool isThreadsafe() const { return true; } // as has no state

	C5ptEssentialMat(const T2dPoints & p1, const T2dPoints & p2, double dUprightThresh=-1) : CImCorrModelHypothesiser(5, p1, p2), dUprightThresh(dUprightThresh) {}

        virtual int modelsPerIteration() const { return 4; /*average*/ }
        virtual int timePerIteration() const { return 2500; /* times time to evaluate a hypothesis */ }
};

class C7ptEssentialMat_GS : public CImCorrModelHypothesiser 
{
protected:
	//Return number of hypotheses
	virtual int getModels_int(const TSubSet & anHypSet, CModels & pModels);
public:
	virtual bool isThreadsafe() const { return true; } // as has no state

	C7ptEssentialMat_GS(const T2dPoints & p1, const T2dPoints & p2) : CImCorrModelHypothesiser(7, p1, p2) {}

        virtual int modelsPerIteration() const { return 1; /*average*/ }
        virtual int timePerIteration() const { return 60; /* times time to evaluate a hypothesis */ }
};

class C2ptMCEssentialMat : public CImCorrModelHypothesiser //Model==essential matrix from >=8 points
{
	const int nHypotheses;
protected:
	//const double dUprightThresh; // If > 0 only upright matrices (within this thresh) returned.

	//Return number of hypotheses
	virtual int getModels_int(const TSubSet & anHypSet, CModels & pModels);
public:
	virtual bool isThreadsafe() const { return true; } // as has no state

	C2ptMCEssentialMat(const T2dPoints & p1, const T2dPoints & p2, const int nHypotheses) : CImCorrModelHypothesiser(2, p1, p2), nHypotheses(nHypotheses) {}

        virtual int modelsPerIteration() const { return nHypotheses; /*average*/ }
        virtual int timePerIteration() const { return 10; /* time time to evaluate a hypothesis */ }
};

//Calculate a basis E from the 5 points with indices in anHypSet, by SVD
void makeBasisForE(const TSubSet & anHypSet, const T2dPoints & m0, const T2dPoints & m1, Eigen::Matrix<double, 9, 4> & EE) HOT;
int calcEssentialMat_7point(const TSubSet & anHypSet, const T2dPoints & m0, const T2dPoints & m1, T3x3MatModels & models, int & FAIL_COUNT);

class CMCEssentialMat : public CImCorrModelHypothesiser 
{
	int nHypotheses;
protected:
	//const double dUprightThresh; // If > 0 only upright matrices (within this thresh) returned.

	//Return number of hypotheses
	virtual int getModels_int(const TSubSet & anHypSet, CModels & pModels);
public:
	virtual bool isThreadsafe() const { return true; } // as has no state

	CMCEssentialMat(const T2dPoints & p1, const T2dPoints & p2, const int nHypotheses) : CImCorrModelHypothesiser(0, p1, p2), nHypotheses(nHypotheses) {}

        virtual int modelsPerIteration() const { return nHypotheses; /*average*/ }
        virtual int timePerIteration() const { return 15; /* time time to evaluate a hypothesis */ }
};

//Calculate E from the 5 points with indices in anHypSet
int calcEssentialMat_5point_Eigen(const TSubSet & anHypSet, const T2dPoints & m0, const T2dPoints & m1, T3x3MatModels & pdEssentialMat, const double dUprightThresh, int & FAIL_COUNT, const bool bUseEngels) HOT;
//void makeClosestE(Eigen::Matrix3d & E);
double getResids(const Eigen::Matrix3d & E, const TSubSet & anHypSet, const T2dPoints & m0, const T2dPoints & m1);
//double getResids(const Eigen::Matrix3d & E, const Eigen::Vector3d * am0, const Eigen::RowVector3d * am1);

#endif /* ESSENTIALMAT_EIGEN_H_ */
