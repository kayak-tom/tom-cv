/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

/*
 * ransac.h
 *
 *  Created on: 6/10/2009
 *      Author: tom
 */

#ifndef RANSAC_H_
#define RANSAC_H_

#include "util/Simple2dPoint.h"

#include "cRansacTerminator.h"
#include "ransacSampler.h" 
#include "ransacIterTerminator.h"
#include "hypothesiser.h"
#include "refiner.h"
#include "models.h"
#include "inlierCounter.h"
#include "findBestMatching.h"

class CRANSACParams;
class CRANSACHomographyParams;


//Returns number of inliers (0 on failure), fills dE with refined matrix (E11, E12, E13, E21....)
int doRansac(CSampler * pSampler,
		CModelHypothesiser * pHypothesise,
		CModelRefiner * pRefine,
		CRansacTerminator * pTerminator,
		CIterTerminator * pIterTerminator,
		CInlierCounter * pCounter,
		CFindBestMatching * pRefineMatches,
		CModel & bestModel, CMask & inliers, int nThreads=1, const bool bVerbose=false);

/* Function to do RANSAC to estimate H */
int getH(const T2dPoints & p1, const T2dPoints & p2, CInlierProbs & adPriorProbs, const CPointIdentifiers & pointIds,
		const CRANSACHomographyParams & PARAMS,
		C3x3MatModel & bestModel, CMask & inliers, const double dFocalLength /*1 if using pixel coordinates*/, const int nThreads);

int getF(const T2dPoints & points0, const T2dPoints & points1, CInlierProbs & adPriorProbs, const CPointIdentifiers & pointIds,
		const CRANSACParams & PARAMS,
		C3x3MatModel & F, CMask & inliers, const int nThreads);

/* Function to do RANSAC to estimate E
 * 
 * Set dUprightThresh = 0.2 or something to force camera to be upright and motion constrained to the x-z (ground) plane.
Normalised direction will have y less than dUprightThresh, and (0,1,0) will be rotated to at most (., 1-dUprightThresh, .)
0.2 works well, 0.1 normally too low.
If camera is exactly upright there are better algorithms for estimating rotation and translation.*/
int getE(const T2dPoints & p1, const T2dPoints & p2, CInlierProbs & adPriorProbs, const CPointIdentifiers & pointIds,
		const CRANSACParams & PARAMS,
		C3x3MatModel & bestModel, CMask & inliers, const double dUprightThresh /*Negative if not necessarily upright*/, const double dFocalLength, const int nThreads);

//Heuristic estimate of what the minimum inlier count from a good model would be.
//dMinPropInliersGood == min overlap between images?
//Maybe dMinPropInliersGood = 0.25 is realistic if there;s only 25% overlap
int getMinInliers(CInlierProbs & adPriorProbs, const double dMinPropInliersGood, const int nNumOfPointsForModel);

//Eigen 2 matrix library required (TIP: increase Eigen's Epsilon for faster eigenvector computation, remove complex EV code)
//Boost is used for threading, and for a fast pool allocator.
//#define _DEBUG to add extra debugging checks
void printModel(const CModel & model_in);

class C3dPoint;
class C3dRotationQuat;
typedef C3dRotationQuat C3dRotation;
//Recover cameras, test if points are in front or behind, update mask
bool RTFromE(CMask & mask, const C3x3MatModel & model, const T2dPoints & points0, const T2dPoints & points1, const CPointIdentifiers & ids, C3dRotation & R, C3dPoint & T, const double THRESH_SQ);
bool RTFromE(CMask & mask, const Eigen::Matrix3d & E, const T2dPoints & points0, const T2dPoints & points1, const CPointIdentifiers & ids, C3dRotation & R, C3dPoint & T, const double THRESH_SQ);
bool RTFromRT(CMask & mask, const T2dPoints & points0, const T2dPoints & points1, const CPointIdentifiers & ids, C3dRotation & R, C3dPoint & T, const double THRESH_SQ);


#endif /* RANSAC_H_ */
