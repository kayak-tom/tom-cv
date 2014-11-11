/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

/*
 * ransacParams.h
 *
 *  Created on: 19/10/2009
 *      Author: tom
 */

#ifndef RANSACPARAMS_H_
#define RANSACPARAMS_H_

#include "params/param.h"
#include "util/convert.h"

PARAMCLASS(RANSAC)
	PARAM(MAX_ITERS, 1, 1000000, 250, "Max number of RANSAC iterations (number of hypothesis sets chosen, may be multiple models per hyp. set). Normally fewer iterations needed")
	PARAM(PROB_SUCCESS, 0.8, 0.999, 0.98, "RANSAC parameter; terminate when this likely to have found inliers.")
	PARAM(E_INLIER_THRESH_PX, 0.001, HUGE, 3, "E inlier threshhold in pixels. Focal length from camera calibration matrix used to translate to image coordinates. Can be set huge to change RANSAC to a simple least-squares fit (a bit hacky)")
	PARAM(E_8PT_CUTOFF, 1, 1000000, 5, "Reject essential matrix candidates from 8pt algorithm where the singular values of F differ by this factor (they would be equal if F was a good estimate for E).")
	PARAM(TOPDOWN_ITERS, 0, 1000, 3, "Remove outliers, then re-fit model to remaining inliers this many times.")
	PARAME(RANSACSampler, BaySAC, "Sampling method")
	PARAME(HypothesiseAlg, 5PtE, "Choose a hypothesis-generation algorithm, also chooses whether to find Essential or Fundamental matrix")
	PARAM(TOPDOWN_EXPAND, .2, 10, 1.0, "Experimental: Scale up inlier threshhold when finding additional inliers (topdown refinement)")
	PARAM(TOPDOWN_SCALEDOWN, 0.1, 1.2, 1.0, "Experimental: Reduce inlier threshhold for each topdown refinement iteration")
	PARAME(RANSACTerminator, SimpleTerminator, "Algorithm to decide when to give up testing data points (when the inlier rate is low). WaldSAC-based terminators are too slow, NoTestTerminator or SimpleTerminator should be used.")
	PARAME(RANSACIterTerminator, ClassicRANSAC, "Algorithm to decide when to stop iterating (because an adequate solution has been found). F+B's paper includes a simple formula (ClassicRANSAC), otherwise terminate only after MAX_ITERS")
	PARAMB(VERBOSE, false, "Output debugging info to cout")
	PARAMB(TERMINATOR_LOG, false, "Output data on when the 'test' stage terminates (e.g. when using WaldSAC)")
	{}

	CNumParam<int> MAX_ITERS;
	CNumParam<double> PROB_SUCCESS;
	//CNumParamDerived<double> E_INLIER_THRESH; //This is now COMPUTED from the E_INLIER_THRESH_PX and the calibration data
	CNumParam<double> E_INLIER_THRESH_PX, E_8PT_CUTOFF;
	CNumParam<int> TOPDOWN_ITERS;
	MAKEENUMPARAM3(RANSACSampler, RANSAC, BaySAC, SimSAC);
	MAKEENUMPARAM11(HypothesiseAlg, 5PtE, 7PtE, 7PtF, 2PtE, MCE, 1PtE, 5Pt_GradientDesc, 4Pt_GradientDesc, 3Pt_GradientDesc, 7PtEFast, 5PtRt);
	CNumParam<double> TOPDOWN_EXPAND, TOPDOWN_SCALEDOWN;

	MAKEENUMPARAM5(RANSACTerminator, NoTestTerminator, SimpleTerminator, WaldSAC, BrownianBridge, BrownianBridgeLinear);
	MAKEENUMPARAM3(RANSACIterTerminator, MaxIters, ClassicRANSAC, TerminateOnPropInliers);

	CNumParam<bool> VERBOSE, TERMINATOR_LOG;
};

PARAMCLASS(RANSACHomography)
	PARAM(MAX_ITERS, 1, 1000000, 250, "Max number of RANSAC iterations for homography estimation (number of hypothesis sets chosen, may be multiple models per hyp. set). Normally fewer iterations needed")
	PARAM(PROB_SUCCESS, 0.8, 0.999, 0.98, "RANSAC parameter; terminate when this likely to have found inliers.")
	PARAM(H_INLIER_THRESH_PX, 0.0001, HUGE, 3, "H inlier threshhold in pixels.")
	PARAM(TOPDOWN_ITERS, 1, 1000, 3, "Remove outliers, then re-fit model to remaining inliers this many times.")
	PARAME(RANSACSampler, BaySAC, "Sampling method")
	PARAM(TOPDOWN_EXPAND, .2, 10, 2.0, "Experimental: Scale up inlier threshhold when finding additional inliers (topdown refinement)")
	PARAM(TOPDOWN_SCALEDOWN, 0.1, 1.2, 0.8, "Experimental: Reduce inlier threshhold for each topdown refinement iteration")
	PARAME(RANSACTerminator, WaldSAC, "Test termination method")
	PARAMB(VERBOSE, false, "Output to cout")
	{}

	CNumParam<int> MAX_ITERS;
	CNumParam<double> PROB_SUCCESS, H_INLIER_THRESH_PX;
	CNumParam<int> TOPDOWN_ITERS;
	MAKEENUMPARAM3(RANSACSampler, RANSAC, BaySAC, SimSAC);
	CNumParam<double> TOPDOWN_EXPAND, TOPDOWN_SCALEDOWN;
	MAKEENUMPARAM3(RANSACTerminator, SimpleTerminator, WaldSAC, BrownianBridge);

	CNumParam<bool> VERBOSE;
};
class CModelHypothesiser;
CModelHypothesiser * getHypothesiser(const CRANSACParams::eHypothesiseAlg alg, const T2dPoints & points0, const T2dPoints & points1, const double dUprightThresh);

#endif /* RANSACPARAMS_H_ */
