#pragma once

#include "params/param.h"
#include "util/stats.h"

PARAMCLASS(EssentialMatGradientDesc)
	PARAM(MAX_ITERS,1,1000,20,"Max iterations before failure")
	PARAM(MAX_FALSE_MINIMA, 1, 1000, 4, "How many times we land in false minima before giving up")
	PARAM(MAX_NUM_SUCCESSFUL_ATTEMPTS,1,100,10, "Num of solutions to attempt to find (less will normally be found as some will be the same)")

	PARAM(DROPOUT_EARLY_ON_FALSE_MINIMA, 0.001, 1000, 1, "Dropout if minima cannot be hit with remaining_iters*this_step_size reduction in error, scaled by this parameter. Higher = dropout later, lower = dropout earlier (low sensitivity).")
	PARAM(DIST_TO_EXISTING_SOLN_THRESH, 0.0, 5.0, 0.25, "Dropout if new solution passes this close to an existing solution. Higher values may miss some minima, but increase rate that solutions are found.")
	//PARAM(FAILURE_PROB_THRESH, 0.0, 1.0, 0.25, "Give-up on this hyp early if the probability of the next run finding a unique solution is below this value")
	PARAM(MAX_SOLUTIONS, 1, 10, 10, "Give-up on this hyp once this number of solutions have been found")
	PARAME(GDAlg,LevenbergMarquardt,"Which gradient descent algorithm to use")
	PARAMB(VERBOSE,false,"Output # iters")
	CHILDCLASS(LM, "Params specific to Levenberg-Marquardt")
	CHILDCLASS(SD, "Params specific to steepest descent and CGD")
	{}
	
	CNumParam<int> MAX_ITERS, MAX_FALSE_MINIMA, MAX_NUM_SUCCESSFUL_ATTEMPTS;
	CNumParam<double> DROPOUT_EARLY_ON_FALSE_MINIMA, DIST_TO_EXISTING_SOLN_THRESH;//, FAILURE_PROB_THRESH;
	CNumParam<int> MAX_SOLUTIONS;
	MAKEENUMPARAM4(GDAlg, LevenbergMarquardt, SteepestDescent, CGD, DownhillSimplex); 
	CNumParam<bool> VERBOSE;

	PARAMCLASS(LM)
		PARAM(LAMBDA,-1000,1000,0.15,"Initial damping factor") //About 0.15 for LEV, 0.01 for MAR
		PARAM(BOOST,1,1000,5,"Damping factor boost")
		PARAM(DROP,0.0000001,1,0.45,"Damping factor drop")
		PARAMB(USE_HESSIAN,false,"Compute full Hessian")
		PARAMB(USE_MARQUARDT,false,"Levenberg-Marquardt damping (diag *= 1+lambda), rather than Levenberg (diag += lambda)")
		{}

		CNumParam<double> LAMBDA, BOOST, DROP;
		CNumParam<bool> USE_HESSIAN, USE_MARQUARDT;
	};

	MAKECHILDCLASS(LM);

	PARAMCLASS(SD)
		PARAM(ALPHA,0.00001,1000,1,"(UNUSED) Initial step size (* size of deriv)") //About 0.15 for LEV, 0.01 for MAR
		PARAM(SMALL_STEP_DROPOUT,1e-20,0.1,1e-13,"We're at a min if steps are this small") //About 0.15 for LEV, 0.01 for MAR
		PARAME(BETA_FN, FletcherReeves, "Function for selecting beta (see Wikipedia)")
		{}

		CNumParam<double> ALPHA, SMALL_STEP_DROPOUT ; 
		MAKEENUMPARAM3(BETA_FN, FletcherReeves, PolakRibiere, HestenesStiefel);
	};

	MAKECHILDCLASS(SD);

};


class CGDStats
{
public:
	enum eTerminationResult { eSuccess, eFalseMin, eRepeatedMin, TERMINTION_TYPES};
	virtual void add(int nIters, eTerminationResult result)
	{

	}

	virtual void reset() {}
        virtual void writeResults(CParam & alg, const char * filename) {}
};

int findE_LMforBasis(const CEssentialMatGradientDescParams & PARAMS, const TSubSet & anHypSet, const T2dPoints & m0, const T2dPoints & m1, TEModels & models, CGDStats & stats) HOT;
int findE_LMforBasis3(const CEssentialMatGradientDescParams & PARAMS, const TSubSet & anHypSet, const T2dPoints & m0, const T2dPoints & m1, T3x3MatModels & models, CGDStats & stats) HOT;
int findE_LMforBasisR3(const CEssentialMatGradientDescParams & PARAMS, const TSubSet & anHypSet, const T2dPoints & m0, const T2dPoints & m1, TEModels & models, CGDStats & stats) HOT;
int findE_LMforBasis_TopRow(const CEssentialMatGradientDescParams & PARAMS, const TSubSet & anHypSet, const T2dPoints & m0, const T2dPoints & m1, T3x3MatModels & models, CGDStats & stats) HOT;

template<typename TM> void EigenNormalise(TM & M) { M /= M.stableNorm(); }