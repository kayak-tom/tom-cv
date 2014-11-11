/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * bowslamParams.h
 *
 *  Created on: 19/10/2009
 *      Author: tom
 */

#ifndef BOWSLAMPARAMS_H_
#define BOWSLAMPARAMS_H_

#include "params/param.h"
#include "util/convert.h"
#include "ransac/ransacParams.h"
#include "bow/bagOfWordsParam.h"
#include "image/imageAccess.h"
#include "description/vectorDescriptor.h"
#include "description/patchParams.h"
#include "featureExtract/cornerDetector.h"
#include "scale/bowSpeedoParams.h"
//Not used any more #include "structSizeParams.h"

//EXPAND_ALL(
WRAPPERPARAMCLASS(BoWSLAM)
	PARAM(START_FRAME, 0, MAX_INT/2, 0, "Delay navigation until this many frames are observed. Use for debugging to get BoW DB initialised or to start somewhere easier.")
	PARAM(MAX_TRACK_LEN, 10, 10000, 10000, "Delete me. If re-implemented us BoW matching")
	PARAM(MIN_TRACK_LEN_RECONSTRUCT, 0, 16, 2, "Not implemented (didn't help)")
	PARAM(MAX_TRIES_FINDING_E, 1, 5, 1, "Give RANSAC a second go at finding E if the first attempt doesn't lead to a camera being selected")
	PARAM(MIN_INLIERS_NEARBY, 8, 50, 20, "Matches with nearby frames need this many inlier correspondences")
	PARAM(MIN_INLIERS_DISTANT, 8, 100, 28, "Matches with far away frames (loop closure) need this many inlier correspondences")
	PARAM(NEARBY_TIME, 0, 250, 2, "Frames within this many timesteps are defined as nearby. See RECURSE_DEPTH.")
	PARAM(RECURSE_DEPTH, 1, 25, 2, "Frames successfully matched (linked) this deeply to the NEARBY_TIME frames are defined to be nearby. Works well even at variable speed because of multiple links.")
	PARAM(TARGET_CORRESPONDENCES, 0, 10000, 500, "Not implemented")
	PARAMB(CORRECT_RD, true, "Turn on/off radial distortion correction (set in cam calibration)")
	PARAM(TOTAL_CORES, 1, 64, 2 DEBUGONLY(-1 /* don't want to step thru MT stuff in Debug */), "Limits number of threads spawned in some circumstances (multithreaded RANSAC and calculating several links at once). Debug build defaults to lower val.")
	PARAM(FRAMERATE_MS, 1, 1, 1, "Not tested. Times==num of frames at the moment.") //Not sure if this will work...
	PARAM(DEPTH_THRESH, 0.01, 200, 0.2, "Threshhold where points become 'at infinity'. Points this far behind chosen cam also used." ) // Scale localisation error (in radians) by this to work out max depth before considered at infinity. High val will alllow many points reconstructed behind cam to be used.
	PARAMB(EXTRAP, false, "Try extrapolating a link or 2 before giving up and starting a new map component. Constant velocity MM, uninformative scale. Sometimes helps.")
	PARAM(READ_AHEAD_LIM, 0, 1000000, 3, "How far ahead can the image-loading thread get. SET LOW=0 FOR OUTPUTTING VIDEO FRAMES")
	PARAMB(MULTI_RUNS, false, "Internal flag turning off output--tuning or something is happening")
	PARAMB(TUNE_PARAMS, false, "Internal flag turning off output for tuning")
    CHILDCLASS(BOW, "Bag-of-Words")
	CHILDCLASS(BOWMatching, "Params for correspondences from Bag-of-Words")
	CHILDCLASS(PatchDescriptor, "Image patch selection params")
	CHILDCLASS(RANSAC, "RANSAC for E params")
	CHILDCLASS(LinkSelection, "How many frames to link to, etc.")
	CHILDCLASS(Im, "Image sizes, channels, calibration etc. Mostly set automatically.")
	CHILDCLASS(Corner, "Corner detection params")
	CHILDCLASS(RefineRT, "Levenberg-Marquardt refinement of rotation and translation params.")
	CHILDCLASS(ResolveScale, "Params for estimating scale from 3d point matches; Grubs test, etc.")
	CHILDCLASS(DescriptorSetClustering, "BoW dictionary generation params")
	CHILDCLASS(BOWSpeedo, "Object recognition params")
	//CHILDCLASS(StructureSize, "Once could apply heuristic constraint on scale of world, not implemented any more?")
	CHILDCLASS(TORO, "Params for outputting optimised TORO files")
	CHILDCLASS(Output, "Params for controlling rendering frequency, EPS output, etc. Rendering is more expensive than everything else combined!")
	CHILDCLASS(Optimise, "Params controlling optimisation around loops (in t-spanner)")
	CHILDCLASS(Mapping, "Params controlling local-global map generation")
	CHILDCLASS(ConstrainScale, "Params controlling constrained scale MM")
	{}

	CNumParam<int> START_FRAME, MAX_TRACK_LEN, MIN_TRACK_LEN_RECONSTRUCT, MAX_TRIES_FINDING_E, MIN_INLIERS_NEARBY, MIN_INLIERS_DISTANT, NEARBY_TIME, RECURSE_DEPTH, TARGET_CORRESPONDENCES;
	CNumParam<bool> CORRECT_RD;
	CNumParam<int> TOTAL_CORES, FRAMERATE_MS;
	CNumParam<double> DEPTH_THRESH;
	CNumParam<bool> EXTRAP;
	CNumParam<int> READ_AHEAD_LIM;
	CNumParamDerived<bool> MULTI_RUNS, TUNE_PARAMS;

	PARAMCLASS(ConstrainScale)
		PARAMB(CONSTRAIN_SCALE, false, "Heuristic constraint on scale of world")
		PARAM(BASELINE_TO_DEPTH_LNPARAM1, 0.01, 100, 4.4, "Scale mean point depth by this LN distn. to give baseline: param 1 (mean of underlying normal)" )
		PARAM(BASELINE_TO_DEPTH_LNPARAM2, 0.01, 100, 1.5, "Scale mean point depth by this LN distn. to give baseline: param 1 (var of underlying normal)" )
		{}

		CNumParam<bool> CONSTRAIN_SCALE;
		CNumParam<double> BASELINE_TO_DEPTH_LNPARAM1, BASELINE_TO_DEPTH_LNPARAM2;
	};

	PARAMCLASS(LinkSelection)
		PARAM(MAX_NUM_TO_TRY_LINKING, 1, 16, 4, "Max number of topologically nearby frames to *attempt to* link to the current frame")
		PARAM(MAX_NUM_TO_TRY_LINKING_LC, 0, 16, 2, "Max number of distant frames to *attempt to* link to the current frame (loop closure)")
		PARAM(NUM_TO_LINK, 1, 16, 3, "Max number of topologically nearby frames to link to the current frame")
		PARAM(NUM_TO_LINK_LC, 0, 16, 1, "Max number of distant frames to link to the current frame (loop closure)")
		PARAM(SIMILARITY_THRESH, 0, 1, 0.85, "Threshhold on distinctiveness needed to be considered LC candidate--lower val will mean more matches will meet threhhold to be considered LC candidates")
		PARAMB(DEBUG_DISABLE_LINKING, false, "For BoW matching debugging")
		PARAMB(ALLOW_ZERO_VELOCITY_LINKS, false, "Add links where camera hasn't moved significantly but has rotated. If added then scale can not be propogated through them. If not added then pure rotation is handled if this frame can be registered to another earlier frame.")
		PARAMB(KEEP_ALL_LINKS, true, "If true then only frames registered successfully to no others will be dropped. False helps keep map sparse, but misses rotations")
		PARAMB(VERBOSE, false, "List strongest BoW matches")
		{}
		CNumParam<int> MAX_NUM_TO_TRY_LINKING, MAX_NUM_TO_TRY_LINKING_LC, NUM_TO_LINK, NUM_TO_LINK_LC;
		CNumParam<double> SIMILARITY_THRESH;
		CNumParam<bool> DEBUG_DISABLE_LINKING, ALLOW_ZERO_VELOCITY_LINKS, KEEP_ALL_LINKS, VERBOSE;
	};

	PARAMCLASS(ResolveScale)
		//PARAME(ALIGNMENT_METHOD, 1dRANSAC, "Not implemented")
		PARAM(ALIGNMENT_THETA_MAX_ERR, 0.01, 0.9, 0.15, "Matches between 3D points considered outliers if the angle between points (in the same frame) differ by more than this amount. Would be 0 in noise-free case.")
		PARAM(MIN_INLIERS_1DALIGN, 3, 1000, 3, "Minimum number of matches between 3D point sets for scale estimate to be used")
		//PARAM(RANSAC_1D_ERR, 0.01, 10, 0.15)
		//PARAM(CONDITION_SCALE, 0.001, 1000, 1)
		PARAMB(MOTION_MODEL, false, "Can turn on a simple constant-velocity motion model, but doesn't seem to help even with auto parameter tuning.")
		PARAM(MOTION_MODEL_MAX_TIME, 2,1000,10, "Links from this long ago are used to estimate future motion")
		PARAM(MOTION_MODEL_SD, 0.00001, 1000, 0.1, "SD of future position estimate as a proportion of current velocity")
		//PARAM(MOTION_MODEL_ZERO_SD, 0.00001, 1000, 5, "")
		PARAM(PROB_NO_MOTION_MODEL, 0.0, 1.0, 0.05, "Probability of any position estimate incompatible with MM occurring (so essentially uniform+Gaussian MM)")
		PARAMB(VERBOSE, false, "MM Debugging output")
		{}

		//MAKEENUMPARAM3(ALIGNMENT_METHOD, 1dRANSAC, 1dMedian, 1dMean)
		CNumParam<double> ALIGNMENT_THETA_MAX_ERR;
		CNumParam<int> MIN_INLIERS_1DALIGN;
		CNumParam<bool> MOTION_MODEL;
		CNumParam<int> MOTION_MODEL_MAX_TIME;
		CNumParam<double> MOTION_MODEL_SD, PROB_NO_MOTION_MODEL;
		CNumParam<bool> VERBOSE;
	};

	PARAMCLASS(TORO)
		PARAMB(SAVE_TORO_MAP, true, "Save a file for input to TORO")
		PARAM(TORO_CONNECTIVITY, 1, 10000, 5, "t-Spanner t value for TORO map output")
		PARAMB(TWO_D, false, "Can output a 2d (X-Z) map to TORO instead--often works better. Full covariances could be used but aren't computed at the moment")
		PARAMB(EXTRA_UNINF_EDGES, false, "Can impose a Motion model by adding extra edges to TORO map. TODO: Investigate further for removing artifacts.")
		{}

		CNumParam<bool> SAVE_TORO_MAP;
		CNumParam<int> TORO_CONNECTIVITY;
		CNumParam<bool> TWO_D, EXTRA_UNINF_EDGES; //save 2d X-Z plane info
	};

	PARAMCLASS(Output)
		PARAMB(SAVE_EPS, true, "Save local map as EPS")
		PARAM(PLOT_INTERVAL, 1, 10000, 1, "Render+show every PLOT_INTERVAL'th frame")
		PARAMB(PRINT_SPEEDS, false, "Print the sequences of edges, and transforms, used to compute")
		PARAMB(PRINT_HEURISTIC_ERRORS, false, "Print warnings when the camera is stationary, fast, moving backwards or rotating out of XZ plane (use to identify failure points)")
		PARAMB(PRINT_POSITION_SOURCES, false, "Display link (yellow line) to most recent frame each frame is registered to")
		PARAMB(OPTIMISE_SCALES, false, "Use SVD to optimise scales around loops")
		PARAM(PRINT_FRAME_NUMS, 0, 10000, 10, "Display every PRINT_FRAME_NUMS'th frame number")
		PARAM(TEST_RD_CORRECTION, -1, MAX_INT, -1, "Use Radial distortion parameters to rectify the TEST_RD_CORRECTION'th frame. Output saved in output dir.")
		PARAM(KILL_FRAME, 0, MAX_INT, MAX_INT, "Force segfault, use with Valgrind to track mem in use")
		PARAMB(OUTPUT_CORR, false, "Save images showing inlier and outlier correspondences after doing RANSAC")
		PARAMB(PRINT_POS_SOURCES, false, "VERY SLOW: Print the sequence of relative poses and scales added to compute each position. Produces a huge amount of output!")
		{}

		CNumParam<bool> SAVE_EPS;
		CNumParam<int> PLOT_INTERVAL;
		CNumParam<bool> PRINT_SPEEDS, PRINT_HEURISTIC_ERRORS, PRINT_POSITION_SOURCES, OPTIMISE_SCALES;
		CNumParam<int> PRINT_FRAME_NUMS, TEST_RD_CORRECTION, KILL_FRAME;
		CNumParam<bool> OUTPUT_CORR, PRINT_POS_SOURCES;
	};

	PARAMCLASS(Mapping)
		PARAME(SET_ORIGIN, AtZero, "Generate either a map relative to the robot's current position, or one realtive to the first frame. Affects where gaps form. Can be moved around in code.")
		PARAMB(VERBOSE, false, "Dijkstra-update output")
		{}

		MAKEENUMPARAM2(SET_ORIGIN, AtZero, AtRobot)
		CNumParam<bool> VERBOSE;
	};

	PARAMCLASS(Optimise)
		PARAM(SPANNER_T, 2, 10000, 50, "t-Spanner t for subgraph used to optimise scale around loops")
		PARAMB(VERBOSE, false, "Print SVD matrix, corrections, etc.")
		{}

		CNumParam<int> SPANNER_T;
		CNumParam<bool> VERBOSE;
	};

	PARAMCLASS(RefineRT)
		PARAMB(ROBUST_COST, false, "Can use a robust cost fn where  we discount cost above a threshhold. Doesn't seem to help")
		PARAM(ROBUST_COST_THRESH, 0, 1, 0.01, "Thresh above which cost is discounted (for robust cost function)")
		PARAM(ROBUST_COST_SCALE, 0, 1, 0.1, "Scale down part of residual greater than thresh by this amount")
		PARAM(ROBUST_COST_CONDITION_THRESH, 0, 1.0, 0.00002, "If av residual exceeds this amount then outliers probably included. TODO: Try outlier removal") 
		PARAM(ROBUST_COST_CONDITION_SCALE, 0, 1000, 5, "If av residual too high then boost G_sq so unlikely to be used.") 
		PARAMB(VERBOSE, false, "Output condition data")
		{}

		CNumParam<bool> ROBUST_COST;
		CNumParam<double> ROBUST_COST_THRESH, ROBUST_COST_SCALE, ROBUST_COST_CONDITION_THRESH, ROBUST_COST_CONDITION_SCALE;
		CNumParam<bool> VERBOSE;
	};


	//BoWSLAM requires these parameterisation so instantiate inside
	MAKECHILDCLASS(BOW)
	MAKECHILDCLASS(BOWMatching)
	MAKECHILDCLASS(PatchDescriptor)
	MAKECHILDCLASS(RANSAC)
	MAKECHILDCLASS(LinkSelection)
	MAKECHILDCLASS(Im)
	MAKECHILDCLASS(Corner)
	MAKECHILDCLASS(RefineRT)
	MAKECHILDCLASS(ResolveScale)
	MAKECHILDCLASS(DescriptorSetClustering)
	MAKECHILDCLASS(BOWSpeedo)
	//MAKECHILDCLASS(StructureSize)
	MAKECHILDCLASS(TORO)
	MAKECHILDCLASS(Output)
	MAKECHILDCLASS(Optimise)
	MAKECHILDCLASS(Mapping)
	MAKECHILDCLASS(ConstrainScale)
};

#endif /* BOWSLAMPARAMS_H_ */
