/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * bagOfWords.h
 *
 *  Created on: 15/10/2009
 *      Author: hilandtom
 */

#ifndef BAGOFWORDSPARAM_H_
#define BAGOFWORDSPARAM_H_

#include "params/param.h"
#include "util/convert.h"
#include "util/Simple2dPoint.h" //for NN_MAX

#define MAX_BOW_LEVELS 10 //32=at least 4 billion words at 2/level
#define LOCATION_STORE_LIM NN_MAX //don't store a list of locations longer than this

//EXPAND_ALL(
PARAMCLASS(BOW)
	PARAME(COMP_METHOD, NisterDist, "Vector distance between word bags. Brute-force compares descriptors directly (slow, no better)")
	PARAME(WB_WEIGHT_METHOD, TF_IDF, "There's a few different TF-IDF weighting functions in the literature. Low sensitivity. Either TF_IDF, from Nister-Stewenius-2006; DF_ITDF (variation with total occurances); TF_IDF_Wikipedia (the one on Wikipedia)")
	PARAM(RWB_THREADS, 1, 64, 2 DEBUGONLY(-1), "Currently blocking, so set num threads high normally.")
	PARAM(QUERY_THREADS, 1, 64, 1, "Can speed-up queries, but fast anyway.")
	CHILDCLASS(BOWClustering, "Parameters for creating hierarchical dictionary (codebook/vocabulary)")
	CHILDCLASS(DescriptorBinning, "Parameters for assigning descriptors to clusters")
	CHILDCLASS(BOWCorrespondenceProb, "Params for assigning prior probabilities to correspondences.")
	CHILDCLASS(FastVectorComparison, "Approximate speeded-up query params")
		{}

	MAKEENUMPARAM3( COMP_METHOD, NisterDist, BruteForce, VectorDistFast ); //Watch dont use more than MAX_NUM_COMP_METHODS
	MAKEENUMPARAM3(WB_WEIGHT_METHOD, TF_IDF, //From Nister-Stewenius-2006
									 DF_ITDF,  //log(1/total occurance count), slightly worse
									 TF_IDF_Wikipedia); //Different TF_IDF calc, from Wikipedia: http://en.wikipedia.org/wiki/Tf-idf
	CNumParam<int> RWB_THREADS;
	CNumParam<int> QUERY_THREADS;

	PARAMCLASS(BOWClustering)
		PARAM(LEVELS, 1, MAX_BOW_LEVELS, 3, "Levels in hierarchical dictionary")
		PARAMB(CLUSTER_IN_SEPERATE_THREAD, false, "Whether to cluster (build a new, better, dictionary) in seperate thread. Necessary for rebuilding dictionary in RT. Probably not wanted if only building dictionary once.")
		PARAM(DESCRIPTORS_PER_WORD, 5, 100, 15, "Controls dictionary size (word count). E.g. 10k images, 300 descriptors/image, DESCRIPTORS_PER_WORD=30 would give 100k words. Sensitivity fairly low, best value unknown (effects matching and recognition differently).")
		PARAM(RECLUSTER_FREQUENCY, -1, 10, -1, "If != -1, a new dictionary will be created (e.g. in background) every time number of descriptors in total grows by this factor.")
		PARAME(BRANCH_METHOD, FixedClusterSizeTarget, "How to choose centre count for sub-clusterings. Sensitivity low")
		PARAM(BF_LEVEL, 0, MAX_BOW_LEVELS, 3, "Level of dictionary used for descriptor partitioning for BoW feature matching") //0==global BF matching. BF_LEVEL <= LEVELS enforced
		{}

		CNumParam<int> LEVELS;
		CNumParam<bool> CLUSTER_IN_SEPERATE_THREAD;
		CNumParam<int> DESCRIPTORS_PER_WORD;
		CNumParam<double> RECLUSTER_FREQUENCY;
		MAKEENUMPARAM2(BRANCH_METHOD, FixedBF, FixedClusterSizeTarget);
	//private: todo: restore--currently hacked to test matching
		CNumParam<int> BF_LEVEL; //Private param only accessible for derived parameters
	public:
		int BRUTEFORCE_MATCHING_LEVEL() const
		{
			static bool bWarned = false;
			if(!bWarned && BF_LEVEL > LEVELS)
			{
				std::cout << "WARNING: BF_LEVEL>LEVELS; setting BF_LEVEL=LEVELS=" << (int)LEVELS << std::endl;
				bWarned = true;
			}
			return std::min<int>(BF_LEVEL, LEVELS);
		}
	};

	MAKECHILDCLASS(BOWClustering)

	PARAMCLASS(DescriptorBinning)
		PARAM(CLUSTER_THREADS, 1, 64, 1, "MT clustering. Normally used if CLUSTER_IN_SEPERATE_THREAD=false, and waiting for new dictionary.")
		PARAM(LOWER_BOUND, 1, 64, 1, "Clusters with fewer than this number of descriptors not used ('stopwords'). Not very sensitive.")
		PARAM(UPPER_BOUND, 1, MAX_INT, MAX_INT, "Clusters bigger than UPPER_BOUND% of average size not used (indistinctive). Not very sensitive.")
		PARAM(RADIUS, 1, MAX_INT, MAX_INT, "Max radius within which descriptors will be mapped to a centre. Inferred and set based on descriptor used (PX_RADIUS) and comparison norm.")
		{}
		CNumParam<int> CLUSTER_THREADS;
		CNumParam<int> LOWER_BOUND;
		CNumParam<int> UPPER_BOUND;
		CNumParamDerived<int> RADIUS; //DO NOT SET! Inferred and set based on descriptor used (PX_RADIUS)
	};
	MAKECHILDCLASS(DescriptorBinning)

	PARAMCLASS(BOWCorrespondenceProb)
		PARAM(MIN_PRIOR, 0, 1, 0.4, "Not used?? Prior prob of poorest 1-1 matches")
		PARAM(MAX_PRIOR, 0.01, 1.0, 0.8, "Prior prob of a 1-1 correspondence. A 2-1 correspondence will have prob MAX_PRIOR/2 etc.")
		{}
		CNumParam<double> MIN_PRIOR, MAX_PRIOR;
	};

	MAKECHILDCLASS(BOWCorrespondenceProb)

	PARAMCLASS(FastVectorComparison)
		PARAM(SUBVEC_CACHE_SIZE, 5, 500, 60, "Todo: higher=more accurate, lower=faster, low sensitivity")
		PARAM(VEC_MIN_LENGTH, 10, 200, 60, "Just look at first VEC_MIN_LENGTH words for first-pass match")
		PARAM(SUBVEC_LENGTH, 10, 200, 60, "Just look at first VEC_MIN_LENGTH words for first-pass match")
		{}
		CNumParam<int> SUBVEC_CACHE_SIZE, VEC_MIN_LENGTH, SUBVEC_LENGTH;
	};

	MAKECHILDCLASS(FastVectorComparison)
};

PARAMCLASS(BOWMatching)
	PARAM(MATCH_NN, 1, LOCATION_STORE_LIM, 3, "Find N-M correspondences with N,M <= MATCH_NN. High (4-12) usually good.")
	PARAM(BF_CORNER_CONDITION, 0.1, 1.0, 0.8, "When finding correspondences by brute-force matching, return matches within this amount of the next closest match.")
	PARAM(OI_NEARBY_ANGLE, 1, 180, 30, "Experimental: only match descriptors with extracted orientations within this thresh when matching nearby")
	PARAME(BF_CORRESPONDENCES, BoWCorrespondences, "Matching method: full brute-force matching of all descriptors, BF matching at course partititioning (BF_LEVEL), BoW matching (all co-occuring words give correspondences)")
	{}
	CNumParam<int> MATCH_NN;
	CNumParam<double> BF_CORNER_CONDITION, OI_NEARBY_ANGLE;
	MAKEENUMPARAM4(BF_CORRESPONDENCES, BoWCorrespondences, BoW_BF_Correspondences, OldBoW_BF_Correspondences, BF_Correspondences);
};

#endif /* BAGOFWORDSPARAM_H_ */
