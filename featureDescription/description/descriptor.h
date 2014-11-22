/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

#include <vector>
#include "util/exception.h"
#include "util/convert.h"
#include "util/dynArray.h"

pragma_warning (push)
pragma_warning (disable: 4127) //const conditional expression(from template params)

class CCluster;
class CClusterSet;

//typedef std::vector<CCluster> TClusterVector;
typedef CDynArrayOwner<CCluster> TClusterVector;

enum eNorm
{
	eCosine, e1Norm, eMaxNorm, eEuclidSquared, eEuclidParallel
}; //for SURF descriptors

enum eKMImplementation
{
	eKMLL, eLloyds, eKMWeighted,//Use my own SAD (??) k-means implementation
	eKMFixedRadius
}; //clusters truncated at a fixed radius

#define MAX_ALLOWED_DIST MAX_INT

#ifdef _DEBUG
#define CAST dynamic_cast
#else
#define CAST static_cast
#endif

#include "util/location.h"

//Abstract class for any descriptor that can be clustered
enum eInvDescriptorType
{
	eSIFT, eSURF, eNone
};

class CCluster;

class CDescriptor
{
	const CCluster * pCluster;
public:
	typedef int TDist;
private:
	TDist closestClusterDist;
public:

	virtual TDist distance(const CDescriptor * pd) const = 0;

	virtual double orientation() const { return 0; }

	//Allow a distance function that is not rotationally invariant, for the last few frames.
	virtual TDist orientedDistance(const CDescriptor * pd, const double OI_NEARBY_ANGLE) const
	{
		if(diffMod(orientation(), pd->orientation(), 360) > OI_NEARBY_ANGLE)
			return MAX_ALLOWED_DIST-1;
		else
			return distance(pd);
	}

	virtual CLocation location() const
	{
		return CLocation(0, 0);
	}

	CDescriptor() : pCluster(0) {}

	virtual ~CDescriptor()
	{
	}

	static int margin()
	{
		return 24;
	}

	virtual bool drawable() const { return false; }
	virtual unsigned char val(int x, int y, int nChannel) const { THROW("This descriptor can't be drawn as an image"); };
	virtual int diameter() const { THROW("This descriptor can't be drawn as an image"); };

	/*virtual*/ inline const CCluster * assignment() const { return pCluster; }
	/*virtual*/ inline const TDist assignmentDist() const { return closestClusterDist; } // MAX_ALLOWED_DIST denotes distance unknown
	/*virtual*/ void assignToCluster(const CCluster * pCentre, TDist closestClusterDist_in);

	//Dodgy k-means access functions
	virtual int size() const = 0;
	virtual int length() const = 0;
	virtual char * ptr() const { return ((char *)(void *)this)+sizeof(void*) /* virtual function ptr */; }; //breaks const

	CDescriptor * clone() const
	{
		CDescriptor * pClone = (CDescriptor *)malloc(size());
		memcpy((void *)pClone, this, size());

		if(IS_DEBUG) CHECK(distance(pClone) != 0, "Clone failed");
		if(IS_DEBUG) CHECK(!(pClone->location() == location()), "Clone failed");

		return pClone;
	}
};

//typedef std::vector<CDescriptor *> TDescriptorVector;
typedef CDynArray<CDescriptor *> TDescriptorVector;
typedef CDynArray<const CDescriptor *> TConstDescriptorVector;

class TMPSet;

class CMatchableDescriptors
{
public:
	class CMatchSettings
	{
		friend class CMatchableDescriptors;
                friend class CGuidedFeatureMatching;
	protected:
		const double MATCH_CONDITION; //Match if better than MATCH_CONDITION * next best match
		const double PP; //Max prior prob
		const int NN; //Max N for N-M
		const int MAX_SEPERATION; //in px. -1 turns off
		const int NEARBYNESS; //Currently -1 for not nearby (use rotation invariance if available, 1 for nearby
		const double OI_NEARBY_ANGLE;
	public:
		CMatchSettings(const double MATCH_CONDITION=0.8, const double PP=0.6, const int NN=1, const int MAX_SEPERATION=-1, const int NEARBYNESS=-1, const double OI_NEARBY_ANGLE=30)
		 : MATCH_CONDITION(MATCH_CONDITION), PP(PP), NN(NN), MAX_SEPERATION(MAX_SEPERATION), NEARBYNESS(NEARBYNESS), OI_NEARBY_ANGLE(OI_NEARBY_ANGLE)
		{}
	};
protected:
    const CBoWCorrespondences * getBruteForceCorrespondenceSet_int(const CMatchableDescriptors * pDS, const CMatchSettings & MS, CBoWCorrespondences * pCorr) const HOT;

	template<bool LIMIT_DISTANCES, bool ROTATION_INVARIANCE>
	void getBFC_matchDescriptors_RI(const int nStart, const int nEnd, const CMatchableDescriptors * pDS, const CMatchSettings & MS, TMPSet * d) const;
	template<bool LIMIT_DISTANCES>
	void getBFC_matchDescriptors(const int nStart, const int nEnd, const CMatchableDescriptors * pDS, const CMatchSettings & MS, TMPSet * d) const;
	void getBFC_matchDescriptorsMT(const int nStart, const int nEnd, const CMatchableDescriptors * pDS, const CMatchSettings & MS, TMPSet * d) const;
public:
	const CBoWCorrespondences * getBruteForceCorrespondenceSet(const CMatchableDescriptors * pDS, const CMatchSettings & MS, CBoWCorrespondences * pCorr) const;
	virtual int Count() const = 0;
	inline const CDescriptor * get_const(int i) const { return (*this)[i]; }
	virtual const CDescriptor * operator[](int i) const { return get_const(i); };

	virtual void Clear() = 0;
};

class CSimpleMatchableDescriptors: public CMatchableDescriptors
{
	TConstDescriptorVector descriptors;
public:
	CSimpleMatchableDescriptors() { descriptors.reserve(32); }
	virtual const CDescriptor * operator[](int i) const { return descriptors[i]; };
	void Push_const(const CDescriptor * pDescriptor);
	void Push_const(const CMatchableDescriptors * pDescriptorSet);
	virtual void Clear() { descriptors.clear(); }
	virtual int Count() const { return descriptors.size(); };
};

class CDescriptorSet : public CMatchableDescriptors
{
private:
	bool Contains(const CDescriptor * pDescriptor) const;

protected:
	//Adds a random subset of this set to pRandSubset, include centres if given
	void RandomSubset(int n, CDescriptorSet * pRandSubset, CClusterSet * pBestClusters = 0, bool bForceSeperation = false);

	void clearClusterAssignments();
public:
	//A function that takes a vector of descriptors and clusters them into k clusters.
	virtual CClusterSet * Cluster(int nClusters, const CCluster * pParentCluster) = 0;

	virtual void Push(CDescriptor * pDescriptor) = 0;
	virtual void Push(CDescriptorSet * pDescriptorSet) = 0;
	//virtual void Clear() = 0;

	virtual CDescriptor * operator[](int i) = 0;
	inline CDescriptor * get(int i) { return (*this)[i]; }

	virtual ~CDescriptorSet()
	{
	}

	const CDescriptor * ClosestDescriptor(const CDescriptor * pDesc) const;
	void Closest2Descriptors(const CDescriptor * pDescCompareTo, CDescriptor const** ppClosest, CDescriptor const** pp2ndClosest) const;
	const CDescriptor * ClosestDescriptor(const CDescriptor * pDesc, double dCondition, int nRadius) const;
	CDescriptor::TDist Closest(const CDescriptor * pDesc) const;
	void descriptorQuality(const CDescriptor * pDescThis, const CDescriptor * pDescOther, double & dQuality, int & nCloser) const;

	virtual CDescriptorSet * makeNewDS(int nEstimatedSize = -1) const = 0;

	//Also deletes all the descriptors
	static void deleteDS(CDescriptorSet ** ppDescriptors)
	{
		deleteDS((const CDescriptorSet **)ppDescriptors);
	}

	static void deleteDS(const CDescriptorSet ** ppDescriptors)
	{
		if(IS_DEBUG) CHECK(!ppDescriptors, "deleteDS: Null parameter");
		const CDescriptorSet * pDescriptors = *ppDescriptors;
		if(pDescriptors)
		{
			for(int i=0; i<pDescriptors->Count(); i++)
				delete pDescriptors->get_const(i); //the descriptor

			delete pDescriptors; *ppDescriptors = 0;
		}
	}
	//Remind these descriptor which cluster they are assigned to
	void assignToCluster(const CCluster * pCentre)
	{
		for(int i=0; i<Count(); i++)
			get(i)->assignToCluster(pCentre, MAX_ALLOWED_DIST);
	}
};

#include "descriptorClustering.h"

pragma_warning (pop)
