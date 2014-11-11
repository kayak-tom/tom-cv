/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

//Todo: Rename to ClusteringDS.h or something

/*pragma_warning (push)
pragma_warning (disable: 4996)
#include "boost/multi_array.hpp"
pragma_warning (pop)*/
#include <iostream>
#include "params/param.h"
#include <boost/smart_ptr.hpp>

PARAMCLASS(DescriptorSetClustering)
	PARAME(ClusteringAlg, CLARA_KMedoids, "Clustering algorithm")
	CHILDCLASS(CLARA_KMedoids, "Parameters for CLARA_KMedoids clustering of descriptors")
	CHILDCLASS(KMeans, "Parameters for kmeans clustering of descriptors")
	CHILDCLASS(Random, "Parameters for kmeans clustering of descriptors")
	{}

	PARAMCLASS(CLARA_KMedoids) //Todo tidy up, add k-means/seperate CLARA
		PARAM(MAX_DESCRIPTORS_KM, 10, 1000, 200, "Max descriptors in each subset that is clustered. Higher=slower and slightly more accurate (cubic complexity)")
		//PARAM(MAX_DESCRIPTORS_KM_TOP)
		PARAM(KMEDOIDS_ITERS, 1, 100, 2, "Iterations in k-medoids alg")
		PARAM(CLARA_ITERS, 1, 100, 1, "Cluster centres are chosen, added to a new subset, refined this many times.")
		PARAM(SPARSE_ASSIGN_COUNT, 1, 100, 1, "Experimental. Speed up intermediate assignments during CLARA clustering.")
		{}

		CNumParam<int> MAX_DESCRIPTORS_KM, /*MAX_DESCRIPTORS_KM_TOP,*/ KMEDOIDS_ITERS, CLARA_ITERS, SPARSE_ASSIGN_COUNT;
	};

	PARAMCLASS(KMeans) //Todo tidy up, add k-means/seperate CLARA
		PARAM(SUBSET_SIZE_CLARA_KMEANS, 10, 10000, 200, "Max descriptors in each subset that is clustered. Higher=slower")
		PARAM(KMEANS_ITERS, 1, 100, 2, "Iterations in k-means alg")
		{}

		CNumParam<int> SUBSET_SIZE_CLARA_KMEANS, KMEANS_ITERS;
	};

	PARAMCLASS(Random) //Todo tidy up, add k-means/seperate CLARA
		PARAM(SUBSET_SIZE_CLARA_RANDOM, 10, 10000, 200, "Max descriptors in each subset that is clustered. Higher=slower")
		{}

		CNumParam<int> SUBSET_SIZE_CLARA_RANDOM;
	};

	MAKEENUMPARAM5(ClusteringAlg, CLARA_KMedoids, CLARA_KMeans, KMeans, RandomCentres, CLARA_RandomCentres);

	MAKECHILDCLASS(CLARA_KMedoids);
	MAKECHILDCLASS(KMeans);
	MAKECHILDCLASS(Random);

	int SUBSET_SIZE(int k) const
	{
		if(ClusteringAlg == eCLARA_KMedoids)
			return CLARA_KMedoids.MAX_DESCRIPTORS_KM;//TODO: Function of branch factor (plot in thesis)
		else if(ClusteringAlg == eCLARA_KMeans)
			return KMeans.SUBSET_SIZE_CLARA_KMEANS;//TODO: Function of branch factor (plot in thesis)
		else if(ClusteringAlg == eCLARA_RandomCentres)
			return k; //Doesn't matter as long as it's higher than the number of centres and lower than the number of
		else
			return MAX_INT;
	}
};

//Macroify template function/class names
#define CTVectorSpaceDescriptor CVectorSpaceDescriptor<elementType, nDescriptorLength, eNormToUse>
#define CTVectorSpaceDescriptorSet CVectorSpaceDescriptorSet<descriptorType, eKMeans>
#define CTLargeVSDescriptorSet CLargeVSDescriptorSet<descriptorType, eKMeans>
#define templateVS template<class elementType, int nDescriptorLength, eNorm eNormToUse>
#define templateVSSet template<class descriptorType, eKMImplementation eKMeans>

#define CTVectorLocationDescriptor CVectorLocationDescriptor<elementType, nDescriptorLength, eNormToUse>

//Class wrapping an array of any suitable type, that can be clustered using k-means
//Still virtual--distance metric should be implemented
templateVS
class CVectorSpaceDescriptor : public CDescriptor
{
	elementType aDescriptorVector[nDescriptorLength];
	//TDist anTwoTimesSumSquare[eNormToUse == eEuclidParallel];

	//Cosine dist for unit vectors
#ifndef __GNUC__
    friend double CTVectorSpaceDescriptor::unitCosineDistance(const CTVectorSpaceDescriptor * d1, const CTVectorSpaceDescriptor * d2);
    friend double CTVectorSpaceDescriptor::euclidDistance(const CTVectorSpaceDescriptor * d1, const CTVectorSpaceDescriptor * d2);
#endif
	static inline CDescriptor::TDist unitCosineDistance(const CTVectorSpaceDescriptor * d1, const CTVectorSpaceDescriptor * d2);
	static inline CDescriptor::TDist euclidDistance(const CTVectorSpaceDescriptor * d1, const CTVectorSpaceDescriptor * d2);

public:
    CDescriptor::TDist distance(const CDescriptor * pd) const;

	CVectorSpaceDescriptor(const elementType * aDescriptor);

	CVectorSpaceDescriptor(const double * aDescriptor, double dScale);
	CVectorSpaceDescriptor(const float * aDescriptor, float dScale);

    ~CVectorSpaceDescriptor() { /*cout << "Deleting a descriptor\n";*/ };

    static int DescriptorLength() { return nDescriptorLength; };
    static int SURFDescriptorLength() { return nDescriptorLength; }; //so will not be overridden
    const elementType * DescriptorVector() const { return aDescriptorVector; };
    typedef elementType elType; //Exposes type to other classes
    static eInvDescriptorType vectorType() { return (((elementType)255) > 0) ? eSIFT : eSURF; }; //SIFT iff unsigned char

	virtual int size() const { return sizeof(this); }
	virtual int length() const { return nDescriptorLength; }
};

/*//Specialisation for char: Template Class Partial Specialization http://www.iis.sinica.edu.tw/~kathy/vcstl/templates.htm
#define templateCharVS template<int nDescriptorLength, eNorm eNormToUse>
#define CCharTVectorSpaceDescriptor CVectorSpaceDescriptor<char, nDescriptorLength, eNormToUse>
templateCharVS
class CCharTVectorSpaceDescriptor : public CCharTVectorSpaceDescriptor
{
    friend static inline double euclidDistance(const CCharTVectorSpaceDescriptor * d1, const CCharTVectorSpaceDescriptor * d2);
	static inline double euclidDistance(const CCharTVectorSpaceDescriptor * d1, const CCharTVectorSpaceDescriptor * d2);
    friend static inline double unitCosineDistance(const CCharTVectorSpaceDescriptor * d1, const CCharTVectorSpaceDescriptor * d2);
	static inline double unitCosineDistance(const CCharTVectorSpaceDescriptor * d1, const CCharTVectorSpaceDescriptor * d2);
};*/

templateVSSet class CVectorSpaceDescriptorSet : public CDescriptorSet
{
	std::vector<const descriptorType * > vDescriptors;
    bool KmeansIterate(CClusterSet * pCentres) const;
public:
    CClusterSet * Cluster(int nClusters) const;
	unsigned int Count() const { return (unsigned int)vDescriptors.size(); };

	CVectorSpaceDescriptorSet(int nDescriptorCountEstimate = 0) { vDescriptors.reserve(nDescriptorCountEstimate); };

    void Push(CDescriptor * pDescriptor);
    void Push(CDescriptorSet * pDescriptorSet);
    void Clear() { vDescriptors.clear(); };

	CDescriptor * operator[](int i) { return vDescriptors[i]; };
	const CDescriptor * operator[](int i) const { return vDescriptors[i]; };

    virtual ~CVectorSpaceDescriptorSet() { vDescriptors.clear(); /*vDescriptors.~vector();*/ };

    //So we can create a new instance of a polymorphic class with the correct type
    virtual CTVectorSpaceDescriptorSet * newDescriptorSet(int nDescriptorCountEstimate = 0) const { return new CTVectorSpaceDescriptorSet(nDescriptorCountEstimate); };
};

//Use CLARA-style k-means
templateVSSet class CLargeVSDescriptorSet : public CTVectorSpaceDescriptorSet
{
    inline CClusterSet * KMCluster(int nClusters) const { return CTVectorSpaceDescriptorSet::Cluster(nClusters); };
public:
    CClusterSet * Cluster(int nClusters) const;
    CLargeVSDescriptorSet(int nDescriptorCountEstimate = 0) : CTVectorSpaceDescriptorSet(nDescriptorCountEstimate) {};

    //So we can create a new instance of a polymorphic class with the correct type
    virtual CTVectorSpaceDescriptorSet * newDescriptorSet(int nDescriptorCountEstimate = 0) const { return new CTLargeVSDescriptorSet(nDescriptorCountEstimate); };
};

//////// Metric Spaces /////////////////////

//class that actually implements clustering. TODO: move to own header
class CMetricSpaceDescriptorSet : public CDescriptorSet, boost::noncopyable
{
protected:
	const CDescriptorSetClusteringParams & DSC_PARAMS;
	CDynArray<CDescriptor * > vDescriptors;

private:
	CClusterSet * kMedoidCluster(int nClusters, const CCluster * pParentCluster) const HOT;
	CClusterSet * kMeansCluster(int nClusters, const CCluster * pParentCluster) HOT;
	CClusterSet * RandomCluster(int nClusters, const CCluster * pParentCluster, bool bDuplicate);

	CClusterSet * cluster_int(unsigned int nClusters, const CCluster * pParentCluster, int nAssignSpeedup);

    void AssignMSDescriptorSets(CClusterSet * pClusters) const;

    typedef int TKMDistType;
    inline static int arrayMinIdx(TKMDistType const * aArray, const int nLength);
    inline static int arrayMinIdx(const boost::scoped_array<TKMDistType> & aArray, int nLength);

public:
	CMetricSpaceDescriptorSet(const CDescriptorSetClusteringParams & DSC_PARAMS, int nDescriptorCountEstimate) : DSC_PARAMS(DSC_PARAMS) { vDescriptors.reserve(nDescriptorCountEstimate); };

	virtual CClusterSet * Cluster(int nClusters, const CCluster * pParentCluster) ; //k-medoids (CLARA)

    void Push(CDescriptor * pDescriptor) ;
    void Push(CDescriptorSet * pDescriptorSet);
    void Clear() { vDescriptors.clear(); };

    int Count() const { return  vDescriptors.size(); }

	CDescriptor * operator[](int i) { return vDescriptors[i]; }
	const CDescriptor * operator[](int i) const { return vDescriptors[i]; }

    virtual ~CMetricSpaceDescriptorSet() { }

	virtual CMetricSpaceDescriptorSet * makeNewDS(int nEstimatedSize = -1) const { return new CMetricSpaceDescriptorSet(DSC_PARAMS, nEstimatedSize==-1 ? Count() : nEstimatedSize); }
};

// Descriptor set that clusters by choosing k random elements. Used to show clustering is worthwhile.
class CRandomDescriptorSet : public CMetricSpaceDescriptorSet
{
    void DownSize(int n);
public:
    CRandomDescriptorSet(const CDescriptorSetClusteringParams & DSC_PARAMS, int nDescriptorCountEstimate) : CMetricSpaceDescriptorSet(DSC_PARAMS, nDescriptorCountEstimate) {};

	virtual CClusterSet * Cluster(int nClusters, const CCluster * pParentCluster) ; //Override k-medoids

	virtual CRandomDescriptorSet * makeNewDS(int nEstimatedSize = -1) const { return new CRandomDescriptorSet(DSC_PARAMS, nEstimatedSize==-1 ? Count() : nEstimatedSize); }
};


/////////// Select the one that's currently in use ///////////

//typedef TOriginalDescriptor TDescriptor;
//typedef CMetricSpaceDescriptorSet TDescriptorSet;
//typedef CLargeVSDescriptorSet<double, 64, eEuclidSquared> TDescriptorSet;
//typedef CVectorSpaceDescriptorSet<char, 64, eCosine> TDescriptorSet;
//typedef CRandomDescriptorSet TDescriptorSet;
//typedef CMetricSpaceDescriptorSet TDescriptorSet;

/////////// Now for the implementation ///////////

//CLARA-style k-means
templateVSSet
CClusterSet * CTLargeVSDescriptorSet::Cluster(int nClusters) const //might need to cast a TDescriptorVector to a TVSDescriptorVector to make this work
{
    const unsigned int MAX_CLARAKMEANS_ITERS = 1;
    const unsigned int MAX_KMEANS_COUNT = 3000;
    const unsigned int KMEANS_SUBSET_COUNT = 2000;

    if(this->Count()<MAX_KMEANS_COUNT)
        return KMCluster(nClusters);

    std::cout << "Clustering: CLARA K-means, " << this->Count() << " descriptors, " << nClusters << " clusters\n";

    unsigned int nIters = this->Count()/(2*KMEANS_SUBSET_COUNT) + 1; //guess
    if(nIters > MAX_CLARAKMEANS_ITERS) nIters = MAX_CLARAKMEANS_ITERS;
    if(nIters < 1) nIters = 1;

	CClusterSet * pClusters=0;
	CClusterSet * pBestClusters=0;
    double dBestClustering = HUGE;

	for(unsigned int nIteration=0; nIteration<nIters; nIteration++)
	{
        std::cout << "CLARA K-means Iteration " << nIteration << '\n';
		//First cluster a random subset (k-medoids)
        CTLargeVSDescriptorSet * pRandSubset = new CTLargeVSDescriptorSet(KMEANS_SUBSET_COUNT);
        RandomSubset(KMEANS_SUBSET_COUNT, pRandSubset); //DO NOT add the best centres--THEY'RE NOT REAL POINTS
                                                        //TODO: Use them to init clustering next time

		pClusters = pRandSubset->KMCluster(nClusters);//, nIteration==CLARA_ITERS-1);

        if(pClusters->Assigned()) pClusters->ResetClusters(); // cos we've only assigned the subset

        double dTotalDist = pClusters->AssignToClusters<true>(this, 1);

        if(dTotalDist < dBestClustering) //if this is the best clustering so far...
        {
            if(pBestClusters)
			    delete pBestClusters;

            pBestClusters = pClusters;
            dBestClustering = dTotalDist;
        }
        else
    	    delete pClusters;

        delete pRandSubset;
	}

    return pBestClusters;
}

//Hacky way of creating new descriptors
CDescriptor * DupDescriptor(const CDescriptor * pCentre, size_t descriptorSize);

//K-means implementation
templateVSSet
CClusterSet * CTVectorSpaceDescriptorSet::Cluster(int nClusters) const //might need to cast a TDescriptorVector to a TVSDescriptorVector to make this work
{
	CHECK(vDescriptors.size()==0 || (unsigned int)nClusters > vDescriptors.size() || nClusters <= 0, "CTVectorSpaceDescriptorSet::Cluster: No descriptors or no images available yet to cluster, or more clusters than descriptors");

    CClusterSet * pClusters = new CClusterSet(nClusters, 0);

    if(eKMeans != eKMLL)
    {
        //k-means implementation
        const unsigned int KMEANS_ITERS = 15; //15 is plenty

        //static descriptorType::elementType * aTemp = new descriptorType::elementType[descriptorType::DescriptorLength()];
        //static descriptorType::elementType * aZero = (descriptorType::elementType *)memset(aTemp, 0, sizeof(descriptorType::elementType)*descriptorType::DescriptorLength());

        //Start with a random set of centres
        CTVectorSpaceDescriptorSet * pRandSubset = this->newDescriptorSet(nClusters);
        RandomSubset(nClusters, pRandSubset);

        // Initiate each cluster with a zero-descriptor and set containing one centre
        for(int nCentre = 0; nCentre < nClusters; nCentre++)
        {
            CDescriptorSet * pClusterMembers = this->newDescriptorSet(0);
            const CDescriptor * pCentre = (*pRandSubset)[nCentre];
            //pClusterMembers->Push(pCentre);

            //CDescriptor * pDescriptor = new descriptorType(aZero);
            CDescriptor * pDescriptor = DupDescriptor(pCentre, sizeof(descriptorType));
        	//CDescriptor * pDescriptor = new TDescriptor((char*)(void*)pCentre, 0,0);

            THROW("Uncomment+Need more args to c'tor below...");
            //pClusters->push_back(new CCluster(pDescriptor, true, &pClusterMembers));
        }
        delete pRandSubset;

        for(unsigned int nIters = KMEANS_ITERS; nIters > 0; nIters--)
        {
            //bool bLastTime = (nIters == KMEANS_ITERS-1);
            if(KmeansIterate(pClusters)) break;
        }

        CHECK(!pClusters, "Cluster set not initialised ?!!");
        return pClusters; //should be assigned
    }
#ifndef __GNUC__
#if 0
    //----------------------------------------------------------------------
    //  Global clustering parameters (some are set in getArgs())
    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    //  Termination conditions
    //	These are explained in the file KMterm.h and KMlocal.h.  Unless
    //	you are into fine tuning, don't worry about changing these.
    //----------------------------------------------------------------------
    // todo: move
    TS_STATIC KMterm TerminationSettings(20, 0, 0, 0,		// run for 100 stages (default)
		0.10,			// min consec RDL
		0.10,			// min accum RDL
		3,			    // max run stages
		0.50,			// init. prob. of acceptance
		10,			    // temp. run length
		0.95);			// temp. reduction factor

    //Apply k-means here (friend needed??)
    //Copy my points into vectors

    KMdata DataPts(descriptorType::DescriptorLength(), Count());		// allocate data storage

    for(unsigned int nDescriptor=0;nDescriptor<vDescriptors.size();nDescriptor++)
    {
        arrayCopy(vDescriptors[nDescriptor]->DescriptorVector(), DataPts[nDescriptor], descriptorType::DescriptorLength());
        //memcpy(DataPts[nDescriptor], vDescriptors[nDescriptor]->aDescriptorVector, descriptorType::DescriptorLength()*sizeof(elementType));
    }

    if(vDescriptors.size() != (unsigned int)DataPts.getNPts()) throw new CException("CTVectorSpaceDescriptorSet::Cluster: Inconsistant descriptor count");
    DataPts.setNPts((int)vDescriptors.size());			// set actual number of pts

    DataPts.buildKcTree(); // build filtering structure

    KMfilterCenters clusters(nClusters, DataPts);		// allocate centers (TODO: Use centres from last time?? How do we do this? Probably not practical below the top level)
    /*clusters.genRandom(); Uncomment to prevent purify access errors
    clusters.getSums();
    clusters.getDists();
    clusters.getSumSqs();
    clusters.getWeights(); Uncomment to prevent purify access errors*/

    cout << "Clustering: Lloyd's, " << DataPts.getNPts() << " descriptors, " << nClusters << " clusters\n";

    //KMlocalEZ_Hybrid kmEZ_Hybrid(clusters, *pTerminationSettings);	// EZ-Hybrid heuristic
    //clusters = kmEZ_Hybrid.execute();
    						// TODO: Try different algorithms
    KMlocalLloyds kmLloyds(clusters, TerminationSettings);	// repeated Lloyd's //todo: get rid of static
    clusters = kmLloyds.execute();

    //KMcenterArray clusterCenterArray = clusters.getAssignments(); //todo: is this useful?

    //Now add clusters to a cluster object to return
    KMcenterArray clusterCenterArray = clusters.getCtrPts();
//    int * clusterCountArray = clusters.getWeights(true); //todo: try false--not recalculating distortion todo: does this leak?

    for(int nWord=0; nWord < nClusters; nWord++)
    {
        CDescriptorSet * pClusterMembers = this->newDescriptorSet(0);
        //pClusters->push_back(new CCluster(new descriptorType(clusterCenterArray[nWord]), true, &pClusterMembers));
        //pClusters->push_back(new CCluster(new descriptorType(clusterCenterArray[nWord], 1), true, &pClusterMembers));
        throw "Don't use this KM";
        //pClusters->push_back(new CCluster(DupDescriptor(clusterCenterArray[nWord], sizeof
    }

    //Now assign all descriptors to a cluster
    /*double dDist =*/ pClusters->AssignToClusters<true>(this, 1);

//    cout << dDist/Count() << "=av dist\n" ;
    /*/Output lengths of cluster centres:
    double d[64];
    memset(d,0,sizeof(double)*64);
    TOriginalDescriptor dZero(d);
    for(int i=0; i < nClusters; i++)
    {
        cout << (*pClusters)[i]->Centre()->distance(&dZero) << " ";
    }
    cout << "\n" ;*/

	return pClusters;
#endif
#endif
}

templateVSSet
void CTVectorSpaceDescriptorSet::Push(CDescriptor * pDescriptor)
{
    if(IS_DEBUG) CHECK(!pDescriptor, "CTVectorSpaceDescriptorSet::Push: Bad parameters");
    const descriptorType * pVSDescriptor = CAST<const descriptorType *>(pDescriptor);
    if(IS_DEBUG) CHECK(!pVSDescriptor, "CTVectorSpaceDescriptorSet::Push: Cast failed");
	vDescriptors.push_back(pVSDescriptor);
}

templateVSSet
void CTVectorSpaceDescriptorSet::Push(CDescriptorSet * pDescriptorSet)
{
    if(IS_DEBUG) CHECK(!pDescriptorSet, "CTVectorSpaceDescriptorSet::Push: Bad parameters");
    const CTVectorSpaceDescriptorSet * pVSDescriptorSet = CAST<const CTVectorSpaceDescriptorSet *>(pDescriptorSet);
    if(IS_DEBUG) CHECK(!pVSDescriptorSet, "CTVectorSpaceDescriptorSet::Push: Cast failed");
	vDescriptors.reserve(vDescriptors.size()+pDescriptorSet->Count());
	vDescriptors.insert(vDescriptors.end(), pVSDescriptorSet->vDescriptors.begin(), pVSDescriptorSet->vDescriptors.end());
}

templateVS
CTVectorSpaceDescriptor::CVectorSpaceDescriptor(const elementType * aDescriptor)
{
	memcpy(aDescriptorVector, aDescriptor, nDescriptorLength*sizeof(elementType));
}

templateVS
CTVectorSpaceDescriptor::CVectorSpaceDescriptor(const double * aDescriptor, double dScale) //todo: if elementType==double ignore
{
    arrayCopyScale<double, elementType>(aDescriptor, aDescriptorVector, dScale, nDescriptorLength);
}

templateVS
CTVectorSpaceDescriptor::CVectorSpaceDescriptor(const float * aDescriptor, float dScale) //todo: if elementType==float ignore
{
    arrayCopyScale(aDescriptor, aDescriptorVector, dScale, nDescriptorLength);
}

//This is in [0,2] and is equivalent to euclidian dist squared for unit vectors
templateVS
inline CDescriptor::TDist CTVectorSpaceDescriptor::unitCosineDistance(const CTVectorSpaceDescriptor * d1, const CTVectorSpaceDescriptor * d2)
{
    elementType d=1;
    for(int i=0;i<nDescriptorLength;i++)
        d -= d1->aDescriptorVector[i] * d2->aDescriptorVector[i]; //todo: pointer arithmatic

    return (CDescriptor::TDist)d;
}

templateVS
CDescriptor::TDist CTVectorSpaceDescriptor::distance(const CDescriptor * d) const
{
    const CTVectorSpaceDescriptor * pd = CAST<const CTVectorSpaceDescriptor *>(d);

    switch(eNormToUse)
    {
        case eEuclidSquared:
            return (CDescriptor::TDist)euclidDist<elementType>(aDescriptorVector, pd->aDescriptorVector, nDescriptorLength);
        case eCosine:
            return (CDescriptor::TDist)intLookup::Sqrt(euclidDist<elementType>(aDescriptorVector, pd->aDescriptorVector, nDescriptorLength));
//            return (CDescriptor::TDist)cosDist(aDescriptorVector, pd->aDescriptorVector, nDescriptorLength);
        case e1Norm:
            return (CDescriptor::TDist)L1Dist<elementType>(aDescriptorVector, pd->aDescriptorVector, nDescriptorLength);
        case eMaxNorm:
            return (CDescriptor::TDist)MaxDist<elementType>(aDescriptorVector, pd->aDescriptorVector, nDescriptorLength);
        /*case eEuclidParallel:
        	return
        anTwoTimesSumSquare*/

        //case eEuclidSquared:
        //    return euclidDistance(this, CAST<const CTVectorSpaceDescriptor *>(d));
        //case eCosine:
        //    return unitCosineDistance(this, CAST<const CTVectorSpaceDescriptor *>(d));
        default:
        	THROW("Norm not handled");
    }
}

templateVSSet
bool CTVectorSpaceDescriptorSet::KmeansIterate(CClusterSet * pCentres) const
{
    //one iteration of k-means
    //assign, then move centres. Return centres without any assignment

    bool bAssignedAlready = pCentres->Assigned();
    double dOldScore = pCentres->TotalDistance();
    if(bAssignedAlready) pCentres->ResetClusters(); //cos we'll reassign

    pCentres->AssignToClusters<true>(this);

    if(bAssignedAlready && pCentres->TotalDistance() >= dOldScore)//check for improving score. Stop if doesn't improve
        return true; //finished converging

    for(CClusterSet::iterator ppCentre = pCentres->begin(); ppCentre != pCentres->end(); ppCentre++)
    {
    	CCluster * pCentre = *ppCentre;
        pCentre->RecalcCentre(eKMeans);
    }

    return false; //haven't finished converging
}


templateVS
class CVectorLocationDescriptor : public CTVectorSpaceDescriptor
{
    CLocation Location;
public:
    CVectorLocationDescriptor(const double * aDescriptor, double dScale, CLocation loc) : CTVectorSpaceDescriptor(aDescriptor, dScale), Location(loc) {}
    CVectorLocationDescriptor(const float * aDescriptor, float dScale, CLocation loc) : CTVectorSpaceDescriptor(aDescriptor, dScale), Location(loc) {}
    virtual CLocation location() const { return Location; }

	virtual int size() const { return sizeof(this); }
};

templateVS
class CVectorLocationOrientationDescriptor : public CTVectorLocationDescriptor
{
    const double orientationAngle;
public:
    CVectorLocationOrientationDescriptor(const double * aDescriptor, double dScale, double orientationAngle, CLocation loc) : CTVectorLocationDescriptor(aDescriptor, dScale, loc), orientationAngle(orientationAngle) {};
    CVectorLocationOrientationDescriptor(const float * aDescriptor, float dScale, double orientationAngle, CLocation loc) : CTVectorLocationDescriptor(aDescriptor, dScale, loc), orientationAngle(orientationAngle) {}
	virtual double orientation() const { return orientationAngle; }

	virtual int size() const { return sizeof(this); }
};
