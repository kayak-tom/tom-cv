/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once
//Templated stuff relating to clustering

#include "description/descriptor.h"
#include <iostream>
#include <boost/noncopyable.hpp>

//One cluster
class CCluster : boost::noncopyable
{
public:
	class CBigClustersFirstPred;
private:
	friend class CCluster::CBigClustersFirstPred;

	CDescriptor * pClusterCentre;
	CDescriptorSet * pvMembers; // this can be clustered!!
	bool bIOwnCentreMem;
	CDescriptor::TDist nHalfMinDist; //Half the distance to the nearest cluster centre--descriptors within this rad are definitely in this cluster
	CDescriptor::TDist nRadiusSD; //SD of dispersion--can be used to determine radius dynamically. Make sure one centre doesn't make one catch-all cluster with massive variance...
	int nDescriptorCountKMed;
	int nMyIdx;
	const CCluster * pParentCluster;
	//const CCluster & operator=(const CCluster & otherCluster) { THROW("DO NOT USE") }
public:
	CCluster(CDescriptor * pCentre_in, bool bIOwnCentreMem_in, CDescriptorSet ** ppDescriptorSet, CDescriptor::TDist nHalfMinDist, CDescriptor::TDist nRadiusSD, int nDescriptorCountKMed, const int nMyIdx, const CCluster * pParentCluster);
	//CCluster & operator=(CCluster & otherCluster);
	//CCluster(CCluster & otherCluster);
	//CCluster(const CCluster & otherCluster) { THROW("DO NOT USE") }

	~CCluster();

	const CCluster * parentCluster() const { return pParentCluster; }

	int idx() const { return nMyIdx; }

	CDescriptor::TDist halfMinDist() const { return nHalfMinDist; };
	int KMCountEst() const { return nDescriptorCountKMed; }

	void AssignDescriptorSet(CDescriptorSet ** ppDescriptorSet);

	const CDescriptor * Centre() const
	{
		if(IS_DEBUG) CHECK(!pClusterCentre, "Missing cluster centre")
		return pClusterCentre;
	}
	CDescriptor * Centre()
	{
		if(IS_DEBUG) CHECK(!pClusterCentre, "Missing cluster centre")
		return pClusterCentre;
	}

	inline int Count() const
	{
		return pvMembers->Count();
	}

	bool IWantCentreMemory()
	{
		bool bTemp = bIOwnCentreMem;
		bIOwnCentreMem = false;
		//pClusterCentre = 0;
		return bTemp;
	}

	inline CDescriptorSet * Members()
	{
		return pvMembers;
	}
	inline const CDescriptorSet * Members() const
	{
		return pvMembers;
	}

	void Push(CDescriptor * pDescriptor)
	{
		if(IS_DEBUG) CHECK(pDescriptor->assignment() != this, "Warning: Assignment has failed")
		pvMembers->Push(pDescriptor);
	}

	enum eWeightMethod
	{
		eNone, eRadius, eDistInv
	};

	void RecalcCentre(eKMImplementation weightMethod); //only relevent for vector descriptors
	void ResetClusters();

	class CBigClustersFirstPred
	{
	public:
		bool operator()(const CCluster * c1, const CCluster * c2) const
		{
			return c1->nDescriptorCountKMed > c2->nDescriptorCountKMed;
		}
	};
};

//A clustering
class CClusterSet: public TClusterVector
{
	static inline const double NOT_ASSIGNED()
	{
		return -1;
	}
	const CDescriptor::TDist * aaDistances;

	double dTotalDistance;
public:
	CClusterSet(unsigned int nSize, CDescriptor::TDist ** paaDistances)
	: aaDistances(paaDistances ? *paaDistances : 0), dTotalDistance ( NOT_ASSIGNED())
	{
		reserve(nSize);
		if(paaDistances)
		{
			*paaDistances = 0;
		}
	}


	~CClusterSet()
	{
		delete [] aaDistances;
	}

	template<bool bQuiet>
	double AssignToClusters(CDescriptorSet * pDescriptors, int nAssignSpeedup) HOT; //Assign all in this set to clusters. NB sometimes the ret val isn't used, and doesn't need to be calculated
	void CheckDisimilarity() const;
	void CheckAssignment() const;

	void RecalcCentres(eKMImplementation weightMethod);
	void ResetClusters();

	inline double TotalDistance() const
	{
		return dTotalDistance;
	}

	inline bool Assigned() const
	{
		return dTotalDistance != NOT_ASSIGNED();
	}

	void remindDescriptorsOfAssignments();
};

// If nAssignSpeedup > 1 we don't actually assign descriptors to clusters, only estimate disp.
template <bool bQuiet>
double CClusterSet::AssignToClusters(CDescriptorSet * pDescriptors, int nAssignSpeedup)
{
    CHECK(dTotalDistance != CClusterSet::NOT_ASSIGNED(), "CClusterSet::AssignToClusters: Already assigned");
    double dTotalDist = 0;
    int nBorderlinePoints = 0;
    int nDescriptors = pDescriptors->Count();

    ARRAY(bool, aTooFarAway, size());
	static int s_nComparisons=0, s_nFarAway=0, s_nAvoidedByBeingClose=0;

    for(int nDescriptor=0; nDescriptor < nDescriptors; nDescriptor += nAssignSpeedup)
    {
        CDescriptor * pDescriptor = (*pDescriptors)[nDescriptor];

        CClusterSet::iterator ppClosestCluster=begin(); //Will be used when there's only 1 option
        CCluster * pClosestCluster = *ppClosestCluster;

        if(pDescriptor->assignment())
        {
        	dTotalDist += pDescriptor->assignmentDist();
        	pClosestCluster = const_cast<CCluster *>( pDescriptor->assignment() ); //Think this is ok...
        }
        else
        {
			CDescriptor::TDist dClosestDist=MAX_ALLOWED_DIST;
			if(size() > 1)
			{
				CDescriptor::TDist dNextClosestDist=MAX_ALLOWED_DIST;

				setConstant(PTR(aTooFarAway), false, size());

				//int nNextIdx = 0;

				int nActualComparisons=0, nFarAway=0;
				for(CClusterSet::iterator ppCluster=begin(); ppCluster != end(); ppCluster++)
				//while(nNextIdx >= 0)
				{
			        CCluster * pCluster = *ppCluster;
					//CCluster * pCluster = pBegin + nNextIdx;

					int nThisCentreIdx = pCluster->idx();

					//cout << "Trying c " << nNextIdx << " idx " << pCluster->idx() << endl;

					if(aTooFarAway[nThisCentreIdx])
					{
						DEBUGONLY(CDescriptor::TDist dDist=pDescriptor->distance(pCluster->Centre());
						if(dDist < dClosestDist)
							cout << "ERROR: Closer centre missed\n");
						nFarAway++;
						continue;
					}

					//std::cout << pCluster->KMCountEst() << ' ';
					CDescriptor::TDist dDist=pDescriptor->distance(pCluster->Centre());
					nActualComparisons++;
					if (dDist<dClosestDist)
					{
						if(!bQuiet)
							dNextClosestDist = dClosestDist;

						dClosestDist=dDist;
						pClosestCluster=pCluster;

						/*if(dDist < pCluster->halfMinDist())
						{
							break;
						}*/
					}

					if(aaDistances)
					{
						//CDescriptor::TDist distFromCentreUB = 2*dDist; //For any metric
						//CDescriptor::TDist distFromCentreUB = 4*dDist; //For Euclidean norm
						//CDescriptor::TDist distFromCentreUB = dDist+dClosestDist; //For any metric
						CDescriptor::TDist distFromCentreUB = 2*(dDist+dClosestDist); //UB1 For Euclidean norm
						//CDescriptor::TDist distFromCentreUBexact = sqr(sqrt(dDist)+sqrt(dClosestDist)); //UB exact
						//CDescriptor::TDist distFromCentreLBexact = sqr(sqrt(dDist)-sqrt(dClosestDist)); //LB exact
						//cout << distFromCentreUB << " approx, " << distFromCentreUBexact << " exact\n";
						//cout << distFromCentreUB << " approx, " << distFromCentreLBexact << " LB exact\n";

						aTooFarAway[nThisCentreIdx] = true;

						/*CDescriptor::TDist dClosestToThis = MAX_ALLOWED_DIST;
						nNextIdx = -1;

						 * Choosing furthest descriptor saves about 10%
						 * 			nearest 					    16%
						 * 			nearest to distance			    16%
						*/

						bool * pTooFarAway = PTR(aTooFarAway);
						for(int nOtherClusterIdx=0; nOtherClusterIdx < size(); nOtherClusterIdx++, pTooFarAway++)
						{
							if(*pTooFarAway)
								continue;

							CDescriptor::TDist dDistToOtherCentre = aaDistances[nThisCentreIdx * size() + nOtherClusterIdx];
							if(distFromCentreUB <= dDistToOtherCentre)
							{
								//then definitely not closer to other centre, by Triangle ineq.
								*pTooFarAway = true;
							}
							/*else if(!aTooFarAway[nOtherClusterIdx] && distFromCentreUBexact <= dDistToOtherCentre)
								cout << "1 in gap\n"; Extremely few in gap
							else if(distFromCentreLBexact >= dDistToOtherCentre) Extremely few can be eliminated this way--1 in 250 descriptors eliminate 1 centre
								cout << "1 within LB\n";*/


							/*else if(fabs(dDistToOtherCentre - dDist) < dClosestToThis)
							{
								nNextIdx = nOtherClusterIdx;
								dClosestToThis = fabs(dDistToOtherCentre - dDist);
							}*/
						}
					}
				}
				//std::cout << "Assigned to " << size() << " centres with " << i << " comparisons, " << nFarAway << " avoided by being too far away\n";
				s_nComparisons += size();
				s_nFarAway += nFarAway;
				s_nAvoidedByBeingClose += size()-nActualComparisons;
				dTotalDist += dClosestDist;
			}
			//else
				//cout << "1 cluster\n";

			//if(nAssignSpeedup==1)
	        pDescriptor->assignToCluster(pClosestCluster, dClosestDist);
        }

        if(nAssignSpeedup==1)
        {
            pClosestCluster->Push(pDescriptor);
        }


        /*if(!bQuiet)
        {
            if(dNextClosestDist-dClosestDist < 5) //count points on the border of 2 clusters
                nBorderlinePoints++;
        }*/
    }
	//std::cout << "RUNNING TOTAL: Assigned to " << s_nComparisons << " centres, saving " << s_nAvoidedByBeingClose << " comparisons by knowing radius; " << s_nFarAway << " avoided by being too far away\n";
//Example RUNNING TOTAL: Assigned to 41799618 centres, saving 11301745 comparisons by knowing radius; 8120828 avoided by being too far away

    if(nAssignSpeedup==1)
        dTotalDistance = dTotalDist; //otherwise not assigned

    if(!bQuiet)
    {
        std::cout << "CLUSTERING: " << (int)size() << " clusters, " << pDescriptors->Count() << " descriptors\n";
        std::cout << nBorderlinePoints << " are borderline, proportion = " << nBorderlinePoints/(double)pDescriptors->Count() << "\n";
        std::cout << dTotalDist << "=total dispersion, average = " << dTotalDist/(double)pDescriptors->Count() << "\n";
        //return nBorderlinePoints;
    }

    return dTotalDist;
}
