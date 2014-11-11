/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#include "descriptorClustering.h"
#include "vectorDescriptor.h"
#include <set>
#include <functional>

using namespace std;

///////////  Metric space clustering /////////////////////////
void CMetricSpaceDescriptorSet::Push(CDescriptor * pDescriptor)
{
	//if(IS_DEBUG) CHECK(vDescriptors.contains(pDescriptor), "Pushing same descriptor to DS twice"); // remove slow check
	vDescriptors.push_back(pDescriptor);
}

void CMetricSpaceDescriptorSet::Push(CDescriptorSet * pDescriptorSet)
{
    const CMetricSpaceDescriptorSet * pMSDescriptorSet = CAST<const CMetricSpaceDescriptorSet * >(pDescriptorSet);
	vDescriptors.copy_back(pMSDescriptorSet->vDescriptors.begin(), pMSDescriptorSet->vDescriptors.end());
}

int pointerSortPredicate(void * p1, void * p2) { return p1>p2; }

//CLARA algorithm, Kaufman and Rousseeuw p144
CClusterSet * CMetricSpaceDescriptorSet::Cluster(int nClusters, const CCluster * pParentCluster)
{
	const CDescriptorSetClusteringParams::CCLARA_KMedoidsParams & CKM_PARAMS = DSC_PARAMS.CLARA_KMedoids;

	nClusters = min<int>(nClusters, CKM_PARAMS.MAX_DESCRIPTORS_KM); //Check not clustering into too many clusters

	//If problem is small use k-medoids
	if((int)Count()<=CKM_PARAMS.MAX_DESCRIPTORS_KM)
		return cluster_int(nClusters, pParentCluster, 1);

    //**cout << "Clustering: CLARA, " << Count() << " descriptors, " << nClusters << " clusters...";

	//Else: CLARA
	CClusterSet * pClusters=0;
	CClusterSet * pBestClusters=0;
    double dBestClustering = HUGE;
    //int nIters = CLARA_ITERS;

    int nIters = 1+Count()/(2*CKM_PARAMS.MAX_DESCRIPTORS_KM);
    if(nIters>CKM_PARAMS.CLARA_ITERS) nIters = CKM_PARAMS.CLARA_ITERS;

    //Use smaller subsets sometimes
    int nSubSetSize=nClusters+10+(Count()/4);
    if (nSubSetSize > CKM_PARAMS.MAX_DESCRIPTORS_KM) nSubSetSize = CKM_PARAMS.MAX_DESCRIPTORS_KM;
    //int nSubSetSize = MAX_DESCRIPTORS_KM;

    int nAssignSpeedup = (CKM_PARAMS.CLARA_ITERS>=2 && (int)Count()>8*CKM_PARAMS.MAX_DESCRIPTORS_KM) ? CKM_PARAMS.SPARSE_ASSIGN_COUNT : 1;

    bool bAssignedToBestClusters = false;

	for(int nIteration=0; nIteration<nIters; nIteration++)
	{
		//First cluster a random subset (k-medoids)

        CMetricSpaceDescriptorSet * pRandSubset = makeNewDS();
        RandomSubset(nSubSetSize, pRandSubset, pBestClusters);
        //**cout << "Random subset selected, size=" << pRandSubset->Count() << "\n";

		clearClusterAssignments();


		pClusters = pRandSubset->cluster_int(nClusters, pParentCluster, nAssignSpeedup);

        //pClusters->ResetClusters();
        //double dTotalDist = pClusters->AssignToClusters<true>(this, nAssignSpeedup);
        //pClusters->CheckAssignment();
        //**cout << "Assigned to cluster, dist=" << dTotalDist << "\n";
		double dTotalDist = pClusters->TotalDistance();

        if(dTotalDist < dBestClustering) //if this is the best clustering so far...
        {
            ////**cout<< "Best clustering so far: saving\n";
            if(pBestClusters)
			    delete pBestClusters;

            pBestClusters = pClusters;
            dBestClustering = dTotalDist;

            bAssignedToBestClusters = true;
        }
        else
        {
    	    delete pClusters;
    	    bAssignedToBestClusters = false;
        }
        delete pRandSubset;
        //**cout << "CLARA iteration complete\n";
	}

	if(!bAssignedToBestClusters)
	{
		clearClusterAssignments();
		if(nAssignSpeedup == 1)
			pBestClusters->remindDescriptorsOfAssignments();
	}

    if(nAssignSpeedup != 1)
    {
    	AssignMSDescriptorSets(pBestClusters);
        pBestClusters->AssignToClusters<true>(this, 1);
    }
    //else if(!bAssignedToBestClusters)
   		//pBestClusters->AssignToClusters<true>(this, 1);

    return pBestClusters;
}

CClusterSet * CMetricSpaceDescriptorSet::cluster_int(unsigned int nClusters, const CCluster * pParentCluster, int nAssignSpeedup)
{
	CClusterSet * pClusters=0;

	if(DSC_PARAMS.ClusteringAlg == CDescriptorSetClusteringParams::eCLARA_KMedoids)
		pClusters=kMedoidCluster(nClusters, pParentCluster);
	else if(DSC_PARAMS.ClusteringAlg == CDescriptorSetClusteringParams::eCLARA_KMeans || DSC_PARAMS.ClusteringAlg == CDescriptorSetClusteringParams::eKMeans)
		pClusters=kMeansCluster(nClusters, pParentCluster);
	//else if(DSC_PARAMS.ClusteringAlg == CDescriptorSetClusteringParams::eCLARA_RandomCentres || DSC_PARAMS.ClusteringAlg == CDescriptorSetClusteringParams::eRandomCentres)
		//pClusters=RandomCluster(nClusters, pParentCluster);

    if(nAssignSpeedup == 1)
    {
        AssignMSDescriptorSets(pClusters);
    }

    pClusters->AssignToClusters<true>(this, nAssignSpeedup);
    return pClusters;
}

//TODO: pointer arith.
inline int CMetricSpaceDescriptorSet::arrayMinIdx(TKMDistType const * aArray, const int nLength)
{
    TKMDistType const * pEl = aArray+1;
    TKMDistType dMin = *aArray;
    int nMin=0;
    for(int i=nLength-1; i>0; i--,pEl++)
    {
        if(*pEl<dMin)
        {
            dMin=*pEl;
            nMin=nLength-i;
        }
    }
    return nMin;
}
inline int CMetricSpaceDescriptorSet::arrayMinIdx(const boost::scoped_array<TKMDistType> & aArray, int nLength)
{
    TKMDistType dMin=aArray[0];
    int nMin=0;
    for(int i=1;i<nLength;i++)
    {
        if(aArray[i]<dMin)
        {
            dMin=aArray[i];
            nMin=i;
        }
    }
    return nMin;
}

CClusterSet * CMetricSpaceDescriptorSet::kMeansCluster(int nClusters, const CCluster * pParentCluster)
{
	THROW("K-means clustering not implemented at the moment. Coming soon...");
	CClusterSet * pClusters = RandomCluster(nClusters, pParentCluster, true);
	return pClusters;
}

CClusterSet * CMetricSpaceDescriptorSet::RandomCluster(int nClusters, const CCluster * pParentCluster, bool bDuplicate)
{
	const CDescriptorSetClusteringParams::CKMeansParams & CKM_PARAMS = DSC_PARAMS.KMeans;

	int nCentreIdx = 0;
	CMetricSpaceDescriptorSet randSubset(DSC_PARAMS, nClusters); //randSubset.reserve(nClusters);
	RandomSubset(nClusters, &randSubset);
	CClusterSet * pClusters = 0;

    for(int nNewCluster=0; nNewCluster<nClusters; nNewCluster++)
    {
    	CDescriptor * pCentre = 0;
    	if(bDuplicate)
    		pCentre = randSubset[nNewCluster]->clone();
    	else
    		pCentre = randSubset[nNewCluster];

    	pClusters->push_back(new CCluster(pCentre, bDuplicate, 0, MAX_ALLOWED_DIST, MAX_ALLOWED_DIST, 0, nNewCluster, pParentCluster));
    }

    return pClusters;
}

//PAM k-medoids algorithm, Kaufman and Rousseeuw p102
//Also makes use of the distances computed to assign clusters
CClusterSet * CMetricSpaceDescriptorSet::kMedoidCluster(int nClusters, const CCluster * pParentCluster) const
{
	const CDescriptorSetClusteringParams::CCLARA_KMedoidsParams & CKM_PARAMS = DSC_PARAMS.CLARA_KMedoids;
	const int MAX_DESCRIPTORS_KM = CKM_PARAMS.MAX_DESCRIPTORS_KM;

    int nDescriptors = (unsigned int)vDescriptors.size();
    if(nDescriptors>MAX_DESCRIPTORS_KM || nClusters>nDescriptors || nClusters<=0)
    CHECK(nDescriptors>MAX_DESCRIPTORS_KM || nClusters>nDescriptors || nClusters<=0, "CMetricSpaceDescriptorSet::kMedoidCluster: Too many clusters for k-medoids");

    //**cout << "Clustering: K-medoids, " << Count() << " descriptors, " << nClusters << " clusters...";

    ARRAY(TKMDistType, cost, MAX_DESCRIPTORS_KM); //Cost of adding centre i (i.e. negative is good)

    //Pre-calculate all distances
#ifdef __GNUC__
    TKMDistType d[MAX_DESCRIPTORS_KM][MAX_DESCRIPTORS_KM];
#else
    CHECK(MAX_DESCRIPTORS_KM>450, "K-Medoids: Static alloc size too small");
    TKMDistType d[450][450]; //Must be constant in Win
#endif

    ARRAY(int, aAssignments, MAX_DESCRIPTORS_KM); // i is assigned to centre aAssignments[i]
    ARRAY(int, aNewCentres, MAX_DESCRIPTORS_KM);
/*
    boost::multi_array<TKMDistType, 1> cost(boost::extents[MAX_DESCRIPTORS_KM]);
    boost::multi_array<TKMDistType, 2> d(boost::extents[MAX_DESCRIPTORS_KM][MAX_DESCRIPTORS_KM]);
    boost::multi_array<unsigned int, 1> aAssignments(boost::extents[MAX_DESCRIPTORS_KM]);
*/
    /*TKMDistType**d=0;//(TKMDistType**)(void*)d_stack;
    if(nDescriptors>200)
    {
        d=new TKMDistType*[nDescriptors];
        for (unsigned int i=0;i<nDescriptors;i++)
        {
            d[i]=new TKMDistType[nDescriptors];
        }
    }*/

    for (int nDesc1=0;nDesc1<nDescriptors;nDesc1++)
    {
        for (int nDesc2=nDesc1+1;nDesc2<nDescriptors;nDesc2++)
        {
            TKMDistType dist = (TKMDistType)(vDescriptors[nDesc1]->distance(vDescriptors[nDesc2]) );
            if(dist<0)
            	dist = (TKMDistType)(vDescriptors[nDesc1]->distance(vDescriptors[nDesc2]) );
            if(IS_DEBUG) CHECK(dist < 0, "CMetricSpaceDescriptorSet::kMedoidCluster: negative distance");
            dist=max<TKMDistType>(dist, 1); //ensure distances are positive
            d[nDesc2][nDesc1] = d[nDesc1][nDesc2] = dist;
        }
        d[nDesc1][nDesc1] = 0;

    //}
    //for (unsigned int nDesc1=0;nDesc1<nDescriptors;nDesc1++)
    //{
        //Pre-calculate sums of distances as well
        TKMDistType dTotalDist=0;
        for (int nDesc2=0;nDesc2<nDescriptors;nDesc2++)
        {
            dTotalDist += d[nDesc1][nDesc2];
        }
        cost[nDesc1] = dTotalDist;
    }


    //Choose k centres BUILD

    //First add the most central cluster center (given previous centres)
    int idxMin = arrayMinIdx(cost, nDescriptors);

    //TODO memset(aAssignments, idxMin, sizeof(unsigned int)*nDescriptors); //Assign all to i
    for(int i=0;i<nDescriptors; i++)
        aAssignments[i] = idxMin;

    //Now the distance to descriptor i's cluster centre is d[i, aAssignments[i]]

    for(int nToAdd=1; nToAdd<nClusters; nToAdd++)
    {
        for(int nNewPotentialCentre=0; nNewPotentialCentre<nDescriptors; nNewPotentialCentre++) // for each descriptor that's not a centre
        {
            cost[nNewPotentialCentre] = 0;
            if(aAssignments[nNewPotentialCentre] != nNewPotentialCentre) //not a centre
            {
                //Work out the decreased sum dist of making this descriptor a centre
                for(int nDescriptor=0; nDescriptor<nDescriptors; nDescriptor++) // for each descriptor that's not a centre
                {
                    TKMDistType dToCurrentCentre=d[nDescriptor][ aAssignments[nDescriptor]];
                    TKMDistType dToNewPotentialCentre=d[nDescriptor][ nNewPotentialCentre];
                    TKMDistType dCost = dToNewPotentialCentre-dToCurrentCentre;
                    if(dCost < 0)
                        cost[nNewPotentialCentre] += dCost;

                    //Don't calc list of new assignments here--will do that later
                }
            }
        }
        int idxMinCost = arrayMinIdx(cost, nDescriptors); //Find the new centre that minimises the cost
        if(aAssignments[idxMinCost] == idxMinCost)
        	idxMinCost = arrayMinIdx(cost, nDescriptors); //Find the new centre that minimises the cost
        CHECK(aAssignments[idxMinCost] == idxMinCost, "CMetricSpaceDescriptorSet::kMedoidCluster: idxMinCost is already a cluster centre");//(this should never be an existing centre);

        //And add it as a centre. Reassign all other descriptors
        for(int nDescriptor=0; nDescriptor<nDescriptors; nDescriptor++) // for each descriptor that's not a centre
        {
            TKMDistType dToCurrentCentre=d[nDescriptor][ aAssignments[nDescriptor]];
            TKMDistType dToNewCentre=d[nDescriptor][ idxMinCost];
            if(dToNewCentre < dToCurrentCentre)
            {
                aAssignments[nDescriptor] = idxMinCost;
            }
        }

        CHECK(aAssignments[idxMinCost] != idxMinCost, "CMetricSpaceDescriptorSet::kMedoidCluster: idxMinCost hasn't been added as a cluster");
    }

    // check we have nClusters cluster centres
    DEBUGONLY(
        int nCentreCount = 0;
        for(int i=0; i<nDescriptors; i++) // for each centre i
        {
            if(aAssignments[i] == i) nCentreCount ++;
        }
        CHECK(nCentreCount != nClusters, "CMetricSpaceDescriptorSet::kMedoidCluster: Wrong number of centres");
    );

    //Try moving each centre in turn to a random non-centre element SWAP
    for(int iter=0; iter<CKM_PARAMS.KMEDOIDS_ITERS; iter++)
    {
        bool bSwapped=false; // have we swapped any in this iteration?

        for(int i=0; i<nDescriptors; i++) // for each centre i
        {
            if(aAssignments[i] == i)
            {
            	int * pNewCentresEnd = PTR(aNewCentres);
                for(int nNewCluster=0; nNewCluster<nDescriptors; nNewCluster++)
                {
                    if(aAssignments[nNewCluster] == nNewCluster && nNewCluster != i) //For each centre that's not i
                    {
                    	*pNewCentresEnd = nNewCluster;
                    	pNewCentresEnd++;
                    }
                }

            	for(int h=0; h<nDescriptors; h++) //for each non-centre h
                {
                    TKMDistType dCost=0; //Calculate the cost (benefit) of doing this swap

                    if(aAssignments[h] != h)
                    {
                        for(int j=0; j<nDescriptors; j++) //for each non-centre that's not h
                        {
                            if(aAssignments[j] != j && j != h)
                            {
                                TKMDistType dToNewCentre=d[j][h];
                                if(aAssignments[j] == i)//if j was in cluster i
                                {
                                    TKMDistType dToCurrentCentre=d[j][i];
                                    TKMDistType dCostImprovment=dToNewCentre-dToCurrentCentre;
                                    if(dCostImprovment < 0) //if j is closer to h than it was to i
                                    {
                                        dCost += dCostImprovment;
                                    }
                                    else
                                    {
                                        //assign j to another cluster (may be h)
/*                                        for(unsigned int nNewCluster=0; nNewCluster<nDescriptors; nNewCluster++)
                                        {
                                            if(aAssignments[nNewCluster] == nNewCluster && nNewCluster != i) //For each centre that's not i*/
                                    	for(int * pNewCentre = PTR(aNewCentres); pNewCentre < pNewCentresEnd; pNewCentre++)
                                    	{
											TKMDistType dToThisCentre=d[j][*pNewCentre];
											if(dToThisCentre < dToNewCentre)
											{
												dToNewCentre = dToThisCentre;
											}
										}
                                        TKMDistType dCostImprovment=dToNewCentre-dToCurrentCentre; //might be positive (bad)
                                        //Don't do anything with assignment
                                        dCost += dCostImprovment;
                                    }
                                }
                                else
                                {
                                    TKMDistType dToCurrentCentre=d[j][aAssignments[j]];
                                    TKMDistType dCostImprovment=dToNewCentre-dToCurrentCentre;

                                    if ( dCostImprovment < 0 )//if j closer to h than its current centre
                                    {
                                        dCost += dCostImprovment;

                                    } //otherwise no effect on cost
                                }
                            }//isn't centre
                        }//for every other descriptor
                    }//isn't centre
                    cost[h] = dCost; //will be 0 if h is a centre
                }//For each possible candidate to swap to

                int idxMinCost = arrayMinIdx(cost, nDescriptors); //Find the new centre that minimises the cost
                if(cost[idxMinCost]<0)
                {
                    if(IS_DEBUG) CHECK(aAssignments[idxMinCost] == idxMinCost, "CMetricSpaceDescriptorSet::kMedoidCluster: idxMinCost is already a cluster centre");//(this should never be an existing centre);
                    //This swap is beneficial so do it
                    bSwapped=true;
                    // add it as a centre. Reassign all other descriptors
                    for(int nDescriptor=0; nDescriptor<nDescriptors; nDescriptor++) // for each descriptor that's not a centre
                    {
                        if(aAssignments[nDescriptor] == i)
                        {
                            //assign nDescriptor to another cluster (may be h)
                            int nNewClosestCluster = idxMinCost;
                            TKMDistType dToNewCentre=d[nDescriptor][idxMinCost];
                            for(int nNewCluster=0; nNewCluster<nDescriptors; nNewCluster++) //for each non-centre h
                            {
                                if(aAssignments[nNewCluster] == nNewCluster && nNewCluster != i) //For each new centre
                                {
                                    TKMDistType dToThisCentre=d[nDescriptor][nNewCluster];
                                    if(dToThisCentre < dToNewCentre)
                                    {
                                        dToNewCentre = dToThisCentre;
                                        nNewClosestCluster = nNewCluster;
                                    }
                                }
                            }
                            aAssignments[nDescriptor] = nNewClosestCluster;
                        }
                        else
                        {
                            TKMDistType dToCurrentCentre=d[nDescriptor][aAssignments[nDescriptor]];
                            TKMDistType dToNewCentre=d[nDescriptor][idxMinCost];
                            if(dToNewCentre < dToCurrentCentre)
                            {
                                aAssignments[nDescriptor] = idxMinCost;
                            }
                        }
                    }

                    if(IS_DEBUG) CHECK(aAssignments[idxMinCost] != idxMinCost, "CMetricSpaceDescriptorSet::kMedoidCluster: idxMinCost hasn't been added as a cluster");
                }//else no swap
            }//is a centre
        }//for each descriptor that...

        if (!bSwapped) break; //we didn't swap any centres this time
    }//iterate

    // check we have nClusters cluster centres
    DEBUGONLY(
        nCentreCount = 0;
        for(int i=0; i<nDescriptors; i++) // for each centre i
        {
            if(aAssignments[i] == i) nCentreCount ++;
        }
        CHECK(nCentreCount != nClusters, "CMetricSpaceDescriptorSet::kMedoidCluster: Wrong number of centres");
    );

    //Now estimate dispersion, dist to other centres, etc.
    ARRAYZ(int, aClusterCount, nDescriptors);
    ARRAY(TKMDistType, aMinDist, nDescriptors);
    setConstant(PTR(aMinDist), MAX_ALLOWED_DIST, nDescriptors);
    ARRAYZ(TKMDistType, aDispersionSD, nDescriptors);

    TKMDistType * aaDistances = new TKMDistType[nClusters*nClusters]; // matrix of distances between centres

    int nCentreIdx = 0;
    for(int nDesc=0; nDesc<nDescriptors; nDesc++)
    {
    	int nThisDescriptorsCluster = aAssignments[nDesc];
    	aClusterCount[nThisDescriptorsCluster]++;
    	TKMDistType nThisDistToCentre = d[nDesc][nThisDescriptorsCluster];
    	aDispersionSD[nThisDescriptorsCluster] += nThisDistToCentre;

    	if(nThisDescriptorsCluster == nDesc) //if this is a centre
    	{
    		aaDistances[nCentreIdx + nClusters*nCentreIdx] = 0;
    		int nNewCentreIdx = 0;

			for(int nNewCluster=0; nNewCluster<nDescriptors; nNewCluster++) //for each other centre
			{
				if(aAssignments[nNewCluster] == nNewCluster)//For each new centre
				{
					if(nNewCluster != nDesc) //not the same...
					{
						TKMDistType nDistToOtherCentre = d[nDesc][nNewCluster];
						if(IS_DEBUG) CHECK(nDistToOtherCentre<=0, "Bad distance");
						if(nDistToOtherCentre < aMinDist[nDesc])
							aMinDist[nDesc] = nDistToOtherCentre;
						if(IS_DEBUG) CHECK(aMinDist[nDesc]<=0, "Bad aMinDist[nDesc] distance");

						aaDistances[nNewCentreIdx + nClusters*nCentreIdx] = nDistToOtherCentre;
					}
					nNewCentreIdx++;
				}
			}
			nCentreIdx++;
    	}
    }

    CClusterSet * pClusters=new CClusterSet(nClusters, &aaDistances);
    ARRAY(int, aCentreNewIndices, nDescriptors);
    //ARRAY(CCluster *, aIndexToSortedCluster, nClusters);

    nCentreIdx = 0;
    for(int nNewCluster=0; nNewCluster<nDescriptors; nNewCluster++)
    {
        if(aAssignments[nNewCluster] == nNewCluster) //For each new centre
        {
        	if(IS_DEBUG) CHECK(aMinDist[nNewCluster] == 0, "Failed to get min dist to other centre"); //Might be big if singleton centre
        	if(aDispersionSD[nNewCluster] == 0)
        		aDispersionSD[nNewCluster] = aMinDist[nNewCluster]; //singleton centre
        	else
        		aDispersionSD[nNewCluster] /= aClusterCount[nNewCluster];

            pClusters->push_back(new CCluster(vDescriptors[nNewCluster], false, 0, aMinDist[nNewCluster]/2, aDispersionSD[nNewCluster], aClusterCount[nNewCluster], nCentreIdx, pParentCluster));
            aCentreNewIndices[nNewCluster] = nCentreIdx;
            nCentreIdx++;

            //Also make use of the assignments:

        }
    }

    for(int nDesc=0; nDesc<nDescriptors; nDesc++)
    {
    	int nThisDescriptorsCluster = aAssignments[nDesc];
    	TKMDistType dist = d[nThisDescriptorsCluster][nDesc];
    	CCluster * pClosestCluster = (*pClusters)[aCentreNewIndices[nThisDescriptorsCluster]];
    	vDescriptors[nDesc]->assignToCluster(pClosestCluster, dist);
    }

    std::sort(pClusters->begin(), pClusters->end(), CCluster::CBigClustersFirstPred() );

	//Also make use of the assignments:
	/*nCentreIdx = 0;
	for(CClusterSet::iterator ppCluster = pClusters->begin(); ppCluster != pClusters->end(); ppCluster++, nCentreIdx++)
	{
		CCluster * pCluster = *ppCluster;
		aIndexToSortedCluster[pCluster->idx()] = pCluster;
	}

    for(int nDesc=0; nDesc<nDescriptors; nDesc++)
    {
    	int nThisDescriptorsCluster = aAssignments[nDesc];
    	TKMDistType dist = d[nThisDescriptorsCluster][nDesc];
    	CCluster * pClosestCluster = aIndexToSortedCluster[aCentreNewIndices[nThisDescriptorsCluster]];
    	vDescriptors[nDesc]->assignToCluster(pClosestCluster, dist);
    }*/

    return pClusters;
}

void CMetricSpaceDescriptorSet::AssignMSDescriptorSets(CClusterSet * pClusters) const
{
    for(CClusterSet::iterator ppCluster=pClusters->begin(); ppCluster != pClusters->end(); ppCluster++)
    {
    	CCluster * pCluster = *ppCluster;
        CMetricSpaceDescriptorSet * pDS = makeNewDS(32);//Todo: size estimates
        CDescriptorSet * pDS2 = (CDescriptorSet *)pDS;
        pCluster->AssignDescriptorSet(&pDS2);
    }
}
