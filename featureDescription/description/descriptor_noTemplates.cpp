/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#include "descriptor.h"
#include "vectorDescriptor.h"
#include <set>
#include <functional>
#include "util/random.h"
#include "util/set2.h"
#include "util/Simple2dPoint.h"
#include <boost/thread.hpp>

using namespace std;
void CDescriptor::assignToCluster(const CCluster * pCentre, TDist closestClusterDist_in) {
    if(IS_DEBUG) CHECK(closestClusterDist_in != 0 && UNINIT(pCentre), "Bad centre set");
    if(IS_DEBUG) CHECK(pCentre != 0 && pCentre->parentCluster() && UNINIT(pCentre->parentCluster()), "Bad centre set");
    //if(IS_DEBUG) CHECK(pCentre && pCluster, "Already assigned; should have been cleared") This is fine if we're just k-medoid clustering
    pCluster = pCentre;
    closestClusterDist = closestClusterDist_in;
}

typedef std::pair<CDescriptor::TDist, int> TMatchPair;
class CMatchSort : binary_function<const TMatchPair &, const TMatchPair &, bool> {
public:
    result_type operator() (first_argument_type a,
            second_argument_type b) {
        return (result_type) (a.first < b.first || (a.first == b.first && a.second < b.second)); //Want to allow multiple equal dists
    }
};
class TMPSet : public set2<TMatchPair, CMatchSort/*, CIndividualPool_NoFree_Allocator< TMatchPair, NN_MAX + 1 >*/ > {
    //Make forward declaration easier
};
template<bool LIMIT_DISTANCES, bool ROTATION_INVARIANCE>
void CMatchableDescriptors::getBFC_matchDescriptors_RI(const int nStart, const int nEnd, const CMatchableDescriptors * pDS, const CMatchableDescriptors::CMatchSettings & MS, TMPSet * d) const {
    //const double dCondition = MS.MATCH_CONDITION;
    //const double dPP = MS.PP;
    const int NN = MS.NN;
    const int MAX_SEPERATION = MS.MAX_SEPERATION;
    const int MAX_SEPERATION_SQ = sqr(MAX_SEPERATION * SUBPIX_RES);
    //const int nCount1 = Count();
    const int nCount2 = pDS->Count();
    for (int i = nStart; i < nEnd; i++) {
        const CDescriptor * pDesc1 = get_const(i);
        CLocation loc1 = pDesc1->location();
        int nX1 = loc1.xAsIntFast();
        int nY1 = loc1.yAsIntFast();

        for (int j = 0; j < nCount2; j++) {
            const CDescriptor * pDesc2 = pDS->get_const(j);
            bool bCloseEnough = true;

            if (LIMIT_DISTANCES) {
                CLocation loc2 = pDesc2->location();
                int nX2 = loc2.xAsIntFast();
                int nY2 = loc2.yAsIntFast();
                bCloseEnough = sqr(nX1 - nX2) + sqr(nY1 - nY2) < MAX_SEPERATION_SQ;
            }

            if (bCloseEnough) {
                CDescriptor::TDist dist = ROTATION_INVARIANCE ? pDesc1->distance(pDesc2) : pDesc1->orientedDistance(pDesc2, MS.OI_NEARBY_ANGLE);
                /*if(!ROTATION_INVARIANCE)
                {
                    cout << pDesc1->orientation() << ' ' << pDesc2->orientation() << endl;
                }*/

                if ((int) d[i].size() < NN + 1) {
                    d[i].insert(TMatchPair(dist, j));
                } else {
                    TMPSet::iterator pBack = d[i].end();
                    pBack--;
                    CDescriptor::TDist lowestDist = pBack->first;
                    if (dist < lowestDist) {
                        d[i].erase(pBack);
                        d[i].insert(TMatchPair(dist, j));
                    }
                }
            }
        }
    }
}
template<bool LIMIT_DISTANCES>
void CMatchableDescriptors::getBFC_matchDescriptors(const int nStart, const int nEnd, const CMatchableDescriptors * pDS, const CMatchableDescriptors::CMatchSettings & MS, TMPSet * d) const {
    if (MS.NEARBYNESS >= 0)
        getBFC_matchDescriptors_RI<LIMIT_DISTANCES, false > (nStart, nEnd, pDS, MS, d);
    else
        getBFC_matchDescriptors_RI<LIMIT_DISTANCES, true > (nStart, nEnd, pDS, MS, d);
}
void CMatchableDescriptors::getBFC_matchDescriptorsMT(const int nStart, const int nEnd, const CMatchableDescriptors * pDS, const CMatchableDescriptors::CMatchSettings & MS, TMPSet * d) const {
    if (MS.MAX_SEPERATION > 0)
        getBFC_matchDescriptors < true > (nStart, nEnd, pDS, MS, d);
    else
        getBFC_matchDescriptors < false > (nStart, nEnd, pDS, MS, d);
}
const CBoWCorrespondences * CMatchableDescriptors::getBruteForceCorrespondenceSet_int(const CMatchableDescriptors * pDS, const CMatchableDescriptors::CMatchSettings & MS, CBoWCorrespondences * pCorrIn) const {
    const int nCount1 = Count();
    const int nCount2 = pDS->Count();

    const double dCondition = MS.MATCH_CONDITION;
    const double dPP = MS.PP;
    const int NN = MS.NN;

    ARRAY(TMPSet, d, nCount1);

    if (nCount1 < 30 || nCount2 < 30) //Single thread. TODO argument?
        getBFC_matchDescriptorsMT(0, nCount1, pDS, MS, PTR(d));
    else {
        int nMid = nCount1 / 2;
        boost::thread firstHalf(boost::bind(&CMatchableDescriptors::getBFC_matchDescriptorsMT, this, 0, nMid, pDS, boost::ref(MS), PTR(d)));
        getBFC_matchDescriptorsMT(nMid, nCount1, pDS, MS, PTR(d));
        firstHalf.join();
    }

    CBoWCorrespondences * pCorr = pCorrIn;
    if (!pCorrIn)
        pCorr = new CBoWCorrespondences(/*nCount1 * NN*/);

    ARRAYZ(int, anNumTimesLeftMatched, nCount1);
    ARRAYZ(int, anNumTimesRightMatched, nCount2);

    for (int i = 0; i < nCount1; i++) {
        int nNumNN = 0;

        if (d[i].size() > 0) {
            double oldDist = d[i].top().first;
            for (TMPSet::const_iterator pMP = d[i].begin(); pMP != d[i].end(); pMP++, nNumNN++) {
                double dNewDist = pMP->first;
                if (oldDist < dNewDist * dCondition) {
                    //We're at the end of the n-n
                    break;
                } else {
                    //anNumTimesRightMatched[pMP->second]++;
                    oldDist = dNewDist;
                }
            }
        }

        anNumTimesLeftMatched[i] = nNumNN;

        if (nNumNN <= NN) {
            for (TMPSet::const_iterator pMP = d[i].begin(); nNumNN > 0; pMP++, nNumNN--) {
                anNumTimesRightMatched[pMP->second]++;
            }
        }
    }

    for (int i = 0; i < nCount1; i++) {
        int nNumMatchedLeft = anNumTimesLeftMatched[i];
        if (nNumMatchedLeft <= NN) {
            for (TMPSet::const_iterator pMP = d[i].begin(); nNumMatchedLeft > 0; pMP++, nNumMatchedLeft--) {
                int j = pMP->second;
                if (anNumTimesRightMatched[j] <= NN) {
                    double dProb = dPP / std::max<int>(anNumTimesLeftMatched[i], anNumTimesRightMatched[j]); //Total prob for first DS is ok, not sure about 2nd. Doesn't matter if not disjoint anyway.
                    pCorr->push_back(CCorrespondence(get_const(i)->location(), pDS->get_const(j)->location(), dProb));
                }
            }
        }
    }

    return pCorr;
}
const CBoWCorrespondences * CMatchableDescriptors::getBruteForceCorrespondenceSet(const CMatchableDescriptors * pDS, const CMatchSettings & MS, CBoWCorrespondences * pCorr) const {
    //cout << "DS Sizes: " << Count() << ',' << pDS->Count() << endl;
    return getBruteForceCorrespondenceSet_int(pDS, MS, pCorr);
}
void CSimpleMatchableDescriptors::Push_const(const CDescriptor * pDescriptor) {
    descriptors.push_back(pDescriptor);
}
void CSimpleMatchableDescriptors::Push_const(const CMatchableDescriptors * pDescriptorSet) {
    descriptors.reserve(descriptors.size() + pDescriptorSet->Count());
    for (int i = 0; i < pDescriptorSet->Count(); i++)
        Push_const(pDescriptorSet->get_const(i));
}

//Only using a set to avoid duplicates so should be fast
class descriptorFastCompare : binary_function<const CDescriptor *, const CDescriptor *, bool> {
public:
    result_type operator() (first_argument_type a,
            second_argument_type b) {
        return (result_type) (a < b);
    };
};

//Adds a random subset of this set to pRandSubset, including the cluster centres
void CDescriptorSet::RandomSubset(int n, CDescriptorSet * pRandSubset, CClusterSet * pBestClusters, bool bForceSeperation) {
    int nCount = (int) Count();
    if(IS_DEBUG) CHECK(pRandSubset->Count() || (int) n > nCount || (pBestClusters && n < (int) pBestClusters->size()), "CDescriptorSet::RandomSubset: Bad params");

    typedef set2_NF<CDescriptor *, descriptorFastCompare> TDescSet;
    TDescSet descSet;

    //Include cluster centres from each try in next random set (5 times)
    if (pBestClusters) {
        int nClusters = (int) pBestClusters->size();

        for (int i = 0; i < nClusters; i++) {
            CDescriptor * pCentre = (*pBestClusters)[i]->Centre();
            descSet.insert(pCentre);
        }
    }

    int insertPos = 0; //Try to add sequentially to avoid random access cache misses
    const int nStepMax = doubleToInt((2 * (int) nCount) / (double) (n - descSet.size()));

    if (bForceSeperation) //force points to be far from any existing points (to seed k-means)
    {
        int nSkipped = 0;
        while ((int) descSet.size() < n) {
            insertPos = (insertPos + CRandom::Uniform(1, nStepMax)) % nCount;

            CDescriptor * pDesc = get(insertPos);
            if (Closest(pDesc) > 25 || nSkipped > 30) //todo: constants
                descSet.insert(pDesc);
            else
                nSkipped++;
        }
    } else {
        while ((int) descSet.size() < n) {
            insertPos = (insertPos + CRandom::Uniform(1, nStepMax)) % nCount;
            descSet.insert(get(insertPos));
        }
    }

    for (TDescSet::iterator ppDesc = descSet.begin(); ppDesc != descSet.end(); ppDesc++)
        pRandSubset->Push(*ppDesc);
}
bool CDescriptorSet::Contains(const CDescriptor * pDescriptor) const {
    for (int i = Count() - 1; i >= 0; i--) //likely to be towards the end
    {
        if (get_const(i) == pDescriptor)
            return true;
    }
    return false;
}
const CDescriptor * CDescriptorSet::ClosestDescriptor(const CDescriptor * pDescCompareTo) const {
    CDescriptor::TDist dClosestDist = MAX_ALLOWED_DIST;
    const CDescriptor * pClosestDesc = 0;
    for (int nDesc = 0; nDesc < Count(); nDesc++) {
        const CDescriptor * pDesc = get_const(nDesc);
        CDescriptor::TDist dDist = pDesc->distance(pDescCompareTo);

        if (dDist < dClosestDist) {
            dClosestDist = dDist;
            pClosestDesc = pDesc;
        }
    }
    return pClosestDesc;
}
void CDescriptorSet::Closest2Descriptors(const CDescriptor * pDescCompareTo, CDescriptor const** ppClosest, CDescriptor const** pp2ndClosest) const {
    CDescriptor::TDist dClosestDist = MAX_ALLOWED_DIST, d2ndClosestDist = MAX_ALLOWED_DIST;
    //const CDescriptor * pClosestDesc = 0;
    for (int nDesc = 0; nDesc < Count(); nDesc++) {
        const CDescriptor * pDesc = get_const(nDesc);
        CDescriptor::TDist dDist = pDesc->distance(pDescCompareTo);

        if (dDist < dClosestDist) {
            d2ndClosestDist = dClosestDist;
            *pp2ndClosest = *ppClosest;
            dClosestDist = dDist;
            *ppClosest = pDesc;
        } else if (dDist < d2ndClosestDist) {
            d2ndClosestDist = dDist;
            *pp2ndClosest = pDesc;
        }
    }
}
const CDescriptor * CDescriptorSet::ClosestDescriptor(const CDescriptor * pDescCompareTo, double dCondition, int nRadius) const {
    CDescriptor::TDist dClosestDist = MAX_ALLOWED_DIST, dNextClosestDist = 0;
    const CDescriptor * pClosestDesc = 0;
    for (int nDesc = 0; nDesc < Count(); nDesc++) {
        const CDescriptor * pDesc = get_const(nDesc);
        CDescriptor::TDist dDist = pDesc->distance(pDescCompareTo);

        if (dDist < dClosestDist) {
            dNextClosestDist = dClosestDist;
            dClosestDist = dDist;
            pClosestDesc = pDesc;
        }
    }
    if (nRadius < dClosestDist) return 0;
    if ((double) dNextClosestDist * dCondition > (double) dClosestDist) return pClosestDesc;
    return 0;
}

//Return the quality of this match: How many are closer, and the ratio of the average match strength to the closest match
void CDescriptorSet::descriptorQuality(const CDescriptor * pDescThis, const CDescriptor * pDescOther, double & dQuality, int & nCloser) const {
    double dAvDist = 0;
    nCloser = 0;
    CDescriptor::TDist dThisDist = pDescThis->distance(pDescOther);

    for (int nDesc = 0; nDesc < Count(); nDesc++) {
        CDescriptor::TDist dDist = pDescOther->distance(get_const(nDesc));
        dAvDist += dDist;

        if (dDist < dThisDist) {
            nCloser++;
        }
    }
    dAvDist /= Count();
    dQuality = dAvDist / dThisDist;
}
void CDescriptorSet::clearClusterAssignments() {
    for (int nDesc = 0; nDesc < Count(); nDesc++) {
        get(nDesc)->assignToCluster(0, 0);
    }
}
CDescriptor::TDist CDescriptorSet::Closest(const CDescriptor * pDesc) const {
    CDescriptor::TDist dClosestDist = MAX_ALLOWED_DIST;
    for (int nDesc = 0; nDesc < Count(); nDesc++) {
        CDescriptor::TDist dDist = pDesc->distance(get_const(nDesc));

        if (dDist < dClosestDist) {
            dClosestDist = dDist;
        }
    }
    return dClosestDist;
}

///// Random clustering ///////
CClusterSet * CRandomDescriptorSet::Cluster(int nClusters, const CCluster * pParentCluster) //Override k-medoids
{

    CRandomDescriptorSet * pRandSubset = this->makeNewDS(0);
    RandomSubset(nClusters, pRandSubset, 0, true);

    CClusterSet * pClusters = new CClusterSet(pRandSubset->Count(), 0);

    for (int nCentre = 0; nCentre < pRandSubset->Count(); nCentre++) {
        cout << "TODO: Uncomment me\n";
        //CDescriptorSet * pClusterMembers = this->makeNewDS(0);pClusters->push_back(CCluster((*pRandSubset)[nCentre], false, &pClusterMembers, 0, 0, 0, nCentre, pParentCluster));
        cout << "TODO: Rand cluster needs counts etc. for fast assignment\n";
    }

    pClusters->AssignToClusters < false > (this, 1);

    return pClusters;
}

/*void CClusterSet::RecalcCentres(eKMImplementation weightMethod)
{
    THROW( "Uncomment next lines")
    //For each cluster...
    / *for(CClusterSet::iterator pCluster = begin(); pCluster != end(); pCluster++)
    {
        pCluster->RecalcCentre(weightMethod);
    }* /
}*/

#define KM_VERBOSE(c)
///*
/*void CCluster::RecalcCentre(eKMImplementation weightMethod) //Only works for type char DO NOT DELETE! Works for char[] descriptors
{
    const TDescriptor * pVCentre = CAST<const TDescriptor *>(pClusterCentre);
    if(IS_DEBUG) CHECK(!bIOwnCentreMem || !pVCentre, "CCluster::RecalcCentre: This method is for k-means and needs vector descriptors and clusters owning their centre");

    if(!Count())
    {
        cout << "ERROR: Cluster with no centres!\n";
        return;
    }

    int nVectorLen = TDescriptor::DescriptorLength();

    TS_STATIC int * aTempVector = new int[nVectorLen];

    memset(aTempVector, 0, nVectorLen*sizeof(int));

    int nSumWeights = 0;

    KM_VERBOSE(cout << Count()  << " descriptors in this cluster\n");

    //Calc weight, sum descriptors * weight
    for(unsigned int nDesc = 0; nDesc < Count(); nDesc++)
    {
        const int SCALE = 10000; //should accomodate up to to 400k descriptors with lengths generally >1
        const CDescriptor::TDist MAX_RAD = 1;

        const CDescriptor * pDescriptor = / *CAST<const TDescriptor *>* /((*pvMembers)[nDesc]);

        CDescriptor::TDist dist;

        int nWeight;
        switch(weightMethod) //todo: template it
        {
        case eLloyds:
            nWeight = 1;
            break;
        case eKMWeighted:
            dist = Centre()->distance(pDescriptor);
            nWeight = (SCALE/(dist+1));
            nWeight = max<int>(nWeight, 1);
            break;
        case eKMFixedRadius:
            dist = Centre()->distance(pDescriptor);
            nWeight = (dist > MAX_RAD) ? 0 : 1;
            break;
        default:
            throw new CException("CCluster::RecalcCentre: Bad comparison method");
        }
        nSumWeights += nWeight;
        KM_VERBOSE(cout << nWeight  << "=w " << dist << "=d " );

        const CDescriptorFactory::elType * pDescriptorVec = pDescriptor->DescriptorVector();

        //Add weighted vector
        int * pTempVector = aTempVector;
        for(int i=nVectorLen; i>0; i--)
        {
 *pTempVector += ((int) *pDescriptorVec) * nWeight;
            pTempVector++;
            pDescriptorVec++;
        }
    }
    KM_VERBOSE(cout << "=weights\n");

    //Divide by sum weight and copy into centre vector
    int * pTempVector = aTempVector;
    TDescriptor::elType * pCentreVector = const_cast<CDescriptorFactory::elType *>(pVCentre->DescriptorVector()); //eeugh
    KM_VERBOSE(cout << "Centre =");
    for(int i=nVectorLen; i>0; i--)
    {
 *pCentreVector = (CDescriptorFactory::elType)(*pTempVector / nSumWeights) ;
        KM_VERBOSE(cout << (int)*pCentreVector << ' ');
        pTempVector++;
        pCentreVector++;
    }
    KM_VERBOSE(cout << '\n');

    KM_VERBOSE(cout << nSumWeights << "=sum weights\n");

    //Now check the centre is still near some of its old members
    //dTotDist = 0;
    //for(unsigned int nDesc = 0; nDesc < Count(); nDesc++)
    //{
    //    const TDescriptor * pDescriptor = CAST<const TDescriptor *>((*pvMembers)[nDesc]);
    //    dTotDist += Centre()->distance(pDescriptor);
    //}
    //cout << dTotDist << "=dist after\n";
    TS_STATIC_DELETE( [] aTempVector);
} //DO NOT DELETE! Works for char[] descriptors*/

////Hacky way of creating new descriptors
//CDescriptor * DupDescriptor(const CDescriptor * pCentre/*, size_t descriptorSize*/)
//{
//	return new TDescriptor(sizeof(void *)+(char*)(void*)pCentre);
//	/*
//    CDescriptor * pDescriptor = (CDescriptor *)malloc(descriptorSize);
//    memcpy(pDescriptor, pCentre, descriptorSize);
//    return pDescriptor;*/
//}
