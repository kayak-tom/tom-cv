/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#include "descriptor.h"
#include <set>
#include <functional>

using namespace std;

CCluster::CCluster(CDescriptor * pCentre_in, bool bIOwnCentreMem_in, CDescriptorSet ** ppDescriptorSet, CDescriptor::TDist nHalfMinDist, CDescriptor::TDist nRadiusSD, int nDescriptorCountKMed, const int nMyIdx, const CCluster * pParentCluster)
   	    : pClusterCentre(pCentre_in), pvMembers(0), bIOwnCentreMem(bIOwnCentreMem_in),
   	      nHalfMinDist(nHalfMinDist), nRadiusSD(nRadiusSD), nDescriptorCountKMed(nDescriptorCountKMed), nMyIdx(nMyIdx), pParentCluster(pParentCluster)
{
	if(IS_DEBUG) CHECK(pParentCluster && UNINIT(pParentCluster), "Bad centre set");
    if(ppDescriptorSet)
    {
        AssignDescriptorSet(ppDescriptorSet);
    }
}
/*CCluster & CCluster::operator=(CCluster & otherCluster)
{
	pClusterCentre = otherCluster.pClusterCentre; otherCluster.pClusterCentre = 0;
	pvMembers = otherCluster.pvMembers; otherCluster.pvMembers = 0;
	bIOwnCentreMem = otherCluster.bIOwnCentreMem;
	nHalfMinDist = otherCluster.nHalfMinDist;
	nRadiusSD = otherCluster.nRadiusSD;
	nDescriptorCountKMed = otherCluster.nDescriptorCountKMed;
	nMyIdx = otherCluster.nMyIdx;
	pParentCluster = otherCluster.pParentCluster;
	return *this;
}
CCluster::CCluster(CCluster & otherCluster)
{
	pClusterCentre = otherCluster.pClusterCentre; otherCluster.pClusterCentre = 0;
	pvMembers = otherCluster.pvMembers; otherCluster.pvMembers = 0;
	bIOwnCentreMem = otherCluster.bIOwnCentreMem;
	nHalfMinDist = otherCluster.nHalfMinDist;
	nRadiusSD = otherCluster.nRadiusSD;
	nDescriptorCountKMed = otherCluster.nDescriptorCountKMed;
	nMyIdx = otherCluster.nMyIdx;
	pParentCluster = otherCluster.pParentCluster;
}*/

void CCluster::AssignDescriptorSet(CDescriptorSet ** ppDescriptorSet)
{
    if(IS_DEBUG) CHECK(!ppDescriptorSet || !*ppDescriptorSet || pvMembers, "CCluster::AssignDescriptorSet: Bad parameter");
    pvMembers = *ppDescriptorSet; *ppDescriptorSet = 0;
}

CCluster::~CCluster()
{
    if(bIOwnCentreMem)
        delete pClusterCentre;

    delete pvMembers;
}

//Output the smallest distance between a pair of centres
void CClusterSet::CheckDisimilarity() const
{
    //For each centre, find the closest centre
    unsigned int nCentres = (unsigned int)size();
    double dClosestDist=MAX_ALLOWED_DIST;

    for(unsigned int nCentre=0;nCentre<nCentres; nCentre++)
    {
        for(unsigned int nCentreComp=nCentre+1;nCentreComp<nCentres; nCentreComp++)
        {
            double dDist = (*this)[nCentreComp]->Centre()->distance((*this)[nCentre]->Centre());
            //cout << dDist << "=disimilarity\n";
            if(dDist < dClosestDist) dClosestDist=dDist;
        }
    }
    cout << dClosestDist << "=least disimilarity\n";
}

void CClusterSet::CheckAssignment() const
{
    //Check each centre has some descriptors assigned!
    unsigned int nCentres = (unsigned int)size();

    for(unsigned int nCentre=0;nCentre<nCentres; nCentre++)
    {
        cout << ((*this)[nCentre])->Count() << ' ';
    }
    cout << "=cluster counts\n";
}

//Clear descriptor sets so we can reassign
void CClusterSet::ResetClusters()
{
    for(CClusterSet::iterator ppCluster=begin(); ppCluster<end(); ppCluster++)
    {
    	CCluster * pCluster = *ppCluster;
        pCluster->Members()->Clear();
    }

    dTotalDistance = NOT_ASSIGNED();
}

void CClusterSet::remindDescriptorsOfAssignments()
{
    for(CClusterSet::iterator ppCluster=begin(); ppCluster<end(); ppCluster++)
    {
    	CCluster * pCluster = *ppCluster;
        pCluster->Members()->assignToCluster(pCluster);
    }
}

#define LogClusterQuality
#ifdef LogClusterQuality

#else
double CClusterSet::AssignToClusters(const CDescriptorSet * pDescriptors) const
{
    double dTotalDist = 0;
    unsigned int nDescriptors = pDescriptors->Count();
    for(unsigned int nDescriptor=0; nDescriptor < nDescriptors; nDescriptor++)
    {
        double dClosestDist=MAX_ALLOWED_DIST;
        int nClosestCluster=-1;
        for(unsigned int nCluster=0; nCluster<size(); nCluster++)
        {
            double dDist=(*pDescriptors)[nDescriptor]->distance((*this)[nCluster]->Centre());
            if (dDist<dClosestDist)
            {
                dClosestDist=dDist;
                nClosestCluster=nCluster;
            }
        }
        if(IS_DEBUG) CHECK(nClosestCluster<0, "CClusterSet::AssignToClusters: Closest cluster not found");

        if(IS_DEBUG) CHECK(pDescriptors->Count() <= nDescriptor, "CClusterSet::AssignToClusters: nDescriptor out of bounds");

        (*this)[nClosestCluster]->Push((*pDescriptors)[nDescriptor]);
        dTotalDist += dClosestDist;
    }
    dTotalDistance = dTotalDist;
    return dTotalDist;
}
#endif

//Doesn't calc dist (+ is quiet)
/*void CClusterSet::QuickAssignToClusters(const CDescriptorSet * pDescriptors) const
{
    unsigned int nDescriptors = pDescriptors->Count();
    for(unsigned int nDescriptor=0; nDescriptor < nDescriptors; nDescriptor++)
    {
        double dClosestDist=MAX_ALLOWED_DIST;
        int nClosestCluster=-1;
        for(unsigned int nCluster=0; nCluster<size(); nCluster++)
        {
            double dDist=(*pDescriptors)[nDescriptor]->distance((*this)[nCluster]->Centre());
            if (dDist<dClosestDist)
            {
                dClosestDist=dDist;
                nClosestCluster=nCluster;
            }
        }
        if(IS_DEBUG) CHECK(nClosestCluster<0, "CClusterSet::AssignToClusters: Closest cluster not found");
        if(IS_DEBUG) CHECK(pDescriptors->Count() <= nDescriptor, "CClusterSet::AssignToClusters: nDescriptor out of bounds");

        (*this)[nClosestCluster]->Push((*pDescriptors)[nDescriptor]);
    }
    dTotalDistance = 0;
}*/
