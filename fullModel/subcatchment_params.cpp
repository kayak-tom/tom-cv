#include "subcatchment_params.h"
#include "catchment_setup_params.h"
#include <iostream>

using std::cout;
using std::endl;

double CSubcatchmentParams::getParam(const char * szParamName) const
{
	return pParent->getParam(szParamName, nrch);
}

CSubcatchmentParams::CSubcatchmentParams(const int nrch, const int rchid, const CCatchmentSetupParams * pParent) : nrch(nrch), rchid(rchid), pParent(pParent), pDownstreamCatchment(0)
{
	
	cout << "Creating subcatchment with rchid " << nrch << endl;
	cout << "Manning's n " << getParam("rchman_n") << endl;
	cout << "Area km2 " << getParam("basarea")/1000000 << endl;
	
}

CSubcatchmentParams::~CSubcatchmentParams()
{
}

void CSubcatchmentParams::setDownstreamCatchment(CSubcatchmentParams * pDownstreamCatchment_in)
{
	pDownstreamCatchment = pDownstreamCatchment_in;
	pDownstreamCatchment->addUpstreamCatchment(this);
}

void CSubcatchmentParams::addUpstreamCatchment(CSubcatchmentParams * pUpstreamCatchment_in)
{
	apUpstreamCatchments.push_back(pUpstreamCatchment_in);
}

double CSubcatchmentParams::getC_canopy() const
{
    return getParam("canscap");
}

double CSubcatchmentParams::getC_rfactor() const
{
    return getParam("canenhf");
}

double CSubcatchmentParams::getHydraulicCon() const
{
    return getParam("hydcon0");
}

double CSubcatchmentParams::getTopnet_m() const
{
    return getParam("topmodm");

}


