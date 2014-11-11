#include "model_base.h"
#include "catchment_setup_params.h"
#include <boost/foreach.hpp>
#include "subcatchment_params.h"
#include "flow_measurements.h"
#include "subcatchment_base.h"
#include "lm_block_model.h"

CModel_base * CModel_base::makeModel(const CCatchmentSetupParams * pCatchmentStructure, const TTimeDuration & timeStep)
{
	CModel_base * pNewCatchment = new CLMBlockModel(pCatchmentStructure, timeStep);	
	
	return pNewCatchment;
}

void CModel_base::setupSubcatchments()
{
	//Find downstream-most catchment:
	const CSubcatchmentParams * pMostDownstream = pCatchmentStructure->getDownstreamCatchment();
	
	setupSubcatchments_recurse(pMostDownstream, 0);
}

/* static */ void CModel_base::setupSubcatchments_recurse(const CSubcatchmentParams * pSubcatchmentParams, CSubcatchment_base * pDownstreamCatchment)
{
	CSubcatchment_base * pNewSubcatchment = setupSubcatchment(pSubcatchmentParams, pDownstreamCatchment);
	aAllSubcatchments.push_back(pNewSubcatchment);
	
	BOOST_FOREACH(const CSubcatchmentParams * pUpstreamSubcatchmentParams, pSubcatchmentParams->getUpstreamCatchments())
	{
		setupSubcatchments_recurse(pUpstreamSubcatchmentParams, pNewSubcatchment);
	}

}

//Setup state only when needed
/*void CModel_base::setupStates(const TTime & timeFrom,const TTime & timeTo)
{
	
	
}*/
	
void CModel_base::setRainfallMeasurements(const std::string strNCFilename, const TTime & timeFrom, const TTime & timeTo)
{
	BOOST_FOREACH(CSubcatchment_base * pSubcatchment, aAllSubcatchments)
	{
		CInputArray rainfall(strNCFilename, "aprecip", pSubcatchment->getParams()->getNrch());
		pSubcatchment->setRainfallMeasurements(rainfall, timeFrom, timeTo);
	}

}

void CModel_base::setPETMeasurements(const std::string strNCFilename, const TTime & timeFrom, const TTime & timeTo)
{
	BOOST_FOREACH(CSubcatchment_base * pSubcatchment, aAllSubcatchments)
	{
		CInputArray PET(strNCFilename, "potevap", pSubcatchment->getParams()->getNrch());
		pSubcatchment->setPETMeasurements(PET, timeFrom, timeTo);
	}
}

void CModel_base::setOutflowMeasurements(const std::string strNCFilename, const TTime & timeFrom, const TTime & timeTo)
{
	CInputArray flow(strNCFilename, "flow"); 

	getDownstreamCatchment()->setOutflowMeasurements(flow, timeFrom, timeTo);
}	
CSubcatchment_base * CModel_base::getDownstreamCatchment()
{
    BOOST_FOREACH(CSubcatchment_base * pSubcatchment, aAllSubcatchments)
	{
		if(!pSubcatchment->getParams()->getDownstreamCatchment())
			return pSubcatchment;
    }	
    throw "No downstream catchment";
}
	