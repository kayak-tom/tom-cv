#ifndef CMODEL_BASE_H
#define CMODEL_BASE_H

#include "main.h"

class CModel_base
{
protected:
	const TTimeDuration timeStep;
	const CCatchmentSetupParams * pCatchmentStructure;
	std::vector<CSubcatchment_base *> aAllSubcatchments;
public:
	CModel_base(const CCatchmentSetupParams * pCatchmentStructure, const TTimeDuration timeStep) : timeStep(timeStep), pCatchmentStructure(pCatchmentStructure)
	{
	}
	virtual ~CModel_base()
	{
	}
	
	virtual CSubcatchment_base * setupSubcatchment(const CSubcatchmentParams * pSubcatchmentParams, CSubcatchment_base * pDownstreamCatchment) = 0;
	
	void setupSubcatchments_recurse(const CSubcatchmentParams * pSubcatchmentParams, CSubcatchment_base * pDownstreamCatchment);
	void setupSubcatchments();
	static CModel_base * makeModel(const CCatchmentSetupParams * pCatchmentStructure, const TTimeDuration & fiveMinsTimeStep); //Also sets up the subcatchment structure

	//The rainfalls should have an uninformative distn when setup.
	void setupStates(const TTime & timeFrom,const TTime & timeTo);
	
	void setRainfallMeasurements(const std::string strNCFilename, const TTime & timeFrom, const TTime & timeTo);
	void setPETMeasurements(const std::string strNCFilename, const TTime & timeFrom, const TTime & timeTo);
	void setOutflowMeasurements(const std::string strNCFilename, const TTime & timeFrom, const TTime & timeTo);
	
	virtual void run() = 0;
	virtual void optimise() = 0;
	virtual void saveStates(const std::string strTSVFilename) const = 0;
	virtual void completeSetup() = 0;
    
    CSubcatchment_base * getDownstreamCatchment();
};

#endif // CMODEL_BASE_H
