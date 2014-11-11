#ifndef CSUBCATCHMENT_BASE_H
#define CSUBCATCHMENT_BASE_H

#include "main.h"

class CSubcatchment_base
{
protected:
	const TTimeDuration timeStep;
	const CSubcatchmentParams * pSubcatchmentParams;
public:
	CSubcatchment_base(const CSubcatchmentParams * pSubcatchmentParams,  const TTimeDuration timeStep) : timeStep(timeStep),pSubcatchmentParams(pSubcatchmentParams)
	{
	}
	
	virtual ~CSubcatchment_base()
	{
	}
	
	
	const CSubcatchmentParams * getParams() const { return pSubcatchmentParams; }
	
	virtual void setRainfallMeasurements(const CInputArray & rainfall, const TTime & timeFrom, const TTime & timeTo) = 0;
	virtual void setPETMeasurements(const CInputArray & PET, const TTime & timeFrom, const TTime & timeTo) = 0;
	virtual void setOutflowMeasurements(const CInputArray & flow, const TTime & timeFrom, const TTime & timeTo) = 0;
    
};

#endif // CSUBCATCHMENT_BASE_H
