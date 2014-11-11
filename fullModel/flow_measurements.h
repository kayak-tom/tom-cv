#ifndef CFLOWMEASUREMENTS_H
#define CFLOWMEASUREMENTS_H

#include "main.h"

/*class CTimeStamp
{
	int nSecondsSinceStart;
public:
	CTimeStamp(const int nSecondsSinceStart) : nSecondsSinceStart(nSecondsSinceStart)
	{
	
	}
	
	int getSecondsSinceStart() const { return nSecondsSinceStart; }
	
	bool operator<(const CTimeStamp & other) const { return getSecondsSinceStart()<other.getSecondsSinceStart(); }
};*/

/*class CFlowMeasurement
{
	double dFlow;
public:
	CFlowMeasurement(const double dFlow) : dFlow(dFlow) 
	{
		
	}
	
	double getMeasurem() const { return dFlow; }
};*/

class CInputArray
{
public:
	typedef std::map<TTime, double> TMeasurements;
private:
	TMeasurements aMeasurements;
	void loadFlows(const std::string & strFullFilename, const std::string & strVariable);
public:
	CInputArray(const std::string strNetCDFDir, const std::string strVariable, const int nReachIdx);
	CInputArray(const std::string strNetCDFDir, const std::string strVariable);
	const TMeasurements & getMeasurements() const { return aMeasurements; }

	virtual ~CInputArray()
	{
	}

};

#endif // CFLOWMEASUREMENTS_H
