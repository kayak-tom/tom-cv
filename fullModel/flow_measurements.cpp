#include "flow_measurements.h"
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <netcdfcpp.h>

CInputArray::CInputArray(const std::string strNetCDFDir, const std::string strVariable, const int nReachIdx)
{
	if(!boost::filesystem::exists(strNetCDFDir))
		throw "strNetCDFDir doesn't exist";

	NcFile inputDataFile(strNetCDFDir.c_str());

	if (!inputDataFile.is_valid()) {
		throw "Couldn't open file!";
	}

	NcVar * pVar = inputDataFile.get_var(strVariable.c_str());
	NcAtt* pTime = inputDataFile.get_att("forecast_reference_time");
//	cout << pTime->as_string(0) << endl;

	const int nDims = pVar->num_dims(); //Some values have ensemble # as a 3rd dimension (1st is non-perturbed)

	const int nNumTimesteps = pVar->get_dim(0)->size();
	const int nNumReaches = pVar->get_dim(1)->size();
	const int nNumEpochs = (nDims > 2) ? pVar->get_dim(2)->size() : 1;
	cout << "nNumTimesteps=" << nNumTimesteps << endl;

	const int nNumVals = nNumEpochs*nNumReaches*nNumTimesteps;

	double * adData = new double[nNumVals];
	pVar->get(adData, nNumTimesteps, nNumReaches, (nDims > 2) ? nNumEpochs : 0);

	const TTime startTime = boost::date_time::parse_delimited_time<boost::posix_time::ptime>( pTime->as_string(0), ' ' );

	double dMinMeasurement = HUGE, dMaxMeasurement = -HUGE;
	int nIdx = nReachIdx * nNumTimesteps * nNumEpochs;
	for(int nTime=0; nTime<nNumTimesteps; nTime++) {
		const double dVal = adData[nIdx];
		//cout << strVariable << ": nTime=" << nTime << " nReach=" << nReachIdx << " nIdx = " << nIdx << " Val = " << dVal << endl;

		const TTimeDuration offset = boost::posix_time::hours(nTime);
		const TTime thisTime = startTime+offset;

		if (dVal < dMinMeasurement)
			dMinMeasurement = dVal;
		if (dVal > dMaxMeasurement)
			dMaxMeasurement = dVal;

		if(dVal>=0)
			aMeasurements[thisTime] = dVal;
		else
			cout << "Warning: negative idx" << endl;
		nIdx+=nNumEpochs;
	}

	cout << "# measurements: " << aMeasurements.size() << endl;
	cout << "dMinMeasurement: " << dMinMeasurement << endl;
	cout << "dMaxMeasurement: " << dMaxMeasurement << endl;

	if(aMeasurements.size() == 0)
		throw "Zero measurements loaded";
}

void CInputArray::loadFlows(const std::string & strFullFilename, const std::string & strVariable)
{
	if(!boost::filesystem::exists(strFullFilename))
		throw "strFullFilename doesn't exist";

	NcFile inputDataFile(strFullFilename.c_str());

	if (!inputDataFile.is_valid()) {
		throw "Couldn't open file!";
	}

	cout << "Load 'flow' variable from " << strFullFilename << endl;
	NcVar * pVar = inputDataFile.get_var(strVariable.c_str());
	//NcAtt* pTime = inputDataFile.get_att();
	//cout << pTime->as_string(0) << endl;

	const int nDims = pVar->num_dims();
	if(nDims != 2)
		throw "Unlikely to work for other flow dims or multiple stations";

	const int nNumTimesteps = pVar->get_dim(0)->size();
	const int nNumStations = pVar->get_dim(1)->size();
	cout << "nNumTimesteps=" << nNumTimesteps << endl;

	const int nNumVals = nNumStations*nNumTimesteps;

	double * adData = new double[nNumVals];
	pVar->get(adData, nNumTimesteps, nNumStations);

	//const TTime startTime = boost::posix_time::from_iso_string(pTime->as_string(0)); //(boost::gregorian::date(2010, boost::gregorian::Sep, 28), boost::posix_time::time_duration(0,0,0));
	//const TTime startTime = boost::date_time::parse_delimited_time<boost::posix_time::ptime>( pTime->as_string(0), ' ' );
	const int startPos = strFullFilename.find("streamq_") + strlen("streamq_");
	const int length = strFullFilename.find("_utc") - startPos;
	const std::string timeString = strFullFilename.substr(startPos, length) + "0000";

	cout << timeString << endl;

	//const TTime startTime = boost::posix_time::from_iso_string(timeString);
	boost::gregorian::date date = boost::gregorian::date(atoi(timeString.substr(0,4).c_str()), atoi(timeString.substr(4,2).c_str()), atoi(timeString.substr(6,2).c_str()));
	TTimeDuration time = boost::posix_time::time_duration(atoi(timeString.substr(8,2).c_str()),0,0);
	const TTime startTime(date, time);

	double dMinMeasurement = HUGE, dMaxMeasurement = -HUGE;
	for(int nTime=0; nTime<nNumTimesteps; nTime++) {
		const double dVal = adData[nTime];

		const TTimeDuration offset = boost::posix_time::hours(nTime);
		const TTime thisTime = startTime+offset;
		cout << strVariable << ": nTime=" << nTime << " time=" << thisTime << " Val = " << dVal << endl;

		if (dVal < dMinMeasurement)
			dMinMeasurement = dVal;
		if (dVal > dMaxMeasurement)
			dMaxMeasurement = dVal;

		if(dVal>=0)
			aMeasurements[thisTime] = dVal;
		else
			cout << "Warning: negative idx" << endl;
	}

	cout << "# measurements: " << aMeasurements.size() << endl;
	cout << "dMinMeasurement: " << dMinMeasurement << endl;
	cout << "dMaxMeasurement: " << dMaxMeasurement << endl;

	/*if(aMeasurements.size() == 0)
		//throw "Zero measurements loaded";*/
}

CInputArray::CInputArray(const std::string strNetCDFDir, const std::string strVariable)
{
	//A bit of a script specific for the flows in each file
	const std::string cacheFlows = "cacheflows.tsv";

	if(!boost::filesystem::exists(cacheFlows)) {
		boost::filesystem::path targetDir(strNetCDFDir);
		boost::filesystem::directory_iterator it(targetDir), eod;
		std::string pattern = "streamq";

		BOOST_FOREACH(boost::filesystem::path const &p, std::make_pair(it, eod)) {
			if(boost::filesystem::is_regular_file(p)) {
				if(p.string().find(pattern) != std::string::npos) {
					loadFlows(p.string(), strVariable);
				}
			}
		}

		std::ofstream cacheFlowsTSV(cacheFlows.c_str());

		for(CInputArray::TMeasurements::const_iterator measurement = getMeasurements().begin(); measurement != getMeasurements().end(); measurement++) {
			cacheFlowsTSV << measurement->first << '\t' << measurement->second << endl;
		}
	}

	std::ifstream cacheFlowsTSV(cacheFlows.c_str());
	while(!cacheFlowsTSV.eof()) {
		TTime time;
		double dFlow = -1;
		cacheFlowsTSV >> time;
		cacheFlowsTSV >> dFlow;

		aMeasurements[time] = dFlow;
	}
}
