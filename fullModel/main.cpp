#include <stdio.h>
#include <iostream>

#include "catchment_setup_params.h"
#include <cstdio>
#include <boost/scoped_ptr.hpp>
#include "model_base.h"
#include "flow_measurements.h"

using boost::scoped_ptr;

int main(int argc, char **argv)
{
	const boost::posix_time::time_duration fiveMinsTimeStep(0,5,0), oneHour(1,0,0);

	// load params from file
	std::string szFilenameForCatchment = "/media/tom/My Passport/Einar_backup/data/mcmillanhk/topnet_hiltom/ancillary_data/topnet/spatial_rec_04017064_strahler3.nc";
	CCatchmentSetupParams params(szFilenameForCatchment);
	
	//Setup model
	scoped_ptr<CModel_base> pModel ( CModel_base::makeModel(&params, oneHour) );
	pModel->setupSubcatchments();
	
	std::string strCombinedFilename = "/media/tom/My Passport/Einar_backup/data/mcmillanhk/topnet_hiltom/output/merged/streamq_2010092800_2012092800_00-06_utc_topnet_04017064_strahler3-RK.nc";
	std::string strFlowDir = "/media/tom/My Passport/Einar_backup/data/mcmillanhk/topnet_hiltom/input/topnet";

	//Load data from file
	const TTime timeFrom(boost::gregorian::date(2010, boost::gregorian::Sep, 28), 
         boost::posix_time::time_duration(0,0,0));
	const TTime timeTo(boost::gregorian::date(2010, boost::gregorian::Sep, 30), 
         boost::posix_time::time_duration(0,0,0));
	
	pModel->setRainfallMeasurements(strCombinedFilename, timeFrom, timeTo);
	pModel->setPETMeasurements(strCombinedFilename, timeFrom, timeTo);
	pModel->setOutflowMeasurements(strFlowDir, timeFrom, timeTo);
	
	pModel->completeSetup();

	//pModel->run();
	pModel->optimise();
	
	pModel->saveStates("output.tsv");
	
	//std::getchar();
	return 0;
}
