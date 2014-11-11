#include "catchment_setup_params.h"
#include "netcdfcpp.h"
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include "subcatchment_params.h"



CCatchmentSetupParams::CCatchmentSetupParams(const std::string & szFilenameForCatchment)
{
	if(!boost::filesystem::exists(szFilenameForCatchment))
		throw "szFilenameForCatchment doesn't exist";
	
	pCatchmentDescriptionFile = new NcFile(szFilenameForCatchment.c_str());
	
	if (!pCatchmentDescriptionFile->is_valid())
    {
        throw "Couldn't open file!";
    }
	
	NcVar * pReachIDs = pCatchmentDescriptionFile->get_var("rchid");
	NcDim* numReaches = pReachIDs->get_dim(0);
	const int nNumReaches = numReaches->size();
	
	cout << "numReaches=" << nNumReaches << endl;
	
	//std::vector<int> aReachIDs(nNumReaches, -1);
	int anReaches[nNumReaches];
	pReachIDs->get(anReaches, nNumReaches); //This is the mapping from rchid to nrch (=idx) 
	
	for(int nrch=0;nrch<nNumReaches;nrch++)
	{
		std::cout << "Reach: " << anReaches[nrch] << std::endl;
		aSubcatchments[nrch] = new CSubcatchmentParams(nrch, anReaches[nrch], this);
	}
	
	int anDownstreamReaches[nNumReaches];
	pCatchmentDescriptionFile->get_var("dsrch_nrch")->get(anDownstreamReaches, nNumReaches); 
	for(int nrch=0;nrch<nNumReaches;nrch++)
	{
		if(anDownstreamReaches[nrch]>=0)
		{
			std::cout << "Catchment " << nrch << " Downstream reach: " << anDownstreamReaches[nrch] << std::endl;
			//aSubcatchments[nrch] = new CSubcatchmentParams(nrch, this);
			aSubcatchments[nrch]->setDownstreamCatchment(aSubcatchments[anDownstreamReaches[nrch]]);
		}
	}
	
	std::cout << std::endl;
	
}

CCatchmentSetupParams::~CCatchmentSetupParams()
{
	delete pCatchmentDescriptionFile;
}

double CCatchmentSetupParams::getParam(const char * szParamName, int nrch) const
{
	NcVar * pParam = pCatchmentDescriptionFile->get_var(szParamName);
	return pParam->as_double(nrch);
}

const CSubcatchmentParams * CCatchmentSetupParams::getDownstreamCatchment() const
{
	const CSubcatchmentParams * pSubcatchment = aSubcatchments.begin()->second;
	while(pSubcatchment->getDownstreamCatchment())
	{
		pSubcatchment = pSubcatchment->getDownstreamCatchment();
	};
	return pSubcatchment;
}