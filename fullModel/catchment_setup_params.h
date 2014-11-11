#ifndef CCATCHMENTSETUPPARAMS_H
#define CCATCHMENTSETUPPARAMS_H

#include "main.h"

class NcFile;
class CSubcatchmentParams;

class CCatchmentSetupParams
{
	NcFile* pCatchmentDescriptionFile;
	
	typedef std::map<int, CSubcatchmentParams *> TSubcatchmentParams;
	TSubcatchmentParams aSubcatchments;
public:
	CCatchmentSetupParams(const std::string & szFilenameForCatchment);
	virtual ~CCatchmentSetupParams();

	double getParam(const char * szParamName, int nrch) const;
	
	const CSubcatchmentParams * getDownstreamCatchment() const;
};

#endif // CCATCHMENTSETUPPARAMS_H
