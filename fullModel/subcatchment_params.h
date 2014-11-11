#ifndef CSUBCATCHMENTPARAMS_H
#define CSUBCATCHMENTPARAMS_H

#include <vector>
#include <string>

class CCatchmentSetupParams;

class CSubcatchmentParams
{
	const int nrch, rchid /* index */;
    const CCatchmentSetupParams * pParent;
	std::vector<CSubcatchmentParams *> apUpstreamCatchments; 
	CSubcatchmentParams * pDownstreamCatchment;
	void addUpstreamCatchment(CSubcatchmentParams * pUpstreamCatchment_in);
public:
	CSubcatchmentParams(const int nrch, const int rchid, const CCatchmentSetupParams * pParent);
	virtual ~CSubcatchmentParams();
	
	void setDownstreamCatchment(CSubcatchmentParams * pDownstreamCatchment_in);

	double getParam(const char * szParamName) const;
	
	const CSubcatchmentParams * getDownstreamCatchment() const { return pDownstreamCatchment; }
	const std::vector<CSubcatchmentParams *> getUpstreamCatchments() const { return apUpstreamCatchments; }
	int getRchid() const { return rchid; }
	int getNrch() const { return nrch; }
    double getC_canopy() const;
    double getC_rfactor() const;
    double getHydraulicCon() const;
    double getTopnet_m() const;
};

#endif // CSUBCATCHMENTPARAMS_H
