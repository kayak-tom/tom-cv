#ifndef CLMBLOCKMODEL_H
#define CLMBLOCKMODEL_H

#include "model_base.h" // Base class: CModel_base

class CLMBlockModelSubcatchment;
class CLMFullModel;

class CLMBlockModel : public CModel_base
{
	friend class CLMFullModel;
public:
	CLMBlockModel(const CCatchmentSetupParams * pCatchmentStructure, const TTimeDuration & timeStep);
	virtual ~CLMBlockModel();

	virtual void run();
	virtual void optimise();
	virtual void saveStates(const std::string strTSVFilename) const;
	virtual CSubcatchment_base* setupSubcatchment(const CSubcatchmentParams* pSubcatchmentParams, CSubcatchment_base* pDownstreamCatchment);
	virtual void completeSetup();
	
protected:
	void recomputeAll();
	void recomputeRiverFlows_UpstreamDown();
    
	CLMBlockModelSubcatchment * getSubcatchment(int nSubcatchment);
	const CLMBlockModelSubcatchment * getSubcatchment(int nSubcatchment) const;

};

#endif // CLMBLOCKMODEL_H
