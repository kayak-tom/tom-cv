#include <set>
#include "lm_block_model.h"
#include "lm_block_model_subcatchment.h"
#include <geom/levMarNumerical.h>

CLMBlockModel::CLMBlockModel(const CCatchmentSetupParams * pCatchmentStructure, const TTimeDuration & timeStep) : CModel_base(pCatchmentStructure, timeStep)
{
}

CLMBlockModel::~CLMBlockModel()
{
}

CLMBlockModelSubcatchment * CLMBlockModel::getSubcatchment(int nSubcatchment)
{
	return dynamic_cast<CLMBlockModelSubcatchment *>(aAllSubcatchments[nSubcatchment]);
}
const CLMBlockModelSubcatchment * CLMBlockModel::getSubcatchment(int nSubcatchment) const
{
	return dynamic_cast<const CLMBlockModelSubcatchment *>(aAllSubcatchments[nSubcatchment]);
}
void CLMBlockModel::recomputeAll()
{
	//for each subcatchment
	for(int i=0; i<(int)aAllSubcatchments.size(); i++)
		getSubcatchment(i)->recomputeAll();
    recomputeRiverFlows_UpstreamDown();
}

void CLMBlockModel::recomputeRiverFlows_UpstreamDown()
{
	//Now need to recompute river flows for each subcatchment's reach, in order, from upstream to downstream
	std::set<CLMBlockModelSubcatchment *> aComputed;
	do
	{
		for(int i=0; i<(int)aAllSubcatchments.size(); i++)
		{
			//If all upstream catchments are computed then can compute this one
			if(!getSubcatchment(i)->getUpstreamCatchment() || aComputed.find(getSubcatchment(i)->getUpstreamCatchment()) != aComputed.end())
			{
                getSubcatchment(i)->recomputeAllRiverFlows();
				
				aComputed.insert(getSubcatchment(i));
			}
		}
	} while (aComputed.size() < aAllSubcatchments.size());
}

void CLMBlockModel::run()
{
	recomputeAll();

}

class CLMFullModel : public CLMFunction 
{
	CLMBlockModel * pBlockModel;
	int nInputs, nRainfallValues, nOutflowValues;
	std::map<int, CLMBlockModelSubcatchment::COneTimestepState *> aIdxToStates;
public:
    virtual int inputs() const { return nInputs; }
    virtual int values() const { return nRainfallValues+nOutflowValues; }

    //If only one parameter has changed since this residual vector was calculated, nParamChanged is set to that parameter index (the function only needs to update relevent residuals). Otherwise nParamChanges = -1 
    virtual eLMSuccessStatus function(const Eigen::VectorXd &x, Eigen::VectorXd &resids, bool bVerbose, const int nParamChanged)
	{
		//First update:
		if(nParamChanged == -1)
		{
			for(int i=0; i<nInputs; i++)
			{
				aIdxToStates[i]->setRainfallEstimate(x(i), false); //Just recomputes this state. They are ordered *within each catchment*, so all later states will be updated by the loop
			}
		}
		else			
		{
			if(nParamChanged > 0)
				aIdxToStates[nParamChanged-1]->setRainfallEstimate(x(nParamChanged-1), false);//Just recomputes this state.
			aIdxToStates[nParamChanged]->setRainfallEstimate(x(nParamChanged), true);// recompute all future states.
		}
        
		//Now recompute river flows:
        pBlockModel->recomputeRiverFlows_UpstreamDown();
        
		//Next compute residuals: [TODO: only have to recompute those which have changed-- 1 or 2 subcatchments]
        int nResidRainfall = 0;
		for(int i=0; i<(int)pBlockModel->aAllSubcatchments.size(); i++)
		{
            const CLMBlockModelSubcatchment::TStates & aStates = pBlockModel->getSubcatchment(i)->aStates;
			for(CLMBlockModelSubcatchment::TStates::const_iterator state = aStates.begin(); state != aStates.end(); state++)
            {
                resids(nResidRainfall++) = state->second->computeRainfallError();
                //resids(nResidFlow++) = state->second->computeFlowError();
            }
        }

        int nResidFlow = nResidRainfall;
        
        //Downstream catchment
        const CLMBlockModelSubcatchment::TStates & aStates = dynamic_cast<const CLMBlockModelSubcatchment *>(pBlockModel->getDownstreamCatchment())->aStates;
        for(CLMBlockModelSubcatchment::TStates::const_iterator state = aStates.begin(); state != aStates.end(); state++)
        {
            resids(nResidFlow++) = state->second->computeFlowError();
        }
        
		if(nResidFlow != values())
            throw "Resid counting error";
            
        return eLMSuccess;
	}

    virtual Eigen::VectorXd init() 
    {
        Eigen::VectorXd initParams = Eigen::VectorXd::Zero(inputs());
        
        for(int i=0; i<nInputs; i++)
        {
            initParams(i) = aIdxToStates[i]->getInitialRainfallEstimate();
        }
        
        return initParams; 
    }
	CLMFullModel(CLMBlockModel * pBlockModel) : pBlockModel(pBlockModel), nInputs(0), nRainfallValues(0), nOutflowValues(0)
	{
        int nParam = 0;
		for(int i=0; i<(int)pBlockModel->aAllSubcatchments.size(); i++)
		{
            CLMBlockModelSubcatchment::TStates & aStates = pBlockModel->getSubcatchment(i)->aStates;
   			const int nNumStates = aStates.size();
			nInputs += nNumStates;
			nRainfallValues += nNumStates;

			for(CLMBlockModelSubcatchment::TStates::const_iterator state = aStates.begin(); state != aStates.end(); state++)
            {
                aIdxToStates[nParam++]= state->second; //Order is important: upstream to downstream
            }
		}
        nOutflowValues += (int)dynamic_cast<const CLMBlockModelSubcatchment *>(pBlockModel->getDownstreamCatchment())->aStates.size();
        
		if(nParam != inputs())
            throw "Param counting error";
	}
	
};

void CLMBlockModel::optimise()
{
    this->run(); //Full recompute of all the 0-states, which will never be recomputed again.
	CLMFullModel fullModelFn(this);
    CLevMar LM(fullModelFn, true, 1e-7);
    Eigen::VectorXd params = fullModelFn.init();
    cout << params.transpose() << endl;
    LM.minimise(params, 5);
    cout << "Optimisation complete" << endl;
}

void CLMBlockModel::completeSetup()
{
	for(int i=0; i<(int)aAllSubcatchments.size(); i++)
	{
		getSubcatchment(i)->completeSetup();
	}
}

void CLMBlockModel::saveStates(const std::string strTSVFilename) const
{
	std::ofstream outputFile(strTSVFilename.c_str());
	for(int i=0; i<(int)aAllSubcatchments.size(); i++)
	{
		if(!getSubcatchment(i)->getDownstreamCatchment())
			getSubcatchment(i)->saveStates(outputFile);
	}
}

CSubcatchment_base* CLMBlockModel::setupSubcatchment(const CSubcatchmentParams* pSubcatchmentParams, CSubcatchment_base* pDownstreamCatchment)
{
	CLMBlockModelSubcatchment * pSubcatchment = new CLMBlockModelSubcatchment(pSubcatchmentParams, dynamic_cast<CLMBlockModelSubcatchment *>(pDownstreamCatchment), timeStep);
	return pSubcatchment;
}
