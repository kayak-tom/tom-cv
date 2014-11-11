#include "lm_block_model_subcatchment.h"
#include "flow_measurements.h"

CLMBlockModelSubcatchment::CLMBlockModelSubcatchment(const CSubcatchmentParams * pSubcatchmentParams, CLMBlockModelSubcatchment * pDownstreamCatchment, const TTimeDuration timeStep) : CSubcatchment_base(pSubcatchmentParams, timeStep), pDownstreamCatchment(pDownstreamCatchment), pUpstreamCatchment(0), streamNetwork(timeStep), riverNetwork(timeStep), bRiverFlowDirty(true), pInitState(0)
{
	if(pDownstreamCatchment)
		pDownstreamCatchment->setUpstreamCatchment(this);
}

CLMBlockModelSubcatchment::~CLMBlockModelSubcatchment()
{
}

void CLMBlockModelSubcatchment::completeSetup()
{
	pInitState = new COneTimestepState(aStates.begin()->first - timeStep, this);
	aStates.begin()->second->setPrevState(pInitState);
	
	for(TStates::const_iterator state = aStates.begin(); state != aStates.end(); state++)
	{
		cout << state->first << endl;
		cout << state->second->getTime() << endl;
		if(!state->second->getPrevState())
		{
			throw "No previous state set";
		}			
	}
}
	
void CLMBlockModelSubcatchment::setUpstreamCatchment(CLMBlockModelSubcatchment * pUpstreamCatchment_in)
{
	pUpstreamCatchment = pUpstreamCatchment_in;
}

void CLMBlockModelSubcatchment::setupState(const TTime& timeFrom, const TTime& timeTo)
{
	for(TTime time = timeFrom; time <= timeTo; time += timeStep)
	{
		 TStates::iterator pTime = aStates.find(time);
		 if(pTime == aStates.end())
		 {
			 aStates[time] = new COneTimestepState(time, this);
			 cout << "Added state at time " << time << endl;

			 TStates::iterator pLastTime = aStates.find(time - timeStep);
			 
			 if(pLastTime != aStates.end()) //last state one exists
			 {
				 aStates[time]->setPrevState(pLastTime->second);
			 }
		 
			 TStates::iterator pNextTime = aStates.find(time + timeStep);
			 
			 if(pNextTime != aStates.end()) //last state one exists
			 {
				 pNextTime->second->setPrevState(aStates[time]);
			 }
		 }
	}
}

void CLMBlockModelSubcatchment::setOutflowMeasurements(const CInputArray& flow, const TTime& timeFrom, const TTime& timeTo)
{
	setupState(timeFrom, timeTo);
    cout << "Time from=" << timeFrom << endl;
	
	for(CInputArray::TMeasurements::const_iterator measurement = flow.getMeasurements().begin(); measurement != flow.getMeasurements().end(); measurement++)
	//BOOST_FOREACH(CInputArray::TMeasurements::value_type & measurement, flow.getMeasurements())
	{
		if(measurement->first >= timeFrom && measurement->first <= timeTo)
		{
			 TStates::iterator pTime = aStates.find(measurement->first);
			 if(pTime == aStates.end())
			 {
				 cout << measurement->first << endl;
				 throw "Timestep not setup";
			 }
			 pTime->second->setFlowMeasurement(measurement->second);
		}
	}
    
    cout << "Time to=" << timeTo << endl;

}

void CLMBlockModelSubcatchment::recomputeAll()
{
    pInitState->recompute(true); // and recurse

	/*for(TStates::const_iterator state = aStates.begin(); state != aStates.end(); state++)
	//BOOST_FOREACH(const typename TStates::value_type & state, aStates)
	{
		state->second->recompute(false);
	}*/
}

void CLMBlockModelSubcatchment::recomputeAllRiverFlows()
{
    if(!bRiverFlowDirty)
        return;
        
    if(pUpstreamCatchment && pUpstreamCatchment->bRiverFlowDirty)
        throw "Upstream catchment not recomputed";
    
	for(TStates::const_iterator state = aStates.begin(); state != aStates.end(); state++)
	//BOOST_FOREACH(const typename TStates::value_type & state, aStates)
	{
		state->second->recomputeRiverFlows();
	}
    
    bRiverFlowDirty = false;
}

void CLMBlockModelSubcatchment::setPETMeasurements(const CInputArray& PET, const TTime& timeFrom, const TTime& timeTo)
{
	setupState(timeFrom, timeTo);
	
	for(CInputArray::TMeasurements::const_iterator measurement = PET.getMeasurements().begin(); measurement != PET.getMeasurements().end(); measurement++)
	//BOOST_FOREACH(const CInputArray::TMeasurements::value_type & measurement, PET.getMeasurements())
	{
		if(measurement->first >= timeFrom && measurement->first <= timeTo)
			aStates[measurement->first]->setPET(measurement->second);
	}
	
}

void CLMBlockModelSubcatchment::setRainfallMeasurements(const CInputArray& rainfall, const TTime& timeFrom, const TTime& timeTo)
{
	setupState(timeFrom, timeTo);
	
	for(CInputArray::TMeasurements::const_iterator measurement = rainfall.getMeasurements().begin(); measurement != rainfall.getMeasurements().end(); measurement++)
	//BOOST_FOREACH(const CInputArray::TMeasurements::value_type & measurement, rainfall.getMeasurements())
	{
		if(measurement->first >= timeFrom && measurement->first <= timeTo)
		{
			if(aStates.find(measurement->first) == aStates.end())
			{
				cout << measurement->first << endl;
				throw "Timestep not set up";
			}
			aStates[measurement->first]->setRainfallMeasurement(measurement->second);
		}
	}
	
}

double CLMBlockModelSubcatchment::(const double S_c)
{
    
}

//A2c fraction of rainfall that throughfalls
double CLMBlockModelSubcatchment::f(const double S_canopy)
{
    const double C_canopy = pSubcatchmentParams->getParam()->getC_canopy();
    
    return (S_canopy/C_canopy)*(2 - S_canopy/C_canopy);
}

void CLMBlockModelSubcatchment::COneTimestepState::recompute(const bool bRecurse) //recursively update all descendent states = later timesteps and flows
{
	//should be all time dependent
    if(pPreviousState) //Not the very first state...
    {
        if(pPreviousState->bStateDirty)
            throw "previous timestep hasn't been recomputed";
            
        const double fS_canopy = pParentSubcatchment->f(pParentSubcatchment->S_c(pPreviousState->dCanopyStore));
        const double p_throughfall = dRainfallEstimate * fS_canopy; //0.2*pPreviousState->dCanopyStore;
        
        const double c_r = pParentSubcatchment->pSubcatchmentParams->getC_rfactor();
        const double e_canopy = dPET*c_r*fS_canopy;
        const double dS_canopy_dt = dRainfallEstimate - p_throughfall - e_canopy; 
        dCanopyStore = pPreviousState->dCanopyStore + dS_canopy_dt;//Check max?
        
        const double K_0 = pParentSubcatchment->pSubcatchmentParams->getHydraulicCon();
        const double m = pParentSubcatchment->pSubcatchmentParams->getTopnet_m();
        const double dMaxInfiltration = todo 
        const double dInfiltration = std::min<double>(dMaxInfiltration, p_throughfall); //This is the rain that actually soaks in
        
        const double dExcessToOverlandFlow = p_throughfall - dInfiltration;
        
        const double drainageToAquifer = dInfiltration*pPreviousState->dRootzoneStore; 
        const double e_rootzone = 0.1*dPET*pPreviousState->dRootzoneStore;
        const double dS_rootzone_dt = dInfiltration - drainageToAquifer - e_rootzone;
        dRootzoneStore = pPreviousState->dRootzoneStore + dS_rootzone_dt;
        
        const double q_baseflow = 0.01*pPreviousState->dAquiferStore;
        const double dS_aquifer_dt = drainageToAquifer - q_baseflow;
        dAquiferStore = pPreviousState->dAquiferStore + dS_aquifer_dt;
        
        /////
        //Streamflow
        const double dInfiltrationExcess = p_throughfall-dInfiltration;
        const double dS_streamflow_dt = q_baseflow + dInfiltrationExcess;
        
        pParentSubcatchment->streamNetwork.distributeWater(this, dS_streamflow_dt);
    }
	
	bStateDirty = false;
    
    pParentSubcatchment->setRiverFlowDirty();
    
    if(bRecurse && pNextState)
        pNextState->recompute(bRecurse);
}

void CLMBlockModelSubcatchment::COneTimestepState::recomputeRiverFlows()
{
    if(bStateDirty)
        throw "Catchment state hasn't been updated yet";
	//Sum over upstream contributions [and upstream flows?]
	dOutflowEstimate = sumStreamFlowContributions() + sumUpstreamRiverFlowContributions();
			
	CLMBlockModelSubcatchment * pDownstream = getParent()->getDownstreamCatchment();

	if(!pDownstream)
		return;
	
	//now distribute over X future timesteps
	//Each river flow state for this subcatchment will be recomputed, and includes a ptr to this state
	//Distribute to a map of packets of water at the next catchment.
	CLMBlockModelSubcatchment::COneTimestepState * pDownstreamSubcatchment_sameTime = pDownstream->getState(getTime());
	
	for(CDistributionFn::TDistributeMap::const_iterator distribution = getParent()->getRiverNetworkDistn().getDistributionFn().getRelDistributions().begin(); distribution != getParent()->getRiverNetworkDistn().getDistributionFn().getRelDistributions().end(); distribution++)
	//BOOST_FOREACH(CDistributionFn::TDistributeMap::value_type & distribution, distnFn.getRelDistributions())
	{
		CLMBlockModelSubcatchment::COneTimestepState * pDownstreamState = pDownstreamSubcatchment_sameTime->getNeighbouringState(distribution->first);
		if(!pDownstreamState)
			break;

		if(getTime() + distribution->first != pDownstreamState->getTime())
			throw "Time stepping mismatch";
			
		pDownstreamState->updateOutflowContributionsFromUpstream(getTime(), distribution->second*dOutflowEstimate);
	}
}

void CLMBlockModelSubcatchment::saveStates(std::ofstream & outputFile) const
{
	for(TStates::const_iterator state = aStates.begin(); state != aStates.end(); state++)
	//BOOST_FOREACH(const typename TStates::value_type & state, aStates)
	{
		state->second->saveState(outputFile);
	}
}

void CLMBlockModelSubcatchment::COneTimestepState::saveState(std::ofstream & outputFile) const
{
	outputFile << getTime() << '\t' << dRainfallMeasurement << '\t' << dRainfallEstimate << '\t' << dOutflowMeasurement << '\t' << dOutflowEstimate << endl;
}
