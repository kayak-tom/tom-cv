#ifndef CLMBLOCKMODELSUBCATCHMENT_H
#define CLMBLOCKMODELSUBCATCHMENT_H

#include "subcatchment_base.h" // Base class: CSubcatchment_base
#include <util/exception.h>

class CLMFullModel;
/*class CFlow
{
	double dFlow;
public:
	

}*/

class CLMBlockModelSubcatchment : public CSubcatchment_base
{
	friend class CLMFullModel;
/////////// Internal classes /////////////////////////////
	class COneTimestepState;
	
	class CDistributionFn
	{
	public:
		typedef std::map<TTimeDuration, double> TDistributeMap;
	private:
		TDistributeMap aDistributeMap;
	public:

		CDistributionFn(const TTimeDuration timeStep) //The number of timesteps water is distributed over should remain unchanged
		{
			double dSum = 0;
			TTimeDuration currentTimeStep;
			for(int i=0; i<10;i++) 
			{
				const double dVal = 1.0/(fabs(i-3) + 3);
				aDistributeMap[currentTimeStep] = dVal;
				dSum += dVal;
				currentTimeStep += timeStep;
			}
			
			currentTimeStep = TTimeDuration();
			for(int i=0; i<10;i++)
			{
				aDistributeMap[currentTimeStep] /= dSum;
				currentTimeStep += timeStep;
			}
		}

		const TDistributeMap & getRelDistributions() const { return aDistributeMap; }
	};
	
	//Not sure if this is one state or all states...
	/*class CRiverState
	{
		double dStoreVal; //value stored at this time
		CDistributionFn distnFn;
	public:
		CRiverState(const TTimeDuration timeStep) : dStoreVal(0), distnFn(timeStep)
		{
			
		}
	};*/
	
	
	/**
	 * @class CStreamNetworkDistribution. One per subcatchment.
	 * @brief 
	 */
	class CStreamNetworkDistribution
	{
		CDistributionFn distnFn;
	public:
		CStreamNetworkDistribution(const TTimeDuration timeStep) : distnFn(timeStep) {}
	
		void distributeWater(CLMBlockModelSubcatchment::COneTimestepState * pSubcatchmentState, const double dVolOfWater)
		{
			//now distribute over X future timesteps
			//Each river flow state for this subcatchment will be recomputed, and includes a ptr to this state
			//Distribute to a map of packets of water at the next catchment.
			
			for(CDistributionFn::TDistributeMap::const_iterator distribution = distnFn.getRelDistributions().begin(); distribution != distnFn.getRelDistributions().end(); distribution++)
			//BOOST_FOREACH(const typename CDistributionFn::TDistributeMap::value_type & distribution, distnFn.getRelDistributions())
			{
				CLMBlockModelSubcatchment::COneTimestepState * pForwardSubcatchment = pSubcatchmentState->getNeighbouringState(distribution->first);
				if(!pForwardSubcatchment)
					break;
				
				if(pSubcatchmentState->getTime() + distribution->first != pForwardSubcatchment->getTime())
					throw "Time stepping mismatch";
					
				pForwardSubcatchment->updateStreamWaterToRiverPacket(pSubcatchmentState->getTime(), distribution->second*dVolOfWater);
			}
			
			pSubcatchmentState->getParent()->setRiverFlowDirty();
		}
	};

	/**
	 * @class CRiverDistribution. One per subcatchment? flows should depend on flow volume
	 * @brief 
	 */
	class CRiverDistribution
	{
		CDistributionFn distnFn;
	public:
		CRiverDistribution(const TTimeDuration timeStep) : distnFn(timeStep) {}
		const CDistributionFn & getDistributionFn() const { return distnFn; }
	};	
	
	class COneTimestepState
	{
		const TTime time;
		double dRainfallEstimate;
		double dRainfallMeasurement, dPET;
		double dCanopyStore, dRootzoneStore, dAquiferStore;
		double dOutflowEstimate, dOutflowMeasurement;
		bool bStateDirty;
		
		COneTimestepState * pPreviousState, * pNextState; 
		CLMBlockModelSubcatchment * pParentSubcatchment;
		
	public:
		typedef std::map<TTime, double> TFlowContributions;
	private:
		TFlowContributions aStreamflowContributions,  //The stream outflow in this reach is the sum of flows out of this subcatchment's stream networks at earlier timesteps
						   aUpstreamReachOutflowContributions;
								
	public:
		COneTimestepState(const TTime & time, CLMBlockModelSubcatchment * pParentSubcatchment) : time(time), dRainfallEstimate(0), dRainfallMeasurement(0), dPET(0), dCanopyStore(0), dRootzoneStore(0), dAquiferStore(0), dOutflowEstimate(0), dOutflowMeasurement(0), bStateDirty(true), pPreviousState(0), pNextState(0), pParentSubcatchment(pParentSubcatchment) {}
		
		void setPET(const double dPET_in) { dPET=dPET_in; bStateDirty = true; } 
		void setRainfallMeasurement(const double dRain_in)
		{
			if(dRainfallEstimate != 0)
				throw "Rainfall estimate already initialised for this state";
			dRainfallMeasurement = dRain_in;
			dRainfallEstimate = dRainfallMeasurement;
		}
		
		void setRainfallEstimate(const double dRain_in, const bool bRecurseAndRecompute)
		{
			if(dRainfallMeasurement == 0)
				return;
			
			dRainfallEstimate = std::max<double>(dRain_in, 0);
			bStateDirty = true;
			
            recompute(bRecurseAndRecompute);
		}
        
        double getInitialRainfallEstimate()
        {
            return dRainfallMeasurement;
        }
        
		void setFlowMeasurement(const double dFlow_in)
		{
			if(dOutflowEstimate != 0)
				throw "Flow estimate already initialised for this state";
			dOutflowMeasurement = dFlow_in;
			dOutflowEstimate = dOutflowMeasurement;
            cout << "Set flow at time " << time << " = " << dOutflowMeasurement << endl;
		}
		
		void recompute(const bool bRecurse); //recursively update all descendent states = later timesteps. Flows are updated seperately
		
		double computeRainfallError() const
		{	
            if(bStateDirty)
                throw "Not recomputed";
                
			if(dRainfallMeasurement == 0)
				return 0;
				
			const double dRainfallSD = 0.5*dRainfallMeasurement;
			const double dNormalisedResid = (dRainfallEstimate - dRainfallMeasurement)/dRainfallSD;
            CHECKBADNUM(dNormalisedResid);
            return dNormalisedResid;
		}
		double computeFlowError() const
		{	
			if(pParentSubcatchment->bRiverFlowDirty)
				throw "River flows not recomputed";
			if(dOutflowMeasurement == 0)
				throw "Outflow not set for this state";//or return 0;
				
			const double dFlowSD = 0.1*dOutflowMeasurement;
			const double dNormalisedResid =  (dOutflowEstimate - dOutflowMeasurement)/dFlowSD;
            CHECKBADNUM(dNormalisedResid);
            return dNormalisedResid;
		}
		
		const TTime & getTime() const { return time; }
		
		void setPrevState(COneTimestepState * pPreviousState_in)
		{
			if(pPreviousState)
				throw "Already have previous state";
			pPreviousState = pPreviousState_in;
            
            pPreviousState->pNextState = this;
		}
		
		COneTimestepState * getNeighbouringState(const TTimeDuration & timeDiff) const 
		{
			return pParentSubcatchment->getState(time+timeDiff);
		}
		const COneTimestepState * getPrevState() const { return pPreviousState; }
		const COneTimestepState * getNextState() const { return pNextState; }
		COneTimestepState * getNextState() { return pNextState; }

		/////////Flow distribution from elsewhere: /////////////
		void updateStreamWaterToRiverPacket(const TTime sourceTime, const double dNewVolume)
		{
			aStreamflowContributions[sourceTime] = dNewVolume;
		}
		/////////Flow distribution from elsewhere: /////////////
		void updateOutflowContributionsFromUpstream(const TTime sourceTime, const double dNewVolume)
		{
			aUpstreamReachOutflowContributions[sourceTime] = dNewVolume;
		}
		
		CLMBlockModelSubcatchment * getParent() { return pParentSubcatchment; }
		
		//const TFlowContributions & getStreamflowContributions() const { return aStreamflowContributions; }  
		//const TFlowContributions & getRiverflowContributions() const { return aUpstreamReachOutflowContributions; }  
		const double sumStreamFlowContributions() const
		{
			double dSumContributions = 0;
			for(TFlowContributions::const_iterator contrib=aStreamflowContributions.begin(); contrib != aStreamflowContributions.end(); contrib++)
			//BOOST_FOREACH(const TFlowContributions::value_type & contrib, aStreamflowContributions)
			{
				dSumContributions += contrib->second;				
			}
			return dSumContributions;
		}
		const double sumUpstreamRiverFlowContributions() const
		{
			double dSumContributions = 0;
			for(TFlowContributions::const_iterator contrib=aStreamflowContributions.begin(); contrib != aStreamflowContributions.end(); contrib++)
			//BOOST_FOREACH(const TFlowContributions::value_type & contrib, aStreamflowContributions)
			{
				dSumContributions += contrib->second;				
			}
			return dSumContributions;
		}
		
		void recomputeRiverFlows();
		
		void saveState(std::ofstream & outputFile) const;
	};
	
///////////////////// Subcatchment state //////////////////////////

	CLMBlockModelSubcatchment * pDownstreamCatchment, * pUpstreamCatchment; //Links to other subcatchments
	
	typedef std::map<const TTime, COneTimestepState *> TStates; // Timesteps within this catchment
	TStates aStates;
	
	CStreamNetworkDistribution streamNetwork; //Stream network distribution for this catchment
	CRiverDistribution riverNetwork; //River outflow distribution for this catchment
	bool bRiverFlowDirty;
	
	COneTimestepState * pInitState; //This one is only for some initial conditions, and is never used

///////////////////// 


/////////////////////
	
	void setUpstreamCatchment(CLMBlockModelSubcatchment * pUpstreamCatchment_in);
	
	COneTimestepState * getState(const TTime & time) const 
	{
		TStates::const_iterator pState = aStates.find(time);
		if(pState == aStates.end())
			return 0;
		return pState->second;
	}
	
	void setupState(const TTime& timeFrom, const TTime& timeTo);
	
public:
	CLMBlockModelSubcatchment(const CSubcatchmentParams* pSubcatchmentParams, CLMBlockModelSubcatchment * pDownstreamCatchment, const TTimeDuration timeStep);
	virtual ~CLMBlockModelSubcatchment();

	CLMBlockModelSubcatchment * getUpstreamCatchment() { return pUpstreamCatchment; }
	const CLMBlockModelSubcatchment * getUpstreamCatchment() const { return pUpstreamCatchment; }
	CLMBlockModelSubcatchment * getDownstreamCatchment() { return pDownstreamCatchment; }
	const CLMBlockModelSubcatchment * getDownstreamCatchment() const { return pDownstreamCatchment; }

	void setRiverFlowDirty() { bRiverFlowDirty = true; if(pDownstreamCatchment) pDownstreamCatchment->setRiverFlowDirty(); }
	
	void completeSetup();
	void recomputeAll();
	void recomputeAllRiverFlows();
	
	const CRiverDistribution & getRiverNetworkDistn() const { return riverNetwork; } //River outflow distribution for this catchment

	
	virtual void setOutflowMeasurements(const CInputArray& flow, const TTime& timeFrom, const TTime& timeTo);
	virtual void setPETMeasurements(const CInputArray& PET, const TTime& timeFrom, const TTime& timeTo);
	virtual void setRainfallMeasurements(const CInputArray& rainfall, const TTime& timeFrom, const TTime& timeTo);
	
	void saveStates(std::ofstream & outputFile) const;
};

#endif // CLMBLOCKMODELSUBCATCHMENT_H
