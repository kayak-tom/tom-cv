/*Each subcatchment is a group of vertices and edges at one timestep
 * */
#ifndef _SUBCATCHMENT_H
#define _SUBCATCHMENT_H

#include "flowsAndStores.h"
#include <boost/optional.hpp>
using boost::optional;

/* A flow out of a store, and the corresponding drainage function/constraint connecting it with its source store
 */
class CFlowBetweenStoresElement : public CFlowVertex
{
    COutFlowConstraint * pOutFlowConstraint;
public:
    CFlowBetweenStoresElement(CStoreElement * pStoreSource, CStoreElement * pStoreDest, const CDrainageFunction * pDrainage) : CFlowVertex(pStoreSource->time()) {
        pOutFlowConstraint = new COutFlowConstraint(pStoreSource, this, pDrainage);
    }
};

/* A flow out of a catchment (ET or final outflow), and the corresponding drainage function/constraint connecting it with the store it is from
 */
class CFlowSinkElement : public CFlowVertex
{
    COutFlowConstraint * pOutFlowConstraint=0;
public:
    CFlowSinkElement(CStoreElement * pStoreSource, const CDrainageFunction * pDrainage) : CFlowVertex(pStoreSource->time()) {
        if(pDrainage)
            pOutFlowConstraint = new COutFlowConstraint(pStoreSource, this, pDrainage);
        else
            pOutFlowConstraint = 0;
    }
};

//need one of these from the canopy, one from the soil. Don't own constraints (so relatively simple). Do we actually need?
class CPETDependentETElement : public CFlowVertex
{
    //COutFlowTemperatureConstraint * pOutFlowConstraint;
public:
    CPETDependentETElement(const int nTime /*CStoreElement * pStoreSource, CTemperatureVertex * pPETEstimate, const CPETDependantCanopyDrainageFunction_base * pDrainage*/) : CFlowVertex(nTime) {
        //pOutFlowConstraint = new COutFlowTemperatureConstraint_TODONEW(pStoreSource, pTempEstimate, pHumidityEstimate, this, pDrainage);

        if(bVeryVerbose) cout << "T" << nTime << ": Added ET flow edge" << endl;
    }
};

class CStore;

//Parent of input and inter-store flows. Owns all of the flows for this particular catchment element
class CFlowBase
{
protected:
    vector<CFlowVertex *> apFlowElements;

public:
    const string label;

    CFlowVertex * getFlow(int nTime) {
        return apFlowElements[nTime];
    }
    const CFlowVertex * getFlow(int nTime) const {
        return apFlowElements[nTime];
    }
    
    virtual ~CFlowBase() {}
    
    virtual void addHeaders(std::ostream & of) const
    {
        of << label << '\t';
    }
    
    /* outflow */
    virtual void addData(const int nTime, std::ostream & of) const
    {
        of << getFlow(nTime)->deltaV() << '\t';
    }
    
    CFlowBase(const string label) : label(label) {}
    virtual void addFlowVertex(const int nTime) = 0;
};

/* Flow out of a store (out of catchment, or ET)
 * 
 * Essentially just subtracted from the store, doesn't go to any sink.
 *
 * */
class CSinkFlow : public CFlowBase
{
protected:
    CStore * pSource;
    const CDrainageFunction * pDrainageFunction;
public:
    CSinkFlow(const string label, CStore * pSource, CDrainageFunction *pDrainageFunction);

    virtual void addFlowVertex(const int nTime);
};

/* Flow between 2 stores, including outflow constraints
 *
 * */
class CFlow : public CFlowBase
{
    CStore * pSource;
    CStore * pDest;
    const CDrainageFunction * pDrainageFunction;
public:
    CFlow(const string label, CStore * pSource, CStore * pDest, CDrainageFunction *pDrainageFunction);

    virtual void addFlowVertex(const int nTime);
};

/* Flow out of a store (a sink), temperature dependent. Over time.
 * 
 * 
class CTemperatureDependentOutFlow : public CFlowBase
{
    CStore * pSource;
    const CTempDependantDrainageFunction *pDrainageFunction;

    // CFlowBase::apFlowElements contains the estimates of the actual flow out by ET
    //TODO: merge these 4 vectors into one timestep struct.
    vector<CTemperatureMeasurement *> apTemperatureElements, apHumidityElements; //Actual temperature measurements
    vector<CTemperatureVertex *> apTemperatureEstimates, apHumidityEstimates; //Estimated temperature measurements
public:
    CTemperatureDependentOutFlow(const string label, CStore * pSource, CTempDependantDrainageFunction *pDrainageFunction);

    virtual void addFlowVertex(const int nTime); //Called by the store (e.g. surface water), actually adds the 3ary constraint on temp,humidity,store,flowOut
    
    CTemperatureVertex * getTemperatureElement(const int nTime)
    {
        return apTemperatureEstimates[nTime];
    }
    CTemperatureVertex * getHumidityElement(const int nTime) //TODO Estimate vs element naming
    {
        return apHumidityEstimates[nTime];
    }
    
    void addTemperatureMeasurement(const int nTime, const double dTemperatureMeasurement, const double dHumidityMeasurement) { 
        / * vertex with temp estimate * /CTemperatureVertex* pTemperatureEstimateOneTimestep = new CTemperatureVertex(nTime);
        apTemperatureEstimates.push_back(pTemperatureEstimateOneTimestep);
        
        const double dTemperatureSD = 0.5, dHumiditySD=0.5;
        //MUST have a measurement for every estimate, otherwise MLE will be unconstrained (could add forward-backward constraints instead?). Each CTemperatureMeasurement constrains a temperature estimate
        / * unary edge constraining temp estimate * / CTemperatureMeasurement * pTemperatureMeasurement = new CTemperatureMeasurement(pTemperatureEstimateOneTimestep, dTemperatureMeasurement, sqr(dTemperatureSD));
        apTemperatureElements.push_back(pTemperatureMeasurement);

        / * vertex with humidity estimate * /CTemperatureVertex* pHumidityEstimateOneTimestep = new CTemperatureVertex(nTime);
        apHumidityEstimates.push_back(pHumidityEstimateOneTimestep);
        / * unary edge constraining humidity estimate * / CTemperatureMeasurement * pHumidityMeasurement = new CTemperatureMeasurement(pHumidityEstimateOneTimestep, dHumidityMeasurement, sqr(dHumiditySD));
        apHumidityElements.push_back(pHumidityMeasurement);
    }

    virtual void addHeaders(std::ostream & of) const
    {
        of << label << "_TempMeasured" << '\t' << label << "_TempMLE" << '\t'; //TODO: humidity
        CFlowBase::addHeaders(of);
    }
    
    / * outflow * /
    virtual void addData(const int nTime, std::ostream & of) const
    {
        of << apTemperatureElements[nTime]->measurement() << '\t' << apTemperatureEstimates[nTime]->T() << '\t'; //TODO: humidity
        CFlowBase::addData(nTime, of);
    }    
};*/

/* Flow out of the canopy and soil, which is dependent on PET from TOPNET. Over time.
 * 
 * 
 * */
class CPETDependentETOutFlow// : public CFlowBase
{
    CStore * pSoilSource, * pCanopySource;
    CSinkFlow * pSoilETFlow, * pCanopyETFlow;
    const CComplexPETDrainageFunction *pETFunction;

    // CFlowBase::apFlowElements contains the estimates of the actual flow out by ET
    vector<CPETMeasurement *> apPETElements; //Actual PET estimates from TOPNET
    vector<CPETVertex *> apPETEstimates; //Estimated PET 
public:
    CPETDependentETOutFlow(const string label, CStore * pSoilSource, CStore * pCanopySource, CComplexPETDrainageFunction *pETFunction);

    //virtual void addETFlowVertex(const int nTime); todo call only once //Called by the store (e.g. surface water), actually adds the 3ary constraint on temp,humidity,store,flowOut
    
    CPETVertex * getPETElement(const int nTime)
    {
        return apPETEstimates[nTime];
    }
    
    void addPETEstimate(const int nTime, const double dPETEstimate) { 
        /* vertex with PET estimate */CPETVertex* pPETEstimateOneTimestep = new CPETVertex(nTime);
        apPETEstimates.push_back(pPETEstimateOneTimestep);
        
        const double dPETSD = 0.01 * dPETEstimate; //very low
        //MUST have a measurement for every estimate, otherwise MLE will be unconstrained (could add forward-backward constraints instead?). Each CTemperatureMeasurement constrains a temperature estimate
        /* unary edge constraining temp estimate */ CPETMeasurement * pPETMeasurement = new CPETMeasurement(pPETEstimateOneTimestep, dPETEstimate, sqr(dPETSD));
        apPETElements.push_back(pPETMeasurement);

    }

    virtual void addHeaders(std::ostream & of) const
    {
        pCanopyETFlow->addHeaders(of);
        pSoilETFlow->addHeaders(of);
        /*of << label << "_PET_fromTOPNET" << '\t' << label << "_PET_MLE" << '\t'; //TODO: humidity
        CFlowBase::addHeaders(of);*/
    }
    
    /* outflow */
    virtual void addData(const int nTime, std::ostream & of) const
    {
        pCanopyETFlow->addData(nTime, of);
        pSoilETFlow->addData(nTime, of);
        //of << apPETElements[nTime]->measurement() << '\t' << apPETEstimates[nTime]->T() << '\t'; //TODO: humidity
        //CFlowBase::addData(nTime, of);
    }   
    void addETComplexFlowVertex(const int nTime);
 
};

/*Represents one store over time,
 */
class CStore
{
public:
    const string label;
private:
    vector<CStoreElement *> apStoreElements;
    vector<CTotalFlowConstraint *> apFlowConservationConstraints;
    CStoreInitialisingAssumption * pInitAssumption;
    //CStore * pDownstreamStore;
    vector<CFlowBase *> apInflows, //rainfall or other stores
                        apOutflows; //To other stores. Last timestep has no outflow
public:
    CStoreElement * getStoreElement(int nTime)
    {
        return apStoreElements[nTime];
    }
    const CStoreElement * getStoreElement(int nTime) const
    {
        return apStoreElements[nTime];
    }
    CTotalFlowConstraint * getConservationConstraint(int nTime)
    {
        return apFlowConservationConstraints[nTime];
    }
    

    void addOutflow(CFlowBase * pOutFlow) {
        if(bVeryVerbose) cout << "Adding outflow to " << label << ": " << pOutFlow->label << endl;
        apOutflows.push_back(pOutFlow);
    }

    void addInflow(CFlowBase * pInFlow) {
        apInflows.push_back(pInFlow);
    }

    CStore(const string label) : label(label), pInitAssumption(0) {
        CStoreElement * pFirstStoreElement = new CStoreElement(0);
        apStoreElements.push_back(pFirstStoreElement);

        if(bVeryVerbose) cout << "Creating first store with weak prior on initial storage" << endl;

        pInitAssumption = new CStoreInitialisingAssumption(pFirstStoreElement);
    }

    //Set up the catchment *before* adding timesteps
    void addTimestep(const int nTime) { //for now, times must start from 0 and be integer steps
        CStoreElement * pNextStoreElement = new CStoreElement(nTime+1); //Actually adds a t+1 timestep, for this time to flow into, plus the constraint on time t->t+1
        apStoreElements.push_back(pNextStoreElement);
    }

    //Adds all the flows out of this store (all the nTime+1 stores should be created now)
    void addOutflows(const int nTime) {
        BOOST_FOREACH(CFlowBase * pOutflow, apOutflows) {
            pOutflow->addFlowVertex(nTime);
        }
    }

    //Adds the total flow constraint (all the nTime -> nTime+1 flows should be created now)
    void addTotalFlowConstraint(const int nTime) {
        CStoreElement * pThisStoreElement = apStoreElements[nTime];
        CStoreElement * pNextStoreElement = apStoreElements[nTime+1];
        CTotalFlowConstraint * pConserveFlow = new CTotalFlowConstraint(pThisStoreElement, pNextStoreElement);

        apFlowConservationConstraints.push_back(pConserveFlow);

        if(nTime > 0) //1st conservation constrain has no inflows
        {
            BOOST_FOREACH(CFlowBase * pInflow, apInflows) {
                CFlowVertex * pInflowVertex = pInflow->getFlow(nTime-1);
                pConserveFlow->addInflow(pInflowVertex);
            }
        }
        
        BOOST_FOREACH(CFlowBase * pOutflow, apOutflows) {
            CFlowVertex * pOutflowVertex = pOutflow->getFlow(nTime);
            pConserveFlow->addOutflow(pOutflowVertex);
        }
    }
    
    CFlowVertex * getOutflow(const int nTime)
    {
        if(apOutflows.size() == 0)
            throw "River has no outflows";
            
        return apOutflows[0]->getFlow(nTime);
    }

    /* write store and outflows */
    virtual void addHeaders(std::ostream & of) const
    {
        of << label << '\t';
        BOOST_FOREACH(CFlowBase * pOutflow, apOutflows) {
            if(bVeryVerbose) cout << "outflow:" << pOutflow->label << " ";
            pOutflow->addHeaders(of);
        }
    }
    
    /* write store and outflows */
    virtual void addData(const int nTime, std::ostream & of) const
    {
        of << getStoreElement(nTime)->volume() << '\t';
        
        BOOST_FOREACH(CFlowBase * pOutflow, apOutflows) {
            pOutflow->addData(nTime, of);
        }        
    }    
};

CFlow::CFlow(const string label, CStore * pSource, CStore * pDest, CDrainageFunction *pDrainageFunction) : CFlowBase(label), pSource(pSource), pDest(pDest), pDrainageFunction(pDrainageFunction)
{
    pSource->addOutflow(this);
    pDest->addInflow(this);
}

void CFlow::addFlowVertex(const int nTime)
{
    CFlowBetweenStoresElement * pElement = new CFlowBetweenStoresElement(pSource->getStoreElement(nTime), pDest->getStoreElement(nTime+1), pDrainageFunction);
    apFlowElements.push_back(pElement);
}

CSinkFlow::CSinkFlow(const string label, CStore * pSource, CDrainageFunction * pDrainageFunction) : CFlowBase(label), pSource(pSource), pDrainageFunction(pDrainageFunction)
{
    pSource->addOutflow(this);
}

void CSinkFlow::addFlowVertex(const int nTime)
{
    CFlowSinkElement * pElement = new CFlowSinkElement(pSource->getStoreElement(nTime), pDrainageFunction);
    apFlowElements.push_back(pElement);
}

CPETDependentETOutFlow::CPETDependentETOutFlow(const string label, CStore * pSoilSource, CStore * pCanopySource, CComplexPETDrainageFunction *pETFunction)
:  pSoilSource(pSoilSource), pCanopySource(pCanopySource), pSoilETFlow(new CSinkFlow(label+"_S2ET", pSoilSource, 0)), pCanopyETFlow(new CSinkFlow(label+"_C2ET", pCanopySource, 0)), pETFunction(pETFunction)
{
    pSoilSource->addOutflow(pSoilETFlow);
    pCanopySource->addOutflow(pCanopyETFlow);
}

//Adds the 5-way constraint on soil store, canopy store, PET, ET soil and canopy outflows 
void CPETDependentETOutFlow::addETComplexFlowVertex(const int nTime)
{
    //CPETDependentETElement * pCanopyETElement = new CPETDependentETElement(nTime);
    //CPETDependentETElement * pSoilETElement = new CPETDependentETElement(nTime);
    pSoilETFlow->addFlowVertex(nTime);
    pCanopyETFlow->addFlowVertex(nTime);
    
    //apPETEstimates.push_back
    
    new CETConstraint(pSoilSource->getStoreElement(nTime), pSoilETFlow->getFlow(nTime), pCanopySource->getStoreElement(nTime), pCanopyETFlow->getFlow(nTime), getPETElement(nTime), pETFunction);
    //Probably only used for logging
    //apFlowElements.push_back(pCanopyETElement); //TODO: what are these used for? We are only tracking one of the 2 ET elements
    //apSoilFlowElements.push_back(pSoilETElement);
    
}

class CSubcatchmentRainfall : public CFlowBase
{
    // apFlowElements contains the estimates of the actual rainfall
    vector<CFlowMeasurement *> apRainGaugeMeasurements;
public:
    CSubcatchmentRainfall(const string label) : CFlowBase(label) {

    }

    void addRainfallMeasurement(const int nTime, const double dMeasurement) { //in cubic metres per timestep
        CFlowVertex* pRainfallEstimateOneTimestep = new CFlowVertex(nTime);
        apFlowElements.push_back(pRainfallEstimateOneTimestep);

        const double dRainfallSD = 0.1*dMeasurement + 1.0/EXACT_INFO; //todo: measurement should be distributed about actual rainfall
        //MUST have a measurement for every estimate, otherwise MLE will be unconstraines
        CFlowMeasurement * pRainGaugeMeasurement = new CFlowMeasurement(pRainfallEstimateOneTimestep, dMeasurement, sqr(dRainfallSD));
        apRainGaugeMeasurements.push_back(pRainGaugeMeasurement);
    }

    virtual void addHeaders(std::ostream & of) const
    {
        of << label << "_measured" << '\t';
        CFlowBase::addHeaders(of);
    }
    
    /* outflow */
    virtual void addData(const int nTime, std::ostream & of) const
    {
        of << apRainGaugeMeasurements[nTime]->measurement() << '\t';
        CFlowBase::addData(nTime, of);
    }    

private:
    virtual void addFlowVertex(const int nTime)
    {
        throw "This is an inflow (unlike ET) so don't call me--use addRainfallMeasurement instead";
    }
};

/**
 * @class CET
 * @brief Represents a set of temperature measurements over time (plus associated ET flows out of a CStore surface water store, so is a CFlowBase)
 */
/*class CET : public CTemperatureDependentOutFlow
{
    CETDrainageFn etDrainageFn;
//	CStoreInitialisingAssumption * pInitAssumption;
public:
    CET(const string label, CStore * pSoilStore) : CTemperatureDependentOutFlow(label, pSoilStore, &etDrainageFn), etDrainageFn(0.03) / * max 2% ET per timestep * / {

    }
    
};
 


double CPETDependantSoilDrainageFunction_base::outflow(const HyperGraph::Vertex * pStore, const HyperGraph::Vertex * pe_potVertex, const HyperGraph::Vertex * pCanopyETOutflowVertex) const 
{
    const double dVolumeInSoilStore = CAST<const CStoreElement*>(pStore)->volume();
    const double e_pot = CAST<const CTemperatureVertex*>(pe_potVertex)->T();
    const double e_c = CAST<const CFlowSinkElement*>(pCanopyETOutflowVertex)->deltaV();
    const double dOutflowETFromSoil = outflow_function(dVolumeInSoilStore, e_pot, e_c);
    
    if(dOutflowETFromSoil < -0.001)
        throw "Warning: negative dOutflowETFromSoil out of store (possible when computing derivatives near to 0)";
    
    return dOutflowETFromSoil;
}*/


class CSubcatchmentParams
{
    string label;
    double dArea_km2;

public:
    double C_c_m3, c_t; //parameters for ET from canopy
    double theta_pa_m3; //ET from soil
    
    const string getLabel() const { return label; }
    //Needed? const double getArea_km2() const { return dArea_km2; }
    
    CSubcatchmentParams(const string label, const double dArea_km2, const double C_c_m3, const double c_t, const double theta_pa_m3) : label(label), dArea_km2(dArea_km2), C_c_m3(C_c_m3), c_t(c_t), theta_pa_m3(theta_pa_m3)
    {
        
        
    }
    
    const double mmPrecipToCubicMetres(double dPrecip_mmPerTimestep) const
    {
        return dPrecip_mmPerTimestep * dArea_km2 * 1000; /* 0.001=mm to m * sqr(1000=km to m) */
    }

};

class CSubCatchment
{
    CSubcatchmentParams params;
    
    CSubcatchmentRainfall rainfall;

    CStore CanopyStore, GWStore, SoilStore, RiverStore;
    
    CPETDependentETOutFlow * pETFlows;
    
    CSubCatchment * pDownstreamCatchment;
    
    vector<CStore *> apStores;
    vector<CFlowBase *> apInputs; //rainfall measurement
public:
    //Create a subcatchment,
    CSubCatchment(const CSubcatchmentParams & params) : params(params), rainfall(label()+"_RainMLE"), CanopyStore(label() + "_Canopy"), GWStore(label()+"_GW"), SoilStore(label()+"_Soil"), RiverStore(label()+"_RW"), /*et(label()+"_ET", &SoilStore),*/ pDownstreamCatchment(0) {
        //Add intra-catchment flows

        /*CETFromCanopyDrainageFn * pCanopyET = new CETFromCanopyDrainageFn(params.C_c_m3, params.c_t);
        new CSinkFlow(label() + "_CanopyET", &CanopyStore, &sink, pCanopyET);*/
        
        //new CETConstraint()

        CDrainageFunction * pCanopyToSurface = new CLinearDrainageFunction(0.1);
        new CFlow(label() + "_C2S", &CanopyStore, &SoilStore, pCanopyToSurface);

        CDrainageFunction * pSurfaceToGroundwater = new CLinearDrainageFunction(0.1);
        new CFlow(label() + "_S2G", &SoilStore, &GWStore, pSurfaceToGroundwater);

        //CETFromSoilDrainageFn * pSurfaceET = new CETFromSoilDrainageFn(params.theta_pa_m3);
        //new CSinkFlow(label() + "_SoilET", &SoilStore, &sink, pSurfaceET);

        CDrainageFunction * pSurfaceToRiver = new CLinearDrainageFunction(0.2);
        new CFlow(label() + "_S2R", &SoilStore, &RiverStore, pSurfaceToRiver);

        CDrainageFunction * pGWToRiver = new CLinearDrainageFunction(0.01);
        new CFlow(label() + "_GW2R", &GWStore, &RiverStore, pGWToRiver);

        //Add a rainfall estimate
        CanopyStore.addInflow(&rainfall);
        apInputs.push_back(&rainfall);
        
        //Make sure all stores are in vector
        apStores.push_back(&CanopyStore);
        apStores.push_back(&SoilStore);
        apStores.push_back(&GWStore);
        apStores.push_back(&RiverStore);
        
        CComplexPETDrainageFunction * pETFn = new CComplexPETDrainageFunction(params.C_c_m3,params.c_t,params.theta_pa_m3);
        pETFlows=new CPETDependentETOutFlow("ET", &SoilStore, &CanopyStore, pETFn);
    }

    const string label() const { return params.getLabel(); }

    //Creates the next store states at nTime+1
    void addTimestep(int nTime) {
        //Add store vertices
        for(CStore * pStore : apStores) {
            pStore->addTimestep(nTime);
        }
    }

    //Creates flows from nTime to nTime+1
    void addFlows(const int nTime, const double dRainfallMeasurement, const double dPETEstimate) {
        //Add flow vertices with a constraint on each
        rainfall.addRainfallMeasurement(nTime, dRainfallMeasurement);
        
        for(CStore * pStore : apStores) {
            pStore->addOutflows(nTime);
        }
        
        pETFlows->addPETEstimate(nTime, dPETEstimate);
        
        pETFlows->addETComplexFlowVertex(nTime);
    }

    //Constrains total flow at nTime
    void addTotalFlowConstraints(int nTime) {
        for(CStore * pStore : apStores) {
            pStore->addTotalFlowConstraint(nTime);
        }
    }

    //void addUpstreamCatchment(CSubCatchment * pUpstream)

    void setDownstreamCatchment(CSubCatchment * pDownstream) {
        pDownstreamCatchment = pDownstream;
        //pDownstream->addUpstreamCatchment(this);

        //Add inter-catchment
        //=Add the flow out of this catchment and into the river of the next catchment
        CStore * pDownstreamRiverStore = pDownstreamCatchment->getRiverStore();
        CDrainageFunction * pRiverToDownstream = new CLinearDrainageFunction(0.9);
        new CFlow(label()+"_to_"+pDownstream->label(), &RiverStore, pDownstreamRiverStore, pRiverToDownstream);
    }
    
    /*This is the bottom of the catchment*/
    void setOutflow() 
    {
        CDrainageFunction * pRiverToDownstream = new CLinearDrainageFunction(0.9); //todo: duplication to remove
        new CSinkFlow(label()+"_out", &RiverStore, pRiverToDownstream);
    }    

    CStore * getRiverStore() {
        return &RiverStore;
    }
    const CStore * getRiverStore() const {
        return &RiverStore;
    }

    virtual void addHeaders(std::ostream & of) const
    {
        BOOST_FOREACH(const CFlowBase * pInputFlow, apInputs) {
            pInputFlow->addHeaders(of);
        }        
        for(const CStore * pStore : apStores) {
            if(bVeryVerbose) cout << "Adding headers for store " << pStore->label << "...";
            pStore->addHeaders(of);
            if(bVeryVerbose) cout << endl;
        }        
    }
    
    virtual void addData(const int nTime, std::ostream & of) const
    {
        BOOST_FOREACH(const CFlowBase * pInputFlow, apInputs) {
            pInputFlow->addData(nTime, of);
        }        
        BOOST_FOREACH(const CStore * pStore, apStores) {
            pStore->addData(nTime, of);
        }        
    }
    
};

class CCatchmentBase
{
public:
    virtual void pp(std::ostream & of) const = 0;
    virtual ~CCatchmentBase() {}
};

/*
 * Timesteps
 * 
 * * Parameters must all be converted to rates-per-timestep when the catchments are setup.
 * * Model will operate on steps of a few minutes, with volume units of cubic metres. 
 * * Measurements will be per hour, or per day, or per few-minutes, and should be upsampled appropriately.
 * 
 * Volumes in stores could be hundreds to humdreds of millions of cubic metres. Flows should usually be a few cubic metres.
 * 
 * If *everything* passed to the model is per-timestep then units are relatively simple. Need to make sure conversion is done elsewhere. Need to know the timestep in order to output cumecs.
 * 
 * 
 */
class CCatchment : public CCatchmentBase
{
    vector<CSubCatchment *> apCatchments;
    
    int nEndTime;//keep track of most recent time
    const double dTimestepLengthSeconds;
    
    std::map<int, CFlowMeasurement*> apFlowMeasurements; //Occasionally measure the flow
    
public:
    CCatchment(const double dTimestepLengthSeconds, const std::vector<CSubcatchmentParams> & aSubcatchmentsParams)
        //: northBranch("northBranch"), southBranch("southBranch"), downstreamBranch("downstreamBranch"), 
        : nEndTime(-1), dTimestepLengthSeconds(dTimestepLengthSeconds) {
        
        BOOST_FOREACH(const CSubcatchmentParams & params, aSubcatchmentsParams)
        {
            apCatchments.push_back(new CSubCatchment(params));
        }
        /*northBranch.setDownstreamCatchment(&downstreamBranch);
		southBranch.setDownstreamCatchment(&downstreamBranch);
		downstreamBranch.setOutflow(); //end of catchment*/
        apCatchments[0]->setDownstreamCatchment(apCatchments[2]);
        apCatchments[1]->setDownstreamCatchment(apCatchments[2]);
        apCatchments[2]->setOutflow();
    }

    void addTimestep(const int nTime, const double dActualTime_TODO, const double dRainfallMeasurement_mmPerTimestep, const double dPET_m3PerTimestep) {
        nEndTime = nTime;
        
        for(CSubCatchment * pCatchment : apCatchments) {
            pCatchment->addTimestep(nTime);
        }
        for(CSubCatchment * pCatchment : apCatchments) {
            pCatchment->addFlows(nTime, dRainfallMeasurement_mmPerTimestep, dPET_m3PerTimestep);
        }
        for(CSubCatchment * pCatchment : apCatchments) {
            pCatchment->addTotalFlowConstraints(nTime);
        }
        
        for(CTotalFlowConstraint * pTFConstraint : g_aTotalFlowConstraints)
            optimizer.addEdge(pTFConstraint);
        g_aTotalFlowConstraints.clear();
    }

    void addFlow(const int nTime, const double dMeasurement) {
        CFlowVertex * pMeasuredFlow = apCatchments[2]->getRiverStore()->getOutflow(nTime);
        CFlowMeasurement * pFlowMeasurement = new CFlowMeasurement(pMeasuredFlow, dMeasurement, 0.0001*dMeasurement);
        apFlowMeasurements[nTime] = pFlowMeasurement;
    }
    
private:
    virtual void addHeaders(std::ostream & of) const
    {
        BOOST_FOREACH(CSubCatchment * pCatchment, apCatchments) {
            pCatchment->addHeaders(of);
        }        
        of << endl;
    }
    
    virtual void addData(const int nTime, std::ostream & of) const
    {
        BOOST_FOREACH(CSubCatchment * pCatchment, apCatchments) {
            pCatchment->addData(nTime, of);
        }
        
        //Add flow measurements, if we have them
        if(apFlowMeasurements.find(nTime) != apFlowMeasurements.end())
            of << (*apFlowMeasurements.find(nTime)).second->measurement();
        
        of << endl;
    }
public:
    void pp(std::ostream & of) const
    {
        //Add timestamp
        time_t currentTime = time(0);
        string strTime = ctime(&currentTime);
        of << strTime.substr(0, strTime.length()-1) << '\t';
        addHeaders(of);
        
        for(int nTime=0; nTime < nEndTime; nTime++)
        {
            of << nTime << '\t';
            addData(nTime, of);
        }
        
    }
};

/*
 * Waihua catchment has the same Y-shape form as original sim catchment. We have [forecast] measurements of outflow, temperature, humidity, precip
 *
class CWaihuaCatchment : public CCatchmentBase
{
	CSubCatchment northBranch, southBranch, downstreamBranch;

	vector<CSubCatchment *> apCatchments;
	
	int nEndTime;//keep track of most recent time
	
	std::map<int, CFlowMeasurement*> apFlowMeasurements; //Occasionally measure the flow
	
public:
	CWaihuaCatchment() : northBranch("northBranch"), southBranch("southBranch"), downstreamBranch("downstreamBranch"), nEndTime(-1) {
		northBranch.setDownstreamCatchment(&downstreamBranch);
		southBranch.setDownstreamCatchment(&downstreamBranch);
		downstreamBranch.setOutflow(); //end of catchment

		apCatchments.push_back(&northBranch);
		apCatchments.push_back(&southBranch);
		apCatchments.push_back(&downstreamBranch);
	}

	void addTimestep(const int nTime, const double dRainfallMeasurementMM, const double dTemperatureMeasurement, const double dHumidityMeasurement/ *, const optional<double> pdFlowMeasurement* /) {
		nEndTime = nTime;
		
		BOOST_FOREACH(CSubCatchment * pCatchment, apCatchments) {
			pCatchment->addTimestep(nTime);
		}
		BOOST_FOREACH(CSubCatchment * pCatchment, apCatchments) {
			pCatchment->addFlows(nTime, dRainfallMeasurement, dTemperatureMeasurement);
		}
		BOOST_FOREACH(CSubCatchment * pCatchment, apCatchments) {
			pCatchment->addTotalFlowConstraints(nTime);
		}
		
		/ *if(pdFlowMeasurement)
			addFlow(nTime, *pdMeasurement);* /
	}

	void addFlow(const int nTime, const double dMeasurement) {
		CFlowVertex * pMeasuredFlow = downstreamBranch.getRiverStore()->getOutflow(nTime);
		CFlowMeasurement * pFlowMeasurement = new CFlowMeasurement(pMeasuredFlow, dMeasurement, 0.0001*dMeasurement);
		apFlowMeasurements[nTime] = pFlowMeasurement;
	}
	
private:
	virtual void addHeaders(std::ostream & of) const
	{
		BOOST_FOREACH(CSubCatchment * pCatchment, apCatchments) {
			pCatchment->addHeaders(of);
		}		
		of << endl;
	}
	
	virtual void addData(const int nTime, std::ostream & of) const
	{
		BOOST_FOREACH(CSubCatchment * pCatchment, apCatchments) {
			pCatchment->addData(nTime, of);
		}
		
		//Add flow measurements, if we have them
		if(apFlowMeasurements.find(nTime) != apFlowMeasurements.end())
			of << (*apFlowMeasurements.find(nTime)).second->measurement();
		
		of << endl;
	}
public:
	void pp(std::ostream & of) const
	{
		//Add timestamp
		time_t currentTime = time(0);
		string strTime = ctime(&currentTime);
		of << strTime.substr(0, strTime.length()-1) << '\t';
		addHeaders(of);
		
		for(int nTime=0; nTime < nEndTime; nTime++)
		{
			of << nTime << '\t';
			addData(nTime, of);
		}
		
	}
};*/

#endif
