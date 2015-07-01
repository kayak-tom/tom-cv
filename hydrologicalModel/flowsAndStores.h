/*Vertices and edges of hypergraph
 * */
#ifndef _FLOWSANDSTORES_H
#define _FLOWSANDSTORES_H

#include "graphbase.h"
#include "flowUtil.h"
#include "boost/foreach.hpp"
#include <boost/math/distributions/normal.hpp>

static const bool bVeryVerbose = false;

class CFlowBase;
class CPETDependentETOutFlow;
class CStore;

/**
 * @class CFlowVertex
 * @brief A flow from a store at time nTimeFrom to nTimeFrom+1
 */
class CFlowVertex : public CVertex1d
{
    const int nTimeFrom;
    const CFlowBase * pParentFlow;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    
    CFlowVertex(const int nTimeFrom=-1, const CFlowBase * pParentFlow=0) : nTimeFrom(nTimeFrom), pParentFlow(pParentFlow) 
    {
        if(bVeryVerbose) cout << "T" << nTimeFrom << ": Added flow vertex to following time" << endl;
    }
    
    double deltaV() const { return pseudoHuber(_estimate); }
    const int timeFrom() const { return nTimeFrom; }
    
    string label() const;
    
    virtual void oplusImpl(const double* v) 
    {
        _estimate += *v;
        if(_estimate < -10 || _estimate > 10000)
        {
            cout << label() << ": _estimate=" << _estimate << " warning large flow" << endl;
            //throw "Flow estimate going negative";
        }
    }

    //! sets the node to the origin (used in the multilevel stuff)
    virtual void setToOriginImpl() 
    {
        setEstimate(stateToEstimate(0.1) /*cu m/timestep*/);
    }    
};


/* A PET measurement at time t */
class CPETVertex : public CVertex1d
{
    const int nTimeFrom;
    const CPETDependentETOutFlow * pParent;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    
    CPETVertex(const int nTimeFrom=-1, const CPETDependentETOutFlow *pParent=0) : nTimeFrom(nTimeFrom), pParent(pParent)
    {
        if(bVeryVerbose) cout << label() << ": Added temperature estimate" << endl;
    }
    
    double T() const { return  estimateToState(_estimate); }
    const int timeFrom() const { return nTimeFrom; }

    string label() const;
    
    virtual void oplusImpl(const double* v)
    {
        _estimate += *v;
        if(_estimate < -100 || _estimate > 100)
        {
            cout << label() << ": _estimate=" << _estimate << " warning large ET (m^3 per timestep)" << endl;
            //throw "Flow estimate going negative";
        }
    }

    //! sets the node to the origin (used in the multilevel stuff)
    virtual void setToOriginImpl() { setEstimate(stateToEstimate(1)/*cu m per timestep*/); }
    
};

class CStoreElement : public CVertex1d
{
    const int nTime;
    const CStore * pParent;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    CStoreElement(const int nTime=-1, const CStore * pParent=0) : nTime(nTime), pParent(pParent)
    {
        
    }

    string label() const;
    
    double volume() const 
    {
        return estimateToState(_estimate);
    }
    
    virtual void oplusImpl(const double* v) 
    {
        _estimate += *v;
        
        if(_estimate < -10 || _estimate > 100)
        {
            cout << label() << ": _estimate=" << _estimate << " warning large ET (m^3 per timestep)" << endl;
            //throw "Flow estimate going negative";
        }
    }


    //! sets the node to the origin (used in the multilevel stuff)
    virtual void setToOriginImpl() 
    {
        setEstimate(stateToEstimate(1) /*cu m*/);
    }
    
    
    const int time() const { return nTime; }
};

//Measure a flow (either a river flow, an input flow (rainfall) or an output flow (a flow gauge)
class CFlowMeasurement : public BaseUnaryEdge<1, double, CFlowVertex>
{
    typedef BaseUnaryEdge<1, double, CFlowVertex> Base;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    CFlowMeasurement(CFlowVertex* pFlowVertex, const double dFlowMeasurement, const double dCovariance)
    {
        _vertices[0] = pFlowVertex;
        setMeasurement(dFlowMeasurement);
        Base::InformationType info;
        info(0) = 1.0/dCovariance;
        setInformation(info);
        
        optimizer.addEdge(this);
    }

    void computeError() {
        const CFlowVertex* pThisState = CAST<const CFlowVertex*>(_vertices[0]);

        _error(0) = (pThisState->deltaV() - _measurement);
        if(bVeryVerbose) cout << "Flow measurement error " << _error << endl;
    }
    
    virtual bool read(std::istream& is)
    {
        /*string label;
		is >> label;
		is >> dFlowMeasurement;
		is >> dCovariance;*/

        throw true;
    }
    
    virtual bool write(std::ostream& os) const {
        throw true;
    }    
};

/**
 * @class CPETMeasurement
 * @brief Measurement of PET (cu m per timestep)
 */
class CPETMeasurement : public BaseUnaryEdge<1, double, CFlowVertex>
{
    typedef BaseUnaryEdge<1, double, CFlowVertex> Base;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    CPETMeasurement(CPETVertex* pPETVertex, const double dPETMeasurement, const double dPETVariance)
    {
        _vertices[0] = pPETVertex;
        setMeasurement(dPETMeasurement);
        Base::InformationType info;
        info(0) = 1.0/dPETVariance;
        setInformation(info);
        
        optimizer.addEdge(this);
        
        if(bVeryVerbose) cout << "T" << pPETVertex->timeFrom() << ": Added temperature measurement of " << dPETMeasurement << endl;
    }

    void computeError() {
        const CPETVertex* pThisState = CAST<const CPETVertex*>(_vertices[0]);

        _error(0) = (pThisState->T() - _measurement);
        if(bVeryVerbose) cout << "PET error = " << _error << endl;
    }
    
    virtual bool read(std::istream& is)
    {
        throw true;
    }
    
    virtual bool write(std::ostream& os) const {
        throw true;
    }    
};

/**
 * @class CStoreInitialisingAssumption
 * @brief Add an initialising assumption for a store. Shouldn't actually be needed once flows are known. Equivalent to fixing an origin. Very important that this is strictly positive, to ensure that all subsequent measurements are > 0
 */
class CStoreInitialisingAssumption : public BaseUnaryEdge<1, double, CStoreElement>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    CStoreInitialisingAssumption(CStoreElement * pStoreVertex)
    {
        _vertices[0] = pStoreVertex;
        
        const double dPrior = 0.1;
        
        setMeasurement(dPrior);
        setInformation(UNINFORMATIVE_INFO_1D);

        optimizer.addEdge(this);
    }

    void computeError() {
        const CStoreElement* pThisState = CAST<const CStoreElement*>(_vertices[0]);

        _error(0) = (pThisState->volume() - _measurement);
        const double dSmallPositiveVal = 1e-4, PENALTY = 1e+10; //Penalty should be big compared to UNINFORMATIVE_INFO_1D*dSmallPositiveVal
        if(pThisState->volume() < dSmallPositiveVal)
        {
            _error(0) += (pThisState->volume()-dSmallPositiveVal)*PENALTY ; 
            if(bVeryVerbose) cout << "Error (weak initial constraint): " << _error  << " value " << pThisState->volume() << endl;
        }
    }
    
    virtual bool read(std::istream& is)
    {
        throw true;
    }
    
    virtual bool write(std::ostream& os) const {
        throw true;
    }    
};

//Add an initialising assumption for a Flow. Shouldn't actually be needed once flows are known. Equivalent to fixing an origin
//Very important that this is strictly positive, to ensure that all subsequent measurements are > 0
/*class CFlowInitialisingAssumption : public BaseUnaryEdge<1, double, CFlowVertex> 
{

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	CFlowInitialisingAssumption(CFlowVertex * pFlowVertex)
	{
		_vertices[0] = pFlowVertex;
		
		const double dPrior = 0.1;
		
		setMeasurement(dPrior);
		setInverseMeasurement(-dPrior);Temperature
		setInformation(EXACT_INFO_1D);

		optimizer.addEdge(this);
	}CPETDependantCanopyDrainageFunction_base

	void computeError() {Temperature
		const CFlowVertex* pThisState = CAST<const CFlowVertex*>(_vertices[0]);

		if(pThisState->deltaV() < 1e-4)
			_error(0) = -HUGE;
		else
			_error(0) = ((pThisState->deltaV()) - _measurement);
		//cout << "Error (weak initial constraint): " << _error << endl;
	}
	
	virtual bool read(std::istream& is)
	{
		throw true;
	}
	
	virtual bool write(std::ostream& os) const {
		throw true;
	}	
};*/

class CDrainageFunction
{
protected:
    virtual double outflow_function(const double dVolume) const = 0;
public:
    double outflow(const HyperGraph::Vertex * pStore) const 
    {
        const double dOutflow = outflow_function(CAST<const CStoreElement*>(pStore)->volume());
        
        if(dOutflow < -0.001)
            throw "Warning: negative flow out of store (possible when computing derivatives near to 0)";
        
        return dOutflow;
    }
};

class CLinearDrainageFunction : public CDrainageFunction
{
protected:
    virtual double outflow_function(const double dVolume) const 
    {
        return dDrainageRate * dVolume; 
        /*const double eps = 0.001;
		if(dDrainageRate>0)
			return dDrainageRate * dVolume //at 0: Gradient dVolume, value 0
		return eps * exp(dDrainageRate*dVolume); //at 0 value eps, gradient dDrainageRate*eps */
    }
    const double dDrainageRate;
public:
    CLinearDrainageFunction(const double dDrainageRate) : dDrainageRate(dDrainageRate) {}
};

/*class CTempDependantDrainageFunction//like a drainage function but temperature-dependent (for snowmelt or ET)
{
protected:
    virtual double outflow_function(const double dVolume, const double dTemperature, const double dHumidity) const = 0;
public:
    double outflow(const HyperGraph::Vertex * pStore, const HyperGraph::Vertex * pTemp, const HyperGraph::Vertex * pHumidity) const 
    {
        const double dVolumeInStore = CAST<const CStoreElement*>(pStore)->volume();
        const double dTemperature = CAST<const CPETVertex*>(pTemp)->T();
        const double dHumidity = CAST<const CPETVertex*>(pHumidity)->T();
        const double dOutflow = outflow_function(dVolumeInStore, dTemperature, dHumidity);
        
        if(dOutflow < -0.001)
            throw "Warning: negative flow out of store (possible when computing derivatives near to 0)";
        
        return dOutflow;
    }

    virtual ~CTempDependantDrainageFunction() {}
};*/

/*class CPETDependantCanopyDrainageFunction_base//like a drainage function but dependent on the PET value from TOPNET
{
protected:
	virtual double outflow_function(const double dVolume, const double dPET) const = 0;
public:
	double outflow(const HyperGraph::Vertex * pStore, const HyperGraph::Vertex * pe_potVertex) const 
	{
		const double dVolumeInCanopyStore = CAST<const CStoreElement*>(pStore)->volume();
		const double dPET = CAST<const CTemperatureVertex*>(pe_potVertex)->T();
		const double dOutflowETFromCanopy = outflow_function(dVolumeInCanopyStore, dPET);
		
		if(dOutflowETFromCanopy < -0.001)
			throw "Warning: negative dOutflowETFromCanopy out of store (possible when computing derivatives near to 0)";
		
		return dOutflowETFromCanopy;
	}

	virtual ~CPETDependantCanopyDrainageFunction_base() {}
};c_t

class CPETDependantSoilDrainageFunction_base//like a drainage function but dependent on the PET value from TOPNET, and the PET used up by evaporation from the canopy
{
protected:
	virtual double outflow_function(const double dVolumeInSoilStore, const double e_pot, const double e_c) const = 0;
public:
	double outflow(const HyperGraph::Vertex * pStore, const HyperGraph::Vertex * pe_potVertex, const HyperGraph::Vertex * pCanopyETOutflowVertex) const;

	virtual ~CPETDependantSoilDrainageFunction_base() {}
};

class CETFromCanopyDrainageFn : public CPETDependantCanopyDrainageFunction_base //Like a drainage function for the ET rate
{
	double f(const double S_c) const
	{
		double dFrac = S_c/C_c_m3;
		return dFrac*(2-dFrac);
	}
protected:
	/ **
	 * @brief Equations from A2
	 * @param dVolume
	 * @param dPET
	 * @return 
	 * /
	virtual double outflow_function(const double dVolume, const double epot) const 
	{
		//const boost::math::normal_distribution<double> tempCDFdistn(20, 5);
		//const double dTemperatureERF = boost::math::cdf(tempCDFdistn, dTemperature);
		//return dMaxDrainageRate * dVolume * dTemperatureERF;
		const double fS_c = f(dVolume);
		const double e_c=epot * c_t * fS_c;
		if(e_c > dVolume)
			cout << "Warning: e_c greater than canopy store volume" << endl;
			
		if(e_c<0)
			throw "e_c is negative";
		
		return std::min<double>(dVolume, e_c);
	}
	const double C_c_m3, c_t;
public:
	CETFromCanopyDrainageFn(const double C_c_m3 / * water holding capacity of the canopy * /, const double c_t / * free parameter in A2d * /) : C_c_m3(C_c_m3), c_t(c_t) {}
};

class CETFromSoilDrainageFn : public CPETDependantSoilDrainageFunction_base //Like a drainage function for the ET rate
{
protected:
	/ **
	 * @brief Equation from A4f
	 * @param dVolume
	 * @param dPET
	 * @return 
	 * /
	virtual double outflow_function(const double dVolumeInSoilStore, const double e_pot, const double e_c) const 
	{
		const double e_r=(e_pot - e_c) * std::min<double>(dVolumeInSoilStore/theta_pa_m3, 1); //TODO: smooth threshold
		if(e_r > dVolumeInSoilStore)
			cout << "Warning: e_r greater than canopy store volume" << endl;
			
		if(e_r<0)
			throw "e_r is negative";
		
		return std::min<double>(dVolumeInSoilStore, e_r);
	}
	const double theta_pa_m3;
public:
	CETFromSoilDrainageFn(const double theta_pa_m3 / * plant-available water content * /) : theta_pa_m3(theta_pa_m3) {}
};*/









//Constrain a flow out of a store to be some function of the value of the store
class COutFlowConstraint : public BaseBinaryEdge<1, double, CStoreElement, CFlowVertex>
{
    const CDrainageFunction * pDrainage;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    COutFlowConstraint(CStoreElement * pStoreVertex, CFlowVertex* pFlowVertex, const CDrainageFunction * pDrainage) : pDrainage(pDrainage)
    {
        _vertices[0] = pStoreVertex;
        _vertices[1] = pFlowVertex;
        setInformation(EXACT_INFO_1D);
        optimizer.addEdge(this);
    }

    void computeError() {
        const CFlowVertex* pFlow = CAST<const CFlowVertex*>(_vertices[1]);

        _error(0) = pDrainage->outflow(_vertices[0]) - pFlow->deltaV();
        if(bVeryVerbose) cout << "Error " << _error << endl;
    }
    
    virtual bool read(std::istream& is)
    {
        throw true;
    }
    
    virtual bool write(std::ostream& os) const {
        throw true;
    }    
};


//Constrain a flow out of a store to be some function of the value of the store and the PET directly from TOPNET. PET can either be fixed or can be optimised. The PET value is represented by a CTemperatureVertex
/*class CPETDependentFlowConstraint : public BaseMultiEdge<1, double>
{
    CTemperatureVertex * pTempEstimate;
    const CPETDependantSoilDrainageFunction_base * pDrainage;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    CPETDependentFlowConstraint(CStoreElement * pStoreVertex, CTemperatureVertex * pPETEstimate, CFlowVertex* pFlowVertex, const CPETDependantSoilDrainageFunction_base * pDrainage) : pTempEstimate(pTempEstimate), pDrainage(pDrainage)
    {
        _vertices.push_back(pStoreVertex);
        _vertices.push_back(pPETEstimate);
        _vertices.push_back(pFlowVertex);
        setInformation(EXACT_INFO_1D);

        resize(_vertices.size());

        optimizer.addEdge(this);
        
        if(bVeryVerbose) cout << "T" << pStoreVertex->time() << ": Added PET-dependent rate controlling edge " << endl;
    }

    void computeError() {
        const CFlowVertex* pFlow = CAST<const CFlowVertex*>(_vertices[3]); //TODO: indices

        _error(0) = pDrainage->outflow(_vertices[0], _vertices[1], _vertices[2]) - pFlow->deltaV();
        if(bVeryVerbose) cout << "ET error " << _error << endl;
    }
    
    virtual bool read(std::istream& is)
    {
        throw true;
    }
    
    virtual bool write(std::ostream& os) const {
        throw true;
    }    
};

//Constrain a flow out of a store to be some function of the value of the store and a temperature
class COutFlowTemperatureConstraint : public BaseMultiEdge<1, double>
{
    CPETVertex * pTempEstimate;
    const CTempDependantDrainageFunction * pDrainage;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    COutFlowTemperatureConstraint(CStoreElement * pStoreVertex, CPETVertex * pTempEstimate, CPETVertex * pHumidityEstimate, CFlowVertex* pFlowVertex, const CTempDependantDrainageFunction * pDrainage) : pTempEstimate(pTempEstimate), pDrainage(pDrainage)
    {
        _vertices.push_back(pStoreVertex);
        _vertices.push_back(pTempEstimate);
        _vertices.push_back(pHumidityEstimate);
        _vertices.push_back(pFlowVertex);
        setInformation(EXACT_INFO_1D);

        resize(_vertices.size());

        optimizer.addEdge(this);
        
        if(bVeryVerbose) cout << "T" << pStoreVertex->time() << ": Added ET rate controlling edge " << endl;
    }

    void computeError() {
        const CFlowVertex* pFlow = CAST<const CFlowVertex*>(_vertices[3]); //TODO: indices

        _error(0) = pDrainage->outflow(_vertices[0], _vertices[1], _vertices[2]) - pFlow->deltaV();
        if(bVeryVerbose) cout << "ET error " << _error << endl;
    }
    
    virtual bool read(std::istream& is)
    {
        throw true;
    }
    
    virtual bool write(std::ostream& os) const {
        throw true;
    }    
};*/

class CComplexPETDrainageFunction
{
    const double C_c_m3, c_t; //Canopy parameters
    const double theta_pa_m3; //Soil parameter
    
    double f(const double S_c) const
    {
        double dFrac = S_c/C_c_m3;
        return dFrac*(2-dFrac);
    }
    
    virtual double canopyET(const double dVolume, const double epot) const 
    {
        //const boost::math::normal_distribution<double> tempCDFdistn(20, 5);
        //const double dTemperatureERF = boost::math::cdf(tempCDFdistn, dTemperature);
        //return dMaxDrainageRate * dVolume * dTemperatureERF;
        const double fS_c = f(dVolume);
        const double evaporationRate = estimateToState(c_t * fS_c, 0.0001);
        
        if(evaporationRate < 0 || evaporationRate > 1)
            throw "evaporationRate is the rate at which water can be evaporated from the canopy";
        const double e_c = epot * evaporationRate;
        if(e_c > dVolume)
            throw "Error: e_c greater than canopy store volume (I don't see how the equations prevent this from happening? Maybe just needs a min(.))";
            
        if(e_c<0)
            throw "e_c is negative";
        
        return std::min<double>(dVolume, e_c);
    }

    /**
	 * @brief Equation from A4f
	 * @param dVolume
	 * @param dPET
	 * @return 
	 */
    virtual double soilET(const double dVolumeInSoilStore, const double e_pot, const double e_c) const 
    {
        const double e_r=(e_pot - e_c) * std::min<double>(dVolumeInSoilStore/theta_pa_m3, 1); //TODO: smooth threshold
        if(e_r > dVolumeInSoilStore)
            cout << "Warning: e_r greater than canopy store volume" << endl;
            
        if(e_r<0)
            throw "e_r is negative";
        
        return std::min<double>(dVolumeInSoilStore, e_r);
    }

public:

    CComplexPETDrainageFunction(const double C_c_m3 /* water holding capacity of the canopy */, const double c_t /* free parameter in A2d */, const double theta_pa_m3) : C_c_m3(C_c_m3), c_t(c_t), theta_pa_m3(theta_pa_m3) {}

    //e_pot is the volume-per-timestep-over-the-catchment which could be evaporated
    Eigen::Vector2d soilCanopyJointETFunction(const double dVolumeInSoilStore, const double dVolumeInCanopyStore, const double e_pot) const
    {
        Eigen::Vector2d ETFlowOutFromSoilCanopy;
        
        ETFlowOutFromSoilCanopy(1) = canopyET(dVolumeInCanopyStore, e_pot);
        ETFlowOutFromSoilCanopy(0) = soilET(dVolumeInSoilStore, e_pot, ETFlowOutFromSoilCanopy(1));
        
        return ETFlowOutFromSoilCanopy;        
    }

};


/*
 * Constrains outflows from soil and canopy to match PET
 */
class CETConstraint : public BaseMultiEdge<2, Eigen::Vector2d>
{
    CStoreElement * pSoilStoreVertex;
    CFlowVertex* pSoilETOutVertex; 
    CStoreElement * pCanopyStoreVertex;
    CFlowVertex* pCanopyETOutVertex;
    CPETVertex * pPETEstimate;

    const CComplexPETDrainageFunction * pETFlowFn;
    enum eComplexETConstraintVertices { eSoilStore,eSoilETOut, eCanopyStore, eCanopyETOut, ePETEstimate, NUM_VERTICES };
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
//pSoilSource->getStoreElement(nTime),  pSoilETElement,pCanopySource->getStoreElement(nTime), pCanopyETElement, getPETElement(nTime), pDrainageFunction
    CETConstraint(CStoreElement * pSoilStoreVertex, CFlowVertex* pSoilETOutVertex, CStoreElement * pCanopyStoreVertex, CFlowVertex* pCanopyETOutVertex, CPETVertex * pPETEstimate, const CComplexPETDrainageFunction * pETFlowFn) : /*pTempEstimate(pTempEstimate),*/
        pSoilStoreVertex(pSoilStoreVertex),pSoilETOutVertex(pSoilETOutVertex), pCanopyStoreVertex(pCanopyStoreVertex), pCanopyETOutVertex(pCanopyETOutVertex), pPETEstimate(pPETEstimate), pETFlowFn(pETFlowFn)
    {
        _vertices.resize(NUM_VERTICES);
        _vertices[eSoilStore] = pSoilStoreVertex;
        _vertices[eSoilETOut] = pSoilETOutVertex;
        _vertices[eCanopyStore] = pCanopyStoreVertex;
        _vertices[eCanopyETOut] = pCanopyETOutVertex;
        _vertices[ePETEstimate] = pPETEstimate;
        
        setInformation(Eigen::Matrix2d::Identity()*EXACT_INFO_1D(0,0));
        
        resize(NUM_VERTICES);
        
        optimizer.addEdge(this);
    }

    void computeError() {
        _error = pETFlowFn->soilCanopyJointETFunction(pSoilStoreVertex->volume(), pCanopyStoreVertex->volume(), pPETEstimate->T()) - Eigen::Vector2d(pSoilETOutVertex->deltaV(), pCanopyETOutVertex->deltaV());
        cout << "Error " << _error << endl;
    }
    
    virtual bool read(std::istream& is)
    {
        throw true;
    }
    
    virtual bool write(std::ostream& os) const {
        throw true;
    }    
};

/* Hyperedge constraining total flow in/out of a store. Ensures that flow is conserved at a store over a timestep. 
 * store + total_flow_in -total_flow_out = nextstore
 * */
class CTotalFlowConstraint : public BaseMultiEdge<1, double>
{
    static double sumFlows(const vector<CFlowVertex*> & aFlows) 
    {
        double dTotalFlow = 0;
        BOOST_FOREACH(const CFlowVertex* pFlow, aFlows)
            dTotalFlow += pFlow->deltaV();
        return dTotalFlow;
    }
    
    vector<CFlowVertex*> _inflow, _outflow;
    CStoreElement * pStore_timeT, * pStore_timeT_plus_1;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    CTotalFlowConstraint(CStoreElement * pStore_timeT, CStoreElement * pStore_timeT_plus_1) : pStore_timeT(pStore_timeT), pStore_timeT_plus_1(pStore_timeT_plus_1)
    {
        _vertices.push_back(pStore_timeT);
        _vertices.push_back(pStore_timeT_plus_1);
        //setMeasurement(dFlowMeasurement);
        //setInverseMeasurement(-dFlowMeasurement);
        setInformation(EXACT_INFO_1D);

        resize(_vertices.size());
        g_aTotalFlowConstraints.push_back(this);
        //optimizer.addEdge(this); not until we've done adding inflows and outflows
    }
    
    void addInflow(CFlowVertex* pFlowVertex)
    {
        _vertices.push_back(pFlowVertex);
        _inflow.push_back(pFlowVertex);
        resize(_vertices.size());
    }
    
    void addOutflow(CFlowVertex* pFlowVertex)
    {
        _vertices.push_back(pFlowVertex);
        _outflow.push_back(pFlowVertex);
        resize(_vertices.size());
    }
    
    void computeError() {
        const double dTotalFlowIn = sumFlows(_inflow);
        const double dTotalFlowOut = sumFlows(_outflow);
        const double dStoreBefore = pStore_timeT->volume();
        const double dStoreAfter = pStore_timeT_plus_1->volume();
        _error(0) = (dStoreBefore + dTotalFlowIn - dTotalFlowOut) - dStoreAfter;
        
        /*#ifdef QT_NO_DEBUG
		cout << "Store before " << dStoreBefore << endl;
		cout << "Inflows " << dTotalFlowIn << endl;
		cout << "Outflows " << dTotalFlowOut << endl;
		cout << "Store after " << dStoreAfter << endl;
		cout << "Error " << _error << endl;
		#endif*/
    }
    
    virtual bool read(std::istream& is)
    {
        throw true;
    }
    
    virtual bool write(std::ostream& os) const {
        throw true;
    }    
};




#endif
