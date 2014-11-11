/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * EdgeScaleObserver.h
 *
 *  Created on: 24/08/2009
 *      Author: tom
 *
 *  Knows about 3d structure from edges
 */
#include "util/dynArray.h"
#include "geom/geom.h"
#include "BoWSpeedo.h"
#include "bow/scales.h"

#ifndef EDGESCALEOBSERVER_H_
#define EDGESCALEOBSERVER_H_

class CEdge
{
protected:
	const imageNum nId1, nId2;
	CScale SCOREscale;
	CBoWSpeedo::CObservedObjectsVec * pObservedObjects;;
public:
	CEdge(imageNum nId1, imageNum nId2) : nId1(nId1), nId2(nId2), pObservedObjects(0) {}
	static inline const double NO_SCALE() { return -HUGE; }

	CBoWSpeedo::CObservedObjectsVec * getObservations() const
	{
		return pObservedObjects;
	}

	void setObservedObjVec(CBoWSpeedo::CObservedObjectsVec * pObservedObjects_in)
	{
		pObservedObjects = pObservedObjects_in;
	}

	virtual const CScale SLAMscale() const = 0;
	bool hasScale() const { return SLAMscale().hasScale(); }
	virtual void notifyUpdatedORScale() = 0;
	virtual void updateScaleWithOR(double dORPredictedVal, double dORPredictedVar)
	{
		SCOREscale = CScale(dORPredictedVal, dORPredictedVar);
		if(SCOREscale.hasScale())
			notifyUpdatedORScale();
	}
	void updateScale(CBoWSpeedo::TObservedObjectsVec & aObservedObjects, int nEdgeIdx = -1);

	//Return length^2, NOT including estimated scale
	virtual double measureObject(CBoWSpeedo::CBoWObjectOccurance * pObjOccurance) = 0;

	//Return log of length
	double measureObject_LN(CBoWSpeedo::CBoWObjectOccurance * pObjOccurance)
	{
		double dMeasurement = measureObject(pObjOccurance);
		if(dMeasurement > 0)
			return 0.5*log(dMeasurement); //0.5 undoes square
		else
			return NO_SCALE();
	}

	bool hasORscale() const { return SCOREscale.hasScale(); }

	virtual void findReconstructedPairs(CBoWSpeedo::TObjectsAndLocations & localObjectFinder) const = 0;

	/*double ORscale() const
	{
		if(IS_DEBUG) CHECK(dORscale<=0, "Scale not available; shouldn't be accessed");
		return dORscale;
	};
	double ORvariance() const
	{
		//variance might be very big though--indicates pos came from high badness observations
		if(IS_DEBUG) CHECK(dORvariance<=0, "Variance not available; shouldn't be accessed");
		return dORvariance;
	};*/

	void findObjects(CBoWSpeedo::CObjectFinder & objectFinder);

	imageNum id1() const { return nId1; }
	imageNum id2() const { return nId2; }

	void reset()
	{
		SCOREscale.setUninit();
	}
};

class CImageSource;
class CImParams;

typedef std::map<const CBoWSpeedo::CBoWObject *, int> TObsCounter;

struct _IplImage;

class CEdgeScaleObserver : public CBoWSpeedo::CScaleObserver
{
	CBoWSpeedo::CBoWSpeedometer * pSpeedo;
	CImageSource * pImageLoader;
	typedef CDynArray<CEdge *> TEdgeVec;
	TEdgeVec vEdges;
	const CImParams & IM_PARAMS;

	void resetScales();
    void markOneObjectObservation(int id, _IplImage *pFrameCol1, const CBoWSpeedo::CBoWObjectOccurance *pObservedObj, _IplImage * pFrameCol2) const;
    void measureObjectOccurancesInEdge(CBoWSpeedo::CObservedObjectsVec & observedObjects, CEdge * pEdge, const int nEdgeIdx, TObsCounter *pCountObservations) const;
public:
	CEdgeScaleObserver(const CImParams & IM_PARAMS) : pSpeedo(0), pImageLoader(0), IM_PARAMS(IM_PARAMS) {}
	virtual ~CEdgeScaleObserver();
	virtual void reobserveScales(bool bDrawObjects);
	virtual void addEdge(CEdge * pEdge);
	virtual void removeEdge(CEdge * pEdge);

	virtual void init(CBoWSpeedo::CBoWSpeedometer * speedo) { pSpeedo = speedo; }
	virtual void initImLoader(CImageSource * pImageLoader_in) { pImageLoader = pImageLoader_in; }
	void drawObjects(const TObsCounter & bestObjects, CBoWSpeedo::CObservedObjectsVec & observedObjects) const;
	void printScales() const;
	virtual void findObjects(const CBoWSpeedo & bow, CBoWSpeedo::CObjectFinder & objectFinder);
	void markObjects(const int nIdLast, const int nIdCurrent, _IplImage * pCurrentFrame) const;
};

//A set of normally distributed observations
class CVariablesFromMultiDistnsObserved
{
	class CMeanInvVarPair
	{
		double dMean, dVar_inv;
	public:
		CMeanInvVarPair(double dMean, double dVar) : dMean(dMean), dVar_inv(1.0/dVar) {}

		double mean() const { return dMean; }
		double var_inv() const { return dVar_inv; }
	};

	typedef CDynArray<CMeanInvVarPair> TNormalDistns;
	TNormalDistns vNormalDistns;
public:
	void addDist(double dMean, double dVar)
	{
		if(dMean > 0 && dVar > 0)
		{
			//breakPoint();
			vNormalDistns.push_back(CMeanInvVarPair(dMean, dVar));
			/*if(dVar<20)
				cout << "Observed length " << dMean << " Var=" << dVar << endl;*/
		}
	}

	//I think MLE of n distn given n observations (repeated below)
	void computeMeanVar(double & dMean, double & dVar)
	{
		double dVar_inv = 0;
		dMean = dVar = 0;
		if(vNormalDistns.size() == 0)
			return;

		for(TNormalDistns::const_iterator pDist = vNormalDistns.begin(); pDist != vNormalDistns.end(); pDist++)
		{
			dVar_inv += pDist->var_inv();
			dMean += pDist->mean() * pDist->var_inv();
			//cout << "Using length " << pDist->mean() << " Var=" << 1.0/pDist->var_inv() << endl;
		}

		dVar = 1.0/dVar_inv;
		dMean *= dVar;

		//cout << "Calculated length " << dMean << " Var=" << dVar << " from " << vNormalDistns.size() << " observations" << endl;

		if(IS_DEBUG) CHECK(!(dMean > 0 && dMean < HUGE && dVar > 0 && dVar < HUGE), "Combining distns failed");
	}

	void reset()
	{
		vNormalDistns.clear();
	}

	//I think MLE of 2 distn given 2 observations (repeated above)
	static void combineTwoObservations(const double dMean1, const double dVar1, const double dMean2, const double dVar2, double & dMean, double & dVar)
	{
		if(dVar1==0 /*|| dMean2 == 0*/)
		{
			dMean = dMean1;
			dVar = dVar1;
			return;
		}
		else if(dVar2==0 /*|| dMean1 == 0*/)
		{
			dMean = dMean2;
			dVar = dVar2;
			return;
		}

		if(IS_DEBUG) CHECK(!(dVar1*dVar2 > 0), "combineTwoObservations: All bad distns");

		double dV1_inv = 1.0/dVar1, dV2_inv = 1.0/dVar2;

		dMean = dMean1*dV1_inv+dMean2*dV2_inv;

		double dVar_inv = dV1_inv + dV2_inv;

		dVar = 1.0/dVar_inv;
		dMean *= dVar;

		if(IS_DEBUG) CHECK(!(dMean > -HUGE && dMean < HUGE && dVar > 0 && dVar < HUGE), "Combining distns failed");
	}
};

//The proper scale estimation from the paper.
class CScaleFromObservingAndMeasuringKnownObjs
{
	class CMeanAndInvSdAndMeasurement
	{
		double dMeasurement, dMean, dVar_inv;
	public:
		CMeanAndInvSdAndMeasurement(double dMeasurement, double dMean, double dVar) : dMeasurement(dMeasurement), dMean(dMean), dVar_inv(1.0/dVar) {}

		double measurement() const { return dMeasurement; } //Baseline 1 new measurement
		double mean() const { return dMean; } //Object size measurement
		//double sd_inv() const { return dSD_inv; }
		double var_inv() const { return sqr(dVar_inv); } //Object class var
	};

	typedef CDynArray<CMeanAndInvSdAndMeasurement> TNormalDistns;
	TNormalDistns vNormalDistns;
public:
	void addObjectMeasurement(double dMeasurement, double dMean, double dVar)
	{
		if(dMean > -HUGE && dVar > 0)
		{
			CHECK(dVar > HUGE, "Log-variance too big in addObjectMeasurement")

			vNormalDistns.push_back(CMeanAndInvSdAndMeasurement(dMeasurement, dMean, dVar));
			//if(dVar<20)
				//cout << "Observed length " << dMean << " Var=" << dVar << endl;
		}
	}

	void computeMeanVar(double & dScale, double & dVar)
	{
		dScale = dVar = 0;
		if(vNormalDistns.size() == 0)
		{
			dScale = dVar = CEdge::NO_SCALE();
			return;
		}

		double dVar_denom = 0 ;
		for(TNormalDistns::const_iterator pDist = vNormalDistns.begin(); pDist != vNormalDistns.end(); pDist++)
		{
			dVar_denom += pDist->var_inv();
			dScale += (pDist->mean() - pDist->measurement()) * pDist->var_inv();

		}

		if(dVar_denom==0)
			dScale = dVar = CEdge::NO_SCALE();
		else
		{
			dVar = 1.0/dVar_denom;
			dScale *= dVar;

			std::cout << "Properly calculated scale " << dScale << " Var=" << dVar << " from " << vNormalDistns.size() << " observations" << std::endl;
			if(IS_DEBUG) CHECK(!(dScale > CEdge::NO_SCALE() && dScale < HUGE && dVar > 0 && dVar < HUGE), "Combining distns failed");
		}
	}
	/*void computeMeanVar(double & dScale, double & dVar)
	{
		double dVar_denom = 0 ;
		dScale = dVar = 0;
		if(vNormalDistns.size() == 0)
		{
			dScale = dVar = CEdge::NO_SCALE();
			return;
		}

		for(TNormalDistns::const_iterator pDist = vNormalDistns.begin(); pDist != vNormalDistns.end(); pDist++)
		{
			const double x_on_sigma2 = pDist->measurement() * pDist->var_inv();
			const double x2_on_sigma2 = pDist->measurement() * x_on_sigma2;
			dVar_denom += x2_on_sigma2;

			dScale += pDist->mean() * x_on_sigma2;

			cout << "Using obj length " << pDist->mean() << " observed length " << pDist->measurement() << " Var=" << 1.0/pDist->var_inv() << endl;
		}

		if(dVar_denom==0)
			dScale = dVar = CEdge::NO_SCALE();
		else
		{
			dVar = 1.0/dVar_denom;
			dScale *= dVar;

			std::cout << "Properly calculated scale " << dScale << " Var=" << dVar << " from " << vNormalDistns.size() << " observations" << std::endl;
			if(IS_DEBUG) CHECK(!(dScale > CEdge::NO_SCALE() && dScale < HUGE && dVar > 0 && dVar < HUGE), "Combining distns failed");
		}
	}*/

	void reset()
	{
		vNormalDistns.clear();
	}

	//I think MLE of 2 distn given 2 observations (repeated above)
	/*static void combineTwoObservations(const double dMean1, const double dVar1, const double dMean2, const double dVar2, double & dMean, double & dVar)
	{
		if(dVar1==0 || dMean2 == 0)
		{
			dMean = dMean1;
			dVar = dVar1;
			return;
		}
		else if(dVar2==0 || dMean1 == 0)
		{
			dMean = dMean2;
			dVar = dVar2;
			return;
		}

		if(IS_DEBUG) CHECK(!(dMean1*dVar1*dMean2*dVar2 > 0), "combineTwoObservations: All bad distns");

		double dV1_inv = 1.0/dVar1, dV2_inv = 1.0/dVar2;

		dMean = dMean1*dV1_inv+dMean2*dV2_inv;

		double dVar_inv = dV1_inv + dV2_inv;

		dVar = 1.0/dVar_inv;
		dMean *= dVar;

		if(IS_DEBUG) CHECK(!(dMean > 0 && dMean < HUGE && dVar > 0 && dVar < HUGE), "Combining distns failed");
	}*/
};

/*class CSCOREScale
{
	static double UNINIT() { return -1; }
	double dLNMean, dLNVar;
public:
	CSCOREScale(double dLNMean, double dLNVar) : dLNMean(dLNMean), dLNVar(dLNVar)
	{
		if(IS_DEBUG) CHECK(dLNMean<=0 || dLNVar <= 0, "CSCOREScale: Bad scale mean/var");
	}
	CSCOREScale() : dMean(UNINIT()), dVar(UNINIT()) {}

	double LNmean() const {
		if(IS_DEBUG) CHECK(dLNMean<=0 || dLNVar <= 0, "CSCOREScale: Bad scale mean/var");
		return dLNMean;
	}

	double LNvar() const
	{
		if(IS_DEBUG) CHECK(dLNMean<=0 || dLNVar <= 0, "CSCOREScale: Bad scale mean/var");
		return dLNVar;
	}

	bool hasScale() const { return dLNVar != UNINIT(); }

	void setUninit()
	{
		dLNMean = UNINIT();
		dLNVar = UNINIT();
	}
};*/

#endif /* EDGESCALEOBSERVER_H_ */
