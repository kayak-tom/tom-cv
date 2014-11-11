#pragma once
#include "slamMap.h"

class CSVDScaleOptimiser : public CScaleOptimiser
{
protected:
	class CEdgeSort2
	{
	public:
		bool operator()(const NSLAMMap::CEdge_relScale * e1, const NSLAMMap::CEdge_relScale * e2) const
		{
			//Todo: Sort to make matrix nice!
			return e1 < e2;
			//return e1->firstNode() < e2->firstNode() || (e1->firstNode() == e2->firstNode() && e1->secondNode() < e2->secondNode() );
		}
	};

	typedef map2< NSLAMMap::CEdge_relScale *, int, CEdgeSort2 > TEdgeIds;
	TEdgeIds aEdgeIds; //Map any edge to an index
	CDynArray<NSLAMMap::CEdge_relScale *> aEdgesInOrder;

	/*class CCoeffEdgePair
	{
		double dCoeff;
		NSLAMMap::CEdge_relScale * pEdge;
	public:
		CCoeffEdgePair(double dCoeff, NSLAMMap::CEdge_relScale * pEdge) : dCoeff(dCoeff), pEdge(pEdge) {}
		inline double coeff() const { return dCoeff; }
		inline NSLAMMap::CEdge_relScale * edge() const { return pEdge; }
	};*/

	typedef map2<int, double> TEquation; //Set of INDICES and coefficients. These are the elements in S
	typedef CDynArray< TEquation > TEquationArray;
	TEquationArray aMatrixRows;
	int nEdges;
	size_t edgeHash;
public:
	virtual void run(const bool bVerbose);
	virtual void addCycle(NSLAMMap::CCycle & c, const bool bVerbose);
	CSVDScaleOptimiser() : nEdges(0), edgeHash(0) {}
	virtual ~CSVDScaleOptimiser();
};

class CCollapsedSVDScaleOptimiser : public CSVDScaleOptimiser
{
public:
	virtual void run(const bool bVerbose) HOT;
};
