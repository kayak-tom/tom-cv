#include "svdScaleOptimiser.h"

#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/LU>
#include "util/dynArray.h"
#include "util/set2.h"
#include "util/opencv.h"

#ifdef USE_SVDLIB
extern "C" {
#include "../../SVDLIBC/svdlib.h"
}
#endif

using namespace std;
using namespace Eigen;

typedef Eigen::Matrix<double, Dynamic, Dynamic, RowMajor + AutoAlign> TSMatrix;
typedef Eigen::Matrix<double, Dynamic, 1, AutoAlign> TDMatrix;

//Also duplicated in scaleOptimiser.cpp
void svd(TSMatrix & A/*, Eigen::Matrix<double, Dynamic, Dynamic, RowMajor + AutoAlign> & U*/, TDMatrix & D, TSMatrix & V)
{
	CvMat cvA = cvMat(A.rows(), A.cols(), CV_64F, A.data());
	//CvMat cvU = cvMat(U.rows(), U.cols(), CV_64F, U.data());
	CvMat cvV = cvMat(V.rows(), V.cols(), CV_64F, V.data());
	CvMat cvD = cvMat(D.rows(), D.cols(), CV_64F, D.data());
	cvSVD(&cvA, &cvD, /*&cvU*/ 0, &cvV, CV_SVD_MODIFY_A);
}

void CSVDScaleOptimiser::run(const bool bVerbose)
{
	//First add an eqn for every edge:
	int nLastId = nEdges;

	//for(TEdgeIds::const_iterator pEdgeIds = aEdgeIds.begin(); pEdgeIds != aEdgeIds.end(); pEdgeIds++)
	int nId = 0;
	for(CDynArray<NSLAMMap::CEdge_relScale *>::iterator ppEdge = aEdgesInOrder.begin(); ppEdge != aEdgesInOrder.end(); ppEdge++, nId++)
	{
		NSLAMMap::CEdge_relScale * pEdge = *ppEdge;
		aMatrixRows.push_back(TEquation());
		TEquation & equation = aMatrixRows.back();
		pEdge->resetOptimisedScale();
		CRelScale relScale = pEdge->getRelScale();
		const double d = relScale.get_d();
		const double g_sq = relScale.get_g_sq();
		const double MAKE_CYCLIC_CONSTRAINTS_IMPORTANT = 0.001; //Also keeps values approx 1
		const double g_inv=MAKE_CYCLIC_CONSTRAINTS_IMPORTANT*1.0/sqrt(g_sq);
		//Want to minimise (d*-d)^2/g_sq
		equation.init(nId, g_inv);
		equation.init(nLastId, -d*g_inv);
	}

	static size_t edgeHashLastTime = 0;

	if(edgeHash == edgeHashLastTime)
	{
		if(bVerbose) cout << "Same optimisation as last time\n";
		return;
	}
	edgeHashLastTime = edgeHash;

	int nRows = aMatrixRows.size();
	int nCols = nLastId+1;

	if(bVerbose) cout << "Starting " << nRows << " by " << nCols << " svd...\n\n";// << S << endl;

	CStopWatch s; s.startTimer();

	const bool USE_EIGEN_OCV = true;
	if(USE_EIGEN_OCV)
	{
		TSMatrix S(nRows, nCols);
		TDMatrix D(nCols, 1);
		TSMatrix V(nCols, nCols);

		S.setZero();
		int nRow = 0;
		for(TEquationArray::iterator pEqn = aMatrixRows.begin(); pEqn != aMatrixRows.end(); pEqn++, nRow++)
		{
			for(TEquation::iterator pEl = pEqn->begin(); pEl != pEqn->end(); pEl++ )
				S(nRow, pEl->first) = pEl->second;
		}

		svd(S, D, V);

		TDMatrix new_d = V.col(nCols-1);

		new_d /= new_d(nCols-1);

		if(bVerbose)
		{
			cout << "New Scales: " << new_d.transpose() << endl;
		}

		nId = 0;
		for(CDynArray<NSLAMMap::CEdge_relScale *>::iterator ppEdge = aEdgesInOrder.begin(); ppEdge != aEdgesInOrder.end(); ppEdge++, nId++)
		{
			NSLAMMap::CEdge_relScale * pEdge = *ppEdge;
			double dNew_d = new_d(nId);
			pEdge->setOptimisedScale(dNew_d, bVerbose);
		}
	}
	else
	{
#ifndef USE_SVDLIB
		THROW("The svdlib dependency has been temporarily removed")
#else
		/*{
			SVDVerbosity = 100;
			smat * pMat = svdNewSMat(2, 2, 2);
			pMat->rowind[0] = 0;
			pMat->rowind[1] = 1;
			pMat->pointr[0] = 0;
			pMat->pointr[1] = 1;
			pMat->pointr[2] = 2;
			pMat->value[0] = 1;
			pMat->value[1] = 1e-4;

			double end[2] = {-1,-1};// {-1.0e-30, 1.0e-30};
			double kappa = 1e-30; //1e-6 default
			int nIterations = 10;
			SVDRec pSvdRes = svdLAS2(pMat, 0, nIterations, end,
	                      				kappa);
	        DMat pVt = pSvdRes->Ut; //want the last row of this

	        cout << pVt->value[0][0] << ',' << pVt->value[0][1] << endl;
	        cout << pVt->value[1][0] << ',' << pVt->value[1][1] << endl;

	        cout << pSvdRes->S[0] << '-' << pSvdRes->S[1] << endl;

		}*/


		int nElements = 0;//nLastId*2;
		for(TEquationArray::iterator pEqn = aMatrixRows.begin(); pEqn != aMatrixRows.end(); pEqn++)
			nElements += pEqn->size();
		cout << nElements << " nz elements\n";

		smat * pMat = svdNewSMat(nCols, nRows, nElements); //Fill-in TRANSPOSED

		int nRow = 0, nEl = 0;
		const double SCALE_UP = 10000;
		for(TEquationArray::iterator pEqn = aMatrixRows.begin(); pEqn != aMatrixRows.end(); pEqn++, nRow++)
		{
			pMat->pointr[nRow] = nEl; //index of this row start in rowind
			for(TEquation::iterator pEl = pEqn->begin(); pEl != pEqn->end(); pEl++, nEl++ )
			{
				pMat->rowind[nEl] = pEl->first;
				pMat->value[nEl] = pEl->second * SCALE_UP; //Scale-up to keep SVs big enough not to be discarded
			}
		}
		pMat->pointr[nRow] = nEl; //last element

		if(IS_DEBUG) CHECK(nElements != nEl, "Element count mismatch");

		{
			DMat pDense = svdConvertStoD(pMat);
			for(int nR=0; nR<nCols; nR++, cout << endl)
				for(int nC=0;nC<nRows; nC++)
					cout << pDense->value[nR][nC] << ' ';
			svdFreeDMat(pDense);
		}

		double end[2] = {-1,-1};// {-1.0e-30, 1.0e-30};
		double kappa = 1e-30; //1e-6 default
		int nIterations = 10;
		SVDRec pSvdRes = svdLAS2(pMat, 0 /*nCols or 0?*/, nIterations, end,
                      				kappa);
        DMat pVt = pSvdRes->Ut; //want the last row of this

		double dLastEl_inv = 1.0 / pVt->value[nCols-1][nCols-1];

		if(bVerbose)
		{
			/*cout << "Vt: \n";
			for(int nR=0; nR<nCols; nR++, cout << endl)
				for(int nC=0;nC<nCols; nC++)
					cout << pVt->value[nR][nC] << ' ';

			cout << "S: \n";
			for(int nC=0;nC<nCols; nC++)
				cout << pSvdRes->S[nC] << ' ';*/

			cout << "\n\nNew Scales: ";
			for(int i=0;i<nCols-1; i++)
				cout << pVt->value[nCols-1][i] * dLastEl_inv << ' ';

			cout << endl;
		}

		nId = 0;
		double dTotalScale = 0;
		for(CDynArray<NSLAMMap::CEdge_relScale *>::iterator ppEdge = aEdgesInOrder.begin(); ppEdge != aEdgesInOrder.end(); ppEdge++, nId++)
		{
			NSLAMMap::CEdge_relScale * pEdge = *ppEdge;
			double dNew_d = pVt->value[nCols-1][nId] * dLastEl_inv;
			pEdge->setOptimisedScale(dNew_d, bVerbose);
			dTotalScale += fabs(dNew_d);
		}

		svdFreeSMat(pMat);
		svdFreeSVDRec(pSvdRes);

		CHECK(dTotalScale == 0, "SVD scale optimisation failed");
#endif
	}

	s.stopTimer();

	if(bVerbose)
		cout << s.getElapsedTime() << " seconds for SVD\n";
}

void CSVDScaleOptimiser::addCycle(NSLAMMap::CCycle & c, const bool bVerbose)
{
	aMatrixRows.push_back(TEquation());
	TEquation & equation = aMatrixRows.back();

	for(NSLAMMap::CCycle::iterator ppEdge = c.begin(); ppEdge != c.end(); ppEdge++)
	{
	 	NSLAMMap::CEdge_relScale * pEdge = ppEdge->second;

	 	edgeHash ^= (size_t)(void *)pEdge;

		int nNewId = 0;
		TEdgeIds::iterator pEdgeId = aEdgeIds.find(pEdge);
		if(pEdgeId == aEdgeIds.end())
		{
			if( aEdgeIds.size()>0 )
				nNewId = nEdges;
			nEdges++;

			aEdgeIds.init(pEdge, nNewId);
			aEdgesInOrder.push_back(pEdge);
			if(IS_DEBUG) CHECK(nEdges != aEdgesInOrder.size(), "Error saving edges in order");

			if(bVerbose) cout << "New id " << nNewId;
		}
		else
		{
			nNewId = pEdgeId->second;
			if(bVerbose) cout << "Existing";
		}

		if(bVerbose) cout << " cycle edge from " << pEdge->firstNode()->id1 << ',' << pEdge->firstNode()->id2 << " to " << pEdge->secondNode()->id1 << ',' << pEdge->secondNode()->id2 << "..." << nNewId << std::endl;

		equation.init(nNewId, ppEdge->first);
	}
}

CSVDScaleOptimiser::~CSVDScaleOptimiser() {}

typedef CDynArray<int> TEdgeDirections; // 0 = not in cycle i, 1 = forwards, -1 = backwards
class CNewEdge
{
	TEdgeDirections dir;
	const int nId;

	class CEdgeDirPair
	{
		NSLAMMap::CEdge_relScale * pEdge;
		int nDir;
	public:
		NSLAMMap::CEdge_relScale * edge() const { return pEdge; }
		int dir() const { return nDir; }

		CEdgeDirPair(NSLAMMap::CEdge_relScale * pEdge, const int nDir) : pEdge(pEdge), nDir(nDir)
		{
			if(IS_DEBUG) CHECK(nDir == 0, "Bad zero direction")
		}
		CEdgeDirPair() : pEdge(0), nDir(0) {}
	};

	typedef CDynArray<CEdgeDirPair> TEdges;
	TEdges edges;
	double dMyG_sq, dMyD;
public:
	CNewEdge(const TEdgeDirections & dir_in, const int nId) : nId(nId), dMyG_sq(0), dMyD(0)
	{
		edges.reserve(100);
		dir.copy_back(dir_in.begin(), dir_in.end());
	}

	static const int NO_ID = -1; //for temporary edge...
	CNewEdge(const TEdgeDirections & dir_in) : nId(NO_ID), dMyG_sq(0), dMyD(0)
	{
		dir.copy_back(dir_in.begin(), dir_in.end());
	}

	bool operator<(const CNewEdge & otherCycle) const
	{
		for(TEdgeDirections::const_iterator pDir = dir.begin(), pDir2=otherCycle.dir.begin(); pDir < dir.end(); pDir++, pDir2++)
		{
			bool b1=*pDir, b2=*pDir2;
			if(b1 < b2)
				return true;
			else if(b1 > b2)
				return false;
		}
		return false;
	}

	void addOldEdge(NSLAMMap::CEdge_relScale * pOldEdge, int nDir)
	{
		edges.push_back(CEdgeDirPair(pOldEdge, nDir));
		dMyG_sq += pOldEdge->getRelScale().get_g_sq();
		dMyD += nDir * pOldEdge->getRelScale().get_d();
	}

	double getD() const { return dMyD; }
	double getG_sq() const { return dMyG_sq; }

	void setD(double dNewD, bool bVerbose = false)
	{
		if(bVerbose)
			cout << "Updating " << dMyD << " to " << dNewD << endl;

		double dUpdatePerG_sq = (dNewD - dMyD) / dMyG_sq;
		double dDTotal = 0;
		for(TEdges::iterator pEdgePair = edges.begin(); pEdgePair != edges.end(); pEdgePair++)
		{
			NSLAMMap::CEdge_relScale * pEdge = pEdgePair->edge();
			double dDOld = pEdge->getRelScale().get_d();
			double dG_sqOld = pEdge->getRelScale().get_g_sq();
			double dDNew = dDOld + pEdgePair->dir() * dUpdatePerG_sq*dG_sqOld;
			pEdge->setOptimisedScale(dDNew, bVerbose);
			dDTotal += pEdgePair->dir() * dDNew; //for debug
			if(bVerbose)
				cout << "Updating scale from " << dDOld << " to " << dDNew << endl;
		}

		if(bVerbose)
			cout << "Updated total D from " << dMyD << " to " << dDTotal << endl;

		if(IS_DEBUG) CHECK(!zero(dDTotal - dNewD), "Error distributing scale error around cycle");
	}

	class CNewEdgePred
	{
	public:
		bool operator()(const CNewEdge * pE1, const CNewEdge * pE2) const
		{
			return (*pE1 < *pE2);
		}
	};
};

#define NOT_FINISHED
#ifdef NOT_FINISHED

void CCollapsedSVDScaleOptimiser::run(const bool bVerbose)
{
	//First add an eqn for every edge:

	static size_t edgeHashLastTime = 0;

	if(edgeHash == edgeHashLastTime)
	{
		if(bVerbose) cout << "Same optimisation as last time\n";
		return;
	}
	edgeHashLastTime = edgeHash;

	CStopWatch s; s.startTimer();

	//First copy every cycle, and add newEdges as they are created.
	const int nCycles = aMatrixRows.size();
	TEquationArray aNewCycles(nCycles);

	//For each edge see which cycles it is in. Add to a newedge (possibly empty)
	typedef set2_NF<CNewEdge *, CNewEdge::CNewEdgePred> TNewEdgeSet;
	TNewEdgeSet aNewEdges;

	CDynArray<CNewEdge *> aNewEdgesInOrder; aNewEdgesInOrder.reserve(nCycles*2);

	int nEdgeId = 0;
	for(TEdgeIds::const_iterator pEdgeIds = aEdgeIds.begin(); pEdgeIds != aEdgeIds.end(); pEdgeIds++)
	{
		pEdgeIds->first->resetOptimisedScale();

		const int nEdge = pEdgeIds->second;
		TEdgeDirections aCycles(nCycles, 0);
		int nCycle=0, nFirstCycle=-1;;
		for(TEquationArray::iterator pEqn = aMatrixRows.begin(); pEqn != aMatrixRows.end(); pEqn++, nCycle++)
		{
			//Does this cycle contain this edge?
			TEquation::const_iterator pEdgeAndDir = pEqn->find(nEdge);
			if(pEdgeAndDir != pEqn->end())
			{
				//cout << "Cycle " << nCycle << " contains edge " << nEdge << endl;
				aCycles[nCycle] = pEdgeAndDir->second;
				if(nFirstCycle == -1)
					nFirstCycle = nCycle;
			}
		}
		CNewEdge wrapper(aCycles);
		TNewEdgeSet::iterator pNewEdgeIter = aNewEdges.find(&wrapper);
		if(pNewEdgeIter == aNewEdges.end()) //A new edge is a new combination of cycles
		{
			for(nCycle = 0; nCycle < nCycles; nCycle++)
				if(aCycles[nCycle] != 0)
					aNewCycles[nCycle].init(nEdgeId, aCycles[nFirstCycle] * aCycles[nCycle]); //Dir may vary because new edges can be traversed in different dirs

			pNewEdgeIter = aNewEdges.insert(new CNewEdge(aCycles, nEdgeId++)).first;
			aNewEdgesInOrder.push_back(*pNewEdgeIter);
		}

		CNewEdge * pNewEdge = *pNewEdgeIter; //This edge is defined to be traversed in a forward direction in its first cycle
		pNewEdge->addOldEdge(pEdgeIds->first, aCycles[nFirstCycle]);

	}

	//First add an eqn for every edge:
	int nLastId = nEdgeId;

	//for(TEdgeIds::const_iterator pEdgeIds = aEdgeIds.begin(); pEdgeIds != aEdgeIds.end(); pEdgeIds++)
	int nId = 0;
	for(CDynArray<CNewEdge *>::iterator ppEdge = aNewEdgesInOrder.begin(); ppEdge != aNewEdgesInOrder.end(); ppEdge++, nId++)
	{
		CNewEdge * pEdge = *ppEdge;
		aNewCycles.push_back(TEquation());
		TEquation & equation = aNewCycles.back();
		const double d = pEdge->getD();
		const double g_sq = pEdge->getG_sq();
		const double MAKE_CYCLIC_CONSTRAINTS_IMPORTANT = 0.05; //Also keeps values approx 1
		const double g_inv=MAKE_CYCLIC_CONSTRAINTS_IMPORTANT*1.0/sqrt(g_sq);
		//Want to minimise (d*-d)^2/g_sq
		equation.init(nId, g_inv);
		equation.init(nLastId, -d*g_inv);
	}

	int nRows = aNewCycles.size();
	int nCols = nLastId+1;

	if(bVerbose)
		cout << "Starting " << nRows << " by " << nCols << " reduced svd...\n\n";// << S << endl;

	TSMatrix S(nRows, nCols);
	TDMatrix D(nCols, 1);
	TSMatrix V(nCols, nCols);

	S.setZero();
	int nRow = 0;
	for(TEquationArray::iterator pEqn = aNewCycles.begin(); pEqn != aNewCycles.end(); pEqn++, nRow++)
	{
		for(TEquation::iterator pEl = pEqn->begin(); pEl != pEqn->end(); pEl++ )
			S(nRow, pEl->first) = pEl->second;
	}

	if(bVerbose)
		cout << S << endl;

	svd(S, D, V);

	TDMatrix new_d = V.col(nCols-1);

	new_d /= new_d(nCols-1);

	if(bVerbose)
	{
		cout << "New Scales: " << new_d.transpose() << endl;
	}

	//residuals around cycles before:
	if(bVerbose)
	{
		int nCycle=0;
		for(TEquationArray::iterator pEqn = aMatrixRows.begin(); pEqn != aMatrixRows.end(); pEqn++, nCycle++)
		{
			//const TEquation & eqn = *pEqn;
			double dTotalLogScale = 0;
			for(TEquation::iterator pEl = pEqn->begin(); pEl != pEqn->end(); pEl++ )
			{
				NSLAMMap::CEdge_relScale * pEdge = aEdgesInOrder[pEl->first];
				dTotalLogScale += pEl->second * pEdge->getRelScale().get_d();
			}
			cout << "Cycle " << nCycle << " has resid " << dTotalLogScale << endl;
		}
	}

	nId = 0;
	for(CDynArray<CNewEdge *>::iterator ppEdge = aNewEdgesInOrder.begin(); ppEdge != aNewEdgesInOrder.end(); ppEdge++, nId++)
	{
		CNewEdge * pNewEdge = *ppEdge;
		double dNew_d = new_d(nId);
		pNewEdge->setD(dNew_d, bVerbose);
	}

	s.stopTimer();

	//Now output residuals around cycles:
	if(true || bVerbose)
	{
		int nCycle=0;
		for(TEquationArray::iterator pEqn = aMatrixRows.begin(); pEqn != aMatrixRows.end(); pEqn++, nCycle++)
		{
			//const TEquation & eqn = *pEqn;
			double dTotalLogScale = 0;
			for(TEquation::iterator pEl = pEqn->begin(); pEl != pEqn->end(); pEl++ )
			{
				NSLAMMap::CEdge_relScale * pEdge = aEdgesInOrder[pEl->first];
				dTotalLogScale += pEl->second * pEdge->getRelScale().get_d();
			}
			cout << "Cycle " << nCycle << " now has resid " << dTotalLogScale << endl;
		}
	}

	//if(bVerbose)
		cout << s.getElapsedTime() << " seconds for collapsed SVD\n";

	cout << "Now correcting rotation around the same cycles...TODO x y z seperately??\n";

}
#endif
