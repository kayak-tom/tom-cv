/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

//============================================================================
// Name        : ransac.cpp
// Author      : Tom Botterill
// Version     :
// Copyright   : Free
// Description : RANSAC/BaySAC/SimSAC Essential Matrix estimation and refinement.
//               Supports a choice of samplers, hypothesis generators, hypothesis refiners, termination conditions.
//               Works with disjoint correspondences (N-M)
//               Also includes PROSAC, Guided-MLESAC sampling, WaldSAC
//============================================================================

#include <iostream>
#include "ransac.h"
#include <boost/thread/mutex.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/bind.hpp>
#include <boost/thread/thread.hpp>
#include <boost/function.hpp>
#include <fstream>

using namespace std;

void printModel(const CModel & model_in) {
    const C3x3MatModel * pModel = dynamic_cast<const C3x3MatModel *> (&model_in);
    if (pModel) {
        const C3x3MatModel & model = *pModel;
        cout << endl;
        for (int r = 0; r < 3; r++) {
            for (int c = 0; c < 3; c++) {
                cout << model[3 * r + c];
                if (c < 2) cout << ", ";
            }
            cout << endl;
        }
        cout << endl;
    } else
        cout << "Can't print this kind of model" << endl;
}

DEBUGONLY(CDynArray<int> vChosen;) //hacky...

typedef boost::function<void( const CModel &, const int, const int, CMask &, double *, int &) > TInlierCounterFn;

int MODELS; //Hack to enable the collection of stats about model counts

void testOneHypothesisSetMT(
        CSampler * pSampler, //Does not need to be TS
        CModelHypothesiser * pHypothesise, //Must be TS
        CIterTerminator * pIterTerminator, // Does NOT need to be TS
        TInlierCounterFn & inlierCounter, //Does need to be TS
        CFindBestMatching * pRefineMatches, //Does NOT need to be TS
        CModel & bestModel, CMask & bestMask, const int nThread, boost::mutex & mxUpdateBGC) {
    const int nHypSetSize = pSampler->hypSetSize();
    const int nNumPoints = pSampler->numPoints();

    ////////START MT

    TSubSet anSample(nHypSetSize);
    //		CHECK(5 != nHypSetSize, "Windows array size mismatch");
    //        int anSample[5];

    {
        boost::unique_lock<boost::mutex> lock(mxUpdateBGC); //Could use a different mutex here
        if (!pSampler->choose(anSample)) return;
    }

#ifdef _DEBUG
    if (nNumPoints > 10 * nHypSetSize) {
        boost::unique_lock<boost::mutex> lock(mxUpdateBGC);

        bool bDumpOut = false;

        for (int i = 0; i < nHypSetSize; i++)
            if (vChosen.count(anSample[i]) > 2) {
                cout << "Warning: RANSAC sampler choosing same data point repeatedly\n";
                bDumpOut = true;
            }

        if (vChosen.size() >= nHypSetSize * 5)
            vChosen.popN(nHypSetSize);

        for (int i = 0; i < nHypSetSize; i++)
            vChosen.push_back(anSample[i]);

        if (bDumpOut) {
            for (int i = 0; i < nHypSetSize; i++)
                cout << anSample[i] << ' ';
            cout << endl;
        }
    }
#endif

    /*bool bOdd = false;
    for(int i=0; i< nHypSetSize; i++)
    {
            cout << anSample[i] << ' ';
            if(anSample[i] % 2 == 1) bOdd = true;
    }
    cout << "=sample\n";
    bOdd || cout << "INLIER SET\n";
     */


    boost::scoped_ptr<const CModels> pModels(pHypothesise->getModels(anSample));

    MODELS += pModels->numModels();

    CMask mask(nNumPoints);

    for (int nModel = 0; nModel < pModels->numModels(); nModel++) {
        if (nModel > 0)
            mask.setZero();

        const CModel & model = pModels->getData(nModel);

        ARRAY(double, adResiduals, nNumPoints);
        double *pdResiduals = 0;
        if (pRefineMatches->supplyResiduals())
            pdResiduals = PTR(adResiduals);

        int nNewGC = 0;
        inlierCounter(model, nThread, pIterTerminator->BGC(), mask, pdResiduals, nNewGC);

        if (pIterTerminator->BGC() < nNewGC) {
            //Prob need to lock here, but not for pIterTerminator. TODO: Reduce num of locks
            boost::unique_lock<boost::mutex> lock(mxUpdateBGC); //Could use a different mutex here
            pRefineMatches->refine(mask, nNewGC, pdResiduals);

            if (pIterTerminator->updateBGC(nNewGC)) {
                model.copyInto(bestModel);
                mask.copyInto(bestMask);
            }
        }
    }
    ////////END MT
}
//void f(int, int, int, int, int, int, int, int, int);

int doRansac(
        CSampler * pSampler,
        CModelHypothesiser * pHypothesise,
        CModelRefiner * pRefine, //Does not have to be TS
        CRansacTerminator * pTerminator,
        CIterTerminator * pIterTerminator, //Does not have to be TS
        CInlierCounter * pCounter,
        CFindBestMatching * pRefineMatches,
        CModel & bestModel, CMask & bestMask, int nThreads, const bool bVerbose) {
    //const int nNumPoints = bestMask.size();
    if(IS_DEBUG) CHECK(!pSampler || !pHypothesise || !pRefine || !pTerminator || !pIterTerminator || !pCounter || !pRefineMatches, "getE: Uninitialised param");
    if(IS_DEBUG) CHECK(nThreads < 1, "getE: Bad number of threads");
    //if(IS_DEBUG) CHECK(p1.size() != p2.size() || nNumPoints < pRefine->minNumPoints(), "getE: Bad number of points");
    if(IS_DEBUG) CHECK(pSampler->numPoints() != bestMask.size(), "getE: Mask size doesn't match points");

    DEBUGONLY(vChosen.clear();) //Todo member...

    if (nThreads > 1) {
        REPEAT(1, cout << "MT RANSAC disabled (is not generally worthwhile anyway).\n");
        nThreads = 1;
    }

    if (nThreads > 1) {
        CHECK(!(pHypothesise->isThreadsafe()), "getE: For multithreaded operation hypothesiser must be threadsafe");
        CHECK(!(pTerminator->isThreadsafe()), "getE: For multithreaded operation terminator must be threadsafe");
    }

    //at most 9 args to boost::bind
    boost::mutex mxUpdateBGC;
    TInlierCounterFn countInliers = boost::bind(&CInlierCounter::countInliers, pCounter, _1, pTerminator, _2, _3, _4, _5, _6, 1.0);

    CDynArrayOwner<boost::thread> apThreads(nThreads - 1);

    do {
        apThreads.clear();

        for (int nThread = 0; nThread < nThreads - 1; nThread++) {
            apThreads.push_back(new boost::thread(
                    boost::bind(testOneHypothesisSetMT,
                    pSampler,
                    pHypothesise,
                    pIterTerminator,
                    boost::ref(countInliers),
                    pRefineMatches,
                    boost::ref(bestModel),
                    boost::ref(bestMask),
                    nThread,
                    boost::ref(mxUpdateBGC))));
        }

        //Use this thread too...
        testOneHypothesisSetMT(
                pSampler,
                pHypothesise,
                pIterTerminator,
                countInliers,
                pRefineMatches,
                bestModel, bestMask, 0, mxUpdateBGC);

        for (int nThread = 0; nThread < nThreads - 1; nThread++) {
            apThreads[nThread]->join();
        }
    } while (!pIterTerminator->terminate(nThreads)); //Tell terminator how many samples we have tried, possibly terminate.

    /*static int s_nItersTotal = 0, s_nInliers = 0, s_nWrongInliers = 0, s_nInliersMissed = 0;
    s_nItersTotal += pIterTerminator->numIters();
    cout << s_nItersTotal << " iterations in total\n";

    s_nInliers += bestMask.countInliers();
    cout << s_nInliers << " inliers in total\n";

    for(int i=0; i<bestMask.size(); i++)
    {
            if(i % 2 == 0)
            {
                    if(!bestMask[i])
                            s_nInliersMissed ++;
            }
            else
            {
                    s_nWrongInliers += bestMask[i];
            }
    }

    cout << s_nInliersMissed << " inliers missed\n";
    cout << s_nWrongInliers << " outliers in solution\n";

    ofstream inlierFile("inFile", ios_base::app);
    inlierFile << "Inliers missed-outliers found\t" << s_nInliersMissed << '\t' << s_nWrongInliers << endl;
    inlierFile.close();*/

    if (bVerbose)
        cout << "Terminated after " << pIterTerminator->numIters() << " of " << pIterTerminator->maxIters() << ", BGC=" << pIterTerminator->BGC() << "/" << bestMask.size() << endl;

    if (pIterTerminator->BGC() < pRefine->minNumPoints())
        return pIterTerminator->BGC();

    if (bVerbose) {
        cout << "Model before refinement:";
        printModel(bestModel);
    }
    
    if (pRefine->fitModel(bestMask, bestModel, bVerbose)) {
        if (bVerbose) {
            cout << "Model after refinement: ";
            printModel(bestModel);
        }
    } else {
        if (bVerbose)
            cout << "Refinement failed\n";
        bestMask.setZero();
    }

    /*for (int nTries=0; nTries<2; nTries++) {
        if (pRefine->fitModel(bestMask, bestModel, bVerbose)) {
            if (bVerbose) {
                cout << "Model after refinement: ";
                printModel(bestModel);
            }
        } else {
            if (bVerbose)
                cout << "Refinement failed\n";
            bestMask.setZero();
        }

        if (bVerbose && bestMask.countInliers() < bestMask.size()) {
            cout << "Trying again from all points: ";
            bestMask.setConst(true);
        } 
        else break;
    }*/

    return bestMask.countInliers();
}

//Return models

const CModels * CImCorrModelHypothesiser::getModels(const TSubSet & anHypSet) {
    CModels * pModels;
    if (modelsPerIteration() <= 10)
        pModels = new T3x3MatModels; //Return a ptr to keep threadsafe
    else
        pModels = new TEModels;

    if(IS_DEBUG) CHECK(pModels->numModels() != 0, "CModelHypothesiser::fitModel: Bad out param");
    if(IS_DEBUG) CHECK(nPoints > (int) p1.size(), "CModelHypothesiser::fitModel: Size mismatch");
    if(IS_DEBUG) CHECK(nPoints > (int) p2.size(), "CModelHypothesiser::fitModel: Size mismatch");
    getModels_int(anHypSet, *pModels);

    if (false && IS_DEBUG) {
        static int s_nFails = 0, s_nSuccesses = 0;
        for (int nModel = 0; nModel < pModels->numModels(); nModel++) {
            const C3x3MatModel & m = static_cast<const C3x3MatModel &> (pModels->getData(nModel));
            for (int nPoint = 0; nPoint < anHypSet.size(); nPoint++) {
                if (!CEssentialMatInlierCounter::isInlier_SE<double>(m.asDouble9(), p1[anHypSet[nPoint]], p2[anHypSet[nPoint]], 0, 0.01)) {
                    s_nFails++;
                    if (s_nFails * 3 > s_nSuccesses) //Occasionally failing
                    {
                        cout << "Model fitting failed\n";
                        if (CEssentialMatInlierCounter::isInlier_SE<double>(m.asDouble9(), p2[anHypSet[nPoint]], p1[anHypSet[nPoint]], 0, 0.01))
                            cout << "[Model fitting succeeded with reversed points]\n";
                    }
                } else {
                    s_nSuccesses++;
                    if (s_nFails * 3 > s_nSuccesses) //Occasionally failing
                        cout << "Model fitting succeeded\n";
                }

            }
        }
    }

    return pModels;
}
