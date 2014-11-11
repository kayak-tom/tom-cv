/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * ransacSampler.cpp
 *
 *  Created on: 12/02/2009
 *      Author: tom
 */

#include "ransacSampler.h"
#include <algorithm>
#include <boost/smart_ptr.hpp>

using namespace std;

bool CRANSACSampler::choose(TSubSet & anSample) {
    const int max_random_iters = 200;

    // choose random <n> (=7) points
    for (int i = 0; i < n; i++) {
        int k = 0;
        for (; k < max_random_iters; k++) {
            int idxi = CRandom::Uniform(N); // cvRandInt(&rng) % count;

            int j = 0;
            for (; j < i; j++) {
                int idxj = anSample[j];
                if (idxj == idxi)
                    break;
            }
            if (j == i) //loop completed
            {
                anSample[i] = idxi;
                break;
            }
        }
        //if(IS_DEBUG) CHECK(k >= max_random_iters , "CRANSACSampler::choose: ERROR: Could not find 7 points in max_random_iters iterations");

        if (k >= max_random_iters)
            return false;
    }
    return true;
}

bool CDisjointRANSACSampler::choose(TSubSet & anSample) {
    const int max_random_iters = 200;

    // choose random <n> (=7) points
    for (int i = 0; i < n; i++) {
        int k = 0;
        for (; k < max_random_iters; k++) {
            int idxi = CRandom::Uniform(N); // cvRandInt(&rng) % count;

            int j = 0;
            for (; j < i; j++) {
                int idxj = anSample[j];
                if (idxj == idxi)
                    break;

                //now check if either of the points selected at (earlier) time j are the same as those now:
                if (pointIds.incompatible(idxi, idxj))
                    break;
            }
            if (j == i) //loop completed
            {
                anSample[i] = idxi;
                break;
            }
        }
        //if(IS_DEBUG) CHECK(k >= max_random_iters, "CRANSACSampler::choose: ERROR: Could not find 7 points in max_random_iters iterations (may happen due to compatable set not being possible)");

        if (k >= max_random_iters)
            return false;
    }
    return true;
}

CGuidedMLESACSampler::CGuidedMLESACSampler(int nSampleSize, const CPointIdentifiers & pointIds, const int * anArrLikelihood, const double dPriorProb, bool bDisjoint) : CDisjointRANSACSampler(nSampleSize, pointIds), adPriorProb_doNotUseDirectly(pointIds.size()), adPriorProb(adPriorProb_doNotUseDirectly), bDisjoint(bDisjoint) {
    for (int i = 0; i < N; i++)
        adPriorProb_doNotUseDirectly[i] = dPriorProb / anArrLikelihood[i];
};

bool CGuidedMLESACSampler::choose(TSubSet & anSample) {
    const int max_random_iters = 200;

    // choose random <n> (=7) points
    for (int i = 0; i < n; i++) {
        int k = 0;
        for (; k < max_random_iters; k++) {
            int idxi = CRandom::Uniform(N); // cvRandInt(&rng) % count;

            if (adPriorProb[idxi] < 1 && adPriorProb[idxi] < CRandom::Uniform())
                continue;

            int j = 0;
            for (; j < i; j++) {
                int idxj = anSample[j];
                if (idxj == idxi)
                    break;

                //now check if either of the points selected at (earlier) time j are the same as those now:
                if (bDisjoint && pointIds.incompatible(idxi, idxj))
                    break;
            }
            if (j == i) //loop completed
            {
                anSample[i] = idxi;
                break;
            }
        }
        //if(IS_DEBUG) CHECK(k >= max_random_iters, "CRANSACSampler::choose: ERROR: Could not find 7 points in max_random_iters iterations (may happen due to compatable set not being possible)");

        if (k >= max_random_iters)
            return false;
    }
    return true;
}

template<class dataType>
struct orderPred {
    static boost::mutex mxSortLock;
    static dataType * aProb;

    bool operator()(int idx1, int idx2) const {
        return aProb[idx1] > aProb[idx2];
    }
};

template<> double * orderPred<double>::aProb = 0;
template<> int * orderPred<int>::aProb = 0;
//Should be template<> boost::mutex orderPred<double>::mxSortLock;
template<class dataType> boost::mutex orderPred<dataType>::mxSortLock;

void CPROSACSampler::init() {
    for (int i = 0; i < N; i++)
        anOrderedIndices[i] = i;

    {
        boost::mutex::scoped_lock scoped_lock(orderPred<double>::mxSortLock); //cos we're using the sort's static
        orderPred<double>::aProb = adPriorProb.begin();
        std::sort(anOrderedIndices.begin(), anOrderedIndices.end(), orderPred<double>());
    }

    if (bSystematic) {
        //Set up systematic ordering--add to a vector, then shuffle, :
        int numEquiv = 0;
        unsigned pow2 = 1 << n;
        for (unsigned i = 0; i < 1000000; i++) {
            if (i >= pow2) {
                numEquiv = 0;
                pow2 *= 2;
            }

            //Count bits:
            int nBits = 0;
            for (int bitId = 0; bitId < 32; bitId++) {
                unsigned bit = 1 << bitId;
                if (i & bit) nBits++;
            }
            if (nBits == n) {
                numEquiv++;
                int * anIndices = new int[n];
                int * pnIndices = anIndices;
                for (int bitId = 0; bitId < 32; bitId++) {
                    unsigned bit = 1 << bitId;
                    if (i & bit) {
                        *pnIndices = bitId;
                        pnIndices++;
                    }
                }
                aSamples.push_back(anIndices);

                //Now swap
                if (numEquiv > 2) {
                    int swapIdx = CRandom::Uniform((int) (aSamples.size() - numEquiv), (int) (aSamples.size() - 1));
                    aSamples[aSamples.size() - 1] = aSamples[swapIdx];
                    aSamples[swapIdx] = anIndices;
                }
            }
        }
    }
}

bool CPROSACSampler::choose(TSubSet & anSample) {
    t++;

    if (bSystematic) {
        const int * aSample = aSamples[t - 1];
        for (int i = 0; i < n; i++) {
            int idxi = anOrderedIndices[aSample[i]];
            anSample[i] = idxi;
        }
    } else {
        const int max_random_iters = 200;
        //calculate g(t);
        //Increment g and T until T>=t
        const int m = n;
        while (T < t) {
            double dT_last = dT;
            dT *= (double) (g + 1) / (double) (g + 1 - m);
            T += (int) ceil(dT - dT_last);
            g++;
        }

        if (g > N) return false;

        if(IS_DEBUG) CHECK(g < n || (int) anOrderedIndices.size() != N, "PROSAC setup error");

        // choose random <n-1> points, plus point g(t)
        int g_idx = anOrderedIndices[g - 1];
        anSample[0] = g_idx;

        int g_t = g - 1; //can be relaxed later for disjointness

        for (int i = 1; i < n; i++) {
            int k = 0;
            for (; k < max_random_iters; k++) {
                int ranki = CRandom::Uniform(g_t); // cvRandInt(&rng) % count;
                int idxi = anOrderedIndices[ranki];

                /*if(idxi>=N || idxi<0)
                        idxi=0;*/

                if(IS_DEBUG) CHECK(idxi >= N || idxi < 0, "PROSAC setup ordered indices error");

                int j = 0;
                for (; j < i; j++) {
                    int idxj = anSample[j];
                    if (idxj == idxi)
                        break;

                    //now check if either of the points selected at (earlier) time j are the same as those now:
                    if (bDisjoint && pointIds.incompatible(idxi, idxj)) {
                        if (k > 10) g_t += 2;
                        if (g_t >= N) return false;
                        break;
                    }
                }
                if (j == i) //loop completed
                {
                    anSample[i] = idxi;
                    break;
                }
            }
            if(IS_DEBUG) CHECK(k >= max_random_iters, "CRANSACSampler::choose: ERROR: Could not find 7 points in max_random_iters iterations (may happen due to compatable set not being possible)");

            if (k >= max_random_iters)
                return false;
        }
    }

    return true;
}

void CSimSACSampler::init() {
    for (int i = 0; i < N; i++)
        anAllIndices[i] = i;

    history.reserve(MAX_ITERS);
    anHistIndices.reserve(MAX_ITERS);
}

bool CSimSACSampler::choose(TSubSet & anSample) {
    if (!chooseML(anSample)) return false;

    int * anSampleHist = new int[N];
    if(IS_DEBUG) CHECK(!anSampleHist, "Alloc failed");

    for (int i = 0; i < n; i++)
        anSampleHist[i] = anSample[i];

    history.push_back(anSampleHist);
    anHistIndices.push_back(history.size() - 1);

    return true;
}

void CSimSACSampler::chooseML_approx(TSubSet & anSampleHist) {
    for (int iter = 0; iter < T; iter++) {
        for (int i = 0; i < N; i++) {
            simStatus.set(i, SS_UNSET);
        }
#ifdef RANDOMIZE
        for (THistory::const_iterator pHistory = history.begin(); pHistory < history.end(); pHistory++) {
            int * anOldHypSet = *pHistory;
#else
        for (int i = 0; i < (int) anHistIndices.size(); i++) {
            int tempIdx = anHistIndices[i];
            int * anOldHypSet = history[anHistIndices[i]];
            //Swap current index with one between 0 and i
            if (i > 1) {
                int swapIdx = CRandom::Uniform(i);
                anHistIndices[i] = anHistIndices[swapIdx];
                anHistIndices[swapIdx] = tempIdx;
            }
#endif
            ///////////////
            bool bContainsOutlier = false;
            bool bOneToSet = false;

            ARRAY(TStatus, aStatuses, n);

            bOneToSet = false;
            for (int i = 0; i < n; i++) //iterate over hyp. set
            {
                //We used to repeat this loop several times (as long as some statuses are 'unset') until we find an outlier (or some outliers)
                //Now: Just choose one of the ones we have set 'true' to be 'false'

                int idx = anOldHypSet[i];
                if (simStatus.get(idx) == SS_UNSET) {
                    TStatus status = CRandom::Uniform() < adPriorProb[i];
                    if (status == SS_FALSE) bContainsOutlier = true;
                    aStatuses[i] = status;
                    bOneToSet = true;
                } else {
                    if (simStatus.get(idx) == SS_FALSE) bContainsOutlier = true;
                }
            }
            if (!bOneToSet && !bContainsOutlier) {
                //No outliers--all of aStatuses are set SS_TRUE
                //Unset one, weighted by *outlier* probability

                //Sum prior probabilities of being outliers
                double dSum = 0;
                for (int i = 0; i < n; i++) {
                    int idx = anOldHypSet[i];
                    if (simStatus.get(idx) == SS_TRUE)
                        dSum += (1 - adPriorProb[idx]);
                }

                //Choose one
                double dSelect = CRandom::Uniform(dSum);
                dSum = 0;
                int i = 0, idx = 0;
                for (; i < n; i++) {
                    idx = anOldHypSet[i];
                    if (simStatus.get(idx) == SS_TRUE) {
                        dSum += (1 - adPriorProb[idx]);
                        if (dSum >= dSelect) break;
                    }
                }
                if(IS_DEBUG) CHECK(i == n, "Failed to select inlier to flip");
                simStatus.set(idx, SS_FALSE);
                bContainsOutlier = true;
            } else if (!bContainsOutlier) {
                //Need to re-choose one of the ones we have just set true and flip it.
                //Sum prior probabilities of being outliers
                double dSum = 0;
                for (int i = 0; i < n; i++) {
                    int idx = anOldHypSet[i];
                    if (simStatus.get(idx) == SS_UNSET && aStatuses[i] == SS_TRUE)
                        dSum += (1 - adPriorProb[idx]);
                }

                //Choose one
                double dSelect = CRandom::Uniform(dSum);
                dSum = 0;
                int i = 0, idx = 0;
                for (; i < n; i++) {
                    idx = anOldHypSet[i];
                    if (simStatus.get(idx) == SS_UNSET && aStatuses[i] == SS_TRUE) {
                        dSum += (1 - adPriorProb[idx]);
                        if (dSum >= dSelect) break;
                    }
                }
                if(IS_DEBUG) CHECK(i == n, "Failed to select inlier to flip");
                aStatuses[i] = SS_FALSE;
                bContainsOutlier = true;
            }

            //} while (!bContainsOutlier);

            if(IS_DEBUG) CHECK(!bContainsOutlier, "We should have set an outlier by now");
            //Now we have observed an outlier so this simulation is good
            //Copy it into simulated set
            for (int i = 0; i < n; i++) //iterate over hyp. set
            {
                int idx = anOldHypSet[i];
                if (simStatus.get(idx) == SS_UNSET)
                    simStatus.set(idx, aStatuses[i]);
            }
        }
        //Now simulate remaining statuses
        for (int i = 0; i < N; i++) //iterate over hyp. set
        {
            if (simStatus.get(i) == SS_UNSET) {
                TStatus status = CRandom::Uniform() < adPriorProb[i];
                simStatus.set(i, status);
            }
        }

        // accumulate histogram
        for (int i = 0; i < N; i++) //iterate over hyp. set
            if (simStatus.get(i) == SS_TRUE)
                anSampleHist[i]++;
    }
}

void CSimSACSampler::chooseML_approx2(TSubSet & anSampleHist) {
    for (int iter = T; iter > 0; iter--) {
        for (int i = 0; i < N; i++) {
            TStatus status = CRandom::Uniform() < adPriorProb[i];
            simStatus.set(i, status);
        }
#ifdef RANDOMIZE
        for (THistory::const_iterator pHistory = history.begin(); pHistory < history.end(); pHistory++) {
            int * anOldHypSet = *pHistory;
#else
        for (int iidx = (int) anHistIndices.size(); iidx > 0; iidx--) {
            int iidx_1 = iidx - 1;
            int tempIdx = anHistIndices[iidx_1];
            int * anOldHypSet = history[tempIdx];
            //Swap current index with one between 0 and i
            if (iidx > 2) {
                int swapIdx = CRandom::UniformFast(iidx_1);
                anHistIndices[iidx_1] = anHistIndices[swapIdx];
                anHistIndices[swapIdx] = tempIdx;
            }
#endif
            ///////////////
            bool bContainsOutlier = false;

            for (int i = n; i > 0; i--) //iterate over hyp. set
            {
                int idx = anOldHypSet[i - 1];
                if (simStatus.get(idx) == SS_FALSE) {
                    bContainsOutlier = true;
                    break;
                }
            }
            if (!bContainsOutlier) {
                //Unset one, weighted by *outlier* probability

                //Sum prior probabilities of being outliers
                double dSum = 0;
                for (int i = 0; i < n; i++) {
                    int idx = anOldHypSet[i];
                    if (simStatus.get(idx) == SS_TRUE)
                        dSum += (1 - adPriorProb[idx]);
                }

                //Choose one
                double dSelect = CRandom::Uniform(dSum);
                dSum = 0;
                int i = 0, idx = 0;
                for (; i < n; i++) {
                    idx = anOldHypSet[i];
                    if (simStatus.get(idx) == SS_TRUE) {
                        dSum += (1 - adPriorProb[idx]);
                        if (dSum >= dSelect) break;
                    }
                }
                if(IS_DEBUG) CHECK(i == n, "Failed to select inlier to flip");
                simStatus.set(idx, SS_FALSE);
            }
        }

        // accumulate histogram
        for (int i = 0; i < N; i++) //iterate over hyp. set
            if (simStatus.get(i) == SS_TRUE)
                anSampleHist[i]++;
    }
}

void CSimSACSampler::chooseML_exact(TSubSet & anSampleHist) {
    for (int iter = 0; iter < T; iter++) {
        for (int i = 0; i < N; i++) //Todo: think about what happens when we select 2 the same...
        {
            simStatus.set(i, CRandom::Uniform() < adPriorProb[i]);
        }

        bool bReject = false;
        for (THistory::const_iterator pHistory = history.begin(); pHistory < history.end(); pHistory++) {
            int * anOldHypSet = *pHistory;
            bool bContainsOutlier = false;
            for (int i = 0; i < n; i++) //iterate over hyp. set
            {
                int idx = anOldHypSet[i];
                if (simStatus.get(idx) == SS_FALSE) {
                    bContainsOutlier = true;
                    break;
                }
            }
            if (!bContainsOutlier) //found an inlier set in history--reject this sim
            {
                bReject = true;
                break;
            }
        }

        if (!bReject) {
            // accumulate histogram
            for (int i = 0; i < N; i++) //iterate over hyp. set
                if (simStatus.get(i) == SS_TRUE)
                    anSampleHist[i]++;
        }
    }
}
//First n+1 must be partially sorted

template <typename orderType>
bool chooseFromPartiallySorted(const CPointIdentifiers & pointIds, TSubSet & anSample, const CDynArray<int> & anPartiallyOrderedIndices, const orderType * aSortVals, bool bDisjoint) {
    const int n = anSample.size();
    const int N = pointIds.size();

    orderType lastSortVal = aSortVals[anPartiallyOrderedIndices[n - 1]];

    if (bDisjoint) //If they're disjoint they prob. have variable probabilities
    {
        //First choose top m that are disjoint, then revise last few.
        cout << 1;

        int idx = 0, i = 0;
        for (; i < n; i++) {
            anSample[i] = -1;
            for (; idx < N;) {
                int candidateIdx = anPartiallyOrderedIndices[idx];
                idx++;
                cout << candidateIdx << "=cand\n";

                //Check not incompat. with previous:
                int j = 0;
                for (; j < i; j++) {
                    //cout << anSample[j] << "=already chosen\n";
                    if (pointIds.incompatible(candidateIdx, anSample[j])) {
                        //cout << "Incompat: " << pointIds[candidateIdx].id1() << ',' << pointIds[candidateIdx].id2() << '-' << pointIds[anSample[j]].id1() << ',' << pointIds[anSample[j]].id2() << '\n';
                        break;
                    }


                    //if(aPoints0[candidateIdx] == aPoints0[anSample[j]] || aPoints1[candidateIdx] == aPoints1[anSample[j]]) break;
                }
                if (j == i) // candidateIdx is compatible with prev chosen indices
                {
                    anSample[i] = candidateIdx;
                    break;
                }
            }
            if (anSample[i] == -1) break;
        }
        if (i != n || idx >= N)
            return false; //could not select disjoint set

        //cout << 4;
        lastSortVal = aSortVals[anPartiallyOrderedIndices[idx]];

        if ((int) anPartiallyOrderedIndices.size() == idx || aSortVals[anSample[n - 1]] != lastSortVal)
            return true;

        //Choose i to n-1 randomly with vals aSortVals[anSample[n-1]]
        while (i > 0 && aSortVals[anSample[i - 1]] == lastSortVal)
            i--;

        //Now choose n-i random el's with prob lastSortVal
        static std::vector<int> anIndices;
        anIndices.clear();
        for (int idx = 0; idx < N; idx++) {
            if (aSortVals[idx] == lastSortVal)
                anIndices.push_back(idx);
        }

        for (int j = i; j < n; j++) {
            int breakOut = 1000;
            for (;; breakOut--) //will never loop infinitely because there must be a possible disjoint set to get this far.
            {
                if (!breakOut) {
                    //std::cout << "Disjoint sampling error\n";
                    return false;
                }
                int nRand = CRandom::UniformFast((int) anIndices.size());
                int idx = anIndices[nRand];
                //Check not already chosen, or incompatible:
                int k = 0;
                for (; k < j; k++) {
                    if (idx == anSample[k]) break;
                    if (pointIds.incompatible(idx, anSample[k])) break;
                }
                if (k == j) {
                    anSample[j] = idx;
                    break;
                }
            }
        }
        return true;
    } else {
        //not disjoint
        if ((int) anPartiallyOrderedIndices.size() == n || lastSortVal != aSortVals[anPartiallyOrderedIndices[n]]) {
            for (int i = 0; i < n; i++)
                anSample[i] = anPartiallyOrderedIndices[i];
            return true;
        }

        int i = 0;
        for (; i < n - 1; i++) {
            int idx = anPartiallyOrderedIndices[i];

            if (aSortVals[idx] == lastSortVal)
                break;

            anSample[i] = idx;
        }

        //Now choose n-i random el's with prob lastSortVal
        static std::vector<int> anIndices;
        anIndices.clear();
        for (int idx = 0; idx < N; idx++) {
            if (aSortVals[idx] == lastSortVal)
                anIndices.push_back(idx);
        }

        int nOptions = anIndices.size();
        if (nOptions == n - i) {
            for (int j = i; j < n; j++) {
                anSample[j] = anIndices[j];
            }
        } else {
            for (int j = i; j < n; j++) {
                for (;;) {
                    int nRand = CRandom::UniformFast(nOptions);
                    int idx = anIndices[nRand];
                    //Check not already chosen:
                    int k = i;
                    for (; k < j; k++) {
                        if (idx == anSample[k]) break;
                    }
                    if (k == j) {
                        anSample[j] = idx;
                        break;
                    }
                }
            }
        }
    }
    return true;
}

bool CSimSACSampler::chooseML(TSubSet & anSample) {
    TSubSet anSampleHist(N, 0);

    if (bExact)
        chooseML_exact(anSampleHist);
    else
        chooseML_approx(anSampleHist);

    //Could just make one iteration, then choose a set of inliers.

    //Now find n peaks in histogram
    {
        boost::mutex::scoped_lock scoped_lock(orderPred<int>::mxSortLock); //cos we're using the sort's static
        orderPred<int>::aProb = anSampleHist.begin();

        CDynArray<int>::iterator endIter = anAllIndices.begin();
        if (bDisjoint) {
            endIter += min<int>(anAllIndices.size(), n * 2 * NN_MAX);
        } else {
            endIter += n;
            if (endIter != anAllIndices.end()) endIter++;
        }

        std::partial_sort(anAllIndices.begin(), endIter, anAllIndices.end(), orderPred<int>());
    }

    //Essentially reverts to normal ransac when few samples have nonzero likelihood
    bool bRes = chooseFromPartiallySorted<int>(pointIds, anSample, anAllIndices, anSampleHist.begin(), bDisjoint);
    cout << endl;
    return bRes;
}

void CNaiveBayesSampler::init() {
    for (int i = 0; i < N; i++) {
        anOrderedIndices[i] = i;
        //adPriorProb[i] += CRandom::Unif0();
    }
}

bool CNaiveBayesSampler::choose(TSubSet & anSample) {
    //Choose n most likely...
    {
        boost::mutex::scoped_lock scoped_lock(orderPred<double>::mxSortLock); //cos we're using the sort's static
        orderPred<double>::aProb = adPriorProb.begin();
        CDynArray<int>::iterator endIter = anOrderedIndices.begin();
        if (bDisjoint) {
            endIter += min<int>(anOrderedIndices.size(), n * 2 * NN_MAX);
        } else {
            endIter += n;
            if (endIter != anOrderedIndices.end()) endIter++;
        }
        //std::partial_sort(anOrderedIndices.begin(), endIter, anOrderedIndices.end(), orderPred<double>());
        //std::sort(anOrderedIndices.begin(), anOrderedIndices.end(), orderPred<double>());
        CDynArray<int>::iterator beginIt = anOrderedIndices.begin();
        CDynArray<int>::iterator endIt = anOrderedIndices.end();
        int * begin = beginIt;
        int * end = endIt;
        int * mid = endIter;
        std::partial_sort(begin, mid, end, orderPred<double>());
    }

    if (!chooseFromPartiallySorted<double>(pointIds, anSample, anOrderedIndices, orderPred<double>::aProb, bDisjoint)) return false;

    //...and update probs
    double dProbAllInliers = 1.0;
    for (int i = 0; i < n; i++)
        dProbAllInliers *= adPriorProb[anSample[i]];

    double dHbadProb_Inv = 1.0 / (1.0 - dProbAllInliers);

#ifdef __GNUC__
    bool abUpdated[N];
#else
    boost::scoped_array<bool> abUpdated(new bool[N]);
#endif

    if (bDisjoint)
        for (int i = 0; i < N; i++)
            abUpdated[i] = false;

    for (int i = 0; i < n; i++) {
        int idx = anSample[i];
        //cout << adPriorProb[idx] << ',';
        adPriorProb[idx] = (adPriorProb[idx] - dProbAllInliers) * dHbadProb_Inv;
        //cout << adPriorProb[idx] << endl;
        abUpdated[idx] = true;

        if (bDisjoint) {
            //Todo: speedup
            for (int j = 0; j < N; j++)
                //TODO: THIS SPEEDUP ONLY SS_TRUE FOR SIM DATA
                /*int lo = max<int>(0, idx-NN_MAX*NN_MAX);
                int hi = min<int>(N, idx+NN_MAX*NN_MAX);
                for(int j=lo;j<hi;j++)*/ {
                if (!abUpdated[j]) //only update each prob once
                {
                    if (pointIds.incompatible(idx, j)) {
                        //if(j != i) //we never choose incompatible sets, so j!=i means j not in set.
                        //{
                        abUpdated[j] = true;

                        //adPriorProb[j] *= dHbadProb_Inv;
                        double newDenom = (adPriorProb[j] + (1 - adPriorProb[j])*(1.0 - dProbAllInliers));
                        adPriorProb[j] /= newDenom;

                        if (adPriorProb[j] > 1)
                            cout << "Update error: adPriorProb[j]==" << adPriorProb[j] << endl;
                        //}
                    }
                }
            }
        }
    }

    return true;
}

void CBaySACSampler::init() {
    for (int i = 0; i < N; i++) {
        sortedProbs.insert(idxPair(i, adPriorProb.begin() + i));
    }
}

#define DEBUG_BAYSAC(...) //__VA_ARGS__

void CBaySACSampler_discretePP::init() {
    typedef map2<int, TFastIntVec, std::less<int>/*, CIndividualPool_NoFree_Allocator<std::pair<int, TFastIntVec>, 512 >*/ > TPointCollisionMaps;
    TPointCollisionMaps samePointLeft, samePointRight;

    ARRAY(TFastIntVec *, avVectorsLeft, N); //Remember locations to prevent a set lookup
    ARRAY(TFastIntVec *, avVectorsRight, N);

    for (int i = 0; i < N; i++) {
        ((TSortedProbs::TMap &)sortedProbVectors)[adPriorProb[i]].insert(i);
        if (bDisjoint) {
            //add i to a list for its left and right points
            TFastIntVec & vLeft = samePointLeft.initOrGet(pointIds[i].id1());
            vLeft.push_back(i);
            TFastIntVec & vRight = samePointRight.initOrGet(pointIds[i].id2());
            vRight.push_back(i);

            //Convert references to pointers...
            //Prob not std c++...
            avVectorsLeft[i] = &vLeft;
            avVectorsRight[i] = &vRight;
        }
    }
    if (bDisjoint) {
        for (int i = 0; i < N; i++) {
            TLongerFastIntVec & aIncompat = aaIncompatable[i];
            /*TFastIntVec & aIncompatLeft = samePointLeft[pointIds[i].id1()],
                                    & aIncompatRight = samePointRight[pointIds[i].id2()];//NB These are DISJOINT except for i*/
            TFastIntVec & aIncompatLeft = *(avVectorsLeft[i]),
                    & aIncompatRight = *(avVectorsRight[i]);

            for (TFastIntVec::const_iterator pLeft = aIncompatLeft.begin(); pLeft != aIncompatLeft.end(); pLeft++) {
                int nLeft = *pLeft;
                if (nLeft != i)
                    aIncompat.push_back(nLeft);
            }
            for (TFastIntVec::const_iterator pRight = aIncompatRight.begin(); pRight
                    != aIncompatRight.end(); pRight++) {
                int nRight = *pRight;
                if (nRight != i)
                    aIncompat.push_back(nRight);
            }
#ifdef _DEBUG
            //check prior prob of a point being an inlier never exceeds 1. This normally means 2 different descriptors in the same location
            double dProbLeft = 0;
            for (TFastIntVec::const_iterator pLeft = aIncompatLeft.begin(); pLeft != aIncompatLeft.end(); pLeft++)
                dProbLeft += adPriorProb[*pLeft];
            if (dProbLeft == 1)
                cout << "Warning: this BaySAC implementation can perform poorly when some points are assigned prior probility 1 of being inliers.\n";
            else if (dProbLeft > 1) {
                cout << "Total prior inlier prob = " << dProbLeft << " (" << aIncompatLeft.size() << " points). This can be caused by the image-motion model when mosaicing\n";
                THROW("CBaySACSampler_discretePP::init(): Total inlier prob for a point exceeds 1");
            }
            //check prior prob of a point being an inlier never exceeds 1
            double dProbRight = 0;
            for (TFastIntVec::const_iterator pRight = aIncompatRight.begin(); pRight != aIncompatRight.end(); pRight++)
                dProbRight += adPriorProb[*pRight];
            if (dProbRight == 1)
                cout << "Warning: this BaySAC implementation can perform poorly when some points are assigned prior probility 1 of being inliers.\n";
            else if (dProbRight > 1) {
                cout << "Total prior inlier prob = " << dProbRight << " (" << aIncompatRight.size() << " points)\n";
                cout << ("CBaySACSampler_discretePP::init(): Total inlier prob for a point exceeds 1\n"); //This can happen with synthetic data that's all tightly clustered
                cout << "This normally means 2 different descriptors in the same location (e.g. different sized SURF blobs or synthetic data)\n";
            }
#endif
        }
    }
    //cout << "INIT COMPLETE" << endl;
}

int CBaySACSampler_discretePP::chooseFromSet(int numAlreadyChosen, TSubSet & anSample, const TIntSet & options, TExcludeSet & excludeSet, int & nNumExcluded) {
    const int numOptions = options.size();
    int numChosen = 0;

    if (numOptions <= n - numAlreadyChosen) {
        for (TIntSet::const_iterator pOption = options.begin(); pOption != options.end(); pOption++) {
            int nCandidate = *pOption;
            if (!bDisjoint) {
                anSample[numAlreadyChosen + numChosen] = nCandidate;
                numChosen++;
            } else {
                int j = 0;

                for (; j < numAlreadyChosen + numChosen; j++) {
                    if (pointIds.incompatible(nCandidate, anSample[j]))
                        break;
                }
                if (j == numAlreadyChosen + numChosen) {
                    anSample[numAlreadyChosen + numChosen] = nCandidate;
                    numChosen++;
                }
            }
        }
        return numChosen;
    }

    ARRAY(int, anOptions, numOptions);

    int * pnOptions = PTR(anOptions);

    TIntSet::const_iterator pOption = options.begin();
    for (int i = numOptions; i > 0; i--) {
        *pnOptions = *pOption;
        pOption++;
        pnOptions++;
    }

    //TIntSet excludeSet;

    //First exclude a few that are incompatible with the indices already chosen:
    if (bDisjoint) {
        for (; nNumExcluded < numAlreadyChosen; nNumExcluded++) {
            int prevIdx = anSample[nNumExcluded];

            const TLongerFastIntVec & aIncompat = aaIncompatable[prevIdx];
            for (TFastIntVec::const_iterator pIncompat = aIncompat.begin(); pIncompat != aIncompat.end(); pIncompat++) {
                int nIncompat = *pIncompat;
                if (options.find(nIncompat) != options.end())
                    excludeSet.insert(nIncompat);
            }
        }
    }

    while (numChosen < n - numAlreadyChosen && (int) excludeSet.size() < numOptions) {
        int nCandidate = anOptions[CRandom::UniformFast(numOptions)];

        while (excludeSet.exists(nCandidate))
            nCandidate = anOptions[CRandom::UniformFast(numOptions)];

        anSample[numAlreadyChosen + numChosen] = nCandidate;
        numChosen++;

        excludeSet.insert(nCandidate);

        if (bDisjoint && numChosen < n - numAlreadyChosen) {
            nNumExcluded++;

            const TLongerFastIntVec & aIncompat = aaIncompatable[nCandidate];
            for (TFastIntVec::const_iterator pIncompat = aIncompat.begin(); pIncompat != aIncompat.end(); pIncompat++) {
                int nIncompat = *pIncompat;
                if (options.find(nIncompat) != options.end())
                    excludeSet.insert(nIncompat);
            }
        }
    }

    return numChosen;
}

bool CBaySACSampler_discretePP::choose(TSubSet & anSample) {
    //Choose n most likely...
    int numChosen = 0, nNumExcluded = 0;

    TExcludeSet excludeSet;
    double dNextProb = 0;
    for (TSortedProbs::iterator pIdxVector = sortedProbVectors.begin(); pIdxVector != sortedProbVectors.end(); pIdxVector++) {
        bool bMoreAtSameProbIfLast = (int) (pIdxVector->second.size()) > (n - numChosen); // In this case any bayes update will change probs
        numChosen += chooseFromSet(numChosen, anSample, pIdxVector->second, excludeSet, nNumExcluded);
        if (numChosen == n) {
            if (!bMoreAtSameProbIfLast) {
                pIdxVector++;
                //could possibly segfault...
                dNextProb = pIdxVector->first;
            }
            break;
        }
    }
    if (numChosen < n)
        return false;

    //Now update these probs
    double dProbAllInliers = 1.0;
    for (int i = 0; i < n; i++) {
        dProbAllInliers *= adPriorProb[anSample[i]];
        DEBUG_BAYSAC(cout << adPriorProb[anSample[i]] << ' ' << anSample[i] << endl);
    }
    DEBUG_BAYSAC(cout << endl);


    /* To ensure we discount probabilities enough to get a 
     * new hyp set (important when n is large and individual 
     * probs are small) we must ensure:
     * dProbAllInliers > (p_n - p_{n+1})/(1 - p_{n+1})
     */
    if (dNextProb) {
        const double EPS = 1e-3;
        const double dMinInlierProb = ((adPriorProb[anSample[n - 1]] - dNextProb) / (1 - dNextProb)) + EPS;
        if (dProbAllInliers < dMinInlierProb) {
            DEBUG_BAYSAC(cout << dProbAllInliers << " updated to " << dMinInlierProb << endl);
            dProbAllInliers = dMinInlierProb;
        }
    }

    double dHbadProb_Inv = 1.0 / (1.0 - dProbAllInliers);
    /*bool abUpdated[N];
    if(bDisjoint)
            for(int i=0;i<N;i++)
                    abUpdated[i] = false;*/

    TIntSet alreadyUpdatedSet;

    for (int i = 0; i < n; i++) {
        int idx = anSample[i];

        const double dOldProb = adPriorProb[idx];

        TSortedProbs::iterator oldLoc = sortedProbVectors.begin(); //Always true if i=0
        if (oldLoc->first != dOldProb)
            oldLoc = sortedProbVectors.find(dOldProb);

        /*if(oldLoc == sortedProbVectors.begin())
                cout << "is begin()\n";
        else
                cout << "is not begin()\n";*/

        if(IS_DEBUG) CHECK(oldLoc == sortedProbVectors.end(), "Old loc doesn't exist")

        oldLoc->second.erase(idx);
        if (oldLoc->second.size() == 0) {
            sortedProbVectors.erase(oldLoc);
        }

        //sortedProbVectors[dOldProb].erase(idx);
        //if(sortedProbVectors[dOldProb].size() == 0) sortedProbVectors.erase(dOldProb);

        const double dNewProb = (dOldProb - dProbAllInliers) * dHbadProb_Inv;
        adPriorProb[idx] = dNewProb;

        ((TSortedProbs::TMap &)sortedProbVectors)[dNewProb].insert(idx);
        //cout<< "Decreasing prob of " << idx << " from " << dOldProb << " to " << dNewProb << endl;

        if (bDisjoint) {
            const TLongerFastIntVec & aIncompat = aaIncompatable[idx];

            for (TFastIntVec::const_iterator pIncompat = aIncompat.begin(); pIncompat != aIncompat.end(); pIncompat++) {
                int j = *pIncompat;

                if (alreadyUpdatedSet.find(j) == alreadyUpdatedSet.end()) //only update each prob once
                {
                    alreadyUpdatedSet.insert(j);

                    double dOldProbDisj = adPriorProb[j];

                    TSortedProbs::iterator pProbSet = sortedProbVectors.find(dOldProbDisj); //Segfault here...
                    if(IS_DEBUG) CHECK(pProbSet == sortedProbVectors.end(), "Prob never existed here");
                    if(IS_DEBUG) CHECK(!(pProbSet->second.exists(j)), "Idx never existed here");

                    pProbSet->second.erase(j);

                    if (pProbSet->second.size() == 0)
                        sortedProbVectors.erase(pProbSet);

                    //adPriorProb[j] *= dHbadProb_Inv;
                    //double newDenom = (dOldProbDisj + (1-dOldProbDisj)*(1.0-dProbAllInliers));
                    double dNewProbDisj = dOldProbDisj * (1 - dNewProb) / (1 - dOldProb);
                    adPriorProb[j] = dNewProbDisj;

                    //cout<< "Increasing prob of " << j << " from " << dOldProbDisj << " to " << dNewProbDisj << endl;

                    if (dNewProbDisj > 1)
                        cout << "adPriorProb[j]==" << adPriorProb[j] << endl;

                    ((TSortedProbs::TMap &)sortedProbVectors)[dNewProbDisj].insert(j);
                }
            }
        }
    }
    /*for(int i=0;i<n;i++)
            cout << anSample[i] << " (" << adPriorProb[anSample[i]] << "), ";
    cout << endl;*/

    return true;
}

bool CBaySACSampler::choose(TSubSet & anSample) {
    //Choose n most likely...
    int i = 0;
    const bool NI_bDisjoint = false; //DISJOINT CASE IS HARDER TO IMPLEMENT--NOT IMPLEMENTED YET

    TSortedProbs::iterator pPoint = sortedProbs.begin();
    for (; pPoint != sortedProbs.end(); pPoint++) {
        int candidateIdx = pPoint->first;

        int j = 0;

        if (NI_bDisjoint) {
            for (; j < i; j++) {
                if (pointIds.incompatible(candidateIdx, anSample[j]))
                    break;
            }
        }

        if (!NI_bDisjoint || j == i) {
            anSample[i] = candidateIdx;
            i++;
            if (i == n) break;
        }
    }
    if(IS_DEBUG) CHECK(!NI_bDisjoint && i < n, "FastBayes: Failed to select n points");
    if (i < n) return false; //could not select n points because of disjointness

    if (i == N) return true;

    //Now randomise. pPoint is a pointer to the last selected
    TSortedProbs::iterator pEndEquiprobRange = pPoint;
    pEndEquiprobRange++;

    double dLatestProb = *(pPoint->second);
    if (*(pEndEquiprobRange->second) == dLatestProb) {
        while (pEndEquiprobRange != sortedProbs.end() && *(pEndEquiprobRange->second) != dLatestProb)
            pEndEquiprobRange++;

        //Now push all candidate indices onto a vector and sample randomly.
        static vector<int> anEquiprobIndices;
        anEquiprobIndices.clear();

        for (;;) {
            pEndEquiprobRange--;

            if (*(pEndEquiprobRange->second) != dLatestProb) break;

            anEquiprobIndices.push_back(pEndEquiprobRange->first);

            if (pEndEquiprobRange == sortedProbs.begin()) break;
        }

        int nSize = (int) anEquiprobIndices.size();

        //Now find which indices we're actually setting
        for (i = 0; i < n; i++)
            if (dLatestProb == adPriorProb[anSample[i]]) break;

        if(IS_DEBUG) CHECK(nSize < n - i - 1, "FastNaiveBayes: should have more equiprobable indices");
        //Now set i..n

        for (; i < n; i++) {
            for (;;) {
                int nCandidateIdx = anEquiprobIndices[CRandom::Uniform(nSize)];
                int j = 0;
                for (; j < i; j++) {
                    if (anSample[j] == nCandidateIdx) break;
                }
                if (j == i) {
                    anSample[i] = nCandidateIdx;
                    break;
                }
            }
        }
    }

    //...and update probs
    double dProbAllInliers = 1.0;
    for (int i = 0; i < n; i++)
        dProbAllInliers *= adPriorProb[anSample[i]];

    double dHbadProb_Inv = 1.0 / (1.0 - dProbAllInliers);

    ARRAY(bool, abUpdated, N);
    if (bDisjoint)
        for (int i = 0; i < N; i++)
            abUpdated[i] = false;

    for (int i = 0; i < n; i++) {
        int idx = anSample[i];
        TSortedProbs::iterator pPoint = sortedProbs.begin();
        while (pPoint->first != idx) pPoint++;
        idxPair moved = *pPoint;

        sortedProbs.erase(pPoint);

        adPriorProb[idx] = (adPriorProb[idx] - dProbAllInliers) * dHbadProb_Inv;

        sortedProbs.insert(moved);

        if (NI_bDisjoint) {
            //Todo: speedup
            //for(int j=0;j<N;j++)
            //TODO: THIS SPEEDUP ONLY SS_TRUE FOR SIM DATA
            int lo = max<int>(0, idx - NN_MAX * NN_MAX);
            int hi = min<int>(N, idx + NN_MAX * NN_MAX);
            for (int j = lo; j < hi; j++) {
                if (!abUpdated[j]) //only update each prob once
                {
                    if (pointIds.incompatible(idx, j)) {
                        if (j != i) //we never choose incompatible sets, so j!=i means j not in set.
                        {
                            abUpdated[j] = true;

                            //adPriorProb[j] *= dHbadProb_Inv;
                            double newDenom = (adPriorProb[j] + (1 - adPriorProb[j])*(1.0 - dProbAllInliers));
                            adPriorProb[j] /= newDenom;

                            if (adPriorProb[j] > 1)
                                cout << "adPriorProb[j]==" << adPriorProb[j] << endl;
                        }
                    }
                }
            }
        }
    }

    return true;
}
