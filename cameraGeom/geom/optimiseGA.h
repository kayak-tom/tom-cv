#include "levMarNumerical.h"

template< int N >
class CGAParamSelector {
    const int DISCRETISE;
    const bool bVerbose;

    int nInserted;

    class CParamSet {
    public:
        double dErr;
        Eigen::VectorXd params;

        CParamSet() : dErr(HUGE) {
        }

        CParamSet(const double dErr, const Eigen::VectorXd & params) : dErr(dErr), params(params) {
        }
    };

    CParamSet paramSet[N];
public:
    // DISCRETISE keeps DISCRETISE-D points (e.g. 2d chunks of the param vec) together.

    CGAParamSelector(const int DISCRETISE = 1, const bool bVerbose=false) : DISCRETISE(DISCRETISE), bVerbose(bVerbose), nInserted(0) {
    }

    void addParams(const Eigen::VectorXd & params, const double dErr) {
        nInserted++;
        //Find insert pos:
        double dMaxErr = 0;
        int nInsertPos = -1, nClosestIdx = -1;

        double dClosestDist = HUGE;

        for (int i = 0; i < N; i++) {
            if (paramSet[i].dErr > dMaxErr) {
                dMaxErr = paramSet[i].dErr;
                nInsertPos = i;
            }

            if (paramSet[i].params.size() > 0) {
                double dDistToParams = (paramSet[i].params - params).squaredNorm();
                if (dDistToParams < dClosestDist) {
                    dClosestDist = dDistToParams;
                    nClosestIdx = i;
                }
            }
        }

        if (nInsertPos >= 0 && dMaxErr > dErr) {
            if (dClosestDist < sqr(0.05)) {
                if (paramSet[nClosestIdx].dErr > dErr) {
                    if(bVerbose) REPEAT(100, cout << "Params too close to existing parameter set, replacing" << endl);
                    nInsertPos = nClosestIdx;
                } else {
                    if(bVerbose) REPEAT(100, cout << "Params too close to existing parameter set, rejecting" << endl);
                    return;
                }
            }

            paramSet[nInsertPos] = CParamSet(dErr, params);
        }
    }

    CParamSet bestParamSet() const {
        if(IS_DEBUG) CHECK(nInserted < 1, "No params inserted yet");
        double dBestErr = HUGE;
        int nBest = -1;
        for (int i = 0; i < N; i++) {
            if (paramSet[i].dErr < dBestErr) {
                dBestErr = paramSet[i].dErr;
                nBest = i;
            }
        }
        return paramSet[nBest];
    }

    Eigen::VectorXd bestParams() const {
        return bestParamSet().params;
    }

    double bestErr() const {
        return bestParamSet().dErr;
    }

    int numSamples() const {
        return nInserted;
    }

    //Bernoulli RV

    bool bernoulli(double p) const {
        return rand() < (double) RAND_MAX*p;
    }

    double scale(const bool bMC) const {
        if(bMC)
            return 0.95 + 0.1 * (rand() / (double) RAND_MAX);
        else
            return 0.75 + 0.5 * (rand() / (double) RAND_MAX);
    }

    int selectParent() const {
        if(IS_DEBUG) CHECK(nInserted < 1, "No params inserted yet");

        for (;;) {
            int nParent = rand() % N;
            if (paramSet[nParent].dErr < HUGE)
                return nParent;
        }
    }

    Eigen::VectorXd select(const bool bMC = false) const {
        const int nParent1 = selectParent();
        const int nParent2 = bMC ? nParent1 : selectParent();

        const int nParams = (int) paramSet[nParent1].params.size();
        if(IS_DEBUG) CHECK(nParams % DISCRETISE, "Params don't split into discrete chunks");

        Eigen::VectorXd newParams(nParams);

        int nChunks = nParams / DISCRETISE;

        for (int nChunk = 0; nChunk < nChunks; nChunk++) {
            int nIdx = nChunk*DISCRETISE;

            const int nParent = (bernoulli(0.5) ? nParent1 : nParent2);
            newParams.segment(nIdx, DISCRETISE) = paramSet[nParent].params.segment(nIdx, DISCRETISE);

            if (bernoulli(0.25)) {
                //Scale each component by some factor from 0.5...1.5
                for (int nComponent = 0; nComponent < DISCRETISE; nComponent++) {
                    newParams(nIdx + nComponent) *= scale(bMC);
                }
            }
        }
        if(bVerbose) { REPEAT(100, cout << "Combined: " << paramSet[nParent1].params.transpose() << endl << "     and: " << paramSet[nParent2].params.transpose() << endl << "  giving: " << newParams.transpose() << endl); }
        return newParams;
    }

    double minimise(CLMFunction & fn, Eigen::VectorXd & params, const int nMaxIters)
    {
        Eigen::VectorXd resids(fn.values());
        const double dInitErr = fn.sumSquare(params, resids);
        addParams(params, dInitErr);

        for(int i=0; i<nMaxIters; i++)
        {
            Eigen::VectorXd paramsGA = select();
            addParams(paramsGA, fn.sumSquare(paramsGA, resids));

            if(bVerbose)
                if(i % (nMaxIters/100) == 0)
                    cout << i << ": " << bestErr() << endl;
        }
        params = bestParams();
        const double dFinalErr = bestErr();

        cout << "GA relative reduction in error of " << ((dInitErr - dFinalErr) / dInitErr) << endl;

        return dFinalErr;
    }
};
