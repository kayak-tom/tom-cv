#pragma once
#include "hypothesiser.h"

class CEssentialMatLM :
public CImCorrModelHypothesiser {
    const int nHypotheses;
protected:

    //Return number of hypotheses
    virtual int getModels_int(const TSubSet & anHypSet, CModels & pModels);
public:

    CEssentialMatLM(const T2dPoints & p1, const T2dPoints & p2, const int nPoints, const int nHypotheses) : CImCorrModelHypothesiser(nPoints, p1, p2), nHypotheses(nHypotheses) {

    }

    ~CEssentialMatLM(void) {
    }

    virtual bool isThreadsafe() const {
        return true;
    } // as has no state

    virtual int modelsPerIteration() const {
        return (nPoints == 5) ? 4 /*about 3.8*/ : nHypotheses;
    }

    virtual int timePerIteration() const {
        return 3000;
    }
};

class CEssentialMatLMbasis :
public CImCorrModelHypothesiser {
protected:

    //Return number of hypotheses
    virtual int getModels_int(const TSubSet & anHypSet, CModels & pModels);
public:

    CEssentialMatLMbasis(const T2dPoints & p1, const T2dPoints & p2, const int nPoints) : CImCorrModelHypothesiser(nPoints, p1, p2) {

    }

    ~CEssentialMatLMbasis(void) {
    }

    virtual bool isThreadsafe() const {
        return true;
    } // as has no state

    virtual int modelsPerIteration() const {
        return 3;
    }//about 3.8

    virtual int timePerIteration() const {
        return 1500;
    } //22 microsecs / 16 ns = 1500
};
