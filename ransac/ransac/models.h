/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

/*
 * models.h
 *
 *  Created on: 8/10/2009
 *      Author: tom
 */

#ifndef MODELS_H_
#define MODELS_H_

#include "util/exception.h"
#include "util/dynArray.h"

//Base for one model (e.g. one essential matrix or homography)
class CModel {
public:

    virtual ~CModel() {
    };
    virtual void copyInto(CModel & dest) const = 0;
};


class C3x3MatModel : public CModel {
    double aaData[9];
public:
    virtual ~C3x3MatModel() {
    };

    double operator[](int n) const {
        return aaData[n];
    }

    double & operator[](int n) {
        return aaData[n];
    }

    double & operator()(int r, int c) {
        return aaData[3 * r + c];
    }

    double operator()(int r, int c) const {
        return aaData[3 * r + c];
    }

    const double * asDouble9() const {
        return reinterpret_cast<const double*> (aaData);
    }

    double * asDouble9() {
        return reinterpret_cast<double*> (aaData);
    }

    virtual void copyInto(CModel & dest) const {
        dynamic_cast<C3x3MatModel &> (dest) = *this;
    }
};

class CModels {
public:
    virtual int numModels() const = 0;
    virtual CModel & addModel() = 0;
    virtual const CModel & getData(int n) const = 0;

    virtual ~CModels() {
    };
};

//A set of models returned by a hypothesis generator

template<typename TModelType, const int MAX_MODELS>
class CModelSet : public CModels {
    TModelType aModels[MAX_MODELS];
    int nNumModels;
public:

    int numModels() const {
        return nNumModels;
    }

    CModelSet() {
        nNumModels = 0;
    }

    void reset() {
        nNumModels = 0;
    }

    virtual CModel & addModel() {
        if(IS_DEBUG) CHECK(nNumModels >= MAX_MODELS, "getData: Too many models");
        return aModels[nNumModels++];
    }

    virtual const CModel & getData(int n) const {
        if(IS_DEBUG) CHECK(n < 0 || n >= nNumModels, "getData: Index OOB");
        return aModels[n];
    }
};

typedef CModelSet<C3x3MatModel, 10 > T3x3MatModels;

template<typename TModelType>
class CModelVector : public CModels {
    std::vector<TModelType> aModels;
public:

    int numModels() const {
        return aModels.size();
    }

    CModelVector() {
        aModels.reserve(100);
    }

    void reset() {
        aModels.clear();
    }

    virtual CModel & addModel() {
        aModels.push_back(C3x3MatModel());
        return aModels[aModels.size() -1];
    }

    virtual const CModel & getData(int n) const {
        return aModels[n];
    }
};

typedef CModelVector<C3x3MatModel> TEModels;

#endif /* MODELS_H_ */
