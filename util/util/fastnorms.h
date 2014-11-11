/*
 * fastnorms.h
 *
 *  Created on: 17/12/2009
 *      Author: hilandtom
 */

#pragma once

#ifndef FASTNORMS_H_
#define FASTNORMS_H_

//TODO: use <mmintrin.h>

#include <iostream>
#include "convert.h"

enum eNormParallel {eL1, eSSD};

template <typename T, eNormParallel norm>
int L1dist(T * ac1, T * ac2, const int LENGTH) HOT /*HARD_INLINE*/;

template <typename T, eNormParallel norm>
int L1dist(T * ac1, T * ac2, const int LENGTH)
{
    int nDist = 0;
    for(int i=LENGTH; i>0; i--)
    {
        //cout << (int)*ac1 << ' ' << (int)*ac2 << ' ';
        volatile int dist = (int)*ac1-(int)*ac2;
        //cout << i << ':' << dist << endl;

        if(norm == eL1)
            nDist += dist > 0 ? dist : -dist;
        else if(norm == eSSD)
            nDist += dist*dist;

        ac1++;
        ac2++;
    }
    return nDist;
}

//template<typename T>
//inline T sqr(T x) {return x*x;}

template<class T>
int sumSquares(T * ac1, const int LENGTH)
{
    int nSS = 0;

    for(int i=LENGTH; i>0; i--)
    {
        nSS += sqr((int)*ac1);
        ac1++;
    }
    return nSS;
}


template<class TInt>
inline int L2distParallel(const unsigned char * ac1, const unsigned char * ac2, const int LENGTH) HOT HARD_INLINE;

template<class TInt>
inline int L2distParallel(const unsigned char * ac1, const unsigned char * ac2, const int LENGTH)
{
    //typedef char TInt_MUST_BE_UNSIGNED[((TInt)(-1))];

    const int NEWLENGTH = LENGTH/sizeof(TInt);

    int nDist = 0, nDist2=0;

    TInt * an1 = (TInt *)(void *)ac1;
    TInt * an2 = (TInt *)(void *)ac2;

    for(int i=NEWLENGTH; i>0; i--)
    {
        //cout << hex << *an1 << ' ' << *an2 << dec << ' ';
        TInt nSums = *an1 + *an2;
        unsigned char * acSums = (unsigned char *)(void *)&nSums;
        for(int j=0; j<(int)sizeof(TInt)/2; j++)
        {
            /*const int sum = (int)*acSums;

            / *int idx = (NEWLENGTH-i)*sizeof(TInt) + j;
            cout << i << ':' << dist << '=' << (int)ac1[idx] << '+' << (int)ac2[idx] << endl;* /
            volatile int distSq = sqr(sum);
            nDist += distSq; //If this is slow then accumulate 2 variables instead and add afterwards
            acSums++;*/
            volatile const int sum = (int)*acSums;
            volatile int distSq = sqr(sum);
            nDist += distSq;

            volatile const int sum2 = (int)*(acSums+1);
            volatile int distSq2 = sqr(sum2);
            nDist2 += distSq2;

            acSums += 2;
        }
        an1++;
        an2++;
    }

    nDist += nDist2;
    unsigned char * pc1=(unsigned char *)(void *)an1;
    unsigned char * pc2=(unsigned char *)(void *)an2;

    for(int i=LENGTH - NEWLENGTH*sizeof(TInt); i>0; i--, pc1++, pc2++)
    {
        const int dist = (int)*pc1+(int)*pc2;
        nDist += sqr(dist);
    }

    return nDist;
}

inline int L2distParallel(const unsigned char * ac1, const unsigned char * ac2, const int LENGTH, const int nac1SS2, const int nac2SS2)  HOT HARD_INLINE;

inline int L2distParallel(const unsigned char * ac1, const unsigned char * ac2, const int LENGTH, const int nac1SS2, const int nac2SS2)
{
    int nSumSquared = 0;
    if(is64bit())
        nSumSquared = L2distParallel<unsigned long long>(ac1, ac2, LENGTH);
    else
        nSumSquared = L2distParallel<unsigned int>(ac1, ac2, LENGTH);

    //L2distParallel returns a.^2 + 2a.b + b.^2
    //We want a.^2 - 2a.b + b.^2 = 2a.^2 + 2b.^2 - L2distCos
    return nac1SS2 + nac2SS2 - nSumSquared;
}

inline void GET_SETUP_MASK(unsigned long long & SETUP)
{
    SETUP = 0x8080808080808080LL;
}

inline void GET_SETUP_MASK(unsigned int & SETUP)
{
    std::cout << sizeof(unsigned long long) << std::endl;
    SETUP = 0x80808080;
}

template<class TInt, eNormParallel norm>
inline int L1distNew(const unsigned char * ac1, const unsigned char * ac2, const int LENGTH)
{
    //typedef char TInt_MUST_BE_UNSIGNED[((TInt)(-1))];

    const int NEWLENGTH = LENGTH/sizeof(TInt);

    int nDist = 0;

    TInt * an1 = (TInt *)(void *)ac1;
    TInt * an2 = (TInt *)(void *)ac2;

    TInt SETUP;
    GET_SETUP_MASK(SETUP);

    for(int i=NEWLENGTH; i>0; i--)
    {
        //cout << hex << *an1 << ' ' << *an2 << dec << ' ';
        TInt diff = SETUP + *an1 - *an2;
        char * acDiffs = (char *)(void *)&diff;
        for(int j=0; j<(int)sizeof(TInt);j++)
        {
            const int dist = (int)*acDiffs; //dist = 128 +- actual_diff (so subtraction doesn't overflow)

            /*int idx = (NEWLENGTH-i)*sizeof(TInt) + j;
            cout << i << ':' << dist << '=' << (int)ac1[idx] << '-' << (int)ac2[idx] << endl;*/
            if(norm == eL1)
                //nDist += dist > 0 ? dist : -dist;
                dist > 0 ? nDist += dist : nDist -= dist;
            else if(norm == eSSD)
            {
                //const int dist2 = dist>0 ? 128-dist : 128+dist;
                //nDist += dist2*dist2;

                nDist += sqr(dist>0 ? 128-dist : 128+dist);

                //Change type of acDiffs to unsigned char to use this:
                //int dist2 = dist>128 ? (int)(char)(dist+128) : (int)(128 - dist);
                //nDist += sqr(dist2);
            }

            /*if(nDist<0)
                cout << nDist << ' ' << dist << endl;*/

            acDiffs++;
        }
        an1++;
        an2++;
    }

    unsigned char * pc1=(unsigned char * const)(void *)an1;
    unsigned char * pc2=(unsigned char * const)(void *)an2;

    for(int i=LENGTH - NEWLENGTH*sizeof(TInt); i>0; i--, pc1++, pc2++)
    {
        const int dist = (int)*pc1-(int)*pc2;
        if(norm == eL1)
            dist > 0 ? nDist -= dist : nDist += dist; //Note sign change
        else
            nDist += sqr(dist);
    }

    if(norm == eL1)
        return NEWLENGTH*sizeof(TInt)*128 - nDist;
    else
        return nDist;
}

inline int L1distParallel(const unsigned char * ac1, const unsigned char * ac2, const int LENGTH) HOT HARD_INLINE;
inline int L1distParallel(const unsigned char * ac1, const unsigned char * ac2, const int LENGTH)
{
    if(is64bit())
        return L1distNew<unsigned long long, eL1>(ac1, ac2, LENGTH);
    else
        return L1distNew<unsigned int, eL1>(ac1, ac2, LENGTH);
}

#endif /* FASTNORMS_H_ */
