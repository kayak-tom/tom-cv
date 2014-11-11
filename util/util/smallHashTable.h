/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

/*
 * fastSet.h
 *
 *  Created on: 18/06/2009
 *      Author: tom
 */

#ifndef SMALLHASH_H_
#define SMALLHASH_H_

#define VERB(x)
VERB(#include "BoWSLAM_Main/geom/geom.h")

template<class THashTable>
class ht_const_iterator;

// # bins is fixed, binsize grows. THash should return a hash between 0 and BINS-1.
template<typename TKey, int BINS, int BINSIZE, class THash>
class CSmallHashTable
{
    struct SBin {
        int hash; TKey val;
        bool operator==(const SBin & s) const { return hash == s.hash && val == s.val; }
    };
    SBin * aDynArr;
    int count, nActualBinsize;
    SBin arr[BINS*BINSIZE];
    int anBinCounts[BINS];
    static const int EMPTY = BINS;
    THash hashFn;

    SBin * binStartPtr(int nBin) const { return aDynArr + (nBin * nActualBinsize); }
    const SBin * endPtr() const { return aDynArr + (BINS * nActualBinsize); }
public:
    typedef ht_const_iterator<CSmallHashTable<TKey, BINS, BINSIZE, THash> > const_iterator;
    friend class ht_const_iterator<CSmallHashTable<TKey, BINS, BINSIZE, THash> >;

    CSmallHashTable(const CSmallHashTable & other) : aDynArr(arr), count(other.count), nActualBinsize(other.nActualBinsize)
    {
        if(nActualBinsize > BINSIZE)
            aDynArr = new SBin[nActualBinsize*BINS];

        SBin * pThisVal = aDynArr;
        const SBin * pEnd = other.endPtr();
        for(SBin * pOtherVal = other.aDynArr; pOtherVal < pEnd; pThisVal++, pOtherVal++)
            *pThisVal = *pOtherVal;

        for(int bin = 0; bin < BINS; bin++)
            anBinCounts[bin] = other.anBinCounts[bin];
    }

    CSmallHashTable() : aDynArr(arr), count(0), nActualBinsize(BINSIZE)
    {
        //set to empty:
        int * pnBinCount = anBinCounts;
        for(int bin=BINS; bin > 0; bin--)
        {
            *pnBinCount = 0;
            pnBinCount++;
        }

        SBin * pDynArr = aDynArr;
        for(int i=BINS*nActualBinsize;i>0;i--)
        {
            pDynArr->hash = EMPTY;
            pDynArr++;
        }
    }

    ~CSmallHashTable()
    {
        if(nActualBinsize > BINSIZE)
            delete [] aDynArr;
    }

    const_iterator begin() const
    {
        return const_iterator(this);
    }
    const_iterator end() const
    {
        return const_iterator(0);
    }

    int size() const { return count; }

    //Returns true if already exists
    bool insert(const TKey & newVal)
    {
        if(count*4 > nActualBinsize*BINS*3)
        {
            //Double size:
            int nNewBinsize = nActualBinsize * 2;
            SBin * aNewDynArr = new SBin[nNewBinsize*BINS];

            SBin * pNewDynArr = aNewDynArr;
            for(int i=nNewBinsize*BINS;i>0;i--)
            {
                pNewDynArr->hash = EMPTY;
                pNewDynArr++;
            }

            //Insert all existing elements into new array
            const SBin * pEnd = endPtr();
            const SBin * pNewBinEnd = aNewDynArr + nNewBinsize*BINS;
            for(SBin * pBin = aDynArr; pBin < pEnd; pBin++)
            {
                if(pBin->hash == BINS)
                    continue; //empty in old

                SBin * pNewBin = aNewDynArr + (pBin->hash * nNewBinsize);
                do
                {
                    if(pNewBin->hash == BINS) //empty in new
                    {
                        *pNewBin = *pBin;
                        break;
                    }
                    pNewBin++;
                    if(pNewBin == pNewBinEnd) pNewBin = aNewDynArr;
                } while(true);
            }
            if(nActualBinsize > BINSIZE)
                delete [] aDynArr;
            aDynArr = aNewDynArr;
            nActualBinsize = nNewBinsize;
        }

        const SBin * pEnd = endPtr();
        const int bin = hashFn(newVal);
        SBin * pBin = binStartPtr(bin), * pEmptyBin = 0;
        int nFound = 0;
        const int nBinCount = anBinCounts[bin];
        do
        {
            if(pBin->hash == bin)
            {
                if(pBin->val == newVal) return true; // already exists
                nFound++;
            } else
                if(pEmptyBin == 0 && pBin->hash == EMPTY)
                    pEmptyBin = pBin;
            pBin++;
            if(pBin == pEnd) pBin = aDynArr;
        } while(nFound < nBinCount || pEmptyBin == 0);

        pEmptyBin->hash = bin, pEmptyBin->val = newVal;
        anBinCounts[bin]++;
        count++;
        return false;
    }

    bool exists(const TKey & val) const
    {
        const int bin = hashFn(val);
        int nFound = 0;
        const int nBinCount = anBinCounts[bin];
        const SBin * pEnd = endPtr();
        for(SBin * pBin = binStartPtr(bin); nFound < nBinCount; pBin++)
        {
            if(pBin == pEnd) pBin = aDynArr;
            if(pBin->hash == bin)
            {
                if(pBin->val == val) return true; // exists
                nFound++;
            }
        }
        return false;
    }

    void erase(const TKey & val)
    {
        const int bin = hashFn(val);
        int nFound = 0;
        const int nBinCount = anBinCounts[bin];
        const SBin * pEnd = endPtr();
        for(SBin * pBin = binStartPtr(bin); nFound < nBinCount; pBin++)
        {
            if(pBin == pEnd) pBin = aDynArr;
            if(pBin->hash == bin)
            {
                if(pBin->val == val)
                {
                    pBin->hash = EMPTY;
                    anBinCounts[bin]--;
                    count--;
                    return;
                }
                nFound++;
            }
        }
    }

    VERB(void printMe()
    {
        std::cout << "# Elements=" << count << ':';
        for(TKey * pElement = aDynArr; pElement < pEnd; pElement++)
            std::cout << *pElement << ',';
        std::cout << std::endl;
    })

    typedef TKey THTKey;
};

template<class THashTable>
class ht_const_iterator
{
    struct THashTable::SBin * p;
    const THashTable * parent;

public:
    ht_const_iterator(const THashTable * pParent = 0) : p(0), parent(pParent)
    {
        if(parent && parent->count > 0)
        {
            p = parent->aDynArr;
            while(p->hash==THashTable::EMPTY) p++;
        }
        else
            parent = 0;
    }

    const typename THashTable::THTKey & operator*() const
    {
        return p->val;
    }

    bool operator!=(const ht_const_iterator<THashTable> & it) const
    {
        return it.p != p || it.parent != parent;
    }

    void operator++()
    {
        const typename THashTable::SBin * pEnd = parent->endPtr();
        do
        {
            p++;
            if(p == pEnd)
            {
                p = 0;
                parent = 0;
                return;
            }
        } while(p->hash == THashTable::EMPTY);
    }
};

template<int BINS>
class intHash
{
public:
    int operator()(int loc) const
    {
        return loc % BINS;
    }
};

#endif /* SMALLHASH_H_ */
