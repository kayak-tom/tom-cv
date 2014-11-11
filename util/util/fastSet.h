/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

/*
 * fastSet.h
 *
 *  Created on: 18/06/2009
 *      Author: tom
 */

#ifndef FASTSET_H_
#define FASTSET_H_

#define VERBOSE(x)
VERBOSE(#include "BoWSLAM_Main/geom/geom.h")

template<typename TKey, int FIXED_SIZE>
class CFastSet
{
    //friend class CFastSet<TKey, FIXED_SIZE>::const_iterator;
    TKey * aDynArr, * pEnd;
    int count, possibleSize;
    const TKey EMPTY;
    TKey arr[FIXED_SIZE];

    void compact()
    {
        //Rearrange
        TKey * pElement = aDynArr;
        for(TKey * pElementToCopyFrom = aDynArr; pElementToCopyFrom < pEnd; pElementToCopyFrom++)
        {
            if(!(*pElementToCopyFrom == EMPTY))
            {
                *pElement = *pElementToCopyFrom;
                pElement++;
            }
        }
        pEnd = pElement;
        VERBOSE(std::cout << "Compacting";)
    }
public:
    CFastSet(const CFastSet & other) : aDynArr(arr), pEnd(arr), count(other.count), possibleSize(other.possibleSize), EMPTY(other.EMPTY)
    {
        if(possibleSize > FIXED_SIZE)
            aDynArr = new TKey[possibleSize];

        TKey * pThisVal = aDynArr;
        for(const TKey * pOtherVal = other.aDynArr; pOtherVal < other.pEnd; pOtherVal++)
        {
            if(!(*pOtherVal == EMPTY))
            {
                *pThisVal = *pOtherVal;
                pThisVal++;
            }
        }
        pEnd = pThisVal;
    }

    CFastSet(const TKey & empty=0) : aDynArr(arr), pEnd(arr), count(0), possibleSize(FIXED_SIZE), EMPTY(empty)
    {
    }

    /*CFastSet() : aDynArr(arr), pEnd(arr), count(0), possibleSize(FIXED_SIZE), EMPTY(0)
    {
    }*/

    ~CFastSet()
    {
        if(possibleSize>FIXED_SIZE)
            delete [] aDynArr;
    }

    /*class const_iterator
    {
        TKey * p;
        CFastSet<TKey, FIXED_SIZE> * parent;

    public:
        const_iterator() : p(0), parent(0) {};
        void operator=

        void operator++()
        {
            do
            {
                p++;
            } while(p < parent.pEnd && *p == parent.END);
        }

        const TKey & operator*() const { return *p; }
    };*/

    //typedef TKey * iterator;
    typedef const TKey * const_iterator;

    //iterator begin() { return aDynArr; }
    //iterator end() { return pEnd; }
    const_iterator begin()
    {
        if(count < pEnd - aDynArr)
        {
            compact();
            VERBOSE(std::cout << "Compacted:"; printMe();)
        }
        return aDynArr;
    }
    const_iterator end() const { return pEnd; }

    int size() const { return count; }

    void insert(const TKey & key)
    {
        if(key == EMPTY)
        {
            VERBOSE(std::cout << "Inserting EMPTY\n";)
            return;
        }
        VERBOSE(std::cout << "Inserting " << key;)
        if(!exists(key))
        {
            if(pEnd == aDynArr+possibleSize)
            {
                if(count == possibleSize)
                {
                    //Realloc
                    int possibleSizeNew = possibleSize * 2;
                    TKey * aDynArrNew = new TKey[possibleSizeNew];
                    TKey * pDynArrNew = aDynArrNew;
                    for(TKey * pElement = aDynArr; pElement < pEnd; pElement++)
                    {
                        *pDynArrNew = *pElement;
                        pDynArrNew++;
                    }
                    if(possibleSize > FIXED_SIZE)
                        delete [] aDynArr;
                    aDynArr = aDynArrNew;
                    pEnd = pDynArrNew;
                    possibleSize = possibleSizeNew;
                }
                else
                {
                    compact(); //alternatively copy into a gap
                }
            }
            *pEnd = key;
            VERBOSE(if(!(*pEnd == key && key == *pEnd) || (EMPTY == *pEnd || *pEnd == EMPTY))
                std::cout << "ERROR" << std::endl;
            if( *((int *)(void*)pEnd) != *((int *)(void*)&key))
                std::cout << "ERROR hack" << std::endl;)
            pEnd++; count++;
        }
        VERBOSE(printMe();)
    }

    bool exists(const TKey & key) const
    {
        for(TKey * pElement = aDynArr; pElement < pEnd; pElement++)
            if(*pElement == key) return true;

        return false;
    }

    void erase(const TKey & key)
    {
        if(key == EMPTY)
        {
            VERBOSE(std::cout << "Erasing EMPTY\n";)
            return;
        }
        for(TKey * pElement = aDynArr; pElement < pEnd; pElement++)
            if(*pElement == key)
            {
                VERBOSE(std::cout << "Erasing " << key << endl;)
                *pElement = EMPTY;
                count--;
                return;
            }

        return;
    }
    VERBOSE(void printMe()
    {
        std::cout << "# Elements=" << count << ':';
        for(TKey * pElement = aDynArr; pElement < pEnd; pElement++)
            std::cout << *pElement << ',';
        std::cout << std::endl;
    })
};

#endif /* FASTSET_H_ */
