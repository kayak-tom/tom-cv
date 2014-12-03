/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

/*
 * CDynArray.h
 *
 *  Created on: 27/08/2009
 *      Author: tom
 */

#ifndef CDynArray_H_
#define CDynArray_H_

#include<vector>
#include "exception.h"
#include "optimisation_attributes.h"
#include<boost/noncopyable.hpp>

/*namespace std
{
    class ostream;
}*/

template<class T>
class CDynArray {
protected:
    T * aData;
    int nCount, nAllocated;
private:
    
    //Disable copying
    CDynArray(const CDynArray &) {
    }

    const CDynArray & operator=(const CDynArray &) {
    }
protected:

    T * alloc(int nNewSize) {
        int size = sizeof (T) * nNewSize;
        if (size == 0) return 0;

        char * pc = new char[size];

        T * mem = reinterpret_cast<T *> (pc);

        for (int i = size; i > 0; i--) {
            *pc = 0;
            pc++;
        }

        /*int * pData = reinterpret_cast<int *>(mem);
        for(int i=size; i>0; i--)
            {*pData=0; pData++;}*/

        return mem;
    };

    void increaseSize() {
        const int nNewSize = nAllocated > 32 ? nAllocated * 2 : 64;
        increaseSize(nNewSize);
    }

    void increaseSize(int nNewSize) {
        T * aNewData = alloc(nNewSize);
        for (int i = 0; i < nCount; i++)
            aNewData[i] = aData[i];

        deleteData();
        aData = aNewData;
        nAllocated = nNewSize;
    }

    void deleteData() {
        for (const_iterator p = begin(); p < end(); p++)
            p->~T(); //Destruct data

        delete [] reinterpret_cast<char *> (aData);
        aData = 0;
    }

public:

    CDynArray(int size) : aData(alloc(size)), nCount(size), nAllocated(size) { //Todo: get rid of this?? Should be ok because calls default c'tor for each member
        //        for(int i=0;i<nCount;i++)
        //            aData[i] = T();
        new(reinterpret_cast<char *> (aData)) T[nCount];

        /*
        #ifdef _GLIBCXX_IOSTREAM
                if(size>0)
                    std::cout << "New dynarray with uninit members...\n";
        #endif*/
    }

    CDynArray(int size, const T & t) : aData(alloc(size)), nCount(size), nAllocated(size) {
        //for(int i=0;i<nCount;i++)
        //  aData[i] = t;
        setConst(t);
    }

    CDynArray() : aData(0), nCount(0), nAllocated(0) {
    }

    ~CDynArray() {
        deleteData();
    }

    int size() const {
        return nCount;
    }

    void setConst(const T & t) {
        for (iterator pd = begin(); pd != end(); pd++)
            *pd = t;
    }

    //Return mean and population var MLEs

    void getMeanVar(T & dMean, T & dVar) const {
        dMean = 0;
        dVar = 0;

        if (size() < 2) {
            if (size() == 1)
                dMean = *begin();
            return;
        }

        for (const_iterator pd = begin(); pd != end(); pd++)
            dMean += *pd;

        dMean /= size();

        for (const_iterator pd = begin(); pd != end(); pd++) {
            T diff = (*pd - dMean);
            dVar += diff*diff;
        }

        dVar /= (size() - 1);
    }

    double sum() const {
        double dSum = 0;

        for (const_iterator pd = begin(); pd != end(); pd++)
            dSum += *pd;

        return dSum;
    }

    //Return min and max. Must be non-empty

    T & minimum() {
        T & dMin = *begin();

        for (const_iterator pd = begin(); pd != end(); pd++) {
            if (*pd < dMin)
                dMin = *pd;
        }
        return dMin;
    }

    T & maximum() {
        T & dMax = *begin();

        for (const_iterator pd = begin(); pd != end(); pd++) {
            if (*pd > dMax)
                dMax = *pd;
        }
        return dMax;
    }

    /*const T & front() const
    {
        if(IS_DEBUG) CHECK(size()==0, "front: Empty CDynArray");
        return *begin();
    }
    const T & back() const
    {
        if(IS_DEBUG) CHECK(size()==0, "back: Empty CDynArray");
        return *(end()-1);
    }*/

    //for debugging (O(n))

    bool contains(const T & el) const {
        for (const_iterator p = begin(); p < end(); p++)
            if (*p == el)
                return true;

        return false;
    }
    //for debugging (O(n))

    int count(const T & el) const {
        int nCount = 0;
        for (const_iterator p = begin(); p < end(); p++)
            if (*p == el)
                nCount++;

        return nCount;
    }

    //    template<class _Cmp=std::less<T> >
    //    void sort() { _Cmp cmp; std::sort(begin(), end(), cmp() ); }

    typedef T const * const_iterator;
#ifndef _DEBUGxx
    typedef T * iterator;
#else

    class iterator {
        T * ptr;

        iterator() {
        };
    public:

        iterator(T * ptr) : ptr(ptr) {
        };

        iterator(const iterator i) : ptr(i.ptr) {
        };

        void operator=(T * ptr_in) {
            ptr = ptr_in;
        };

        void operator=(const iterator i) {
            ptr = i.ptr;
        };

        T & operator*() {
            if(IS_DEBUG) CHECK(!ptr, "Dereferencing uninit iterator");
            return *ptr;
        }

        const T & operator*() const {
            if(IS_DEBUG) CHECK(!ptr, "Dereferencing uninit iterator");
            return *ptr;
        }

        T & operator->() {
            if(IS_DEBUG) CHECK(!ptr, "Dereferencing uninit iterator");
            return *ptr;
        }

        const T & operator->() const {
            if(IS_DEBUG) CHECK(!ptr, "Dereferencing uninit iterator");
            return *ptr;
        }
    };
#endif

    inline iterator begin() HARD_INLINE {
        return aData;
    }

    inline iterator end() HARD_INLINE {
        return begin() + size();
    }

    inline const_iterator begin() const HARD_INLINE {
        return aData;
    }

    inline const_iterator end() const HARD_INLINE {
        return begin() + size();
    }

    inline T & operator[](int i) HARD_INLINE;
    inline const T & operator[](int i) const HARD_INLINE;
    inline const T & get(int i) const HARD_INLINE;

    const T & back() const {
        if(IS_DEBUG) CHECK(size() == 0, "Empty dynArray");
        return *(end() - 1);
    }

    T & back() {
        if(IS_DEBUG) CHECK(size() == 0, "Empty dynArray");
        return *(end() - 1);
    }

    const T & top() const {
        if(IS_DEBUG) CHECK(size() == 0, "Empty dynArray");
        return *(begin());
    }

    T & top() {
        if(IS_DEBUG) CHECK(size() == 0, "Empty dynArray");
        return *(begin());
    }

    inline void push_back(const T & val) HARD_INLINE {
        if (nCount >= nAllocated)
            increaseSize();

        aData[nCount] = val;
        nCount++;
    }

    void push_back(T & val) //no const, for when copying destroys original
    {
        if (nCount >= nAllocated)
            increaseSize();

        aData[nCount] = val;
        nCount++;
    }

    template<typename TIter>
    void copy_back(TIter p, const TIter p_end)
    //    void push_back(std::iterator<std::forward_iterator_tag, T> & p, const std::iterator<std::forward_iterator_tag, T> & p_end)
    {
        reserve(nCount + (int)(p_end - p));

        for (; p != p_end; p++, nCount++)
            aData[nCount] = *p;
    }

    void resize(int nNewSize) {
        if (nNewSize > nAllocated) {
            increaseSize(nNewSize);
            /*#ifdef _GLIBCXX_IOSTREAM
                        std::cout << "Resized dynarray with uninit members...\n";
            #endif*/
        }
        for (int i = nCount; i < nNewSize; i++) //Call default c'tor for each new element
            aData[i] = T();

        for (int i = nNewSize; i < nCount; i++)
            aData[i].~T(); //Destruct data

        nCount = nNewSize;
    }

    void erase(iterator pEl) {
        nCount--;
        for (iterator p = pEl; p < end(); p++)
            *p = *(p + 1);

        end()->~T();
    }
    //Inefficient, erase all elements val

    void erase(const T & val) {
        for (iterator p = begin(); p < end(); p++)
            if (*p == val) {
                erase(p); //This will reduce end(),
                p--; //in case the next element is also p;
            }

    }

    //Remove zeros

    void cleanup() {
        iterator p = begin(), pNew = begin();
        for (; p < end(); p++, pNew++)
            while (*p == 0) {
                p++;
                if (p == end())
                    goto endLine;
            }
        *pNew = *p;

endLine:
        ;
        int nNewCount = pNew - begin();

        for (; pNew < end(); pNew++)
            pNew->~T();

        nCount = nNewCount;
    }

    void pop(const int nRemove=1) {
        for (iterator p = end() - nRemove; p < end(); p++)
            p->~T();

        nCount -= nRemove;
    }

    //Remove nRemove off front. Could just maintain 2 pointers

    void popN(const int nRemove) {
        if (nRemove > 0) {
            if (nRemove >= size())
                clear();
            else {
                for (iterator p = begin(), pOld = begin() + nRemove; pOld < end(); p++, pOld++)
                    *p = *pOld;

                nCount -= nRemove;
            }
        }
    }

    void clear() {
        for (const_iterator p = begin(); p < end(); p++)
            p->~T(); //Destruct data

        nCount = 0;
    }

    void reserve(int nNewSize) {
        if (nNewSize > nAllocated) {
            increaseSize(nNewSize);
        }
    }

    /* Deprecate these--causing issues in Windows
	void apply(T f(T)) {
        for (iterator p = begin(); p < end(); p++)
            *p = f(*p);
    }

    void apply(T f(const T &)) {
        for (iterator p = begin(); p < end(); p++)
            *p = f(*p);
    }*/

    void copyInto(CDynArray<T> & dest) const {
        dest.resize(size());

        iterator pCopy = dest.begin();
        for (const_iterator pThis = begin(); pThis < end(); pThis++, pCopy++)
            *pCopy = *pThis;

        //if(IS_DEBUG) CHECK(dest.countInliers() != countInliers(), "Copy error");
        //if(IS_DEBUG) CHECK(mask.checksum() != checksum(), "Copy error");
    }

        /*    void operator=(CDynArray<T> & source)
        {
            delete aData;
            aData = source.aData;
            nCount = source.nCount;
            nAllocated = source.nAllocated;

            source.aData = 0;
            source.nCount = source.nAllocated = 0;

    #ifdef _GLIBCXX_IOSTREAM
            std::cout << "Warning, vector copy has wiped original. Use pointers instead\n";
    #endif
        }*/

    template<typename U>
    operator CDynArray<U> & () {
        //Check vectors of U and T are compatible:
        CHECK_SIZES_EQUAL_RT(T, U, CDynArray_typecast);
#ifndef __clang__
        U & TYPES_MUST_BE_COMPATIBLE(T()); //Compiler should get rid of this
#endif
        return reinterpret_cast<CDynArray<U> &> (*this);
    }

    template<typename U>
    operator const CDynArray<U> & () const {
        //Check vectors of U and T are compatible:
        CHECK_SIZES_EQUAL_RT(T, U, CDynArray_typecast);
#ifndef __clang__
        U & TYPES_MUST_BE_COMPATIBLE(T()); //Compiler should get rid of this
#endif
        return reinterpret_cast<const CDynArray<U> &> (*this);
    }

    void pp(const char * seperator = ", ", std::ostream & s = std::cout) const {
        for (const_iterator pEl = begin(); pEl != end(); pEl++) {
            if (pEl != begin())
                s << seperator;
            s << *pEl;
        }
        s << std::endl;
    }

    bool operator==(const CDynArray<T> & other) const {
        if (size() != other.size())
            return false;

        for (const_iterator pEl = begin(), pEl2 = other.begin(); pEl != end(); pEl++, pEl2++) {
            if (*pEl != *pEl2)
                return false;
        }
        return true;
    }
};

template<typename T>
inline T & CDynArray<T>::operator[](int i) {
    if(IS_DEBUG) CHECK(i < 0 || i >= size(), "CDynArray::operator[]: Idx OOB")
    return aData[i];
}

template<typename T>
inline const T & CDynArray<T>::operator[](int i) const {
    if(IS_DEBUG) CHECK(i < 0 || i >= size(), "CDynArray::operator[]: Idx OOB")
    return aData[i];
}

template<typename T>
inline const T & CDynArray<T>::get(int i) const {
    if(IS_DEBUG) CHECK(i < 0 || i >= size(), "CDynArray::get: Idx OOB")
    return aData[i];
}

//An array of pointers OWNING THEIR MEMORY!

template<class T>
class CDynArrayOwner : public CDynArray<T *>, boost::noncopyable {
    typedef CDynArray<T *> TDynArray;
protected:

    void deleteAll(int nFromHere = 0) {
        for (typename TDynArray::iterator ppX = TDynArray::begin() + nFromHere; ppX < TDynArray::end(); ppX++) {
            delete *ppX;
            *ppX = 0;
        }
    }
public:

    CDynArrayOwner(int size) : TDynArray(size, 0) {
    };
    //CDynArrayOwner(int size, const T & t) Do not init this way!

    CDynArrayOwner() : TDynArray() {
    }

    ~CDynArrayOwner() {
        deleteAll();
        //Todo: TDynArray::~CDynArray();
        TDynArray::deleteData();
    }

    void resize(int nNewSize) {
        if (nNewSize > TDynArray::nAllocated) {
            TDynArray::increaseSize(nNewSize);
        } else
            deleteAll(nNewSize);

        for (int i = TDynArray::nCount; i < nNewSize; i++)
            TDynArray::aData[i] = 0;

        TDynArray::nCount = nNewSize;
    }

    void clear() {
        deleteAll();
        TDynArray::clear();
    }
};

template<class T, int MAX_LEN>
class CFixedArray {
    int nCount;
    T aData[MAX_LEN];

    template<class TPred>
    int minmaxIdx() const {
        if(IS_DEBUG) CHECK(nCount == 0, "maxIdx: Array empty");

        TPred pred;

        const_iterator p = begin();
        const_iterator max_p = begin();
        p++;
        for (; p != end(); p++) {
            const T & el = *p;
            if (pred(el, *max_p)) {
                max_p = p;
            }
        }
        int maxIdx = max_p - begin();

        return maxIdx;
    }
public:

    CFixedArray() : nCount(0) {
    };

    CFixedArray(int nCount, const T & t) : nCount(nCount) {
        for (iterator p = begin(); p != end(); p++)
            *p = t;
    };

    typedef T const * const_iterator;
    typedef T * iterator;

    int maxIdx() const //index of max element
    {
        return minmaxIdx<std::greater<T> >();
    }

    int minIdx() const //index of min element
    {
        return minmaxIdx<std::less<T> >();
    }

    int size() const {
        return nCount;
    }

    iterator begin() HARD_INLINE {
        return aData;
    }

    iterator end() HARD_INLINE {
        return begin() + size();
    }

    const_iterator begin() const HARD_INLINE {
        return aData;
    }

    const_iterator end() const HARD_INLINE {
        return begin() + size();
    }

    T & operator[](int i) HARD_INLINE {
        return aData[i];
    }

    const T & operator[](int i) const HARD_INLINE {
        return aData[i];
    }

    const T & back() const {
        if(IS_DEBUG) CHECK(size() == 0, "Empty array");
        return *(end() - 1);
    }

    T & back() {
        if(IS_DEBUG) CHECK(size() == 0, "Empty array");
        return *(end() - 1);
    }

    const T & top() const {
        if(IS_DEBUG) CHECK(size() == 0, "Empty array");
        return *(begin());
    }

    T & top() {
        if(IS_DEBUG) CHECK(size() == 0, "Empty array");
        return *(begin());
    }

    void push_back(const T & val) HARD_INLINE {
        if(IS_DEBUG) CHECK(nCount >= MAX_LEN, "No more room in CFixedArray")

        aData[nCount] = val;
        nCount++;
    }

    void clear() {
        nCount = 0;
    }
    
    void resize(int nNewSize) 
    {
        nCount = nNewSize; 
        if(IS_DEBUG) CHECK(nCount >= MAX_LEN || nCount < 0, "No more room in CFixedArray")
    }
};


#endif /* CDynArray_H_ */
