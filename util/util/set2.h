/*
 * set2.h
 *
 *  Created on: 20 Jan 2010
 *      Author: tom
 */

#ifndef SET2_H_
#define SET2_H_

#include <set>
#include <map>
#include "exception.h"
//For MS intrinsics... #include <Intrin.h>

#define DEBUG_LOG(...) // DEBUGONLY(__VA_ARGS__)

class logUse
{
public:
    static void log(const char * szUser, int nCountUses)
    {
        if(!nCountUses)
            std::cout << "Warning: " << szUser << " Pool never used\n";
        else if(nCountUses<0)
            log(szUser, nCountUses+1);
        else //if(nCountUses < N/10)
            std::cout << szUser << " pool used " << nCountUses << " times\n";
    }
};

template<typename T, int N>
class CMemPool {

    class CMemBlock
    {
        char pMem[sizeof(T)];
        CMemPool<T, N> * pParentIfOccupied; //Indicated both occupies status and a pointer to the parent (allows very fast erasing)
    public:
        inline CMemBlock() : pParentIfOccupied(0) {}

        inline bool occupied() const { return pParentIfOccupied != 0; }
        inline T * allocMem(CMemPool<T, N> * pParent)
        {
            if(IS_DEBUG) CHECK(pParentIfOccupied, "Allocating alloc'd mem");
            pParentIfOccupied = pParent;
            return static_cast<T*>((void *)pMem);
        }
        inline CMemPool<T, N> * freeMem() //Returns ptr to parent iff parent was full but now has space
        {
            if(IS_DEBUG) CHECK(!pParentIfOccupied, "Freeing free'd mem");
            if(IS_DEBUG) CHECK(pParentIfOccupied->nFree < 0 || pParentIfOccupied->nFree >= N, "Block count error");
            CMemPool<T, N> * pMemPoolNowWithSpace = (0==pParentIfOccupied->nFree) ? pParentIfOccupied : 0;
            pParentIfOccupied->nFree++;
            pParentIfOccupied->pErasedBlock = this;
            pParentIfOccupied = 0;
            return pMemPoolNowWithSpace;
        }
    };

    friend class CMemBlock;

    int nFree;
    CMemBlock * pNextBlock, * pErasedBlock; //Track recently erased block for fast relocation
    CMemPool<T, N> * pChildBlock;
    const int nBlockId;

    CMemBlock aBlocks[N];

    const CMemBlock * memBlockEnd() const { return aBlocks+N; }

public:
    inline CMemPool(int nId = 0) : nFree(N), pNextBlock(aBlocks), pErasedBlock(0), pChildBlock(0), /*pChildBlockWithFreeSpace(0), pParentBlock(pParentBlock),*/ nBlockId(nId) {};
    inline ~CMemPool() { delete pChildBlock; };

    int id() { return nBlockId; };

    inline T * alloc(CMemPool<T, N> ** ppLowestPoolWithFreeSpace)
    {
        if(nFree == 0)
        {
            if(!pChildBlock)
                pChildBlock = new CMemPool<T, N>(nBlockId+1);

            return pChildBlock->alloc(ppLowestPoolWithFreeSpace);
        }
        else
        {
            *ppLowestPoolWithFreeSpace = this;

            nFree--;

            if(pErasedBlock)
            {
                CMemBlock * pBlockUse = pErasedBlock;
                pErasedBlock=0;
                return pBlockUse->allocMem(this);
            }

            const CMemBlock * pEndBlock = memBlockEnd();
            while(pNextBlock->occupied())
            {
                pNextBlock++;
                if(pNextBlock == pEndBlock)
                    pNextBlock = aBlocks;
            }
            return pNextBlock->allocMem(this);
        }
    }

    inline CMemPool<T, N> * free(T * pMem)
    {
        CMemBlock * pBlock = (CMemBlock*)(void *)pMem;
        return pBlock->freeMem();
    }
};

template<typename T, int N>
class CMemPool_NoFree {
public:
    class CMemBlock
    {
        char pMem[sizeof(T)];
    public:
        inline CMemBlock() {}

        inline T * allocMem()
        {
            return static_cast<T*>((void *)pMem);
        }
    };
private:

    CMemBlock * pNextBlock;//, * pErasedBlock; //Track recently erased block for fast relocation
    CMemPool_NoFree<T, N> * pChildBlock;

    CMemBlock aBlocks[N];

    const CMemBlock * memBlockEnd() const { return aBlocks+N; }

public:
    inline CMemPool_NoFree() : pNextBlock(aBlocks), pChildBlock(0) {};
    inline ~CMemPool_NoFree() { delete pChildBlock; };

    inline T * alloc(CMemPool_NoFree<T, N> ** ppLowestPoolWithFreeSpace)
    {
        if(pNextBlock == memBlockEnd())
        {
            pNextBlock = aBlocks;

            if(!pChildBlock)
                pChildBlock = new CMemPool_NoFree<T, N>;

            return pChildBlock->alloc(ppLowestPoolWithFreeSpace);
        }
        else
        {
            *ppLowestPoolWithFreeSpace = this;

            T * mem = pNextBlock->allocMem();
            pNextBlock++;
            return mem;
        }
    }
};

template<int N>
class CDefineWordT;

template<>
class CDefineWordT<4>
{
    //typedef char FAIL_ON_64_BIT[4-(int)sizeof(void *)];
public:
    typedef unsigned int word_t;
    static inline word_t MSB() { return 0x80000000; };
};

template<>
class CDefineWordT<8>
{
    //typedef char FAIL_ON_32_BIT[ (int)sizeof(void *) - 8 ];
public:
    typedef unsigned long long int word_t;
    static inline word_t MSB() { return 0x8000000000000000LL; }
};

template<typename T, int N>
class CBitfieldMemPool : public CDefineWordT<sizeof(void *)> 
{
    static const int NUM_BITS = sizeof(word_t) * 8;
    static const int NUM_BINS = N/NUM_BITS;

    word_t aBitfield[NUM_BINS];

    //template <typename t> inline static t MSB();
    //template <typename t> inline static t
    //template <> inline static unsigned long long int MSB<unsigned long long int>()
    void print(word_t w)
    {
        for(int i=0;i<(int)sizeof(void*)*8; i++)
            std::cout << (((MSB() >> i) & w) ? '1' : '0');
        std::cout << std::endl;
    }

    inline word_t getBit(int nBit) const
    {
        int nBin = nBit / NUM_BITS;
        int nBitIdx = nBit % NUM_BITS;
        return aBitfield[nBin] & (MSB() >> nBitIdx);
    }

    inline void setBit(int nBit) 
    {
        //std::cout << "[Set bit " << nBit << ']' << std::endl;
        if(IS_DEBUG) CHECK(getBit(nBit), "Already set");
        int nBin = nBit / NUM_BITS;
        int nBitIdx = nBit % NUM_BITS;
        //print(aBitfield[nBin]);
        aBitfield[nBin] |= MSB() >> nBitIdx;
        //print(aBitfield[nBin]);
        if(IS_DEBUG) CHECK(!getBit(nBit), "Not set");
    }

    inline void clearBit(int nBit) 
    {
        //std::cout << "[Clear bit " << nBit << ']' << std::endl;
        if(IS_DEBUG) CHECK(!getBit(nBit), "Already clear");
        int nBin = nBit / NUM_BITS;
        int nBitIdx = nBit % NUM_BITS;
        //print(aBitfield[nBin]);
        aBitfield[nBin] ^= MSB() >> nBitIdx;
        //print(aBitfield[nBin]);
        if(IS_DEBUG) CHECK(getBit(nBit), "Not clear");
    }

    static inline int firstFreeBit(unsigned char c) 
    {
        if(IS_DEBUG) CHECK(c == 0, "No free bit");

        unsigned int n = c;
        for (int i=0;i<7;i++)
            if(n & (128 >> i))
                return i;

        return 7;
    }

    static inline int firstFreeBit(unsigned short s) 
    {
        if(IS_DEBUG) CHECK(s == 0, "No free bit");
        unsigned char * ac = reinterpret_cast<unsigned char *>(&s);
        if(!ac[1])
            return 8 + firstFreeBit(ac[0]);
        return firstFreeBit(ac[1]);
    }
    static inline int firstFreeBit(unsigned int n) 
    {
        if(IS_DEBUG) CHECK(n == 0, "No free bit");
        union { unsigned int n; unsigned short as[2]; } as; as.n = n;
        if(!as.as[1])
            return 16 + firstFreeBit(as.as[0]);
        return firstFreeBit(as.as[1]);
    }
    static inline int firstFreeBit(unsigned long long int l)
    {
        if(IS_DEBUG) CHECK(l == 0, "No free bit");
        unsigned int * an = reinterpret_cast<unsigned int *>(&l);
        if(!an[1])
            return 32 + firstFreeBit(an[0]);
        return firstFreeBit(an[1]);
    }
    /*static inline int firstFreeBit(unsigned int n)
    {
        if(IS_DEBUG) CHECK(n == 0, "No free bit");
        unsigned short * as = reinterpret_cast<unsigned short *>(&n);
        if(!as[1])
            return 16 + firstFreeBit(as[0]);
        return firstFreeBit(as[1]);
    }*/
    //Doesn't matter which bit these return, as long as they do.

    class CBitfieldMemBlock
    {
        char pMem[sizeof(T)];
        //static const int PAD_SIZE = (4 - sizeof(T) % 4) == 4 ? 0 : (4 - sizeof(T) % 4);
        //char padding[PAD_SIZE];
        CBitfieldMemPool<T, N> * pParent; // a pointer to the parent (allows very fast erasing)
        int idx() const 
        {
            return (this - pParent->aBlocks);///sizeof(CBitfieldMemBlock);  
        } // my idx
    public:
        inline CBitfieldMemBlock() : pParent(0) {}

        inline bool occupied() const { return pParent != 0; }
        inline T * allocMem(CBitfieldMemPool<T, N> * pParent_in)
        {
            //if(IS_DEBUG) CHECK(pParent, "Allocating alloc'd mem");
            pParent = pParent_in;
            //std::cout << "Alloc bit " << idx() << std::endl;
            pParent->clearBit(idx());
            return static_cast<T*>((void *)pMem);
        }
        inline CBitfieldMemPool<T, N> * freeMem() //Returns ptr to parent iff parent was full but now has space
        {
            if(IS_DEBUG) CHECK(!pParent, "Freeing mem that was never alloc'd");
            if(IS_DEBUG) CHECK(pParent->nFree < 0 || pParent->nFree >= N, "Block count error");
            CBitfieldMemPool<T, N> * pMemPoolNowWithSpace = (0==pParent->nFree) ? pParent : 0;
            pParent->nFree++;
            pParent->pErasedBlock = this;
            //std::cout << "Free bit " << idx() << std::endl;
            pParent->setBit(idx());
            return pMemPoolNowWithSpace;
        }
    };

    friend class CBitfieldMemBlock;

    int nFree;
    word_t * pNextBlock;
    CBitfieldMemBlock * pErasedBlock; //Track recently erased block for fast relocation
    CBitfieldMemPool<T, N> * pChildBlock;
    const int nBlockId;

    CBitfieldMemBlock aBlocks[N];

    const word_t * memBlockEnd() const { return aBitfield+NUM_BINS; }
public:
    inline CBitfieldMemPool(int nId = 0) : nFree(N), pNextBlock(aBitfield), pErasedBlock(0), pChildBlock(0), nBlockId(nId) 
    {
        /*for(unsigned int i=1; i<2000; i++)
        {
            cout << i << ' ' << firstFreeBit(i) << endl;
            // _mm_cmpestri
        }*/
        int i=NUM_BINS;
        for(word_t * pBF = aBitfield; i>0; i--, pBF++)
            *pBF = (word_t)-1;
        //setConstant<word_t>(aBitfield, -1, NUM_BINS);
    };
    inline ~CBitfieldMemPool() { delete pChildBlock; };

    int id() { return nBlockId; };

    inline T * alloc(CBitfieldMemPool<T, N> ** ppLowestPoolWithFreeSpace)
    {
        if(nFree == 0)
        {
            if(!pChildBlock)
                pChildBlock = new CBitfieldMemPool<T, N>(nBlockId+1);

            return pChildBlock->alloc(ppLowestPoolWithFreeSpace);
        }
        else
        {
            *ppLowestPoolWithFreeSpace = this;

            nFree--;

            if(pErasedBlock)
            {
                CBitfieldMemBlock * pBlockUse = pErasedBlock;
                pErasedBlock=0;
                //std::cout << "Alloc erased block\n";
                return pBlockUse->allocMem(this);
            }

            const word_t * pEndBlock = memBlockEnd();
            while(*pNextBlock == 0)
            {
                pNextBlock++;
                if(pNextBlock == pEndBlock)
                    pNextBlock = aBitfield;
            }

            int nBlock = firstFreeBit(*pNextBlock) + NUM_BITS*(pNextBlock-aBitfield);
            //std::cout << "Alloc new block " << nBlock << std::endl;
            return aBlocks[nBlock].allocMem(this);
        }
    }

    inline CBitfieldMemPool<T, N> * free(T * pMem)
    {
        CBitfieldMemBlock * pBlock = (CBitfieldMemBlock*)(void *)pMem;
        return pBlock->freeMem();
    }
};

template<typename T, int N>
class CMemPoolWrapper //A ref-counted pointer that can be passed about
{
    CMemPool<T, N> myPool, * pLowestPoolWithFreeSpace;
    int nRefs;
    DEBUG_LOG(int nCountUses);
public:
    CMemPoolWrapper() : pLowestPoolWithFreeSpace(&myPool), nRefs(1) DEBUG_LOG(, nCountUses(0))  {}

    T * alloc()
    {
        DEBUG_LOG(nCountUses++);
        return pLowestPoolWithFreeSpace->alloc(&pLowestPoolWithFreeSpace);
    }

    void free(T * p)
    {
        CMemPool<T, N> * pPoolNowWithSpace = myPool.free(p);
        if(pPoolNowWithSpace)
        {
            if(pLowestPoolWithFreeSpace->id() > pPoolNowWithSpace->id())
                pLowestPoolWithFreeSpace = pPoolNowWithSpace;
        }
    }

    void addRef()
    {
        nRefs++;
    }
    void removeRef()
    {
        nRefs--;
        if(nRefs==0)
        {
            DEBUG_LOG(logUse::log("Original", nCountUses));
            delete this;
        }
    }
};

template<typename T, int N>
class CMemPool_NoFree_Wrapper //A ref-counted pointer that can be passed about
{
    CMemPool_NoFree<T, N> myPool_NoFree_, * pLowestPool_NoFree_WithFreeSpace;
    typename CMemPool_NoFree<T, N>::CMemBlock * pErasedBlock;
    int nRefs;
    DEBUG_LOG(int nCountUses);
public:
    CMemPool_NoFree_Wrapper() : pLowestPool_NoFree_WithFreeSpace(&myPool_NoFree_), pErasedBlock(0), nRefs(1) DEBUG_LOG(, nCountUses(0)) {}

    T * alloc()
    {
        DEBUG_LOG(nCountUses++);
        if(pErasedBlock)
        {
            T * mem = pErasedBlock->allocMem();
            pErasedBlock = 0;
            return mem;
        }
        else
            return pLowestPool_NoFree_WithFreeSpace->alloc(&pLowestPool_NoFree_WithFreeSpace);
    }

    void free(T * p)
    {
        pErasedBlock = reinterpret_cast<typename CMemPool_NoFree<T,N>::CMemBlock *>(p);
    }

    void addRef()
    {
        nRefs++;
    }

    void removeRef()
    {
        nRefs--;
        if(nRefs==0)
        {
            DEBUG_LOG(logUse::log("NoFree", nCountUses));
            delete this;
        }
    }
};

template<typename T, int N>
class CBitfieldMemPoolWrapper //A ref-counted pointer that can be passed about
{
    CBitfieldMemPool<T, N> myPool, * pLowestPoolWithFreeSpace;
    int nRefs;
    DEBUG_LOG(int nCountUses);
public:
    CBitfieldMemPoolWrapper() : pLowestPoolWithFreeSpace(&myPool), nRefs(1) DEBUG_LOG(, nCountUses(0)) {}

    T * alloc()
    {
        DEBUG_LOG(nCountUses ++);
        return pLowestPoolWithFreeSpace->alloc(&pLowestPoolWithFreeSpace);
    }

    void free(T * p)
    {
        CBitfieldMemPool<T, N> * pPoolNowWithSpace = myPool.free(p);
        if(pPoolNowWithSpace)
        {
            if(pLowestPoolWithFreeSpace->id() > pPoolNowWithSpace->id())
                pLowestPoolWithFreeSpace = pPoolNowWithSpace;
        }
    }

    void addRef()
    {
        nRefs++;
    }

    void removeRef()
    {
        nRefs--;
        if(nRefs==0)
        {
            DEBUG_LOG(logUse::log("Bitfield", nCountUses));
            delete this;
        }
    }
};

template<typename T>
class CAllocatorBase
{
public:
    typedef T value_type;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;

    inline pointer address(reference r) { return &r; }
    inline const_pointer address(const_reference r) { return &r; }

    inline size_type max_size() const 
    {
        //return std::numeric_limits<size_type>::max() / sizeof(T);
        //typedef char[((size_type)(-1)>0)-1] size_type_must_be_unsigned;
        return (size_type)(-1) / sizeof(T);
    }

    inline void construct(pointer p, const T& t) { new(p) T(t); }
    inline void destroy(pointer p) { p->~T(); }

    /*inline bool operator==(CIndividualPoolAllocator<T, N> const&) {
        return true;
    }
    inline bool operator!=(CIndividualPoolAllocator<T, N> const& a)
    {
        return !operator==(a);
    }*/
};

template<typename T, int N>
class CIndividualPoolAllocator : public CAllocatorBase<T> 
{
    CMemPoolWrapper<T, N> * pMemPool;
public :
    //    convert an allocator<T> to allocator<U>
    template<typename U>
    struct rebind {
        typedef CIndividualPoolAllocator<U, N> other;
    };

    CIndividualPoolAllocator() : pMemPool(0 /*new CMemPoolWrapper<T, N>*/)
    {}
    ~CIndividualPoolAllocator() { if(pMemPool) pMemPool->removeRef();  }
    CIndividualPoolAllocator(CIndividualPoolAllocator<T, N> const& otherPool) : pMemPool(otherPool.pMemPool)
    {
        //std::cout << "Copying pool allocator\n";
        pMemPool->addRef();
    }

#ifdef __GNUC__
    template<typename U>
    CIndividualPoolAllocator(CIndividualPoolAllocator<U, N> const&) : pMemPool(0)
    {
        //This is called all the time to get to the c'tor
        //std::cout << sizeof(U) << ' ' << sizeof(T) << " UT copy\n";
    }
#else
    template<typename U>
    CIndividualPoolAllocator(CIndividualPoolAllocator<U, N> const & otherAlloc) : pMemPool(new CMemPoolWrapper<T, N>)
    {
        //std::cout << "Creating different type of pool allocator\n"; //Called once per set creation in windows. should check the derived one is actually used more than once...
    }
#endif

    inline typename CAllocatorBase<T>::pointer allocate(typename CAllocatorBase<T>::size_type cnt,
       typename std::allocator<void>::const_pointer = 0) {
        //Could return reinterpret_cast<pointer>(::operator new(cnt * sizeof (T)));
        if(IS_DEBUG) CHECK(cnt!=1, "This allocator only supports single blocks!");
        if(!pMemPool)
            pMemPool = new CMemPoolWrapper<T, N>;
        if(IS_DEBUG) CHECK(!pMemPool, "Mem pool is null");
        return pMemPool->alloc();
    }

    inline void deallocate(typename CAllocatorBase<T>::pointer p, typename CAllocatorBase<T>::size_type) {
        //::operator delete(p);
        if(IS_DEBUG) CHECK(!pMemPool, "Mem pool is null")
        pMemPool->free(p);
    }
};   

template<typename T, int N>
class CIndividualPool_NoFree_Allocator : public CAllocatorBase<T> 
{
    CMemPool_NoFree_Wrapper<T, N> * pMemPool_NoFree_;
public :
    //    convert an allocator<T> to allocator<U>
    template<typename U>
    struct rebind {
        typedef CIndividualPool_NoFree_Allocator<U, N> other;
    };

    CIndividualPool_NoFree_Allocator() : pMemPool_NoFree_(0 /*new CMemPool_NoFree_Wrapper<T, N>*/)
    {}
    ~CIndividualPool_NoFree_Allocator() { if(pMemPool_NoFree_) pMemPool_NoFree_->removeRef();  }
    CIndividualPool_NoFree_Allocator(CIndividualPool_NoFree_Allocator<T, N> const& otherPool_NoFree_) : pMemPool_NoFree_(otherPool_NoFree_.pMemPool_NoFree_)
    {
        //std::cout << "Copying Pool_NoFree allocator\n";
        if(pMemPool_NoFree_)
            pMemPool_NoFree_->addRef();
    }

#ifdef __GNUC__
    template<typename U>
    CIndividualPool_NoFree_Allocator(CIndividualPool_NoFree_Allocator<U, N> const&) : pMemPool_NoFree_(0)
    {
        //This is called all the time to get to the c'tor
        //std::cout << sizeof(U) << ' ' << sizeof(T) << " UT copy\n";
    }
#else
    template<typename U>
    CIndividualPool_NoFree_Allocator(CIndividualPool_NoFree_Allocator<U, N> const & otherAlloc) : pMemPool_NoFree_(new CMemPool_NoFree_Wrapper<T, N>)
    {
        //std::cout << "Creating different type of Pool_NoFree allocator\n"; //Called once per set creation in windows. should check the derived one is actually used more than once...
    }
#endif

    inline typename CAllocatorBase<T>::pointer allocate(typename CAllocatorBase<T>::size_type cnt,
       typename std::allocator<void>::const_pointer = 0) {
        //Could return reinterpret_cast<pointer>(::operator new(cnt * sizeof (T)));
        if(IS_DEBUG) CHECK(cnt!=1, "This allocator only supports single blocks!")
        if(!pMemPool_NoFree_)
            pMemPool_NoFree_ = new CMemPool_NoFree_Wrapper<T, N>;
//        if(IS_DEBUG) CHECK(!pMemPool_NoFree_, "Mem Pool_NoFree is null")
        return pMemPool_NoFree_->alloc();
    }

    inline void deallocate(typename CAllocatorBase<T>::pointer p, typename CAllocatorBase<T>::size_type) {
        //::operator delete(p);
        if(IS_DEBUG) CHECK(!pMemPool_NoFree_, "Mem Pool_NoFree is null")
        pMemPool_NoFree_->free(p);
    }
};  

//#define CBitfieldMemPoolWrapper CMemPoolWrapper
template<typename T, int N>
class CBitfieldIndividualPoolAllocator : public CAllocatorBase<T>//, protected CDefineWordT< sizeof(void *) >
{
    typedef char THIRTYTWO_OR_SIXTYFOUR_MUST_DIVIDE_N[-(int)(N % (sizeof(void *)*8))]; //Should use large N, like 1024, anyway
    CBitfieldMemPoolWrapper<T, N> * pMemPool;
public :
    // convert an allocator<T> to allocator<U>
    template<typename U>
    struct rebind {
        typedef CBitfieldIndividualPoolAllocator<U, N> other;
    };

    CBitfieldIndividualPoolAllocator() : pMemPool(0 /*new CBitfieldMemPoolWrapper<T, N>*/)
    {}
    ~CBitfieldIndividualPoolAllocator() { if(pMemPool) pMemPool->removeRef();  }
    CBitfieldIndividualPoolAllocator(CBitfieldIndividualPoolAllocator<T, N> const& otherPool) : pMemPool(otherPool.pMemPool)
    {
        //std::cout << "Copying pool allocator\n";
        if(pMemPool)
            pMemPool->addRef();
    }

#ifdef __GNUC__
    template<typename U>
    CBitfieldIndividualPoolAllocator(CBitfieldIndividualPoolAllocator<U, N> const&) : pMemPool(0)
    {
        //This is called all the time to get to the c'tor
        //std::cout << sizeof(U) << ' ' << sizeof(T) << " UT copy\n";
    }
#else
    template<typename U>
    CBitfieldIndividualPoolAllocator(CBitfieldIndividualPoolAllocator<U, N> const & otherAlloc) : pMemPool(new CBitfieldMemPoolWrapper<T, N>)
    {
        //std::cout << "Creating different type of pool allocator\n"; //Called once per set creation in windows. should check the derived one is actually used more than once...
    }
#endif

    inline typename CAllocatorBase<T>::pointer allocate(typename CAllocatorBase<T>::size_type cnt,
       typename std::allocator<void>::const_pointer = 0) {
        //Could return reinterpret_cast<pointer>(::operator new(cnt * sizeof (T)));
        if(IS_DEBUG) CHECK(cnt!=1, "This allocator only supports single blocks!")
        //if(IS_DEBUG) CHECK(!pMemPool, "Mem pool is null")
        if(!pMemPool)
            pMemPool = new CBitfieldMemPoolWrapper<T, N>;

        return pMemPool->alloc();
    }

    inline void deallocate(typename CAllocatorBase<T>::pointer p, typename CAllocatorBase<T>::size_type) {
        //::operator delete(p);
        if(IS_DEBUG) CHECK(!pMemPool, "Mem pool is null")
        pMemPool->free(p);
    }
};

//TB: This isn't working in Windows (17-6-12)
//template<typename T, typename _Compare = std::less<T>,
//typename _Alloc = CBitfieldIndividualPoolAllocator< T, 128 > >

template<typename T, typename _Compare = std::less<T> >
class set2 : public std::set<T, _Compare>
{
    typedef std::set<T, _Compare> TSet;

//    _Alloc myAllocator;
    _Compare myComp;
public:
    inline bool exists(const T & el) const
    {
        return TSet::find(el) != TSet::end();
    }

    void pp(const char * seperator = ", ") const
    {
        for(typename TSet::const_iterator pEl = TSet::begin(); pEl != TSet::end(); pEl++)
        {
            if(pEl != TSet::begin())
                std::cout << seperator;
            std::cout << *pEl;
        }
        std::cout << std::endl;
    }

    T pop()
    {
        if(IS_DEBUG) CHECK(TSet::size() == 0, "pop: Nothing in map");
        typename TSet::iterator iter=TSet::begin();
        T val = *iter; //Not a valid ref anymore once popped
        TSet::erase(iter);
        return val;
    }

    const T & top() const
    {
        if(IS_DEBUG) CHECK(TSet::size() == 0, "top: Nothing in set");
        return *TSet::begin();
    }
    const T & back() const
    {
        if(IS_DEBUG) CHECK(TSet::size() == 0, "top: Nothing in set");
        typename TSet::const_iterator end = TSet::end();
        end--;
        return *end;
    }

    //set2() : TSet(myComp, myAllocator) {}
};

template<typename Key, typename Val, typename _Compare = std::less<Key> >
// TB: This isn't working in Windows        typename _Alloc = std::allocator<std::pair< Key, Val> > >
class map2 : public std::map<Key, Val, _Compare/*, _Alloc*/>
{
    //_Alloc myAllocator;
    //_Compare myComp;

public:
    typedef std::map<Key, Val, _Compare/*, _Alloc*/> TMap;

    inline bool exists(const Key & el) const
    {
        return TMap::find(el) != TMap::end();
    }

    inline const Val & ifExists(const Key & el, const Val & otherwise) const
    {
        typename TMap::const_iterator iter =  TMap::find(el);
        return iter != TMap::end() ? iter->second : otherwise;
    }

    /*inline const Val & ifExists(const Key & el, const Val otherwise) const
    {
        typename TMap::const_iterator iter =  TMap::find(el);
        return iter != TMap::end() ? iter->second : otherwise;
    }*/

    void pp(const char * seperator = ", ") const
    {
        for(typename TMap::const_iterator pEl = TMap::begin(); pEl != TMap::end(); pEl++)
        {
            if(pEl != TMap::begin())
                std::cout << seperator;
            std::cout << pEl->first << "->" << pEl->second;
        }
        std::cout << std::endl;
    }
    const Val & top() const
    {
        if(IS_DEBUG) CHECK(TMap::size() == 0, "top: Nothing in map");
        return TMap::begin()->second;
    }

    /*Val & top()
    {
        if(IS_DEBUG) CHECK(TMap::size() == 0, "top: Nothing in map");
        return TMap::begin()->second;
    }*/

    const Key & topKey() const
    {
        if(IS_DEBUG) CHECK(TMap::size() == 0, "topKey: Nothing in map");
        return TMap::begin()->first;
    }

    typename TMap::const_iterator backIter() const
    {
        if(IS_DEBUG) CHECK(TMap::size() == 0, "back: Nothing in map");
        typename TMap::const_iterator pBack = TMap::end();
        pBack--;
        return pBack;
    }

    const Val & back() const
    {
        if(IS_DEBUG) CHECK(TMap::size() == 0, "back: Nothing in map");
        return backIter()->second;
    }

    /*Val & back()
    {
        if(IS_DEBUG) CHECK(TMap::size() == 0, "back: Nothing in map");
        return backIter()->second;
    }*/

    Val & pop()
    {
        if(IS_DEBUG) CHECK(TMap::size() == 0, "pop: Nothing in map");
        typename TMap::iterator iter=TMap::begin();
        Val & val = iter->second;
        TMap::erase(iter);
        return val;
    }

    const Key & backKey() const
    {
        if(IS_DEBUG) CHECK(TMap::size() == 0, "backKey: Nothing in map");
        return backIter()->first;
    }

    const Val & operator[](const Key & key) const
    {
        if(IS_DEBUG) CHECK(!exists(key), "Key not found in map");
        return TMap::find(key)->second;
    }

    Val & operator[](const Key & key)
    {
        if(IS_DEBUG) CHECK(!exists(key), "Key not found in map");
        return TMap::operator[](key);
    }

    //Explicit initiation to prevent accidental initiation
    Val & init(const Key & key, const Val & val)
    {
        if(IS_DEBUG) CHECK(exists(key), "Key already exists. Use initOrSet?");
        Val & newVal = TMap::operator[](key);
        newVal = val;
        return newVal;
    }

    //Explicit initiation to prevent accidental initiation
    //Overwrites value if already exists
    Val & initOrSet(const Key & key, const Val & val)
    {
        Val & newVal = TMap::operator[](key);
        newVal = val;
        return newVal;
    }

    //Explicit initiation to prevent accidental initiation
    Val & initOrGet(const Key & key)
    {
        Val & newVal = TMap::operator[](key);
        return newVal;
    }

    //Explicit initiation to prevent accidental initiation
    //Returns old value if already exists
    Val initOrGet(const Key & key, const Val & val)
    {
        typename TMap::iterator loc = TMap::find(key);
        if(loc == TMap::end())
        {
            TMap::operator[](key) = val;
            return val;
        }
        else
        {
            return loc->second;
        }
    }

    //Explicit initiation to prevent accidental initiation
    Val & init(const Key & key)
    {
        if(IS_DEBUG) CHECK(exists(key), "Key already exists");
        Val & newVal = TMap::operator[](key);
        return newVal;
    }

    //map2() : TMap(myComp, myAllocator) {}
};

//Fast insertion and relocation but memory not freed until goes out of scope
template<typename Key, typename Val, typename _Compare = std::less<Key>, int N=256 >
class map2_NF : public map2<Key, Val, _Compare /*, CIndividualPool_NoFree_Allocator<std::pair< Key, Val>, N >*/ >
{
    typedef typename map2<Key, Val /*, _Compare, CIndividualPool_NoFree_Allocator<std::pair< Key, Val>, N >*/ >::TMap TMap;
public:
    typename TMap::iterator move(const Key & oldKey, const Key & newKey)
    {
        typename TMap::iterator pOld = TMap::find(oldKey);

        if(oldKey == newKey)
            return pOld;

        if(IS_DEBUG) CHECK(pOld == TMap::end(), "move: Old key doesn't exist")

        const Val & val = pOld->second;

        TMap::erase(pOld);

        const std::pair<typename TMap::iterator, bool> & insRes = insert(std::pair< Key, Val >(newKey, val) );

        if(IS_DEBUG) CHECK(!insRes.second, "Error re-inserting--is the comparison working?")

        return insRes.first;
    }
};

//Fast insertion and relocation but memory not freed until goes out of scope
template<typename T, typename _Compare = std::less<T>, int N=256 >
class set2_NF : public set2<T, _Compare /*CIndividualPool_NoFree_Allocator< T, N >*/ >
{
public:

};

#endif /* SET2_H_ */
