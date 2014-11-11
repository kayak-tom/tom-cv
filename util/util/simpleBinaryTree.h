/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/* Code by Tom Botterill, documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * simpleBinaryTree.h
 *
 *  Created on: 2/11/2009
 *      Author: tom
 *
 *      Intended to be a much faster alternative to std::set/std::map by having its own memory pool, and by avoiding all the balancing ops
 *      Not guaranteed O(log(n)) but should be close **if keys arrive in random order**
 *      Poor performance if keys arrive sorted
 *
 *      TODO: Tree balance by iterating and inserting into a new tree?? Maybe access elements in array at random and insert.
 */

#ifndef SIMPLEBINARYTREE_H_
#define SIMPLEBINARYTREE_H_

#include "exception.h"
//#include <vector>
#define VERBOSE(x)

template<class T, class CComparison = std::less<T> >
class CSimpleSet;

namespace _hidden
{
    template<class T, class CComparison>
    class const_iterator;

    template<class T, class CComparison>
    class const_iterator //Todo: ordered iter too
    {
        const /*CSimpleSet<T, CComparison>::CStoredSimpleNode*/ void * pStoredNode;
        const CSimpleSet<T, CComparison> * pSet;

        const_iterator() {};
    public:
        const_iterator(const /*CSimpleSet<T, CComparison>::CStoredSimpleNode*/ void * pStoredNode, const CSimpleSet<T, CComparison> * pSet) : pStoredNode(pStoredNode), pSet(pSet) {}
        const T & operator*() const;

        //const CSimpleSet &
        void operator++();

        bool operator!=(const const_iterator & pOther) const { return pStoredNode != pOther.pStoredNode; }
        bool operator==(const const_iterator & pOther) const { return pStoredNode == pOther.pStoredNode; }
    };

}

template<class T, class CComparison>
class CSimpleSet
{
    typedef CSimpleSet<T, CComparison> TSimpleSet;
    CComparison compOp;

    static const int UNSET = -1;

    class CStoredSimpleNode;
    friend class _hidden::const_iterator<T, CComparison>;

public:
    typedef _hidden::const_iterator<T, CComparison> const_iterator;

private:
    class CSimpleNode
    {
        friend class CSimpleSet;
        CComparison compOp;
        int nOffsetLeft, nOffsetRight, nOffsetParent;
        const T nodeVal;
    public:
        CSimpleNode(const T & val, int nOffsetParent) : nOffsetLeft(CSimpleSet::UNSET), nOffsetRight(CSimpleSet::UNSET), nOffsetParent(nOffsetParent), nodeVal(val) {}
        const T & val() const { return nodeVal; }
    };

    class CStoredSimpleNode
    {
        char data[sizeof(CSimpleNode)];
        bool bOccupied;
    public:
        CStoredSimpleNode() : bOccupied(false) {}
        bool occupied() const { return bOccupied; }
        void * getMemory() {
            if(IS_DEBUG) CHECK(bOccupied, "Mem in use");
            bOccupied = true;
            return data;
        }
        CSimpleNode * getNode() {
            if(IS_DEBUG) CHECK(!bOccupied, "Mem not initialised");
            return (CSimpleNode *)data;
        }
        const CSimpleNode * getNode() const {
            if(IS_DEBUG) CHECK(!bOccupied, "Mem not initialised");
            return (const CSimpleNode *)data;
        }
        void erase()
        {
            if(IS_DEBUG) CHECK(!bOccupied, "Erasing unoccupied node");
            bOccupied = false;
        }
    };

    enum eDirection { eFirst, eLast };

    template<eDirection dir>
    bool insert_int(CStoredSimpleNode * pNodeStored, const T & val)
    {
        CSimpleNode * pNode = pNodeStored->getNode();
        int * pnOffset = dir==eFirst ? &(pNode->nOffsetLeft) : &(pNode->nOffsetRight);
        if(*pnOffset != UNSET)
            return insertAtNode(elements + *pnOffset, val);

        //New element
        int nParent = pNodeStored - elements;

        if(nCount >= nAllocated) //Increase size
        {
            // This will invalidate ptr so RESET them
            increaseSize(nAllocated * 4);
            pNodeStored = elements + nParent;
            pNode = pNodeStored->getNode();
            pnOffset = dir==eFirst ? &(pNode->nOffsetLeft) : &(pNode->nOffsetRight);
        }

        int nNewLoc = insertAtMem(val, nParent);
        *pnOffset = nNewLoc;
        return true;
    }

    int insertAtMem(const T & val, int nParent)
    {
        while(pFindNext->occupied())
        {
            pFindNext++;
            if(pFindNext >= ptr_end())
                pFindNext = elements;
        }
        int nOffset = pFindNext - elements;
        VERBOSE(std::cout << "inserted " << val << " at " << nOffset << std::endl;)
        new(pFindNext->getMemory()) CSimpleNode(val, nParent);
        nCount++;

        pFindNext++;
        if(pFindNext >= ptr_end())
            pFindNext = elements;

        return nOffset;
    }

    inline CStoredSimpleNode * ptr_end() const { return elements + nAllocated; }

    //True if new element
    bool insertAtNode(CStoredSimpleNode * pNodeStored, const T & val)
    {
        CSimpleNode * pNode = pNodeStored->getNode();
        if (compOp(val, pNode->nodeVal))
        {
            return insert_int<eFirst>(pNodeStored, val);
        }
        else if (compOp(pNode->nodeVal, val))
        {
            return insert_int<eLast>(pNodeStored, val);
        }
        //Insert failed because already exists (here)
        return false;
    }

    int nCount, nAllocated, nRootOffset;
    CStoredSimpleNode * elements, * pFindNext; //Root node must be 1st in elements

    void increaseSize(int nNewSize)
    {
        VERBOSE(std::cout <<"Increasing size to " << nNewSize << " (" << sizeof(CStoredSimpleNode) << " bytes each)" << std::endl;)
        CStoredSimpleNode * aNewStore = new CStoredSimpleNode[nNewSize];
        //Copy:
        for(int i=0;i<nAllocated;i++)
            aNewStore[i] = elements[i]; // check is binary-copying

        delete [] elements; elements = aNewStore; aNewStore = 0; pFindNext=elements + nAllocated; nAllocated = nNewSize;
    }
    CStoredSimpleNode * findNode(CStoredSimpleNode * pThisNodeStored, const T & t) const
    {
        CSimpleNode * pThisNode = pThisNodeStored->getNode();

        if (compOp(t, pThisNode->nodeVal))
        {
            if(pThisNode->nOffsetLeft == UNSET)
                return 0;

            return findNode(elements+pThisNode->nOffsetLeft, t);
        }
        else if (compOp(pThisNode->nodeVal, t))
        {
            if(pThisNode->nOffsetRight == UNSET)
                return 0;

            return findNode(elements+pThisNode->nOffsetRight, t);
        }
        return pThisNodeStored;
    }
    CStoredSimpleNode * findNode(const T & t) const
    {
        if(nCount==0)
            return 0;

        return findNode(elements+nRootOffset, t);
    }

    template<eDirection dir>
    const CStoredSimpleNode * endNode(const CStoredSimpleNode * pThisNodeStored) const
    {
        const CSimpleNode * pThisNode = pThisNodeStored->getNode();
        int nNextNode = dir==eFirst ? pThisNode->nOffsetLeft : pThisNode->nOffsetRight;
        if(nNextNode == UNSET)
            return pThisNodeStored;
        else
            return endNode<dir>(elements + nNextNode);
    }
    template<eDirection dir>
    CStoredSimpleNode * endNode(CStoredSimpleNode * pThisNodeStored)
    {
        const CSimpleNode * pThisNode = pThisNodeStored->getNode();
        int nNextNode = dir==eFirst ? pThisNode->nOffsetLeft : pThisNode->nOffsetRight;
        if(nNextNode == UNSET)
            return pThisNodeStored;
        else
            return endNode<dir>(elements + nNextNode);
    }

    const CStoredSimpleNode * firstNode() const
    {
        if(IS_DEBUG) CHECK(nCount <= 0, "firstNode: No nodes");
        return endNode<eFirst>(elements + nRootOffset);
    }

    const CStoredSimpleNode * lastNode() const
    {
        if(IS_DEBUG) CHECK(nCount <= 0, "lastNode: No nodes");
        return endNode<eLast>(elements+nRootOffset);
    }

    const_iterator getEndIterator(eDirection dir) const
    {
        //Descend left until 0
        if(nCount == 0)
            return end();

        const CStoredSimpleNode * pFirst = 0;

        if (dir==eFirst)
            pFirst = endNode<eFirst>(elements + nRootOffset);
        else
            pFirst = endNode<eLast>(elements + nRootOffset);

        if(IS_DEBUG) CHECK(!pFirst, "begin(): Failed to find first node");

        return const_iterator(pFirst, this);
    }

    inline const CStoredSimpleNode * nextNodeUp(const CStoredSimpleNode * pPrevNodeStored) const
    {
        if(IS_DEBUG) CHECK(!pPrevNodeStored, "nextNodeUp: Uninit param");
        const CSimpleNode * pLastNodeOnPath = pPrevNodeStored->getNode();

        if(pLastNodeOnPath->nOffsetParent == UNSET)
            return 0;

        const int nLastNodeIdx = pPrevNodeStored - elements;

        const CStoredSimpleNode * pThisNodeStored = elements + pLastNodeOnPath->nOffsetParent;
        const CSimpleNode * pThisNode = pThisNodeStored->getNode();

        //Is this node a successor to pPrevNodeStored?
        if(pThisNode->nOffsetLeft == nLastNodeIdx)
            return pThisNodeStored;
        else
            return nextNodeUp(pThisNodeStored);
    }

    //To traverse in sorted order:
    inline const CStoredSimpleNode * nextNode(const CStoredSimpleNode * pThisNodeStored) const
    {
        if(IS_DEBUG) CHECK(!pThisNodeStored, "nextNode: Uninit param");
        //If we have a right node the next will be somewhere down there
        const CSimpleNode * pThisNode = pThisNodeStored->getNode();
        if(pThisNode->nOffsetRight != UNSET)
            return endNode<eFirst>(elements + pThisNode->nOffsetRight);

        //At the bottom of the tree. Go up until I'm a left node or we find a right node
        return nextNodeUp(pThisNodeStored);
    }


public:
    CSimpleSet() : nCount(0), nAllocated(0), nRootOffset(UNSET), elements(0), pFindNext(0)
    {
        increaseSize(64);
    }

    CSimpleSet(int nSizeEstimate) : nCount(0), nAllocated(0), nRootOffset(UNSET), elements(0), pFindNext(0)
    {
        increaseSize(nSizeEstimate);
    }

    ~CSimpleSet()
    {
        delete [] elements;
    }

    bool exists(const T & t) const { return findNode(t) != 0; }

    //return true if new element
    bool insert(const T & t)
    {
        bool bRes = false;
        const int nOldCount = nCount;
        if(nCount>0)
            bRes = insertAtNode(elements + nRootOffset, t);
        else
        {
            nRootOffset = insertAtMem(t, UNSET);
            bRes = true;
        }

        if(IS_DEBUG) CHECK(nOldCount+(bRes?1:0) != nCount, "insert failed");
        if(IS_DEBUG) CHECK(!exists(t), "insert failed");

        return bRes;
    }

    bool erase(const T & t)
    {
        //Make sure nRootOffset still correct.
        const int nOldCount = nCount;

        CStoredSimpleNode * pStoredNode = findNode(t);
        if(pStoredNode)
        {
            int nStoredNodeIdx = pStoredNode - elements;

            VERBOSE(std::cout << "Erasing  " << t << " at " << nStoredNodeIdx << std::endl;)
            //if parent is bigger insert L into R
            //else R into L
            CSimpleNode * pNode = pStoredNode->getNode();
            int * pnOldOffsetLocation = 0, nNewParentNodeIdx = pNode->nOffsetParent;
            if(pNode->nOffsetParent != UNSET)
            {
                CSimpleNode * pParentNode = elements[pNode->nOffsetParent].getNode();

                if(nStoredNodeIdx == pParentNode->nOffsetLeft)
                    pnOldOffsetLocation = &(pParentNode->nOffsetLeft);
                else
                    pnOldOffsetLocation = &(pParentNode->nOffsetRight);
            }
            else
                pnOldOffsetLocation = &nRootOffset; //Erasing root node

            //Simple case--join one descendent subtree to parent
            if(pNode->nOffsetLeft == UNSET && pNode->nOffsetRight == UNSET)
            {
                VERBOSE(std::cout << "Leaf\n";)
                *pnOldOffsetLocation = UNSET;
            }
            else if(pNode->nOffsetLeft == UNSET)
            {
                VERBOSE(std::cout << "hasRightChildren\n";)
                *pnOldOffsetLocation = pNode->nOffsetRight;
                CSimpleNode * pRightChildNode = elements[pNode->nOffsetRight].getNode();
                pRightChildNode->nOffsetParent = nNewParentNodeIdx;
            }
            else if(pNode->nOffsetRight == UNSET)
            {
                VERBOSE(std::cout << "hasLeftChildren\n";)
                *pnOldOffsetLocation = pNode->nOffsetLeft;
                CSimpleNode * pLeftChildNode = elements[pNode->nOffsetLeft].getNode();
                pLeftChildNode->nOffsetParent = nNewParentNodeIdx;
            }
            else
            {
                //More complex case: 2 descendents
                VERBOSE(std::cout << "has2Children\n";)

                //Attach L to Lmost of R descendents

                CStoredSimpleNode * pInsertPosStored = endNode<eFirst>(elements + pNode->nOffsetRight);
                CSimpleNode * pInsertPos = pInsertPosStored->getNode();
                if(IS_DEBUG) CHECK(!pInsertPosStored || !pInsertPos || pInsertPos->nOffsetLeft != UNSET, "erase: bad insert pos returned");
                pInsertPos->nOffsetLeft = pNode->nOffsetLeft;
                *pnOldOffsetLocation = pNode->nOffsetRight;

                //Update parents
                CSimpleNode * pLeftChildNode = elements[pNode->nOffsetLeft].getNode();
                pLeftChildNode->nOffsetParent = pInsertPosStored - elements;
                CSimpleNode * pRightChildNode = elements[pNode->nOffsetRight].getNode();
                pRightChildNode->nOffsetParent = nNewParentNodeIdx;
            }

            //Erase this node
            pNode->nodeVal.~T(); //should probably call destructor:
            pStoredNode->erase();
            nCount--;

             if(IS_DEBUG) CHECK(nCount != nOldCount -1, "erase: Erase failed");
        }
        else
        {
            if(IS_DEBUG) CHECK(nCount != nOldCount, "erase: Erase failed");
        }

        return (bool)pStoredNode;
    }

    int size() const { return nCount; }

    const_iterator end() const
    {
        return const_iterator(0, this);
    }

    const_iterator begin() const
    {
        return getEndIterator(eFirst);
    }

    const_iterator last() const
    {
        return getEndIterator(eLast);
    }
};

/*template<class T, class CComparison>
void CSimpleSet::const_iterator::operator++(int)
{
    if(IS_DEBUG) CHECK(!pStoredNode, "Incrementing 'end' iterator");
    pStoredNode = pSet->nextNode(pStoredNode);
}*/
namespace _hidden
{

template<class T, class CComparison>
const T & const_iterator<T, CComparison>::operator*() const
{
    if(IS_DEBUG) CHECK(!pStoredNode, "Dereferencing 'end' iterator");
    return reinterpret_cast<const typename CSimpleSet<T, CComparison>::CStoredSimpleNode *>(pStoredNode)->getNode()->val();
}

template<class T, class CComparison>
void const_iterator<T, CComparison>::operator++()
{
    if(IS_DEBUG) CHECK(!pStoredNode, "Incrementing 'end' iterator");
    pStoredNode = pSet->nextNode(reinterpret_cast<const typename CSimpleSet<T, CComparison>::CStoredSimpleNode *>(pStoredNode));
}

}

/*template<eDirection dir>
SimpleSet<T, CComparison>::const_iterator TSimpleSet::getEndIterator() const
{
    //Descend left until 0
    if(nCount == 0)
        return end();

    const CStoredSimpleNode * pFirst = endNode<dir>(elements + nRootOffset);

    if(IS_DEBUG) CHECK(!pFirst, "begin(): Failed to find first node");

    return const_iterator(pFirst, this);
}*/


#endif /* SIMPLEBINARYTREE_H_ */
