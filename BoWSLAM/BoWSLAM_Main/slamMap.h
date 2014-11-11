/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * slamMap.h
 *
 *  Created on: 2/11/2009
 *      Author: tom
 */

#ifndef SLAMMAP_H_
#define SLAMMAP_H_

#include "geom/geom.h"
#include <boost/noncopyable.hpp>

#include "util/set2.h"
#include "util/dynArray.h"
#include "scale/EdgeScaleObserver.h"
#include "fullRelPose.h"
#include "time/SpeedTest.h"
#include "util/fastSet.h"
#include "scaleOptimiser.h"

#include <fstream>

namespace NSLAMMap {

    //A di-edge represents the scale change between 2 nodes

    class CIds {
        int nId1, nId2;
    public:

        CIds(int nId1, int nId2) : nId1(nId1), nId2(nId2) {
            if(IS_DEBUG) CHECK(nId2 <= nId1, "Bad edge ids");
            if(IS_DEBUG) CHECK(nId2 < 0 || nId1 < 0, "Bad edge ids");
        }

        CIds() : nId1(-1), nId2(-1) {
        }

        bool init() const {
            return nId1 >= 0;
        }

        int id1() const {
            return nId1;
        }

        int id2() const {
            return nId2;
        }

        bool operator<(const CIds & ids2) const {
            return id1() < ids2.id1() || (id1() == ids2.id1() && id2() < ids2.id2());
        }

        bool operator==(const CIds & ids2) const {
            return (id1() == ids2.id1() && id2() == ids2.id2());
        }

        bool operator!=(const CIds & ids2) const {
            return !(*this == ids2);
        }
    };

    class CNode_relPos;

    class CEdge_relScale {
        //const int id1, id2, id3; //joins id1,id2 to id2,id3. Belongs to id1,id2. 
        CNode_relPos * pNode1;
        CNode_relPos * pNode2;

        CRelScale relScale;
        double dOptimisedRelScale;

        static inline const double UNINIT_RELSCALE() {
            return 1e+10;
        }
    public:
        void setOptimisedScale(const double optimisedD_in, const bool bVerbose);

        inline void resetOptimisedScale() {
            dOptimisedRelScale = UNINIT_RELSCALE();
        }

        CEdge_relScale(CNode_relPos * pNode1, CNode_relPos * pNode2, const CRelScale relScale);

        CNode_relPos * firstNode() {
            return pNode1;
        }

        CNode_relPos * secondNode() {
            return pNode2;
        }

        const CNode_relPos * firstNode() const {
            return pNode1;
        }

        const CNode_relPos * secondNode() const {
            return pNode2;
        }

        inline const CNode_relPos * otherNode(const CNode_relPos * pThis) const {
            CHECK(pThis != pNode1 && pThis != pNode2, "otherNode: pThis isn't one of the nodes");
            return pThis != pNode1 ? pNode1 : pNode2;
        }

        inline CNode_relPos * otherNode(const CNode_relPos * pThis) {
            CHECK(pThis != pNode1 && pThis != pNode2, "otherNode: pThis isn't one of the nodes");
            return pThis != pNode1 ? pNode1 : pNode2;
        }

        bool linkedTo(int nId) const;

        //double getScaleChange() const { std::cout << "Todo, shouldn't get scale change from here, should apply scale updates to relative poses"; return 1 /*dRelScale*/; }

        inline const CRelScale getRelScale() const {
            if (dOptimisedRelScale == UNINIT_RELSCALE())
                return relScale;
            else {
                REPEAT(1000,
                if (fabs(1 - dOptimisedRelScale) > relScale.get_d()*0.25)
                        std::cout << "Returning rel scale " << dOptimisedRelScale << " rather than " << relScale.get_d() << std::endl);
                return CRelScale(dOptimisedRelScale, relScale.get_g_sq());
            }
        }

        void pp(const CNode_relPos * pLast) const;

        class CSortByBadness {
        public:

            bool operator()(const CEdge_relScale * pEdge1, const CEdge_relScale * pEdge2) {
                const double g1 = pEdge1->getRelScale().get_g_sq();
                const double g2 = pEdge2->getRelScale().get_g_sq();
                if(IS_DEBUG) CHECK(g1 < 0 || g2 < 0, "Should we have uninit rel scales here???")
                return g1 < g2 || (g1 == g2 && pEdge1 < pEdge2);
            }
        };
    };

    class CCycle : public CDynArray<std::pair<double, CEdge_relScale *> > {
    };

    //A node represents the relative position of id1, id2

    class CNode_relPos : boost::noncopyable {
        CDynArray<CEdge_relScale *> vEdges; //This is where edges live...
        CEdge_relScale * pPrevEdge; //Locations with no previous edge are positioned at origin
        int nShortestPath; //For optimisation, not for Dijkstra SPT
        CEdge_relScale * pPrevEdgePosSource; //The edge that we adjust.
    public:
        const int id1, id2; //Doesn't actually need these...

        int shortestPathLength() const {
            return nShortestPath;
        }

        CEdge_relScale * prevEdgePosSource() const {
            return pPrevEdgePosSource;
        }

        CEdge_relScale * parent() const {
            return pPrevEdge;
        }

        void setShortestPathLengthAndPosSource(int nShortestPath_in, CEdge_relScale * pPrevEdgePosSource_in) {
            nShortestPath = nShortestPath_in;
            pPrevEdgePosSource = pPrevEdgePosSource_in;
            if(IS_DEBUG) CHECK(nShortestPath < MAX_INT && !pPrevEdgePosSource, "No pPrevEdgePosSource set")
        }

        inline const int otherId(const int id) const {
            if(IS_DEBUG) CHECK(id != id1 && id != id2, "otherId: Id given is neither of the ids");
            return id != id1 ? id1 : id2;
        }

        //Remove all edges pointing to this id

        void removeEdges(int nId) {
            for (outEdgeIterator ppEdge = begin(); ppEdge != end(); ppEdge++) {
                CEdge_relScale * pEdge = *ppEdge;
                if (pEdge->linkedTo(nId)) {
                    vEdges.erase(ppEdge);
                    ppEdge--;
                }
            }
        }

    private:
        CFullRelPose relPose;
        C3dPose absPoseOfFarId; //Several of these absolute positions exist for every one of the original nodes. Only the first will be linked-in to the map pos
        int nSecondId;
        static const int ID_NOT_SET = -1;
    public:

        CNode_relPos(const int id1, const int id2, const C3dNormalisedPoseWithSD & normalisedPose)
        : pPrevEdge(0), id1(id1), id2(id2), relPose(normalisedPose), nSecondId(ID_NOT_SET) {
            //Doesn't necessarily need any edges.
        }
        ~CNode_relPos();

        int secondId() const {
            return nSecondId;
        }

        class CNodeSortByBadness {
        public:

            bool operator()(const CNode_relPos * p1, const CNode_relPos * p2) const {
                return p1->getBadness() < p2->getBadness() || (p1->getBadness() == p2->getBadness() && p1 < p2);
            }
        };

        inline void setParent(CEdge_relScale * pParentEdge) {
            if (pParentEdge->otherNode(this)->id1 == id1
                    || pParentEdge->otherNode(this)->id2 == id1)
                nSecondId = id2;
            else if (pParentEdge->otherNode(this)->id1 == id2
                    || pParentEdge->otherNode(this)->id2 == id2)
                nSecondId = id1;
            else
                THROW("Bad parent edge ids");

            if (pParentEdge->secondNode() == this) {
                relPose.updateScale(pParentEdge->firstNode()->getScale(), pParentEdge->firstNode()->getSLAMscale(), pParentEdge->getRelScale());
            } else if (pParentEdge->firstNode() == this) {
                relPose.updateScale(pParentEdge->secondNode()->getScale(), pParentEdge->secondNode()->getSLAMscale(), pParentEdge->getRelScale().inverse());
            } else
                THROW("Bad parent edge");

            pPrevEdge = pParentEdge;
        }

        inline void setRootNode(bool bIsFirstComponent) {
            relPose.setRoot(bIsFirstComponent);
            pPrevEdge = 0;
            nSecondId = id2;
        }

        inline void setUnused() {
            relPose.setUnused();
            pPrevEdge = 0;
            nSecondId = ID_NOT_SET;
        }

        typedef CDynArray<CEdge_relScale *>::iterator outEdgeIterator;

        inline outEdgeIterator begin() {
            return vEdges.begin();
        }

        inline outEdgeIterator end() {
            return vEdges.end();
        }

        inline const CScale getSLAMscale() const {
            return relPose.getSLAMscale();
        }

        inline const CScale & getScale() const {
            return relPose.getScale();
        }

        inline const double getBadness() const {
            return relPose.getScale().getG_sq();
        }

        inline void assignNewScale(const CScale & newScale) {
            relPose.setSCOREScale(newScale);
        }

        inline double length() const {
            if(IS_DEBUG) CHECK(!relPose.hasPosition(), "Node is too bad to have a length");
            return relPose.length(); //Includes SCORE update. dAssignedAbsScale > 0 ? dAssignedAbsScale : dThisEdgeLength;
        }

        //Only if the framerate is constant

        double speed() const {
            double dSpeed = length() / fabs(id1 - id2);
            return dSpeed;
        }

        void updateAbsolutePos() {
            //updateScale();
            if (pPrevEdge) {
                CNode_relPos * pPrevNode = pPrevEdge->otherNode(this);
                if (pPrevNode->id1 == id1 || pPrevNode->id2 == id1) {

                    /*TB: Remove as infinite loop can start here... if (pPrevNode->pose(id1).t.length() > 1000 || relPose.position(false).t.length() > 1000) {
                        REPEAT(10000,
                                cout << "Prev: " << pPrevNode->pose(id1) << endl;
                                cout << "Next: " << relPose.position(false) << endl;
                                );
                    }*/
                    int depth=0;
                    absPoseOfFarId = pPrevNode->pose(id1,depth) + relPose.position(false);
                    if(IS_DEBUG) CHECK(nSecondId != id2, "nSecondId appears wrong")
                } else //we're positioning the near id from the far one
                {
                    int depth=0;
                    absPoseOfFarId = pPrevNode->pose(id2,depth) + relPose.position(true);
                    if(IS_DEBUG) CHECK(nSecondId != id1, "nSecondId appears wrong")
                }
            } else {
                //Reversed when accessed if needed...
                absPoseOfFarId = relPose.position(false); //Reset to origin, length 1
                nSecondId = id2;
            }
        }

        inline const C3dPose & pose(int id, int & depth) const {
            if(IS_DEBUG) CHECK(!relPose.hasPosition(), "Node is too bad to be used, should not have been assigned to a map pos");
            if (!pPrevEdge && id == id1) {
                static const C3dPose origin;
                return origin;
            }

            if(IS_DEBUG) CHECK(!pPrevEdge && nSecondId != id, "Requesting pose of the wrong id");

            if (nSecondId == id)
                return absPoseOfFarId;
            else
            {
                const CNode_relPos * pOther = pPrevEdge->otherNode(this);
                depth++;
                
                if(depth>3)
                {
                    cout << pPrevEdge->firstNode() << endl;
                    cout << pPrevEdge->secondNode() << endl;
                    cout << pPrevEdge->otherNode(this) << endl;
                    cout << pPrevEdge->getRelScale() << endl;
                    cout << "Returning incorrect pose (for time " << nSecondId << " rather than " << id << ")" << endl;
                    return absPoseOfFarId;
                }
                
                return pOther->pose(id, depth);
            }
        }

        void addEdge(CEdge_relScale * pEdge) {
            if(IS_DEBUG) CHECK(pEdge->firstNode() != this && pEdge->secondNode() != this, "Node mismatch");
            vEdges.push_back(pEdge); // CEdge_relScale(pThisNode, pOtherNode, relScale));
        }

        void writeEdge(std::ostream & toroGraphFile, const bool b2d) const {
            if(IS_DEBUG) CHECK(!relPose.hasPosition(), "Node is disconnected?? Or too bad");

            toroGraphFile << "EDGE" << (b2d ? " " : "3 ");
            toroGraphFile << id1 << ' ' << id2 << ' ';
            relPose.writeEdge(toroGraphFile, b2d);
            toroGraphFile << '\n';
        }

        void pp() const {
            std::cout << "Edge " << id1 << ' ' << id2 << ' ' << relPose << "\n";

            if (pPrevEdge)
                pPrevEdge->pp(this);
        }

    };

    class CMapPosition {
        const CNode_relPos * pPositionSource; //0 if does not have position OR IS ORIGIN
        bool bOrigin; //Don't strictly need...
        const int nId; //Not strictly needed...
    public:
        //CMapPosition() : pPositionSource(0) {}

        CMapPosition(bool bOrigin, int nId) : pPositionSource(0), bOrigin(bOrigin), nId(nId) {
        }

        bool hasPosition() const {
            return bOrigin || pPositionSource > 0;
        }

        const CNode_relPos * positionSource() const {
            return pPositionSource;
        }

        void setPosition(const CNode_relPos * pPos) {
            if(IS_DEBUG) CHECK(hasPosition(), "This map position hasn't been reset");
            if(IS_DEBUG) CHECK(bOrigin, "This map position is the origin");
            pPositionSource = pPos;
        }

        void reset() {
            pPositionSource = 0;
        }

        const C3dPose & pose() const {
            if(IS_DEBUG) CHECK((bool)pPositionSource == bOrigin, "Map position has no pPositionSource or does and is origin");
            if (bOrigin) {
                static const C3dPose origin;
                return origin;
            } else
            {
                int depth=0;
                return pPositionSource->pose(nId, depth);
            }
        }

        const CScale badness(bool bIsFirstComponent) const {
            if (bOrigin) {
                CScale s;
                s.setOrigin(bIsFirstComponent);
                return s;
            } else if (!pPositionSource)
                return CScale();
            else
                return pPositionSource->getScale();
        }

        double speed() const {
            return bOrigin ? 1 : pPositionSource->speed();
        }

        void write(std::ostream & toroGraphFile, int nMyId, bool bWriteEdge, const bool b2d) {
            if (hasPosition()) {
                toroGraphFile << "VERTEX";
                if (b2d)
                    toroGraphFile << " ";
                else
                    toroGraphFile << "3 ";
                toroGraphFile << nMyId << ' ';
                pose().write(toroGraphFile, b2d);
                toroGraphFile << '\n';

                if (!bOrigin && bWriteEdge)
                    pPositionSource->writeEdge(toroGraphFile, b2d); //Add the edge I got my position from, Will add more edges later
            }
        }

        void pp(int n) const {
            std::cout << "Position " << n << ":\n";
            if (pPositionSource) {
                std::cout << pose() << " from:\n";
                pPositionSource->pp();
            }
        }
    };

    class CSLAMMap {
        typedef map2<CIds, CNode_relPos *> TNodes;
        TNodes allNodes;

        typedef map2<int, CMapPosition> TMapPositions;
        //map2<int, TMapPositions > aAllMapPositions;
        typedef std::vector<CNode_relPos *> TSortedPositions;
        //map2<int, TSortedPositions > aUpdateScalesInOrder, aSortedEdgesNotInMST;
        typedef set2<CNode_relPos *, CNode_relPos::CNodeSortByBadness> TNodesWithPaths;

        const bool bSetNearbyOrigin;
        const bool bVerbose;

        class CComponentData {
            TMapPositions allMapPositions;
            TSortedPositions updateScalesInOrder, sortedEdgesNotInMST;
            bool bDirty, bDirtyScales, bSetOrigin, bFirstComp;
            TNodesWithPaths nodesWithPaths;
            CIds rootIds;
        public:

            CComponentData() : bDirty(true), bDirtyScales(true), bSetOrigin(false), bFirstComp(false) {
            }

            TMapPositions & getAllMapPositions() {
                return allMapPositions;
            }

            const TMapPositions & getAllMapPositions() const {
                return allMapPositions;
            }

            TSortedPositions & getUpdateScalesInOrder() {
                return updateScalesInOrder;
            }

            TSortedPositions & getSortedEdgesNotInMST() {
                return sortedEdgesNotInMST;
            }

            bool isDirty() const {
                return bDirty;
            }

            bool hasOrigin() const {
                return bSetOrigin;
            }

            bool isScalesDirty() const {
                return bDirtyScales;
            }
            //bool isDirty() const { return bDirty; }

            TNodesWithPaths & getNodesWithPaths() {
                return nodesWithPaths;
            }

            const CIds & getRootIds() const {
                if(IS_DEBUG) CHECK(!hasOrigin(), "No root set");
                return rootIds;
            }

            void setClean() {
                bDirty = bDirtyScales = false;
            }

            void setDirty() {
                bDirty = true;
            }

            void removeRoot() {
                cout << "Erasing rootIds " << rootIds.id1() << " to " << rootIds.id2() << endl;
                bSetOrigin = false;
                rootIds = CIds();
            }

            void setScaleDirty() {
                bDirtyScales = true;
            }

            void setOrigin(const CIds & origin/*, bool bNewComponent*/) {
                bSetOrigin = true;
                if (origin != rootIds) {
                    rootIds = origin;
                    bDirty = true;
                }
            }

            void setFirstComponent() {
                bFirstComp = true;
            }

            bool isFirstComponent() const {
                return bFirstComp;
            }

        };

    protected:
        map2<int, CComponentData> aComponentData;
    private:
        CDynArrayOwner<CEdge_relScale> aEdgeOwner;

        CEdge_relScale * addEdge(CNode_relPos * pThisNode, CNode_relPos * pOtherNode, const CRelScale & relScale) {
            if(IS_DEBUG) CHECK(!relScale.notTooBad(), "Scale OOB");

            CEdge_relScale * pEdge = new CEdge_relScale(pThisNode, pOtherNode, relScale);
            pThisNode->addEdge(pEdge);
            pOtherNode->addEdge(pEdge);

            aEdgeOwner.push_back(pEdge);

            return pEdge;
        }

        void checkPositionsUptoDate(CComponentData & comp) {
            if (comp.isDirty())
                updateSPTree(comp);
            else if (comp.isScalesDirty())
                updateScales(comp);
        }

        void updateSPTree(CComponentData & comp) HOT;

        void updateScales(CComponentData & comp) {
            if(IS_DEBUG) CHECK(comp.isDirty() || !comp.isScalesDirty(), "Unnecessary/wrong update");
            if(IS_DEBUG) CHECK(!comp.hasOrigin(), "No origin yet");

            TSortedPositions & updateScalesInOrder = comp.getUpdateScalesInOrder();
            for (TSortedPositions::iterator ppNode = updateScalesInOrder.begin(); ppNode != updateScalesInOrder.end(); ppNode++) {
                (*ppNode)->updateAbsolutePos();
            }
            comp.setClean();

            std::cout << "Finished scale update, " << updateScalesInOrder.size() << " nodes have scales" << std::endl;
        }

        CNode_relPos * tryAddPose(int nId1, int nId2, const C3dNormalisedPoseWithSD & pose12) {
            const CIds pose12id(nId1, nId2);
            TNodes::iterator itPose12 = allNodes.find(pose12id);
            if (itPose12 == allNodes.end()) {
                CNode_relPos * pRelPos12 = new CNode_relPos(nId1, nId2, pose12);
                allNodes.init(pose12id, pRelPos12);
                return pRelPos12;
            } else
                return itPose12->second;
        }
        int nCurrentComponent;
        static const int UNINIT_COMPONENT = -1;
    public:

        bool haveMap() const {
            return nCurrentComponent != UNINIT_COMPONENT;
        }

        int currentComponent() const {
            if(IS_DEBUG) CHECK(nCurrentComponent == UNINIT_COMPONENT, "No edges added yet")
            return nCurrentComponent;
        }

        int currentPos() const {
            if(IS_DEBUG) CHECK(nCurrentComponent == UNINIT_COMPONENT, "No edges added yet")
            return aComponentData[nCurrentComponent].getAllMapPositions().backKey();
        }

        int firstPosThisComponent() const {
            if(IS_DEBUG) CHECK(nCurrentComponent == UNINIT_COMPONENT, "No edges added yet")
            return aComponentData[nCurrentComponent].getAllMapPositions().topKey();
        }

        CSLAMMap(const bool bSetNearbyOrigin, const bool bVerbose) : bSetNearbyOrigin(bSetNearbyOrigin), bVerbose(bVerbose), nCurrentComponent(UNINIT_COMPONENT) {
        }

        virtual ~CSLAMMap() {
            for (TNodes::iterator pNode = allNodes.begin(); pNode != allNodes.end(); pNode++)
                delete pNode->second;
        };

        int getComponentId(const int id1) {
            return aComponents[id1];
        }

        CComponentData & getComponent(const int id1, const int id2) {
            int nComp1 = aComponents[id1];
            int nComp2 = aComponents[id2];
            if(IS_DEBUG) CHECK(nComp1 != nComp2, "Setting component root with bad link");
            return aComponentData[nComp1];
        }


        //Set the edge defined to have length 1 pointing North.

        void setOrigin(const int id1, const int id2) {
            CComponentData & comp = getComponent(id1, id2);
            comp.setOrigin(CIds(id1, id2));
        }

        /*//Get the badness of the scale calculated for this node (rel pos)
        const CScale getRelPosBadness(const int id1, const int id2)
        {
                CComponentData & comp = getComponent(id1, id2);
                checkPositionsUptoDate(comp);

                TNodes::const_iterator pMP = allNodes.find(CIds(id1, id2));
                if(pMP != allNodes.end())
                        return pMP->second->getScale();
                else
                        return CScale();
        }*/
        inline const CScale getSLAMscale(const int id1, const int id2) {
            TComponents::iterator pComp = aComponents.find(id1);

            if (pComp != aComponents.end()) {
                CComponentData & comp = aComponentData[pComp->second];
                checkPositionsUptoDate(comp);

                TNodes::const_iterator pMP = allNodes.find(CIds(id1, id2));
                if (pMP != allNodes.end())
                    return pMP->second->getScale();
            }

            return CScale();
        }

        class CSortNodesByDist {
        public:

            bool operator()(const CNode_relPos * p1, const CNode_relPos * p2) const {
                return (p1->shortestPathLength() < p2->shortestPathLength()) || (p1->shortestPathLength() == p2->shortestPathLength() && p1 < p2);
            }
        };

        void optimiseScaleAroundLoops(CScaleOptimiser & scaleOptimiser, /*CScaleOptimiser & scaleOptimiser2,*/ const int nComponent, const int t, const bool bVerbose) HOT;

        void printSpeeds(int nComponent) {
            CComponentData & comp = aComponentData[nComponent];
            checkPositionsUptoDate(comp);

            TMapPositions & allMapPositions = comp.getAllMapPositions();

            for (TMapPositions::iterator pNode = allMapPositions.begin(); pNode != allMapPositions.end(); pNode++) {
                CMapPosition & mapPos = pNode->second;
                if (mapPos.hasPosition()) {
                    std::cout << pNode->first << ',' << mapPos.speed() << std::endl;
                }
            }

        }

        void adjustScale(int id1, int id2, const CScale & SCOREScale) {
            CComponentData & comp = getComponent(id1, id2);

            if (comp.isDirty())
                updateSPTree(comp); //Don't update scales each time they are invalidated, but check all possible nodes have a SLAM scale.

            comp.setScaleDirty();

            TNodes::const_iterator pMP = allNodes.find(CIds(id1, id2));
            if(IS_DEBUG) CHECK(pMP == allNodes.end(), "adjustScale: ids not found");

            CNode_relPos * pNode = pMP->second;
            pNode->assignNewScale(SCOREScale);
        }

        //Return the id of the other position from which this position is calculated

        int positionSource(const int id) {
            const int nComponent = aComponents[id];
            CComponentData & comp = aComponentData[nComponent];
            checkPositionsUptoDate(comp);

            TMapPositions & allMapPositions = comp.getAllMapPositions();

            TMapPositions::const_iterator pMP = allMapPositions.find(id);
            if (pMP != allMapPositions.end() && pMP->second.hasPosition()) {
                const CNode_relPos * pN = pMP->second.positionSource();
                if (pN) {
                    //if(IS_DEBUG) CHECK(id != pN->secondId(), "This position is clearly not the source of the position for this node")
                    if (id != pN->secondId()) {
                        //pN->pp();
                        cout << "TODO Restore check here..?\n";
                    }
                    return pN->id1;
                } else
                    //This is the origin
                    return -1;
            } else
                THROW("Trying to get the source of a node with no position (maybe its the root node or something?)");
        }

        void pp(const int nComponent) {
            int numPrinted = 0;
            const int NUM_TO_PRINT = 10;
            CComponentData & comp = aComponentData[nComponent];
            checkPositionsUptoDate(comp);


            TMapPositions & allMapPositions = comp.getAllMapPositions();

            for (TMapPositions::const_iterator pMP = allMapPositions.begin(); pMP != allMapPositions.end(); pMP++) {
                if (pMP->second.speed() > 100 || pMP->second.speed() < 0.01 || NUM_TO_PRINT > (int) allMapPositions.size()) {
                    pMP->second.pp(pMP->first);
                    numPrinted++;
                    if (numPrinted > NUM_TO_PRINT)
                        break;
                }
            }
        }

        /*int firstPosition(const int nComponent)
        {
                CComponentData & comp = aComponentData[nComponent];
                checkPositionsUptoDate(comp);


                return aAllMapPositions[nComponent].topKey();
        }*/

        const C3dPose * position(const int id, const int nComp_in) {
            const int nComponent = aComponents.ifExists(id, (int) UNINIT_COMPONENT);
            if (nComp_in != nComponent)
                return 0;

            CComponentData & comp = aComponentData[nComponent];
            checkPositionsUptoDate(comp);

            TMapPositions & allMapPositions = comp.getAllMapPositions();

            TMapPositions::const_iterator pMP = allMapPositions.find(id);
            if (pMP != allMapPositions.end() && pMP->second.hasPosition()) {
                const C3dPose * pPose = &(pMP->second.pose());
                /*double D = pMP->second.positionSource() ? pMP->second.positionSource()->getScale().getD() : 0;
                if(fabs(D) > 3)
                        return 0;
                else*/
                return pPose;
            } else
                return 0;
        }

        const CScale positionQuality(const int id) {
            const int nComponent = aComponents.ifExists(id, (int) UNINIT_COMPONENT);
            if (nComponent == UNINIT_COMPONENT)
                return CScale();

            CComponentData & comp = aComponentData[nComponent];
            checkPositionsUptoDate(comp);

            TMapPositions & allMapPositions = comp.getAllMapPositions();

            TMapPositions::const_iterator pMP = allMapPositions.find(id);
            if (pMP != allMapPositions.end() && pMP->second.hasPosition()) {
                return pMP->second.badness(comp.isFirstComponent());
            } else
                return CScale();
        }

        typedef map2<int, int> TComponents;
        TComponents aComponents;

        //Adds a scale-edge (in both directions). Both kinds of nodes added automatically,

        void addBiDiScaleLink(int nId1, int nId2, int nId3, const C3dNormalisedPoseWithSD & pose12, const C3dNormalisedPoseWithSD & pose23, const CRelScale & relscale123) {
            if (fabs(relscale123.get_d()) > 3)
                cout << "BADSCALE Adding " << nId1 << " to " << nId2 << " to " << nId3 << " with rel scale " << relscale123 << endl;

            if(IS_DEBUG) CHECK(nId1 == nId3, "Circular link not allowed");
            //if(IS_DEBUG) CHECK(!badness.notTooBad(), "Link too bad to be used");
            if(IS_DEBUG) CHECK(!relscale123.notTooBad(), "Scale OOB");
            //Only the connected component this link extends from is 'dirty'. May be 2 of these (if joining components!)

            const int local_UNINIT_COMPONENT = UNINIT_COMPONENT;
            const int nComp1 = aComponents.ifExists(nId1, local_UNINIT_COMPONENT);
            const int nComp2 = aComponents.ifExists(nId2, local_UNINIT_COMPONENT);
            const int nComp3 = aComponents.ifExists(nId3, local_UNINIT_COMPONENT);

            //Should now have at most 2 components. Possibly all uninit.
            bool bMergingComponents = false;
            int anComps[3] = {nComp1, nComp2, nComp3};
            std::sort(anComps, anComps + 3);

            int nMinComp = UNINIT_COMPONENT, nMaxComp = UNINIT_COMPONENT;
            for (int nCompIdx = 0; nCompIdx < 3; nCompIdx++) {
                if (nMinComp == UNINIT_COMPONENT)
                    nMinComp = anComps[nCompIdx];

                nMaxComp = anComps[nCompIdx];
            }

            if (nMaxComp == UNINIT_COMPONENT) //All uninit.
            {
                //start a new component
                nMinComp = aComponentData.size() ? aComponentData.backKey() + 1 : 0;
                std::cout << "Initialising component " << nMinComp << std::endl;
                CComponentData & comp = aComponentData.init(nMinComp);
                CIds newOriginIds = ((nId1 < nId2) ? CIds(nId1, nId2) : CIds(nId2, nId1));

                if (nMinComp == 0)
                    comp.setFirstComponent();

                comp.setOrigin(newOriginIds);
            } else if (nMaxComp == nMinComp) {
                //Within existing componend
            } else {
                std::cout << "Merging components " << nMaxComp << " into " << nMinComp << std::endl;
                bMergingComponents = true;
                if(IS_DEBUG) CHECK(anComps[0] == UNINIT_COMPONENT, "Don't expect to also have uninit nodes here???");
                aComponentData.erase(nMaxComp);

                for (TComponents::iterator pCompId = aComponents.begin(); pCompId != aComponents.end(); pCompId++)
                    if (pCompId->second == nMaxComp)
                        pCompId->second = nMinComp;
            }


            if (nComp1 != nMinComp)
                aComponents.initOrSet(nId1, nMinComp);
            if (nComp2 != nMinComp)
                aComponents.initOrSet(nId2, nMinComp);
            if (nComp3 != nMinComp)
                aComponents.initOrSet(nId3, nMinComp);


            CNode_relPos * pPose12 = 0, * pPose23 = 0;
            if (nId1 < nId2)
                pPose12 = tryAddPose(nId1, nId2, pose12);
            else
                pPose12 = tryAddPose(nId2, nId1, pose12.reverse());


            if (nId2 < nId3)
                pPose23 = tryAddPose(nId2, nId3, pose23);
            else
                pPose23 = tryAddPose(nId3, nId2, pose23.reverse());

            NSLAMMap::CEdge_relScale * pNewEdge = addEdge(pPose12, pPose23, relscale123);

            CComponentData & comp = aComponentData[nMinComp];

            int nMaxTime = max3(nId1, nId2, nId3);
            int nMinTime = min3(nId1, nId2, nId3);

            //Try to avoid setting dirty
            if (bMergingComponents || //merging
                    nMaxComp == UNINIT_COMPONENT || // new component
                    nMaxTime - nMinTime > 25 // loop closure
                    ) {
                comp.setDirty();
                updateSPTree(comp);
            } else {
                if (!pPose23->parent() && pPose12->parent()) {
                    pPose23->setParent(pNewEdge);
                    pPose23->updateAbsolutePos();
                } else if (!pPose12->parent() && pPose23->parent()) {
                    pPose12->setParent(pNewEdge);
                    pPose12->updateAbsolutePos();
                } else if (pPose12->parent() && pPose23->parent()) { //Both have positions, now adding edge between poses. HACK, doesn't update the actual position
                    if (pPose23->getSLAMscale().isMoreAccurateThan(pPose12->getSLAMscale() + pNewEdge->getRelScale()) && pPose12->parent()) {
                        pPose23->setUnused();
                        pPose23->setParent(pNewEdge);
                        pPose23->updateAbsolutePos();
                    } else if ((pPose23->getSLAMscale() + pNewEdge->getRelScale()).isMoreAccurateThan(pPose12->getSLAMscale()) && pPose23->parent()) {
                        pPose12->setUnused();
                        pPose12->setParent(pNewEdge);
                        pPose12->updateAbsolutePos();
                    }
                } else //Neither have parents (one may or may not be the origin)
                {
                    comp.setDirty();
                }
            }

            nCurrentComponent = nMinComp;
        }

        void deleteId(int nId) {
            //Find ALL edges start/ending at this id and erase
            CDynArrayOwner<CNode_relPos> vNodesToDelete;
            for (TNodes::iterator pEdges = allNodes.begin(); pEdges != allNodes.end(); pEdges++) {
                pEdges->second->removeEdges(nId);

                if (pEdges->first.id1() == nId || pEdges->first.id2() == nId)
                    vNodesToDelete.push_back(pEdges->second);
            }
            for (CDynArrayOwner<CNode_relPos>::const_iterator ppNode = vNodesToDelete.begin(); ppNode != vNodesToDelete.end(); ppNode++)
                allNodes.erase(CIds((*ppNode)->id1, (*ppNode)->id2));

            const int nComponent = aComponents.ifExists(nId, (int) UNINIT_COMPONENT);
            if (nComponent != UNINIT_COMPONENT) {
                aComponents.erase(nId);

                CComponentData & comp = aComponentData[nComponent];
                if (comp.hasOrigin() && !allNodes.exists(comp.getRootIds()))
                    comp.removeRoot();
                comp.setDirty();
            }
        }

        typedef std::vector<const CNode_relPos *> TOnePosSources;
        typedef map2<int, TOnePosSources > TPosSources;
        typedef set2<const CNode_relPos *> TGraphEdges;
        typedef map2<int, int> TVisAtDepthMap;

        bool tryFindShortPath(const int nDepthToBottom, const int nId1, const int nId2, const TGraphEdges & augmentedSPTree, TPosSources & positionSources, TVisAtDepthMap & aVisitedAtDepth) {
            if (nId1 == nId2)
                return true;
            else if (nDepthToBottom == 0)
                return false;
            else {
                int nLastVisited = -1;
                TVisAtDepthMap::iterator pVD = aVisitedAtDepth.find(nId1);
                if (pVD == aVisitedAtDepth.end())
                    pVD = aVisitedAtDepth.insert(std::pair<const int, int>(nId1, nDepthToBottom)).first;
                else {
                    //nLastVisited = pVD->second; //Spurious warning apparently...
                    const int & n = pVD->second;
                    nLastVisited = n;
                }

                //int nVisited = aVisitedAtDepth.ifExists(nId1, -1); //inefficient repeated find...
                /*if(nVisited == -1)
                {
                        aVisitedAtDepth.init(nId1, nTemporalConnectivity);
                }*/
                //else
                //{
                if (nLastVisited < nDepthToBottom) {
                    pVD->second = nDepthToBottom; //Increase visited level and recurse from here again
                    if(IS_DEBUG) CHECK(aVisitedAtDepth[nId1] != nDepthToBottom, "Failed to set depth")

                    TPosSources::const_iterator pPosSource = positionSources.find(nId1);
                    if (positionSources.end() != pPosSource)
                        for (std::vector<const CNode_relPos *>::const_iterator pOutEdge = pPosSource->second.begin(); pOutEdge != pPosSource->second.end(); pOutEdge++)
                            if (tryFindShortPath(nDepthToBottom - 1, (*pOutEdge)->otherId(nId1), nId2, augmentedSPTree, positionSources, aVisitedAtDepth))
                                return true;
                }
                //}
            }
            return false;
        }

        void saveToroMap(const char * szFilename, const int nTemporalConnectivity, const int nComponent, const bool bAddConsecutiveEdges, const bool b2d) {
            CComponentData & comp = aComponentData[nComponent];
            checkPositionsUptoDate(comp);


            TMapPositions & allMapPositions = comp.getAllMapPositions();

            int nNumPositions = (int) allMapPositions.size();
            if (nNumPositions < 2)
                return;
            /*
             *  A set of simple text messages to represent nodes and edges of the graph. Note that examples files are in the repository, see folder data.

                    Format of the 2D graph files:

                    Every line in the file specifies either one vertex or one edge

                    The vertices are specified as follws: VERTEX2 id x y orientation (A 2D node in the graph)

                    EDGE2 observed_vertex_id observing_vertex_id forward sideward rotate inf_ff inf_fs inf_ss inf_rr inf_fr inf_sr (A 2D-edge in the graph. inf_xx are the information matrix entries of the constraint)

                    EQUIV id1 id2 (Equivalence constraints between nodes. It merges the node id1 and id2 wrt to the constraint between both vertices.)

                    Format of the 3D graph files:

                    Every line in the file specifies either one vertex or one edge

                    The vertices are specified as follws: VETREX3 x y z phi theta psi

                    The edges are specified as follows: EDGE3 observed_vertex_id observing_vertex_id x y z roll pitch yaw inf_11 inf_12 .. inf_16 inf_22 .. inf_66 (the information matrix is specified via its upper triangular block that means 21 values).
             */
            const bool bSPANNER = true;

            std::ofstream toroGraphFile(szFilename);

            //where does the node at int get its position from, can be many places
            TPosSources positionSources;
            TGraphEdges augmentedSPTree;

            int nLastTime = -1;

            for (TMapPositions::iterator pNode = allMapPositions.begin(); pNode != allMapPositions.end(); pNode++) {
                CMapPosition & mapPos = pNode->second;

                if (mapPos.hasPosition()) //This node and its corresponding edge should be added
                {
                    mapPos.write(toroGraphFile, pNode->first, bSPANNER || nTemporalConnectivity > 0, b2d);

                    if (bAddConsecutiveEdges && nLastTime >= 0 && !b2d) {
                        if (b2d) {
                            toroGraphFile << "EDGE3 " << nLastTime << ' ' << pNode->first << " 0 0 0 0 0 0 ";
                            for (int i = 0; i < 6; i++) {
                                toroGraphFile << " " << (((i < 3) ? 0.1 : 10) / (pNode->first - nLastTime));
                                for (int j = i + 1; j < 6; j++) {
                                    toroGraphFile << " 0 ";
                                }
                            }
                        } else {
                            toroGraphFile << "EDGE " << nLastTime << ' ' << pNode->first << " 0 0 0 ";
                            for (int i = 0; i < 3; i++) {
                                toroGraphFile << " " << (0.1 / (pNode->first - nLastTime));
                                for (int j = i + 1; j < 3; j++) {
                                    toroGraphFile << " 0 ";
                                }
                            }
                        }
                        toroGraphFile << std::endl;
                    }

                    if (mapPos.positionSource()) {
                        augmentedSPTree.insert(mapPos.positionSource());
                        positionSources.initOrGet(mapPos.positionSource()->id1).push_back(mapPos.positionSource());
                        positionSources.initOrGet(mapPos.positionSource()->id2).push_back(mapPos.positionSource());
                    }

                    nLastTime = pNode->first;
                }
            }

            if (bSPANNER) //SPANNER algorithm. Althofer-etal-1993
            {
                //Edge set augmentedSPTree includes SP tree and all connected nodes.
                //foreach edge (sorted by weight):
                //If edge not in augmentedSPTree:
                //Find s.p. in augmentedSPTree, if length > nTemporalConnectivity add it
                TSortedPositions & sortedEdgesNotInMST = comp.getSortedEdgesNotInMST();
                for (TSortedPositions::const_iterator ppEdge = sortedEdgesNotInMST.begin(); ppEdge != sortedEdgesNotInMST.end(); ppEdge++) {
                    CNode_relPos * pEdge = *ppEdge;
                    map2<int, int> aVisitedAtDepth;
                    bool bShortPath = tryFindShortPath(nTemporalConnectivity, pEdge->id1, pEdge->id2, augmentedSPTree, positionSources, aVisitedAtDepth);

                    if (!bShortPath) {
                        pEdge->writeEdge(toroGraphFile, b2d);
                        augmentedSPTree.insert(pEdge);

                        positionSources.initOrGet(pEdge->id1).push_back(pEdge);
                        positionSources.initOrGet(pEdge->id2).push_back(pEdge);
                    }
                }
            } else {

                /*
                 * OLD METHOD:
                 *
                 * Want a graph that uses a minimal number of links to position all nodes.
                 *
                 * Option 1) Every node has a link forwards and backwards in time
                 *
                 * Option 2) Every node has a sequence of n<K edges to both forwards and backwards
                 *
                 * Option 3) Find the shortest path (in terms of additional edges) linking each node forwards/backwards in time
                 *
                 */

                if (nTemporalConnectivity <= 0) // Add every link that's not too bad
                {
                    //Write every edge
                    for (TNodes::iterator ppRelPos = allNodes.begin(); ppRelPos != allNodes.end(); ppRelPos++) {
                        //If the future position gets its position from the future then add this node to candidates to allow it to get its position from the past
                        const CNode_relPos * pRelPos = ppRelPos->second;
                        if (pRelPos->getScale().hasScale())
                            pRelPos->writeEdge(toroGraphFile, b2d);
                    }
                } else // For every pose that is not linked forward/backward to a node within nTemporalConnectivity time steps, add the best link from that node forward/backwards
                    // Not ideal, as leaves with length 1 etc. will gain a link forwards/backwards
                    // TODO: group sets of nTemporalConnectivity nodes together and ensure the *group* has links forwards and backwards
                {

                    //Now find disconnected temporally successive nodes nodes. Need to cope with start+end
                    std::set<int> nodesWithNoLinkToFuture, nodesWithNoLinkToPast;

                    TMapPositions::iterator pBackNode = allMapPositions.end();
                    pBackNode--;
                    for (TMapPositions::iterator pNode = allMapPositions.begin(); pNode != pBackNode; pNode++) {
                        int nTime = pNode->first;
                        bool bHaveFutureLink = false, bHavePastLink = false;

                        std::vector<const CNode_relPos *> & vMyPosSources = positionSources[nTime];
                        for (std::vector<const CNode_relPos *>::const_iterator pPS = vMyPosSources.begin(); pPS != vMyPosSources.end(); pPS++) {
                            const CNode_relPos * pRelPos = *pPS;

                            int otherTime = (nTime != pRelPos->id1) ? pRelPos->id1 : pRelPos->id2;
                            if (otherTime - nTime > 0 && otherTime - nTime <= nTemporalConnectivity)
                                bHaveFutureLink = true;
                            else if (nTime - otherTime >= nTemporalConnectivity)
                                bHavePastLink = true;
                        }

                        if (!bHaveFutureLink)
                            nodesWithNoLinkToFuture.insert(nTime);

                        if (!bHavePastLink)
                            nodesWithNoLinkToPast.insert(nTime);
                    }

                    typedef std::vector<const CNode_relPos *> TLinkCandidates;
                    typedef std::map<int, TLinkCandidates > TExtraLinkCandidates;
                    TExtraLinkCandidates candidatesForFutureLinks, candidatesForPastLinks;

                    for (TNodes::iterator ppRelPos = allNodes.begin(); ppRelPos != allNodes.end(); ppRelPos++) {
                        //If the future position gets its position from the future then add this node to candidates to allow it to get its position from the past
                        const CNode_relPos * pRelpos = ppRelPos->second;
                        if (pRelpos->getScale().hasScale()) {
                            const int id_lo = std::min<int>(pRelpos->id1, pRelpos->id2);
                            const int id_hi = std::max<int>(pRelpos->id1, pRelpos->id2);

                            if (nodesWithNoLinkToFuture.find(id_lo) != nodesWithNoLinkToFuture.end())
                                candidatesForFutureLinks[id_lo].push_back(pRelpos);

                            if (nodesWithNoLinkToPast.find(id_hi) != nodesWithNoLinkToPast.end())
                                candidatesForPastLinks[id_hi].push_back(pRelpos);

                            /*if(nodesWithNoLinkToFuture.find(id1) != nodesWithNoLinkToFuture.end()
                                            || nodesWithNoLinkToPast.find(id2) != nodesWithNoLinkToPast.end())
                            {
                                    pRelpos->writeEdge(toroGraphFile);
                            }*/
                        }
                    }

                    //Now for each node choose the best future/past link found
                    for (TExtraLinkCandidates::iterator pvLinkCandidate = candidatesForFutureLinks.begin(); pvLinkCandidate != candidatesForFutureLinks.end(); pvLinkCandidate++) {
                        TLinkCandidates & vLinkCandidates = pvLinkCandidate->second;
                        if(IS_DEBUG) CHECK(vLinkCandidates.size() <= 0, "Vector but no link candidates found");
                        const CNode_relPos * pBestLC = *(std::max_element(vLinkCandidates.begin(), vLinkCandidates.end(), CNode_relPos::CNodeSortByBadness()));
                        pBestLC->writeEdge(toroGraphFile, b2d);

                        //and remove the second id from list of nodes needing links to past:
                        candidatesForPastLinks.erase(std::max<int>(pBestLC->id1, pBestLC->id2));
                    }
                    for (TExtraLinkCandidates::iterator pvLinkCandidate = candidatesForPastLinks.begin(); pvLinkCandidate != candidatesForPastLinks.end(); pvLinkCandidate++) {
                        TLinkCandidates & vLinkCandidates = pvLinkCandidate->second;
                        if(IS_DEBUG) CHECK(vLinkCandidates.size() <= 0, "Vector but no link candidates found");
                        const CNode_relPos * pBestLC = *(std::max_element(vLinkCandidates.begin(), vLinkCandidates.end(), CNode_relPos::CNodeSortByBadness()));
                        pBestLC->writeEdge(toroGraphFile, b2d);
                    }
                }
            }
        }
    };
};
#endif /* SLAMMAP_H_ */
