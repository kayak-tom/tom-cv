#include "slamMap.h"


namespace NSLAMMap {
    CEdge_relScale::CEdge_relScale(CNode_relPos * pNode1, CNode_relPos * pNode2, const CRelScale relScale) : /*id1(id1),id2(id2),id3(id3)*/pNode1(pNode1), pNode2(pNode2), relScale(relScale), dOptimisedRelScale(UNINIT_RELSCALE()) {
        if(IS_DEBUG) CHECK(!relScale.notTooBad(), "CEdge_relScale: Edge added is too bad to ever be used");
        if(IS_DEBUG) CHECK(!pNode1 || !pNode2, "CEdge_relScale: Uninit node *");
        //if(IS_DEBUG) CHECK(pNode1->id2 != pNode2->id1, "Nodes not successive");
    }
    bool CEdge_relScale::linkedTo(int nId) const {
        return pNode1->id1 == nId
                || pNode1->id2 == nId
                || pNode2->id1 == nId
                || pNode2->id2 == nId;
    }
    void CEdge_relScale::pp(const CNode_relPos * pLast) const {
        std::cout << relScale << ' ' << pNode1->id1 << '-' << pNode1->id2 << '-' << pNode2->id2 << "\n";
        otherNode(pLast)->pp();
    }
    CNode_relPos::~CNode_relPos() {
        /*for(outEdgeIterator ppEdge = begin(); ppEdge != end(); ppEdge++)
        {
                CEdge_relScale * pEdge = *ppEdge;
                if(pEdge->firstNode() == this)
                        delete pEdge;
        } owned by aEdgeOwner instead */
    }
    void CEdge_relScale::setOptimisedScale(const double optimisedD_in, const bool bVerbose) {
        if (bVerbose)
            std::cout << "Set scale for edge " << firstNode()->id1 << ',' << firstNode()->id2 << " to " << secondNode()->id1 << ',' << secondNode()->id2 << " from " << relScale.get_d() << " to " << optimisedD_in << endl;
        dOptimisedRelScale = optimisedD_in;
    }
    void CSLAMMap::updateSPTree(CComponentData & comp) {
        if(IS_DEBUG) CHECK(!comp.isDirty(), "Unnecessary update");

        if (!comp.hasOrigin()) {
            std::cout << "No origin yet\n";
            return;
        }

        TMapPositions & allMapPositions = comp.getAllMapPositions();

        if (bSetNearbyOrigin && allMapPositions.size() > 0) {
            //Choose the edge used to position the last node as origin
            /*
                                    const CNode_relPos * pPositionSource = allMapPositions.back().positionSource();
                                    if(pPositionSource)
                                    {
                                            comp.setOrigin(CIds(pPositionSource->id1, pPositionSource->id2));
                                    }*/


            TMapPositions::const_iterator pMapPos = allMapPositions.end();
            for (pMapPos--;; pMapPos--) {
                const CNode_relPos * pPositionSource = pMapPos->second.positionSource();
                if (pPositionSource) {
                    CIds candidateIds(pPositionSource->id1, pPositionSource->id2);
                    if (allNodes.exists(candidateIds)) {
                        comp.setOrigin(candidateIds);
                        break;
                    }
                }

                if (pMapPos == allMapPositions.begin())
                    THROW("No possible root ids in this component"); //I think it may be possible to get here if a node in a small component is erased...
            }
        }

        CStopWatch s;
        s.startTimer();

        //Perform Dijkstra's algorithm from id1,id2.
        //For all nodes set 'too bad' and clear pPrevEdge
        for (TNodes::iterator pNode = allNodes.begin(); pNode != allNodes.end(); pNode++)
            pNode->second->setUnused();


        for (TMapPositions::iterator pNode = allMapPositions.begin(); pNode != allMapPositions.end(); pNode++)
            pNode->second.reset();

        allMapPositions.clear(); //Todo: Pretty sure this should be here to fully reset map...

        TNodesWithPaths & nodesWithPaths = comp.getNodesWithPaths();
        nodesWithPaths.clear();

        TSortedPositions & updateScalesInOrder = comp.getUpdateScalesInOrder();
        updateScalesInOrder.clear();

        TSortedPositions & sortedEdgesNotInMST = comp.getSortedEdgesNotInMST();
        sortedEdgesNotInMST.clear();

        const CIds & rootIds = comp.getRootIds();
        if (bVerbose) cout << "Root ids: " << rootIds.id1() << ',' << rootIds.id2() << endl;
        CNode_relPos * pRootNode = allNodes[rootIds];

        //Make origin for this component
        allMapPositions.insert(std::pair<int, CMapPosition > (rootIds.id1(), CMapPosition(true, rootIds.id1())));

        pRootNode->setRootNode(comp.isFirstComponent());
        nodesWithPaths.insert(pRootNode);

        while (nodesWithPaths.size() > 0) {
            CNode_relPos * pActiveNode = nodesWithPaths.pop(); //Node with shortest dist
            if (bVerbose) std::cout << "Active node " << pActiveNode->id1 << ' ' << pActiveNode->id2 << " near ";

            //Not actually needed: nodesWithShortestPaths.insert(pActiveNode);
            pActiveNode->updateAbsolutePos(); //Compute position and scale

            TMapPositions::iterator pMP = allMapPositions.find(pActiveNode->secondId());
            if (pMP == allMapPositions.end())
                pMP = allMapPositions.insert(std::pair<int, CMapPosition > (pActiveNode->secondId(), CMapPosition(false, pActiveNode->secondId()))).first;

            if (!pMP->second.hasPosition()) {
                pMP->second.setPosition(pActiveNode); //This map position now knows its immediate descendant, and that descendent has a position
                updateScalesInOrder.push_back(pActiveNode);
                //std::cout << "Positioning node " << pActiveNode->secondId() << "\n";
            } else {
                sortedEdgesNotInMST.push_back(pActiveNode);
                //if(bVerbose) std::cout << pActiveNode->getScale() << "=scale not in map\n";
            }

            //Update distances of outgoing links
            //Shortest paths depend only on SLAM scales, not on SCORE scales
            for (CNode_relPos::outEdgeIterator ppOutEdge = pActiveNode->begin(); ppOutEdge != pActiveNode->end(); ppOutEdge++) {
                CEdge_relScale * pOutEdge = *ppOutEdge;
                CScale newBadness = pActiveNode->getSLAMscale() + pOutEdge->getRelScale(); //Could reduce amount of computation here... but NB direction not necessarily computed correctly
                CNode_relPos * pOtherNode = pOutEdge->otherNode(pActiveNode);

                if (bVerbose) std::cout << '(' << pOtherNode->id1 << ',' << pOtherNode->id2 << ") ";

                CScale currentBadness = pOtherNode->getSLAMscale(); //Might be too bad...

                if (bVerbose) std::cout << currentBadness << " New=" << newBadness << ',';
                if (newBadness.isMoreAccurateThan(currentBadness)) {
                    if (currentBadness.hasScale()) {
                        if (bVerbose) std::cout << "RELOC ";
                        nodesWithPaths.erase(pOtherNode); //Todo implement an insertOrMove method
                    } else if (bVerbose)
                        std::cout << "INSERT ";

                    pOtherNode->setParent(pOutEdge); //BEFORE insert

                    DEBUGONLY(if (pOtherNode->getSLAMscale() != newBadness) {
                        cout << pOtherNode->getSLAMscale() << endl;
                                cout << newBadness << endl;
                    });
                    if(IS_DEBUG) CHECK(pOtherNode->getSLAMscale() != newBadness, "Dijkstra update failed to consistently calculate new scale")

                    nodesWithPaths.insert(pOtherNode);
                }
            }
            if (bVerbose) std::cout << std::endl;
        }
        //Now every node has a badness and parent (except the origin)
        //Every node in nodesWithShortestPaths has a parent and a path back to the root.

        //Now we want to walk this tree to calculate the best position for each real node

        //Each edge has a sequence of predecessor edges eventually leading back to origin
        //Each CNode_relPos has a badness from a particular CEdge_relScale,

        //Each map node has an absolute position from its inbound CNode_relPos with the least badness
        s.stopTimer();

        if (bVerbose) std::cout << "Finished Dijkstra update, " << updateScalesInOrder.size() << " map nodes from " << allNodes.size() << " relative positions have positions, time=" << s.getElapsedTime() << std::endl;
        comp.setClean();
    }
    void CSLAMMap::optimiseScaleAroundLoops(CScaleOptimiser & scaleOptimiser, /*CScaleOptimiser & scaleOptimiser2,*/ const int nComponent, const int t, const bool bVerbose) {
        CComponentData & comp = aComponentData[nComponent];
        checkPositionsUptoDate(comp);

        int nCycles = 0;

        //typedef std::pair<const CNode_relPos *, const CNode_relPos *> TNodePair;

        typedef std::vector<CEdge_relScale *> TAnotherEdgeVec;
        typedef map2<CNode_relPos *, TAnotherEdgeVec > TEdgeMap;
        TEdgeMap edgeSetInSpanner; //, edgeSetNotInSpanner;
        typedef set2<CEdge_relScale *, CEdge_relScale::CSortByBadness > TSortedEdgeSet;
        TSortedEdgeSet sortedEdgeSet;

        TSortedPositions & updateScalesInOrder = comp.getUpdateScalesInOrder();
        TSortedPositions & sortedEdgesNotInMST = comp.getSortedEdgesNotInMST();

        typedef CDynArray<CNode_relPos *> TNodeVec;
        TNodeVec allNodes;
        allNodes.reserve(sortedEdgesNotInMST.size() + updateScalesInOrder.size());
        allNodes.copy_back<TSortedPositions::iterator > (updateScalesInOrder.begin(), updateScalesInOrder.end());
        allNodes.copy_back<TSortedPositions::iterator > (sortedEdgesNotInMST.begin(), sortedEdgesNotInMST.end());

        for (TNodeVec::iterator ppNode = allNodes.begin(); ppNode != allNodes.end(); ppNode++) {
            CNode_relPos * pNode = *ppNode;

            TAnotherEdgeVec & edgeVec = edgeSetInSpanner.initOrGet(pNode);

            if (pNode && pNode->parent()) {
                //pNode->setShortestPathLength(MAX_INT);

                CEdge_relScale * pSourceEdge = pNode->parent();

                edgeVec.push_back(pSourceEdge);

                CNode_relPos * pOtherNode = pSourceEdge->otherNode(pNode);
                TAnotherEdgeVec & otherEdgeVec = edgeSetInSpanner.initOrGet(pOtherNode);
                otherEdgeVec.push_back(pSourceEdge);

                if (bVerbose) std::cout << "Adding edge from " << pNode->id1 << ',' << pNode->id2 << " to " << pOtherNode->id1 << ',' << pOtherNode->id2 << " to map\n";
            }
        }
        for (TNodeVec::iterator ppNode = allNodes.begin(); ppNode != allNodes.end(); ppNode++) {
            CNode_relPos * pNode = *ppNode;

            TAnotherEdgeVec & edgeVec = edgeSetInSpanner[pNode];

            if (pNode && pNode->parent()) {
                CEdge_relScale * pParentEdge = pNode->parent();
                for (CNode_relPos::outEdgeIterator ppOutEdge = pNode->begin(); ppOutEdge != pNode->end(); ppOutEdge++) {
                    CEdge_relScale * pOutEdge = *ppOutEdge;
                    if (bVerbose) std::cout << "Considering edge from " << pOutEdge->firstNode()->id1 << ',' << pOutEdge->firstNode()->id2 << " to " << pOutEdge->secondNode()->id1 << ',' << pOutEdge->secondNode()->id2 << "...";
                    bool bInUse = false;
                    for (TAnotherEdgeVec::const_iterator ppEdge = edgeVec.begin(); ppEdge < edgeVec.end(); ppEdge++)
                        if (*ppEdge == pOutEdge) {
                            //Already in use
                            if (bVerbose) std::cout << "already here\n";
                            bInUse = true;
                            break;
                        }

                    if (!bInUse) {
                        CNode_relPos * pOtherNode = pOutEdge->otherNode(pNode);
                        if (pParentEdge == pOutEdge || pOtherNode->parent() == pOutEdge) {
                            THROW("Node should already have been added")
                                    /*if(bVerbose) std::cout << "sortedEdgeSet add\n";
                                    edgeVec.push_back(pOutEdge);
                                    TAnotherEdgeVec & otherEdgeVec = edgeSetInSpanner.initOrGet(pOtherNode);
                                    otherEdgeVec.push_back(pOutEdge);*/
                        } else {
                            if (bVerbose) std::cout << "sortedEdgeSet add\n";
                            sortedEdgeSet.insert(pOutEdge);
                        }
                    }
                }
            }
        }

        //typedef std::pair<int, CNode_relPos *> TDistNodePair;
        /*class CSortNodesByDist
        {
        public:
                bool operator()(const TDistNodePair & p1, const TDistNodePair & p2) const
                {
                        return (p1.first < p2.first) || (p1.first == p2.first && p1.second < p2.second);
                }
        };*/

        //find SP e1,e2,e3,e4.....
        //TSortedPositions & sortedEdgesNotInMST = comp.getSortedEdgesNotInMST();
        for (TSortedEdgeSet::const_iterator ppEdge = sortedEdgeSet.begin(); ppEdge != sortedEdgeSet.end(); ppEdge++) {
            CEdge_relScale * pEdge = *ppEdge;
            CNode_relPos * pNode1 = pEdge->firstNode();
            CNode_relPos * pNode2 = pEdge->secondNode();

            if (bVerbose) std::cout << "First edge is from " << pNode1->id1 << ',' << pNode1->id2 << " to " << pNode2->id1 << ',' << pNode2->id2 << "...\n";

            //First set all node distances to Inf
            /*for(TEdgeMap::iterator pNodeEdges = edgeSetInSpanner.begin(); pNodeEdges != edgeSetInSpanner.end(); pNodeEdges++)
            {
                    pNodeEdges->first->setShortestPathLengthAndPosSource(MAX_INT, 0);
            }*/
            for (TNodeVec::iterator ppNode = allNodes.begin(); ppNode != allNodes.end(); ppNode++) {
                CNode_relPos * pNode = *ppNode;
                pNode->setShortestPathLengthAndPosSource(MAX_INT, 0);
            }

            //Now find shortest path from pNode1 to pNode2;
            //This time its the same Dijkstra SP as updateSPTree

            typedef set2< CNode_relPos *, CSortNodesByDist > TNodesWithPaths_SortByDist;
            TNodesWithPaths_SortByDist nodesWithPaths;

            CNode_relPos * pRootNode = pNode1;
            pRootNode->setShortestPathLengthAndPosSource(0, pEdge);
            nodesWithPaths.insert(pRootNode);

            int nPathLengthToActiveNode = 0;
            bool bSuccess = false;

            while (nodesWithPaths.size() > 0) {
                CNode_relPos * pActiveNode = nodesWithPaths.pop(); //Node with shortest dist
                if (bVerbose) std::cout << "Active node " << pActiveNode->id1 << ' ' << pActiveNode->id2 << " near ";

                nPathLengthToActiveNode = pActiveNode->shortestPathLength();

                if (pActiveNode == pNode2) {
                    bSuccess = true;
                    break; //we're there
                }

                //Update distances of outgoing links
                TAnotherEdgeVec & edgeVec = edgeSetInSpanner[pActiveNode];

                for (TAnotherEdgeVec::iterator ppOutEdge = edgeVec.begin(); ppOutEdge != edgeVec.end(); ppOutEdge++) {
                    CEdge_relScale * pOutEdge = *ppOutEdge;
                    CNode_relPos * pOtherNode = pOutEdge->otherNode(pActiveNode);
                    if (bVerbose) {
                        std::cout << '(' << pOtherNode->id1 << ',' << pOtherNode->id2 << ": ";
                        std::flush(std::cout);
                    }

                    if (pOtherNode == pNode2 && nPathLengthToActiveNode + 1 <= t) {
                        //SHORT ENOUGH path to ignore this edge
                        bSuccess = true;
                        break; //we're there
                    }

                    if (bVerbose) std::cout << pOtherNode->shortestPathLength() << ") ";
                    if (pOtherNode->shortestPathLength() > nPathLengthToActiveNode + 1) {
                        TNodesWithPaths_SortByDist::iterator pOtherNodeIter = nodesWithPaths.find(pOtherNode);
                        if (pOtherNodeIter != nodesWithPaths.end())
                            nodesWithPaths.erase(pOtherNodeIter);

                        pOtherNode->setShortestPathLengthAndPosSource(nPathLengthToActiveNode + 1, pOutEdge);
                        nodesWithPaths.insert(pOtherNode);
                    }
                }
                if (bVerbose) std::cout << std::endl;
                if (bSuccess)
                    break;
            }

            if (!bSuccess) {
                THROW("Error: No path from node 1 to node 2\n");
            } else {
                //If SP > t add to t-spanner:
                if (nPathLengthToActiveNode > t) {
                    const CNode_relPos * pCycleNode = pNode2;

                    CCycle aCycle;

                    for (;;) {
                        //push cycle to optimiser
                        CEdge_relScale * pCycleEdge = pCycleNode->prevEdgePosSource();

                        double dDir = (pCycleEdge->firstNode() == pCycleNode) ? 1 : -1;

                        aCycle.push_back(std::pair<double, CEdge_relScale *>(dDir, pCycleEdge));
                        pCycleNode = pCycleEdge->otherNode(pCycleNode);

                        if (bVerbose) std::cout << "Cycle edge from " << pCycleEdge->firstNode()->id1 << ',' << pCycleEdge->firstNode()->id2 << " to " << pCycleEdge->secondNode()->id1 << ',' << pCycleEdge->secondNode()->id2 << "..." << dDir << std::endl;

                        if (pCycleNode == pNode1)
                            break;
                    }

                    scaleOptimiser.addCycle(aCycle, bVerbose);
                    //scaleOptimiser2.addCycle(aCycle, bVerbose);

                    edgeSetInSpanner[pEdge->firstNode()].push_back(pEdge);
                    edgeSetInSpanner[pEdge->secondNode()].push_back(pEdge);

                    nCycles++;
                }
            }
        }

        if (nCycles == 0) {
            if (bVerbose) std::cout << "No " << t << "-cycles found\n";
            return;
        }
        if (bVerbose) std::cout << nCycles << " " << t << "-cycles found\n";

        scaleOptimiser.run(bVerbose);
        //scaleOptimiser2.run(true);

        comp.setScaleDirty();
    }

}
