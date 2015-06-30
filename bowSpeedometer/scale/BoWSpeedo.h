/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * BoWSpeedo.h
 *
 *  Created on: 21/08/2009
 *      Author: tom
 */

#ifndef BOWSPEEDO_H_
#define BOWSPEEDO_H_

#include "bow/scales.h"
#include "bow/bagOfWords.h"
#include "bowSpeedoParams.h"

class CBoWSpeedo : public CBoW {

    virtual void resetBeforeRecreateWB() {
        speedo.reset(); //delete existing objects
    }

    virtual void doOR() {
        if (SpeedoPARAMS.NUM_OBJECTS > 0)
            speedo.makeNewObjectDB(SpeedoPARAMS.DRAW_OBJECTS_AFTER_RETRAIN, true); //Find new objects etc.
    }

public:
    const CBOWSpeedoParams & SpeedoPARAMS;

    class CScaleObserver;

    CBoWSpeedo(const CBOWParams & PARAMS, const CBOWSpeedoParams & SpeedoPARAMS, CScaleObserver * pSpeedoScaleObserver) : CBoW(PARAMS), SpeedoPARAMS(SpeedoPARAMS), speedo(*this, pSpeedoScaleObserver) {
    }

    virtual ~CBoWSpeedo() {
    };

    class CBoWObject : boost::noncopyable {
    public:

        class CLengthAndWeight //These are now LOG lengths and weights
        {
            double dLength, dPrecision;
            int nEdgeId;
            CScale SLAMScaleObservedAt;
        public:

            CLengthAndWeight(double dLength, const CScale & SLAMScaleObservedAt, int nEdgeId) : dLength(dLength), dPrecision(-1 /*Todo double check...*/), nEdgeId(nEdgeId), SLAMScaleObservedAt(SLAMScaleObservedAt) {
                if (SLAMScaleObservedAt.hasScale()) {
                    dPrecision = 1.0 / (SLAMScaleObservedAt.getG_sq() + sqr(dLength + SLAMScaleObservedAt.getD()));
                    if(IS_DEBUG) CHECK(!hasSLAMLength(), "Failed to set SLAM length");
                    getPrecision(); //checks val
                }
                if(IS_DEBUG) CHECK(fabs(dLength) >= HUGE, "CLengthAndWeight: bad length/badness observation");
            }

            /*CLengthAndWeight(const CLengthAndWeight & lw) : dLength(lw.dLength), dPrecision(lw.dPrecision), nEdgeId(lw.nEdgeId), SLAMScaleObservedAt(lw.SLAMScaleObservedAt) {}
            void operator=(const CLengthAndWeight & lw) {
                    dLength = (lw.dLength), dPrecision = (lw.dPrecision), nEdgeId = (lw.nEdgeId), SLAMScaleObservedAt = (lw.SLAMScaleObservedAt);
            }*/

            double getBaseline1Length() const {
                if(IS_DEBUG) CHECK(fabs(dLength) >= HUGE, "Invalid length");
                return dLength;
            }

            double getSLAMLength() const {
                if(IS_DEBUG) CHECK(!SLAMScaleObservedAt.hasScale() || fabs(dLength) >= HUGE, "Invalid length");
                return dLength + SLAMScaleObservedAt.getD();
            }

            bool hasSLAMLength() const {
                return SLAMScaleObservedAt.hasScale();
            }

            double getPrecision() const {
                //double dPrecision = 1.0/dVariance;
                CHECKBADNUM(dPrecision);
                if(IS_DEBUG) CHECK(!SLAMScaleObservedAt.hasScale(), "Missing scale");
                //if(dPrecision==0) std::cout << "Warning, zero precision\n";
                return dPrecision;
            }

            double getVar() const {
                if(IS_DEBUG) CHECK(!SLAMScaleObservedAt.hasScale(), "No scale--wasn't observed where a SLAM scale was available");
                if(IS_DEBUG) CHECK(dPrecision <= 0, "Invalid var");
                return 1.0 / dPrecision;
            }

            int getId() const {
                return nEdgeId;
            }
        };

        typedef CDynArray<CLengthAndWeight> TLengths;
    private:
        const CBoWWord * pW1, * pW2;

        TLengths scaleObservations;
        double dMean, dVariance, dSumWeights;
        bool bInit;

    public:

        CBoWObject(const CBoWWord * pW1, const CBoWWord * pW2) : pW1(pW1), pW2(pW2), dMean(0), dVariance(0), dSumWeights(0), bInit(false) {
        }

        void observeOccurance(double dLength, const CScale & SLAMScaleObservedAt, int nEdgeId, bool bVerbose) {
            if(IS_DEBUG) CHECK(scaleObservations.size() > 0 && scaleObservations[scaleObservations.size() - 1].getId() >= nEdgeId, "Error: edges should arrive in strictly increasing order of id");
            if (bVerbose)
                std::cout << "Observed obj " << getObjectId() << " with length " << dLength << " and SLAM scale " << SLAMScaleObservedAt << std::endl;

            scaleObservations.push_back(CLengthAndWeight(dLength, SLAMScaleObservedAt, nEdgeId));

            if (isInit()) {
                if (bVerbose)
                    cout << "Object already initialised; recomputing parameters...\n";
                calculateScale();
            }
        }

        bool isInit() const {
            return bInit;
        }

        double estimatedObjLogSize() const {
            if(IS_DEBUG) CHECK(!bInit, "Not initialised");
            return dMean;
        }

        double estimatedObjLogVar() const {
            if(IS_DEBUG) CHECK(!bInit, "Not initialised");
            return dVariance;
        }

        int count() const {
            return scaleObservations.size();
        }

        void pp() const {
            cout << "Log-lengths: ";
            for (TLengths::const_iterator pLength = scaleObservations.begin(); pLength != scaleObservations.end(); pLength++) {
                if (pLength->hasSLAMLength()) {
                    double dVar = pLength->getVar();
                    std::cout << '[' << pLength->getSLAMLength() << ',' << dVar << "], ";
                } else {
                    std::cout << '[' << pLength->getBaseline1Length() << "], ";
                }
            }
            std::cout << std::endl;

        }

        //Calculate obj size distribution
        void calculateScale();
    private:

        double getRememberedScale(int nEdgeId, int nBegin, int nEnd) const {
            if (nBegin == nEnd) return 0;
            int nMid = nBegin + (nEnd - nBegin) / 2;
            int nMiddleId = scaleObservations[nMid].getId();
            if (nMiddleId > nEdgeId)
                return getRememberedScale(nEdgeId, nBegin, nMid);
            else if (nMiddleId < nEdgeId)
                return getRememberedScale(nEdgeId, nMid + 1, nEnd);
            else
                return scaleObservations[nMid].getBaseline1Length();
        }

    public:

        double rememberScaleObservation(int nEdgeId) const {
            //Binary chop to find if this edge had a scale. 0 on fail.
            return getRememberedScale(nEdgeId, 0, scaleObservations.size());
        }

        const CBoWWord * word1() const {
            return pW1;
        }

        const CBoWWord * word2() const {
            return pW2;
        }

        int getObjectId() const {
            return abs(idFromPtr(this)) % 971;
        }
    };

    //When both words occur in both frames make one of these. Return to BoWSLAM to observe scales.

    class CBoWObjectOccurance : boost::noncopyable {
        CBoWObject * pOb;
        TImPointVec loc1Im1, loc2Im1, //This object's 2 words occur at these locations in Im 1
        loc1Im2, loc2Im2; //and these locations in Im 2

        CLocation l11, l12, l21, l22;
    public:

        CLocation l1Im1() const {
            return l11;
        }

        CLocation l2Im1() const {
            return l21;
        }

        CLocation l1Im2() const {
            return l12;
        }

        CLocation l2Im2() const {
            return l22;
        }

        void setLocs(CLocation lo11, CLocation lo12, CLocation lo21, CLocation lo22) {
            l11 = lo11, l12 = lo12, l21 = lo21, l22 = lo22;
            if(IS_DEBUG) CHECK(l11.zero() || l12.zero() || l21.zero() || l22.zero(), "setLocs: Bad location");
        }

        CBoWObjectOccurance(imageNum nId1, imageNum nId2, CBoWObject * pOb, const CBoWSpeedo & bow)
        : pOb(pOb) {
            bow.get2dObjectFeatures(pOb, nId1, loc1Im1, loc2Im1);
            //should be unneeded: bow.get2dObjectFeatures(pOb, nId2, loc1Im2, loc2Im2);

            //if(IS_DEBUG) CHECK(loc1Im1.size()>0, "Shouldn't be empty") //may be empty if there's more than 8 of each or whatever
            //if(IS_DEBUG) CHECK(loc2Im1.size()>0, "Shouldn't be empty");
            //if(IS_DEBUG) CHECK(loc1Im2.size()>0, "Shouldn't be empty");
            //if(IS_DEBUG) CHECK(loc2Im2.size()>0, "Shouldn't be empty");
        }

        const TImPointVec & getLocOfFeature1InIm1() const {
            return loc1Im1;
        }

        const TImPointVec & getLocOfFeature2InIm1() const {
            return loc2Im1;
        }

        const CBoWObject * object() const {
            return pOb;
        }

        //Points are reconstructed with a baseline 1. EdgeId should be increasing

        void observeOccurance(double dMinLogDist, const CScale & SLAMscale, int nEdgeId) {
            const bool bVerbose = false;
            if (bVerbose)
                std::cout << "Edge: " << nEdgeId;

            if(IS_DEBUG) CHECK(dMinLogDist <= -HUGE || dMinLogDist >= HUGE, "Neg/0 distance--should never occur because never detect 2 features at same point");

            /*Actually move these calculations to observeOccurance
            const double dLogLength = SLAMscale.getD() + dMinLogDist;
            //Something like: const double dLogVar =log(   exp(SLAMscale.getG_sq())*sqr(dMinSqDist);  )
            const double dLogVar = SLAMscale.getG_sq()+sqr(dLogLength); //Todo double check...*/

            pOb->observeOccurance(dMinLogDist, SLAMscale, nEdgeId, bVerbose);
        }
    };

    typedef std::vector<CBoWObjectOccurance*> TObservedObjectsVec;

    class CObservedObjectsVec : public TObservedObjectsVec {
        imageNum nId1, nId2;
    public:

        CObservedObjectsVec(imageNum nId1, imageNum nId2) : nId1(nId1), nId2(nId2) {
        }

        CObservedObjectsVec() : nId1(-1), nId2(-1) {
        }

        void init(imageNum nId1_in, imageNum nId2_in) {
            nId1 = nId1_in, nId2 = nId2_in;
        }

        ~CObservedObjectsVec() {
            for (TObservedObjectsVec::iterator ppOb = begin(); ppOb != end(); ppOb++)
                delete *ppOb;
        }

        void push_back(CBoWObjectOccurance * p) {
            if(IS_DEBUG) CHECK(nId1 < 0 || nId2 < 0, "Not initialised!");
            if(IS_DEBUG) CHECK(!p, "Bad obj occurance");
            TObservedObjectsVec::push_back(p);
        }

        imageNum id1() const {
            if(IS_DEBUG) CHECK(nId1 < 0 || nId2 < 0, "Not initialised!");
            return nId1;
        }

        imageNum id2() const {
            if(IS_DEBUG) CHECK(nId1 < 0 || nId2 < 0, "Not initialised!");
            return nId2;
        }
    };

    typedef CDynArray<CBoWObject *> TObjectsVec;

    class CObjectFinder {

        class CCoOccurance {
            const CBoW::CBoWWord * pW1, * pW2;
            int nCount;
        public:

            CCoOccurance(const CBoW::CBoWWord * pW1, const CBoW::CBoWWord * pW2) : pW1(pW1), pW2(pW2), nCount(0) {
                if(IS_DEBUG) CHECK(pW1 >= pW2, "Cooccurance order constraint broken");
            }

            bool operator<(const CCoOccurance & otherOc) const {
                return pW1 < otherOc.pW1 || (pW1 == otherOc.pW1 && pW2 < otherOc.pW2);
            }

            class CSortCOByWeight {
            public:

                bool operator() (const CCoOccurance & co1, const CCoOccurance & co2) const {
                    return co1.count() > co2.count();
                }
            };

            int count() const {
                return nCount;
            }

            void increment(int nWeight) {
                nCount += nWeight; /* tf-idf or whatever */
                if(IS_DEBUG) CHECK(nWeight <= 0 || nCount <= 0, "Overflow/bad weight");
            }

            CCoOccurance * nonConstPtr() const {
                return const_cast<CCoOccurance *> (this);
            }

            const CBoW::CBoWWord * W1() const {
                return pW1;
            }

            const CBoW::CBoWWord * W2() const {
                return pW2;
            }
        };

        typedef std::set<CCoOccurance> TCoOccurances;
        //typedef std::multiset<CCoOccurance, CCoOccurance::CSortCOByWeight > TBestCoOccurances;
        typedef CDynArray<CCoOccurance> TBestCoOccurances;
        TCoOccurances cooccurances;

    public:

        void addCoOccurance(const CBoWWord * pW1, const CBoWWord * pW2, int nCount1, int nCount2) {
            CCoOccurance co(pW1, pW2);
            TCoOccurances::iterator pCO = cooccurances.find(co);

            if (pCO == cooccurances.end())
                pCO = cooccurances.insert(co).first;

            //const_cast because modifying item in a set isn't allowed
            //pCO->nonConstPtr()->increment(pWord1->FrequencyInImage()*pWord1->Weight() + pWord2->FrequencyInImage()*pWord1->Weight()); //Sum of TF-IDF (dubious)
            //pCO->nonConstPtr()->increment(10000/(pWord1->FrequencyInImage() * pWord2->FrequencyInImage()));
            //Was div-by-DF last time... pCO->nonConstPtr()->increment(1+(10000*nCount1*nCount2*(pW1->TotalOccurances()*pW2->TotalOccurances()))); //Product of TF-IDF (this is same as LSA when freq1*weight1 = freq2*weight2 (if this is always true then the weights are equal so frequencies are equal too)
            pCO->nonConstPtr()->increment(1);
        }

        void countCoOccurances(const CBoW::CBoWWordBag::TWordBag & WB) {
            for (CBoW::CBoWWordBag::TWordBag::const_iterator pWord1 = WB.begin(); pWord1 != WB.end(); pWord1++) {
                for (CBoW::CBoWWordBag::TWordBag::const_iterator pWord2 = pWord1 + 1; pWord2 != WB.end(); pWord2++) {
                    addCoOccurance(pWord1->Word(), pWord2->Word(), pWord1->FrequencyInImage(), pWord2->FrequencyInImage());
                }
            }
        }

        void reset() {
            cooccurances.clear();
        }

        void getObjects(const CBOWSpeedoParams & SpeedoPARAMS, TObjectsVec & objects) {
            TBestCoOccurances bestCOs;
            bestCOs.reserve(cooccurances.size());
            for (TCoOccurances::const_iterator pCOs = cooccurances.begin(); pCOs != cooccurances.end(); pCOs++) {
                bestCOs.push_back(*pCOs);
            }
            int nNumObjects = std::min<int>(SpeedoPARAMS.NUM_OBJECTS, (int) ((double) bestCOs.size() * (double)SpeedoPARAMS.MAX_PROP_COOCCURANCES)); //Assume 99.99% of co-occurances rubbish

            if (nNumObjects > 0) {
                std::partial_sort(bestCOs.begin(), bestCOs.begin() + nNumObjects, bestCOs.end(), CCoOccurance::CSortCOByWeight());
                std::cout << "Best and worst scores selected: " << bestCOs.begin()->count() << ',' << (bestCOs.begin() + nNumObjects - 1)->count() << std::endl;
                for (TBestCoOccurances::iterator pCO = bestCOs.begin(); pCO != bestCOs.begin() + nNumObjects; pCO++)
                    if (pCO->count() > 0)
                        objects.push_back(new CBoWObject(pCO->W1(), pCO->W2()));

                std::cout << nNumObjects << " objects found\n";
            } else
                std::cout << "No objects found yet (" << bestCOs.size() << " cooccurances)" << std::endl;
        }
    };

    class CBoWSpeedometer;

    class CScaleObserver {
    public:
        virtual void reobserveScales(bool bDrawObjects) = 0;
        virtual void addEdge(CEdge * pEdge) = 0;
        virtual void removeEdge(CEdge *) = 0;
        virtual void init(CBoWSpeedo::CBoWSpeedometer * speedo) = 0;

        virtual ~CScaleObserver() {
        };
        virtual void findObjects(const CBoWSpeedo & bow, CObjectFinder & objectFinder) = 0;
    };

    class CNullScaleObserver : public CScaleObserver {
    public:

        virtual void reobserveScales(bool) {
        }

        virtual void addEdge(CEdge *) {
        }

        virtual void removeEdge(CEdge *) {
        }

        virtual void init(CBoWSpeedo::CBoWSpeedometer *) {
        }

        virtual void findObjects(const CBoWSpeedo &, CObjectFinder &) {
        };
    };

    class CObjAndLocations {
        friend class CEdgeScaleObserver;
        const CBoW::CBoWWord * pW1;
        CDynArray<CLocation> aLocations;
        int nReconCount;
    public:

        CObjAndLocations(const CBoW::CBoWWord * pW1) : pW1(pW1), nReconCount(0) {
        }

        CObjAndLocations(const CObjAndLocations & ol) : pW1(ol.pW1), nReconCount(ol.nReconCount) {
            ol.aLocations.copyInto(aLocations);
            //if(IS_DEBUG) CHECK(ol.aLocations.size(), "Should copy nonempty vector");
        }

        void setReconCounts(int n1) {
            nReconCount = n1;
        }

        bool isReconstructed() const {
            return nReconCount > 0;
        }

        const CDynArray<CLocation> & locations() const {
            return aLocations;
        }

        void addLoc(CLocation l) {
            aLocations.push_back(l);
        }
    };

    typedef CDynArrayOwner<CObjAndLocations> TObjectsAndLocations;

    class CBoWSpeedometer //Also needs a poly. class to tell it about edges...
    {
        TObjectsVec objects;
        const CBoWSpeedo & bow;
        CScaleObserver * pScaleObserver;
    public:

        CBoWSpeedometer(const CBoWSpeedo & bow, CScaleObserver * pScaleObserver) : bow(bow), pScaleObserver(pScaleObserver) {
            if (!pScaleObserver) pScaleObserver = new CNullScaleObserver;
            pScaleObserver->init(this);
        }

        void makeNewObjectDB(bool bDrawObjects, bool bUseEdges) {
            //iterate over images counting co-occurances
            //objectFinder.reset();
            CObjectFinder objectFinder;

            if (bUseEdges) {
                pScaleObserver->findObjects(bow, objectFinder);
            } else {
                for (CBoW::constImIt ppIm = bow.vImages.begin(); ppIm != bow.vImages.end(); ppIm++) {
                    const CBoWWordBag * pIm = *ppIm;
                    const CBoW::CBoWWordBag::TWordBag & aWordsBottomLevel = pIm->bottomLevelWords();
                    objectFinder.countCoOccurances(aWordsBottomLevel);
                }
            }

            objectFinder.getObjects(bow.SpeedoPARAMS, objects);

            pScaleObserver->reobserveScales(bDrawObjects); //give them a size distn.
        }

        void updateObjectSizeDistns() {
            for (TObjectsVec::iterator pOb = objects.begin(); pOb != objects.end(); pOb++)
                (*pOb)->calculateScale();
        }

        virtual ~CBoWSpeedometer() {
            reset();
            delete pScaleObserver;
        }

        virtual void reset() {
            cout << "reset: Deleting " << objects.size() << " objects\n";
            for (TObjectsVec::iterator ppOb = objects.begin(); ppOb != objects.end(); ppOb++)
                delete *ppOb;

            objects.clear();
        }

        virtual void getObjectsInBothFrames(CObservedObjectsVec & objectsObserved) {
            imageNum nId1 = objectsObserved.id1();
            imageNum nId2 = objectsObserved.id2();

            for (TObjectsVec::iterator pOb = objects.begin(); pOb != objects.end(); pOb++) {
                try
                {
                    if (bow.vImages[nId1]->containsWordLocations((*pOb)->word1()) && bow.vImages[nId1]->containsWordLocations((*pOb)->word2()) && bow.vImages[nId2]->containsWordLocations((*pOb)->word1()) && bow.vImages[nId2]->containsWordLocations((*pOb)->word2())) {
                        CBoWObjectOccurance * pObOccurance = new CBoWObjectOccurance(nId1, nId2, *pOb, bow);
                        objectsObserved.push_back(pObOccurance);
                    }
                } catch(const char * szErr)
                {
                    cout << "Error finding object occurrances: " << szErr << endl;
                    cout << "ID1: " << nId1 << endl;
                    cout << "ID2: " << nId2 << endl;
                }
            }
        }

        void addSpeedoEdge(CEdge * pEdge) {
            pScaleObserver->addEdge(pEdge);
        }

        void removeSpeedoEdge(CEdge * pEdge) {
            pScaleObserver->addEdge(pEdge);
        }
    };

    void addSpeedoEdge(CEdge * pEdge) {
        speedo.addSpeedoEdge(pEdge);
    }

    void removeSpeedoEdge(CEdge * pEdge) //Not actually needed right now
    {
        speedo.removeSpeedoEdge(pEdge);
    }

    void get2dObjectFeatures(const CBoWObject * pOb, const imageNum nId, TImPointVec & loc1Im1, TImPointVec &loc2Im1) const {
        vImages[nId]->getWordLocations(pOb->word1(), loc1Im1);
        vImages[nId]->getWordLocations(pOb->word2(), loc2Im1);
    }

    void getObjectsAndLocationsInIntersection(TObjectsAndLocations & localObjectFinder, const imageNum nId1, const imageNum nId2) const {
        const CBoWWordBag::TWordBag & pWB1 = vImages[nId1]->bottomLevelWords();
        const CBoWWordBag::TWordBag & pWB2 = vImages[nId2]->bottomLevelWords();
        CBoWWordBag::TWordBag::const_iterator pWord1 = pWB1.begin();
        CBoWWordBag::TWordBag::const_iterator pWord2 = pWB2.begin();
        const CBoWWordBag::TWordBag::const_iterator pEnd1 = pWB1.end();
        const CBoWWordBag::TWordBag::const_iterator pEnd2 = pWB2.end();

        while (pWord1 != pEnd1 && pWord2 != pEnd2) {
            if(IS_DEBUG) CHECK(!pWord1 || !pWord1, "Uninitialised word in bottomLevelWords array")
            while (pWord2 != pEnd2 && pWord2->Word() < pWord1->Word()) {
                pWord2++;
            }

            if (pWord2 != pEnd2 && pWord2->Word() == pWord1->Word()) {
                if (!pWord1->Location(0).zero() && !pWord2->Location(0).zero()) {
                    CObjAndLocations * pObjLoc = new CObjAndLocations(pWord1->Word());
                    for (int i = 0; i < pWord1->FrequencyInImage(); i++)
                        pObjLoc->addLoc(pWord1->Location(i)); //Only need to know about locs in im 1

                    localObjectFinder.push_back(pObjLoc);
                }
            }

            pWord1++;
        }
    }
private:
    CBoWSpeedometer speedo;
};

#endif /* BOWSPEEDO_H_ */
