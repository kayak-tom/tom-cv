/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
    Use proximity constraint to seed clusters.

    Add new images to appropriate clusters as needed. Could cluster by ins data in first instance.

    Now do we match new images to clusters, or to images (but only in some clusters)?
 */

#ifndef _BoW
#define _BoW

#include "bagOfWordsParam.h"

#define templateCompMethod template<CBOWParams::eCOMP_METHOD eCompMethod>

#include "util/dynArray.h"

/*#if CLUSTER_THREADS>0
#define MT
#endif
#if RWB_THREADS>0
#define MT
#endif
#if QUERY_THREADS>0
#define MT
#endif
#ifdef CLUSTER_IN_SEPERATE_THREAD
#define MT
#define MX_NO_DELETE_WHILE_CLUSTER_IN_SEPERATE_THREAD //boost::mutex::scoped_lock scoped_lock(mxClusteringCantDelete);
#else
#define MX_NO_DELETE_WHILE_CLUSTER_IN_SEPERATE_THREAD
#endif*/
#define MX_NO_DELETE_WHILE_CLUSTER_IN_SEPERATE_THREAD boost::mutex::scoped_lock scoped_lock(mxClusteringCantDelete);

#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>
#include <boost/thread/barrier.hpp>
#include <boost/thread/shared_mutex.hpp>
#include <boost/pool/singleton_pool.hpp>
#include <boost/pool/pool_alloc.hpp>

#include <functional>
#include <set>
#include <map>
#include "util/smallHashTable.h"
#include "util/exception.h"
#include "description/descriptor.h"

#define SCOPED_WRITE_LOCK checkClusteringException(); boost::unique_lock<boost::shared_mutex> write_lock(mxClusteringRW);
#define WRITE_LOCK SCOPED_WRITE_LOCK
#define RELEASE_WRITE_LOCK
//#define WRITE_LOCK checkClusteringException(); pthread_mutex_init(&mxClusteringRW_pt, NULL);
//#define RELEASE_WRITE_LOCK pthread_mutex_destroy(&mxClusteringRW_pt);
/*
#define READ_LOCK checkClusteringException(); boost::shared_lock<boost::shared_mutex> reader_lock(mxClusteringRW);
#define UPGRADEABLE_LOCK checkClusteringException(); boost::upgrade_lock<boost::shared_mutex> upg_lock(mxClusteringRW);
#define UPGRADE_LOCK checkClusteringException(); boost::unique_lock<boost::shared_mutex> write_lock(boost::move(upg_lock));
 */

// For Brute force:
#define READ_LOCK WRITE_LOCK
#define UPGRADEABLE_LOCK WRITE_LOCK
#define UPGRADE_LOCK

#define LOCK_WHILE_COUNTING_OCCURANCES boost::mutex::scoped_lock scoped_lock(mxOccuranceCounting);
//************************************************//

typedef unsigned int wordNum; //Index words with an integer
typedef unsigned int imageNum; //Index vImages with an integer

#define DEFAULT_RECLUSTER_FREQUENCY 1.25 //Re-cluster when num of features increases by this factor
#define INT_LOG_CACHE_LIM 1000 //Cache int logs up to this limit
#define INTINT_LOG_CACHE_LIM 10250 //Cache intint logs up to this limit
#define INTINTLOG_SCALE 100 // higher values give more precision, but too big will cause overflow
#define NISTERVECTOR_LENGTH 10000 //should be >> INTINTLOG_SCALE*num words in each image, but too big will cause overflow

#include "BoW.h"

class CCluster;
class CClusterSet;

//Defined elsewhere, for debugging clustering:

class CClusterDrawer {
public:
    virtual void drawCS(const CClusterSet * pCS, int nLevel) = 0;
};

class CEdge;

//Consists of a dictionary of 'words' and a DB of vImages ((sparse??) vectors of word frequencies)

class CBoW {
    boost::shared_mutex mxClusteringRW; //Lock pretty much everything while we rebuild word bags after clustering
    boost::mutex mxOccuranceCounting; //Lock while incrementing word occurance counts
    boost::mutex mxClusteringCantDelete;

protected:
    friend class CBoWWordBag;

    const CBOWParams & PARAMS;

    class CBoWWord {
        const CDescriptor * pDescriptor;
        int nTotalFrequency,
        nTotalOccurances; // how many images does it occur in

        bool bIOwnMyDescriptor;

    public:
        CBoWWord(CCluster * pCluster);
        /*virtual*/ ~CBoWWord();

        inline const CDescriptor * Descriptor() const {
            return pDescriptor;
        };

        inline int TotalFrequency() const {
            return nTotalFrequency;
        };

        inline int TotalOccurances() const {
            return nTotalOccurances;
        };

        inline void IncrementOccurances() {
            nTotalOccurances++;
        };

        inline void DecrementOccurances() {
            nTotalOccurances--;
            if(IS_DEBUG) CHECK(nTotalOccurances < 0, "DecrementOccurances fell below zero");
        };

        inline void DecrementOccurances(int i) {
            nTotalOccurances -= i;
            if(IS_DEBUG) CHECK(nTotalOccurances < 0, "DecrementOccurances fell below zero");
        };

        static bool sortWordsByDescPtr(const CBoWWord * pWord1, const CBoWWord * pWord2) {
            return pWord1->Descriptor() < pWord2->Descriptor();
        }
    };
    static CBoW::CBoWWord * findBinaryChop(const CDescriptor * pDesc, CBoW::CBoWWord * const* apWords, const int nLength);

public:

    class CBoWDictionary {
        CDynArray<CBoWWord *> Words;

        int nLevel, //< 1=bottom level
        nFirstWordIdx; //< Words in this dictionary are indexed starting from this value

        int * anWordLevelCounts; // number of words at all levels below *this* dictionary.
        // Use this from the top dictionary to count words in each level
        // Also use to construct offsets

        CClusterSet * pClusters;

        CException pDictionaryException;

        void checkDictionaryException() {
            if (pDictionaryException) {
                throw pDictionaryException;
            }
        };

        CBoWWord * ClosestWord(const CDescriptor * pDescriptor, CDescriptor::TDist nRadius) const;
    public:

        class CGetClusterNum {
            const int nBranchFactor, nTotalDescriptors, nImages;
        public:

            CGetClusterNum(int nBranchFactor, int nTotalDescriptors, int nImages) : nBranchFactor(nBranchFactor), nTotalDescriptors(nTotalDescriptors), nImages(nImages) {
                if(IS_DEBUG) CHECK(nImages < 1, "Clustering before we have images");
                if(IS_DEBUG) CHECK(nTotalDescriptors < 1, "Clustering before we have descriptors");
                if(IS_DEBUG) CHECK(nBranchFactor < 1, "Bad branch factor");
                /*if(nBranchFactor == 1)
                        nBranchFactor=2;*/
            }

            int numImages() const {
                return nImages;
            }

            int getClusterNum(const int nLevel, const int nDescriptorCount, const CBOWParams::CBOWClusteringParams & BOWCLUSTERPARAMS) const;
        };
    private:

        void ClusterToWords_recurse(CBOWParams::CDescriptorBinningParams const &, const CGetClusterNum & getClusterCount, CCluster *pCluster, const CBOWParams::CBOWClusteringParams & BOWCLUSTERPARAMS, CClusterDrawer ** ppDrawCS);
        void ClusterRangeToWords_recurse(const CBOWParams::CDescriptorBinningParams &, CClusterSet *pClusters, const CGetClusterNum & getClusterCount, int nClusterStart, int nClusterEnd, const CBOWParams::CBOWClusteringParams & BOWCLUSTERPARAMS, CClusterDrawer ** ppDrawCS);
        void ClusterRangeToWords_recurse_int(CBOWParams::CDescriptorBinningParams const &, CClusterSet *pClusters, const CGetClusterNum & getClusterCount, int nClusterStart, int nClusterEnd);
        //int ClusterNum(int nDescriptors, const CGetClusterNum & getClusterCount, const CBOWParams::CBOWClusteringParams & BOWCLUSTERPARAMS);
    public:
        CBoWDictionary(const CCluster * pCluster, //0 at the top
                CDescriptorSet * pDescriptors, //<All the points from which to build this dictionary
                const CGetClusterNum &,
                int nLevel_in,
                const int * anLevelWordCounts_in,
                const CBOWParams::CDescriptorBinningParams & pBinningParams, const CBOWParams::CBOWClusteringParams & BOWCLUSTERPARAMS, CClusterDrawer ** ppDrawCS);

        void LookupWordAllLevels(const CDescriptor * pDescriptor, int nLevel, CBoWWord ** apWords, const CDescriptor ** apClosestDescriptors, const CBOWParams::CDescriptorBinningParams & pBinningParams) const; //Look up word in all levels at once--for adding to dictionary

        ~CBoWDictionary();

        inline int Level() const {
            return nLevel;
        };

        inline const int * WordCountArray() const {
            return anWordLevelCounts;
        };
    };
    typedef CDynArray<CLocation> TImPointVec;

protected:

    class CBoWNodeWord : public CBoWWord {
        CBoWDictionary * pSubDictionary; //This 'word' is a big cluster of descriptors. This dictionary clusters this cluster.
    public:
        CBoWNodeWord(CCluster * pCluster, const CBoW::CBoWDictionary::CGetClusterNum & getClusterCount, int nLevel, const int * anLevelWordCounts_in, const CBOWParams::CDescriptorBinningParams & pBinningParams, const CBOWParams::CBOWClusteringParams & BOWCLUSTERPARAMS, CClusterDrawer ** ppDrawCS);

        ~CBoWNodeWord() {
            delete pSubDictionary;
            ((CBoWWord *)this)->~CBoWWord();
            //CBoW::CBoWWord::~CBoWWord();
        };

        inline const CBoWDictionary * SubDictionary() const {
            return pSubDictionary;
        };
    };

    //Rev. 125: Moved to tree structure for word bag. Leave word-counting in dictionary (anWordLevelCounts) as need per-level word counts for bayes classifier
public:
    class CBoWSpeedo;

protected:

#ifndef __GNUC__
#pragma warning(disable: 4522) //Multiple assignment operators
#endif

    class CBoWWordBag //Todo--replace pParent with CBoWClusteringParams &?
    {
        //    	friend class CBoWSpeedo::CBoWSpeedometer;

        //Represents an image as a bag of words! (at 1 level only)
        //Consists of an array of sorted vectors of CBoWImageWord objects, representing
        friend void LoadSettingsFromCfg(const char * szFN);

        class CBoWImageWord //essentially a struct
        {
            CBoWWord * pWord; //not const cos may change occurance count

            int nFrequency;
            int nWeight;

            union {
                int aLoc[sizeof (CLocation *) / sizeof (int) ];
                CLocation * pLoc;

                CLocation loc() const {
                    return aLoc[0];
                }

                void operator=(CLocation l) {
                    aLoc[0] = l.id();
                }
            } aLocations;
            //CBoWImageWord(const CBoW::CBoWWordBag::CBoWImageWord &word); //Stop copying without properly moving pointer
        public:
            //DEBUGONLY(static int COPY_OK;);

            CBoWImageWord() : pWord(0), nFrequency(0), nWeight(0) {
                aLocations.pLoc = 0;
                //std::cout << "Calling CBoWImageWord()\n";
            };

            CBoWImageWord(const CBoW::CBoWWordBag::CBoWImageWord &word) //must support copying in words with lists of locations
            {
                //if(IS_DEBUG) CHECK(!COPY_OK && nFrequency>1 && word.Location(0).CLocationPtr != 0, "CBoW::CBoWWordBag::CBoWImageWord: Copy constructor should not be used");
                //DEBUGONLY(COPY_OK = 0;);
                if(IS_DEBUG) CHECK(word.FrequencyInImage() < 1, "CBoW::CBoWWordBag::CBoWImageWord: Copy should only be used with real words with locations");
                pWord = word.Word();
                nFrequency = word.FrequencyInImage(); //=1
                //std::cout << "Calling CBoWImageWord(CBoWImageWord &)\n";

                //Acquire memory:
                aLocations = word.aLocations;
                const_cast<CBoW::CBoWWordBag::CBoWImageWord &> (word).aLocations.pLoc = 0;

                nWeight = word.Weight();
            }

            void operator=(CBoW::CBoWWordBag::CBoWImageWord &word) //Stop copying without properly moving pointer
            {
                if(IS_DEBUG) CHECK(word.FrequencyInImage() < 1, "CBoW::CBoWWordBag::CBoWImageWord: Copy should only be used with real words with locations");
                //std::cout << "Calling operator=\n";
                pWord = word.Word();
                nFrequency = word.FrequencyInImage(); //=1
                aLocations = word.aLocations;
                word.aLocations.pLoc = 0;
                nWeight = word.Weight();
            }

            void operator=(const CBoW::CBoWWordBag::CBoWImageWord &word) {
                //std::cout << "Calling operator=(const)\n";
                operator=(const_cast<CBoW::CBoWWordBag::CBoWImageWord &> (word)); //Because STL insert requires const
            }

            CBoWImageWord(CBoW::CBoWWord * pWord_in, const CLocation & Loc);
            void IncrementFrequency(const CLocation & Loc);
            ~CBoWImageWord();
            void DestroyLoc();

            inline void SetWeight(int nWeight_in) {
                if(IS_DEBUG) CHECK(nWeight_in < 0, "Weight less than zero, overflow?")
                nWeight = nWeight_in;
            }

            inline CBoWWord * Word() {
                return pWord;
            }; //not const cos may change occurance count

            inline CBoWWord * Word() const {
                return pWord;
            }; // const one too

            inline int FrequencyInImage() const {
                return nFrequency;
            };

            inline int Weight() const {
                return nWeight;
            };

            inline const CLocation Location(int n) const {
                if (aLocations.loc().zero()) return aLocations.loc(); //equiv return 0

                if (nFrequency == 1) {
                    //if(IS_DEBUG) CHECK(!((aLocations.x() > 0) && (aLocations.x() < 1024) && (aLocations.y() > 0) && (aLocations.y() < 768)), "Dodgy location (uninit?)");

                    if(IS_DEBUG) CHECK(n > 0, "CBoW::CBoWWordBag::CBoWLocation: n OOB");
                    return aLocations.loc();
                }
                //else
                if(IS_DEBUG) CHECK(n >= LOCATION_STORE_LIM, "CBoW::CBoWWordBag::CBoWLocation: n OOB");
                return aLocations.pLoc[n];
            };
        };

        class WordPointerSort : std::binary_function<const CBoWImageWord &, const CBoWImageWord &, bool> {
        public:

            inline result_type operator() (first_argument_type a,
                    second_argument_type b) const {
                return (result_type) (a.Word() < b.Word());
            }
        };


        //static TWordSet * GetWordSet(int nLevels);
        //static void DeleteWordSet(TWordSet * apWordSet, int nLevels);

        const CBoW * pParent;
    public:
        typedef CDynArray<CBoWImageWord> TWordBag;
    private:
        typedef CBoWImageWord * TImWordIt;
        TWordBag * aWords;
        CDescriptorSet * pDescriptors; //does own the memory

        class CBoWWordDescMatch {
            const CBoWWord * pImWord;
            int nDescIdx;
        public:

            const CBoWWord * Word() const {
                return pImWord;
            }

            int DescriptorSetIdx() const {
                return nDescIdx;
            }

            CBoWWordDescMatch(const CBoWWord * pImWord, const int nDescIdx) : pImWord(pImWord), nDescIdx(nDescIdx) {
            }
            //CBoWWordDescMatch() : pImWord(0), nDescIdx(0) {} //for array initialisation

            void operator=(const CBoWWordDescMatch & m) {
                pImWord = m.Word();
                nDescIdx = m.DescriptorSetIdx();
            } //for array initialisation
            //and default copy c'tor is fine

            class CBoWWordDescMatchSort : std::binary_function<const CBoWWordDescMatch &, const CBoWWordDescMatch &, bool> {
            public:

                inline result_type operator() (first_argument_type a,
                        second_argument_type b) const {
                    return (result_type) (a.Word() < b.Word());
                }
            };

        };
        CDynArray<CBoWWordDescMatch> aDescriptorWordsBF; //For smart BF matching

        templateCompMethod
        int Compare(TWordBag *, TWordBag *, int nWordsInCat) const;

        templateCompMethod
        int CompareWithAll(TWordBag * pWordsInImage, int nWords) const;

        boost::mutex mxQuery;

        CException pWBException;

        void checkWBException() {
            if (pWBException) {
                CException pWBExceptionTemp = pWBException;
                throw pWBExceptionTemp;
            }
        };

        templateCompMethod
        void getBoWMatches_Loop(TBoWMatchVector *pvMatches, int nReturnMax, const int * anScoreAgainstBackground, CDynArray<CBoWWordBag *>::const_iterator imageIterBegin, CDynArray<CBoWWordBag *>::const_iterator imageIterEnd) HOT;

        templateCompMethod
        void getBoWMatches_Loop_int(TBoWMatchVector *pvMatches, int nReturnMax, const int * anScoreAgainstBackground, CDynArray<CBoWWordBag *>::const_iterator imageIterBegin, CDynArray<CBoWWordBag *>::const_iterator imageIterEnd) HOT;

        typedef std::multiset<int, std::less<int> > TSortedIntSet;
        inline int VectorCompare(const CBoWWordBag * pWB, int nLevel) const;
        inline int VectorCompare_int(TWordBag::iterator catWordList, TWordBag::iterator catWordListEnd, TWordBag::iterator imageWordList, TWordBag::iterator imageWordListEnd) const;

        template<class TImWordIterator>
        inline int NormalisedVectorCompare_int(TImWordIterator catWordList, const TImWordIterator catWordListEnd, TImWordIterator imageWordList, const TImWordIterator imageWordListEnd) const;
#define fnVectorCompare NormalisedVectorCompare_int
        inline int VectorCompareFast(const CBoWWordBag * pWB, int nLevel, TSortedIntSet * pSubvecScores, int nVecReturnMax) const;

        const int nId;

        //static CBoWCorrGroup * new_cBoWCorrGroup(TWordBag::iterator word1, TWordBag::iterator word2);
        //static int BF_LEVEL; //What level do we use for smart BF matching?
        //static int BF_NN = 3; //return 2-2 and 3-3 corr's

        class CLocationSort : std::binary_function<const CLocation &, const CLocation &, bool> {
        public:

            inline result_type operator() (first_argument_type a,
                    second_argument_type b) const {
                return (result_type) (a < b);
            }
        };

        //typedef set<CLocation, CLocationSort> TLocationSet;

        typedef std::pair<CLocation, CDescriptor::TDist> TLocDistPair;
        typedef boost::details::pool::null_mutex boost_pool_mutex;
        //typedef boost::details::pool::pthread_mutex boost_pool_mutex;
        typedef CSmallHashTable<CLocation, 64, 1, locationHash < 64 > > TSimpleLocSet;

        class CLocDistPairSort : std::binary_function<const TLocDistPair &, const TLocDistPair &, bool> {
        public:

            inline result_type operator() (first_argument_type a,
                    second_argument_type b) const {
                return (result_type) (a.second < b.second);
            }
        };

        class CMatchMap {
            //friend void LoadSettingsFromCfg(const char * szFN);
            double dMaxPrior;
            const double dMinPrior;

            TLocDistPair * aaLDPairsLR; //row of R matches to 1st L match, row of R matches to 2nd L match, etc
            TLocDistPair * aaLDPairsRL; //row of L matches to 1st R match, row of L matches to 2nd R match, etc
            int nLocsLeft, nLocsRight, nMaxSize;
            double dClosest;
        public:

            CMatchMap(const CBOWParams::CBOWCorrespondenceProbParams & PROB_PARAMS) : dMaxPrior(PROB_PARAMS.MAX_PRIOR), dMinPrior(PROB_PARAMS.MIN_PRIOR), aaLDPairsLR(0), aaLDPairsRL(0), nLocsLeft(0), nLocsRight(0), nMaxSize(0), dClosest(MAX_INT) {
                if (dMaxPrior < dMinPrior) dMaxPrior = dMinPrior;
            }

            ~CMatchMap() {
                delete [] aaLDPairsLR;
                delete [] aaLDPairsRL;
            }

            void init(int nLL, int nLR) {
                nLocsLeft = nLL;
                nLocsRight = nLR;

                if (nMaxSize < nLocsLeft * nLocsRight) {
                    nMaxSize = nLocsLeft * nLocsRight;
                    delete [] aaLDPairsLR;
                    delete [] aaLDPairsRL;
                    aaLDPairsLR = new TLocDistPair[nMaxSize];
                    aaLDPairsRL = new TLocDistPair[nMaxSize];
                }
            }

            void add(int nLL, int nLR, CLocation locL, CLocation locR, CDescriptor::TDist dist) {
                aaLDPairsLR[nLL * nLocsRight + nLR] = TLocDistPair(locR, dist);
                aaLDPairsRL[nLR * nLocsLeft + nLL] = TLocDistPair(locL, dist);
            }

            void sort() {
                //Todo: prob only need to sort top 9
                CDescriptor::TDist nClosest = MAX_INT;
                for (int i = 0; i < nLocsLeft; i++) {
                    std::sort(aaLDPairsLR + (i * nLocsRight), aaLDPairsLR + ((i + 1) * nLocsRight), CLocDistPairSort());
                    if (aaLDPairsLR[i * nLocsRight].second < nClosest) nClosest = aaLDPairsLR[i * nLocsRight].second;
                }
                dClosest = nClosest;

                for (int i = 0; i < nLocsRight; i++)
                    std::sort(aaLDPairsRL + (i * nLocsLeft), aaLDPairsRL + ((i + 1) * nLocsLeft), CLocDistPairSort());
            }

            double getPP(int nLeftIdx) {
                //Prior prob is prop. to ratio of match strength to best match strength (in this cluster)
                double dScale = 1.0;
                if (dClosest > 0) {
                    CDescriptor::TDist nThisBestDist = aaLDPairsLR[nLeftIdx * nLocsRight].second;
                    if (nThisBestDist > 0)
                        dScale = dClosest / (double) nThisBestDist;
                }
                double dNewProb = dMinPrior + (dMaxPrior - dMinPrior) * dScale;
                if(IS_DEBUG) CHECK(dNewProb <= 0 || dNewProb >= 1, "Bad prob update");
                return dNewProb;
            }

            void getLRMatchesWC(TSimpleLocSet & locs, int idx, const double dCondition) {
                CDescriptor::TDist nLastScore = MAX_INT;
                TLocDistPair * pPair = aaLDPairsLR + (idx * nLocsRight);
                for (int i = nLocsRight; i > 0; i--) {
                    CDescriptor::TDist nThisScore = pPair->second;

                    if (nLastScore < nThisScore * dCondition)
                        break;

                    locs.insert(pPair->first);

                    pPair++;
                    nLastScore = nThisScore;
                }
            }

            void getRLMatchesWC(TSimpleLocSet & locs, int idx, const double dCondition) {
                CDescriptor::TDist nLastScore = MAX_INT;
                TLocDistPair * pPair = aaLDPairsRL + (idx * nLocsLeft);
                for (int i = nLocsLeft; i > 0; i--) {
                    CDescriptor::TDist nThisScore = pPair->second;

                    if (nLastScore < nThisScore * dCondition)
                        break;

                    locs.insert(pPair->first);

                    pPair++;
                    nLastScore = nThisScore;
                }
            }
        };

        // ## ONLY USE fast_pool_allocator here to avoid sharing pool with BaySAC ####
        typedef std::set<TLocDistPair, CLocDistPairSort, boost::fast_pool_allocator<TLocDistPair, boost::default_user_allocator_new_delete, boost_pool_mutex > > TMatchSet;
        typedef std::map<CLocation, TMatchSet, std::less<CLocation>, boost::fast_pool_allocator<std::pair<CLocation, TMatchSet>, boost::default_user_allocator_new_delete, boost_pool_mutex > > TMatchMap;
        typedef std::map<CLocation, TSimpleLocSet, std::less<CLocation>, boost::fast_pool_allocator<std::pair<CLocation, TSimpleLocSet>, boost::default_user_allocator_new_delete, boost_pool_mutex > > TSimpleMatchMap;

        static void BFCorrFromRange(const CDescriptorSet * pDS1, const CDescriptorSet * pDS2,
                const CDynArray<CBoWWordDescMatch>::const_iterator WB1WordList,
                const CDynArray<CBoWWordDescMatch>::const_iterator WB1WordListRangeEnd,
                const CDynArray<CBoWWordDescMatch>::const_iterator WB2WordList,
                const CDynArray<CBoWWordDescMatch>::const_iterator WB2WordListRangeEnd,
                double dCondition, CBoWCorrespondences * pCorrespondences);
        static void BFCorrFromRangeNN(const CDescriptorSet * pDS1, const CDescriptorSet * pDS2,
                const CDynArray<CBoWWordDescMatch>::const_iterator WB1WordList,
                const CDynArray<CBoWWordDescMatch>::const_iterator WB1WordListRangeEnd,
                const CDynArray<CBoWWordDescMatch>::const_iterator WB2WordList,
                const CDynArray<CBoWWordDescMatch>::const_iterator WB2WordListRangeEnd,
                double dCondition, const int BF_NN, CBoWCorrespondences * pCorrespondences, const CBOWParams::CBOWCorrespondenceProbParams & CORR_PROB_PARAMS);
        static void addMatchesToLeftPoint(const TSimpleMatchMap & aStrongestMatchesLeft, const TSimpleMatchMap & aStrongestMatchesRight, const CLocation locLeft, TSimpleLocSet & matchesLeft, TSimpleLocSet & matchesRight);
        static void addMatchesToRightPoint(const TSimpleMatchMap & aStrongestMatchesLeft, const TSimpleMatchMap & aStrongestMatchesRight, const CLocation locLeft, TSimpleLocSet & matchesLeft, TSimpleLocSet & matchesRight);
        static void cullMatches(const TMatchMap & aStrongestMatchesLeft, TSimpleMatchMap & aStrongestMatchesLeftGood, double dCondition);
        static void ensureBiDiMatches(TSimpleMatchMap & aStrongestMatchesA, TSimpleMatchMap & aStrongestMatchesB);
        static void ensureBiDiMatches_remove(TSimpleMatchMap & aStrongestMatchesA, TSimpleMatchMap & aStrongestMatchesB);
        static void printMap(TSimpleMatchMap & aStrongestMatches);
    public:
        CBoWWordBag(CDescriptorSet * pDescriptors_in, const CBoW * pParent_in, int nId);
        ~CBoWWordBag();

        int BruteForceMatch(CBoW::CBoWWordBag * pWB) const;

        inline unsigned int TotalWordCount() const {
            return (unsigned int) pDescriptors->Count();
        }; //total words in this image=>same at every level

        void RecreateWordBag();

        template<bool DECREMENT>
        void CountWordOccurances() const;

        templateCompMethod
        void WeightWordBag();

        inline CDescriptorSet * DescriptorSet() const {
            return pDescriptors;
        };
        //inline CDescriptorSet * DescriptorSet() { return pDescriptors; };

        inline int id() const {
            return nId;
        };

        templateCompMethod
        int getBoWMatch(const CBoWWordBag * pwbExisting, const int * anScoreAgainstBackground, TSortedIntSet * pSubvecScores, int nVecReturnMax) const HOT;

        templateCompMethod
        TBoWMatchVector * getBoWMatches(int nReturnMax);

        void merge(CBoWWordBag * pOtherWB);

        const CBoWCorrespondences * getCorrespondences(const CBoW::CBoWWordBag * pWB2, const CBOWMatchingParams & MATCHING_PARAMS) const;
        const CBoWCorrespondences * getBruteForceCorrespondences(const CBoW::CBoWWordBag * pWB, const double dCondition, const int BF_NN) const HOT;
        const CBoWCorrespondences * getNewBruteForceCorrespondences(const CBoW::CBoWWordBag * pWB, const CBOWMatchingParams & MATCHING_PARAMS) const HOT;
    private:
        static const CBoWCorrespondences * getCorrespondences(const TWordBag & WB1, const TWordBag & WB2, int nN, int nM);

        //SPEEDO

        static bool containsWordLocations(const CBoWWord * pWord, const TWordBag & aWordsBottomLevel, int nBegin, int nEnd) {
            if (nBegin == nEnd) return false;
            int nMid = nBegin + (nEnd - nBegin) / 2;
            const CBoWWord * pMiddleWord = aWordsBottomLevel[nMid].Word();
            if (pMiddleWord > pWord)
                return containsWordLocations(pWord, aWordsBottomLevel, nBegin, nMid);
            else if (pMiddleWord < pWord)
                return containsWordLocations(pWord, aWordsBottomLevel, nMid + 1, nEnd);
            else
                return !(aWordsBottomLevel[nMid].Location(0).zero());
        }

        static void getWordLocations(const CBoWWord * pWord, const TWordBag & aWordsBottomLevel, int nBegin, int nEnd, TImPointVec & locations) {
            if(IS_DEBUG) CHECK(nBegin == nEnd, "getWordLocations: Should have checked that word exists here first");

            int nMid = nBegin + (nEnd - nBegin) / 2;
            const CBoWWord * pMiddleWord = aWordsBottomLevel[nMid].Word();
            if (pMiddleWord > pWord) {
                getWordLocations(pWord, aWordsBottomLevel, nBegin, nMid, locations);
                return;
            } else if (pMiddleWord < pWord) {
                getWordLocations(pWord, aWordsBottomLevel, nMid + 1, nEnd, locations);
                return;
            }

            if(IS_DEBUG) CHECK(pMiddleWord != pWord, "Chop failed to find word");

            int nOccurances = aWordsBottomLevel[nMid].FrequencyInImage();
            for (int i = 0; i < nOccurances; i++) {
                CLocation loc = aWordsBottomLevel[nMid].Location(i);
                if (loc.zero()) break;
                locations.push_back(loc);
            }
        }

    public:

        const TWordBag & bottomLevelWords() const {
            return aWords[pParent->PARAMS.BOWClustering.LEVELS - 1];
        }

        void getWordLocations(const CBoWWord * pWord, TImPointVec & locations) const {
            const TWordBag & aWordsBottomLevel = bottomLevelWords();
            int nWords = aWordsBottomLevel.size();
            getWordLocations(pWord, aWordsBottomLevel, 0, nWords, locations);
        }

        bool containsWordLocations(const CBoWWord * pWord) const {
            //Sorted, so do binary chop
            const TWordBag & aWordsBottomLevel = bottomLevelWords();
            int nWords = aWordsBottomLevel.size();

            return containsWordLocations(pWord, aWordsBottomLevel, 0, nWords);
        }
    };

    typedef CDynArray<CBoWWordBag *> TImageDB; //Todo tidy up map misuse risk (check existance first)
    typedef TImageDB::const_iterator constImIt;

    class CBoWImageDB : private TImageDB {
        int nAllDescriptorsCount;

        typedef std::pair<imageNum, iterator> WBIdPair;

        class ImageIdSort : std::binary_function<const WBIdPair &, const WBIdPair &, bool> {
        public:

            inline result_type operator() (first_argument_type a,
                    second_argument_type b) const {
                return (result_type) (a.first < b.first);
            };
        };
        std::map<imageNum, int> imIdLookupMap;
        int nErased;

        //Queue images to erase (avoid deleting descriptors that are being clustered)
        TImageDB vImagesToErase;

        inline iterator Find_int(imageNum nID) {
            std::map<imageNum, int>::iterator pMap = imIdLookupMap.find(nID);
            if (pMap != imIdLookupMap.end())
                return TImageDB::begin() + pMap->second;
            else
            {
                breakInCpp();
                //THROW("Find_int: Accessing image that doesnt exist");
                throw "Find_int: Accessing image that doesnt exist"; //so can catch
            }
        }

        inline const_iterator Find_int(imageNum nID) const {
            std::map<imageNum, int>::const_iterator pMap = imIdLookupMap.find(nID);
            if (pMap != imIdLookupMap.end())
                return TImageDB::begin() + pMap->second;
            else
            {
                breakInCpp();
                //THROW("Find_int: Accessing image that doesnt exist");
                throw "Find_int: Accessing image that doesnt exist";//so can catch
            }
        }
    public:

        int totalDescriptorCount() const {
            return nAllDescriptorsCount;
        }

        CBoWImageDB() : nAllDescriptorsCount(0), nErased(0) {
            TImageDB::reserve(1000);
        };

        inline unsigned int size_th() const {
            return (unsigned int) TImageDB::size();
        }; //used ONLY to split into 4 fo threading

        inline int Count() const {
            return TImageDB::size() - nErased;
        };

        inline const_iterator end() const {
            return TImageDB::end();
        };

        inline const_iterator begin() const {
            return TImageDB::begin();
        };

        inline void Push(CBoWWordBag * pWB);
        inline void Erase(imageNum nID, const bool);

        inline CBoWWordBag * Find(imageNum nID) {
            CBoWWordBag * pWB = *Find_int(nID);
            CHECK(!pWB, "Find: No WB found, does id exist?");
            return pWB;
        }

        bool exists(imageNum id) const {
            return imIdLookupMap.find(id) != imIdLookupMap.end();
        }

        CBoWWordBag * operator[](imageNum nId) {
            return *Find_int(nId);
        }

        const CBoWWordBag * operator[](imageNum nId) const {
            return *Find_int(nId);
        }
        //const CBoWWordBag * operator[](imageNum nID) { return *Find(nId); }
        ~CBoWImageDB();

        void cleanUpErasedImages(); //call when *not clustering*

        CDescriptorSet * makeNewDS() const {
            if(IS_DEBUG) CHECK(Count() == 0, "Making DS when there's no descriptors");

            for (constImIt ppIm = begin(); ppIm != end(); ppIm++)
                if (*ppIm)
                    return (*ppIm)->DescriptorSet()->makeNewDS();

            THROW("Error counting images")
        }
    };


    static const int MAX_NUM_COMP_METHODS = 20;
    TBoWMatchVector * (CBoWWordBag::*apfn_getBoWMatches[MAX_NUM_COMP_METHODS])(int); //Force compiler to make all the functions we need
    TBoWMatchVector * (CBoWWordBag::*pfn_getBoWMatches)(int);
    void (CBoWWordBag::*apfn_WeightWordBag[MAX_NUM_COMP_METHODS])(void);
    void (CBoWWordBag::*pfn_WeightWordBag)(void);
    void (CBoW::*apfn_RecreateWordBags[MAX_NUM_COMP_METHODS])(constImIt, constImIt, boost::barrier &);
    void (CBoW::*pfn_RecreateWordBags)(constImIt, constImIt, boost::barrier &);

    class WordBagSort : std::binary_function<const CBoWWordBag *, const CBoWWordBag *, bool> {
    public:

        inline result_type operator() (first_argument_type a,
                second_argument_type b) const {
            return (result_type) (a->id() < b->id());
        };
    };


    void IDToWBIterator(int nId1, int nId2, constImIt &ppWb1, constImIt &ppWb2) const;
    void IDToWBIterator(int nId1, int nId2, constImIt &ppWb1, constImIt &ppWb2);

    //int nLevels;
    //static int s_nLevels;

    //double dPriorProbOfAnyMatch, dProportionBackground;
    int * anMinMatchStrength; //need a different min strength for each level. Todo: can we calculate bounds on this?

    //frame vs. feature.
    CBoWImageDB vImages; //Some sort of tree per frame? We want to ask "What frames have word n?" Stick with arrays for now--we have a fairly small word list.

    CBoWDictionary * pDictionary; //All words, including mid-points

    //Do clustering
    void ClusterDescriptorsIntoWords();

    //May be dispatched in a seperate thread:
    void doClustering();
    void copyDescriptorsTS(CDescriptorSet & myDescriptorSetCopy, unsigned int & nTotalWordsTarget);
    void copyAllDescriptors(CDescriptorSet & myDescriptorSetCopy) const;

    templateCompMethod
    void RecreateWordBags(constImIt ppIm, constImIt ppImEnd, boost::barrier & wordWeightBarrier);
    templateCompMethod
    void RecreateWordBags_int(constImIt ppIm, constImIt ppImEnd, boost::barrier & wordWeightBarrier);
    void ReplaceDictionary_int(CBoW::CBoWDictionary ** ppNewDictionary);
    void ReplaceDictionary(CBoW::CBoWDictionary ** ppNewDictionary);

    void addImage(CBoWWordBag * pWB);

    double getBayesPosteriorFromPriorWithLambda(double dPrior, int nImageWords, int nWordFrequencyTotal, int nWordFrequencyImage, int nWordFrequencyNewImage);
    double getBayesPosteriorFromPriorMEstimate(double dPrior, int nImageWords, int nWordFrequencyTotal, int nWordFrequencyImage, int nWordFrequencyNewImage) const;
    inline static double getBayesScoreM(int nTimesWordAppearsInCat, int nTotalWordsInCat);
    inline static double getLogBayesScoreM(int nTimesWordAppearsInCat, int nTotalWordsInCat);

    void CreateTemplateFns();
    void SetupMinMatchThreshholds(CBOWParams::eCOMP_METHOD eCompMethod);

    TBoWMatchVector * getMatches_int(CBoWWordBag * pwb, int nReturnMax);

    static double intPow(double d, int i);
    static double intLog(int i);
    static double adIntLogCache[INT_LOG_CACHE_LIM]; //TODO: Libraryify
    static void setupIntLogCache();
    static int intintLog(int i);
    static int anIntLogCache[INTINT_LOG_CACHE_LIM];
    static void setupIntIntLogCache();
    static int nCache_i;
    static int nLog_i;

    //CDescriptorSet * pAllDescriptors;

    int nNextClusterCount; //After what descriptor count has been reached we will next cluster

    bool bClustering; // Don't start clustering until the last clustering attempt is done.
    boost::thread * pClusterThread;

    //double dReclusterFrequency;
    //unsigned int nDescriptorsPerWord;
    bool bDeleting;
    CException pClusteringException;
    CClusterDrawer * pDrawCS, * pDrawCS_temp; //For debugging (if non-zero)

    void checkClusteringException() {
        if (pClusteringException) {
            CException pClusteringExceptionTemp = pClusteringException;
            throw pClusteringExceptionTemp;
        }
    };

protected:

    virtual void resetBeforeRecreateWB() {
    };
    void RecreateWB();

    virtual void doOR() {
    };
    void testCorrespondences() const;

    static double geomAdd(double stat);
    static double geomScale(double stat);

    void ensureClusteredOnce(); //should NOT have a write lock when call this
public:
    //static int getNumLevels() { return s_nLevels; };

    COUNTHITS(static int nCache; static int nHit; static int nMiss); //check lookup is efficient

    //Apply a geometric compat. update to top nProcess matches
    //TBoWStatMatchVector * getGeomMatches(TBoWMatchVector * pMatches, int id, int nProcess);
    //class CScaleObserver;

    //CBoW(unsigned int nDescriptorsPerWord, int nLevels_in, eComparisonMethod eCompMethod, double dMinMatchStrength_in, double dPriorProbOfAnyMatch_in, double dProportionBackground_in, double dReclusterFrequency_in = DEFAULT_RECLUSTER_FREQUENCY, CBoWParams::CDescriptorBinningParams &* ppBinningParams = 0, bool bIOwnAddedDSs_in = false, CScaleObserver * pSpeedoScaleObserver = 0);
    CBoW(const CBOWParams & PARAMS);
    virtual ~CBoW();

    static const int DONT_ADD = -1;
    static const int RETURN_ALL = -1;

    //Does NOT add image
    TBoWMatchVector * getMatches(CDescriptorSet * pDescriptors_in, int nReturnMax = RETURN_ALL);

    //Adds image and returns matches. Memory now owned by CBoW

    TBoWMatchVector * getMatches(CDescriptorSet ** ppDescriptors_in, int nId_in, int nReturnMax = RETURN_ALL) {
        addImage(ppDescriptors_in, nId_in);
        return getMatches(nId_in, nReturnMax);
    }

    TBoWMatchVector * getMatches(int nId_in, int nReturnMax = RETURN_ALL);

    //Adds image. Memory now owned by CBoW
    void addImage(CDescriptorSet ** ppDescriptors, int nNewImageId);

    void recreateDictionary(); //Use vAllDescriptors to rebuild dictionary

    void remove(int nId);

    bool contains(int nId) const {
        return vImages.exists(nId);
    }

    void setClusterDebugger(CClusterDrawer * pDrawCS_in) {
        pDrawCS = pDrawCS_in;
    }

    // Get correspondences between 2 images. Includes N-M correspondences (see CBOWMatchingParams)
    const CBoWCorrespondences * getCorrespondences(CDescriptorSet * pDS1, int nId2, const CBOWMatchingParams & MATCHING_PARAMS);
    const CBoWCorrespondences * getCorrespondences(int nId1, CDescriptorSet * pDS2, const CBOWMatchingParams & MATCHING_PARAMS);
    const CBoWCorrespondences * getCorrespondences(CDescriptorSet * pDS1, CDescriptorSet * pDS2, const CBOWMatchingParams & MATCHING_PARAMS);
    const CBoWCorrespondences * getCorrespondences(int nId1, int nId2, const CBOWMatchingParams & MATCHING_PARAMS);

};

#endif
