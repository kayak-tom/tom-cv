/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#include "bagOfWords.h"
//DEBUGONLY(int CBoW::CBoWWordBag::CBoWImageWord::COPY_OK = 0;);
#include <cstdlib>			// C standard includes
#include <string>			// C++ strings
#include <algorithm>
#include "util/convert.h"
#include <boost/smart_ptr.hpp>
#include "time/SpeedTest.h"
#include "util/cout_TS.h"

#ifndef __GNUC__
#include <windows.h>
#else
#include <sched.h>
#endif

#define DEPRECATED(...) THROW( "Calling deprecated code")

pragma_warning(push)
pragma_warning(disable : 4127) //const conditional expression(from template params)

using namespace std;

#define INVALID_WORD (unsigned int)-1

#define DEBUG_CHECKIDX(i) //{if(IS_DEBUG) CHECK(i>MAX_WORDS || i<0, "CBoW::CBoW::CBoWDictionary: Index check failed");}

//Setup statics:

void CBoW::setupIntLogCache() {
    adIntLogCache[0] = 0;
    for (int i = 1; i < INT_LOG_CACHE_LIM; i++)
        adIntLogCache[i] = log((double) i);
}

void CBoW::setupIntIntLogCache() {
    anIntLogCache[0] = 0;
    for (int i = 1; i < INTINT_LOG_CACHE_LIM; i++)
        anIntLogCache[i] = doubleToInt(INTINTLOG_SCALE * log((double) i));
}

double CBoW::adIntLogCache[INT_LOG_CACHE_LIM] = {};
int CBoW::anIntLogCache[INTINT_LOG_CACHE_LIM] = {};
int CBoW::nCache_i = -1;
int CBoW::nLog_i = -1;
COUNTHITS(int CBoW::nHit = 0; int CBoW::nMiss = 0; int CBoW::nCache = 0);
//int CBoW::CBoWWordBag::BF_LEVEL = 4; //What level do we use for smart BF matching?
//int CBoW::CBoWWordBag::BF_NN = 3; //return 2-2 and 3-3 corr's

DEBUGONLY(int g_nClustersTooSmall, g_nClustersJustRight, g_nClustersTooBig, g_nWordsInsideRad, g_nWordsOutsideRad);

void CBoW::CBoWDictionary::ClusterRangeToWords_recurse(const CBOWParams::CDescriptorBinningParams & BINNING_PARAMS, CClusterSet *pClusters, const CGetClusterNum & getClusterCount, int nClusterStart, int nClusterEnd, const CBOWParams::CBOWClusteringParams & BOWCLUSTERPARAMS, CClusterDrawer ** ppDrawCS) {
    if(IS_DEBUG) CHECK(nClusterEnd > pClusters->size(), "Bad cluster centre range\n");
    if(IS_DEBUG) CHECK(nClusterEnd < nClusterStart, "Bad cluster centre range\n");
    if(IS_DEBUG) CHECK(nClusterStart < 0, "Bad cluster centre range\n");
    if(IS_DEBUG) CHECK(!pClusters, "No cluster set given\n");

    try {
        if (/*(nLevel <= 3) &&*/ ppDrawCS && *ppDrawCS) {
            CClusterDrawer * pDrawCS = *ppDrawCS;

            if (nLevel == 1)
                *ppDrawCS = 0;

            pDrawCS->drawCS(pClusters, nLevel);
        }

        for (wordNum nWord = (unsigned int) nClusterStart; nWord < (unsigned int) nClusterEnd; nWord++) {
            CCluster * pCluster = (*pClusters)[nWord];
            ClusterToWords_recurse(BINNING_PARAMS, getClusterCount, pCluster, BOWCLUSTERPARAMS, ppDrawCS);
        }
    } catch (CException pEx) {
        pDictionaryException = pEx;
    } catch (...) {
        pDictionaryException = CException("Unknown exception type caught in clustering sub-thread");
    }
}

//This is where we convert clusters into words, and recurse (cluster each cluster): Use a couple of threads, everything below here should be independent

void CBoW::CBoWDictionary::ClusterToWords_recurse(const CBOWParams::CDescriptorBinningParams & BINNING_PARAMS, const CGetClusterNum & getClusterCount, CCluster *pCluster, const CBOWParams::CBOWClusteringParams & BOWCLUSTERPARAMS, CClusterDrawer ** ppDrawCS) {
    //cout << "Bounds: " << BINNING_PARAMS.LOWER_BOUND() << '-' << BINNING_PARAMS.UPPER_BOUND() << '\n';
    const int LOWER_BOUND = min<int>(getClusterCount.numImages() / 2, BINNING_PARAMS.LOWER_BOUND);
    if (pCluster->Count() > LOWER_BOUND && (nLevel != 1 || pCluster->Count() < BINNING_PARAMS.UPPER_BOUND)) //otherwise discard this cluster and the words it contains -- too big/small
        // Some words may be within the radius of another cluster
    {
        if(IS_DEBUG) CHECK(pCluster->Count() <= 0, "CBoW::CBoWDictionary::CBoWDictionary: Cluster has no descriptors assigned");

        if (nLevel == 1) //Bottom level ('real' 'words'!)
        {
            Words.push_back(new CBoWWord(pCluster));
            DEBUGONLY(g_nClustersJustRight++);
        } else {
            //If this is a distinctive cluster there may only be 1 word, but still keep going to bottom level so that each word has an index at the bottom
            /*int nWordsAtThisLevel = (pCluster->Count()/3)+1;//Todo--const or something
            if(nWordsAtThisLevel>nWordsBelowEachWord) nWordsAtThisLevel=nWordsBelowEachWord; Not done here anymore*/

            CBoWNodeWord * pWordBelow = new CBoWNodeWord(pCluster, getClusterCount, nLevel - 1, anWordLevelCounts + 1, BINNING_PARAMS, BOWCLUSTERPARAMS, ppDrawCS);

            Words.push_back(pWordBelow);

            const int * anWordCountsBelow = pWordBelow->SubDictionary()->WordCountArray();

            for (int nEachLevel = 1; nEachLevel < nLevel; nEachLevel++) {
                DEBUG_CHECKIDX(nEachLevel);
                anWordLevelCounts[nEachLevel] += anWordCountsBelow[nEachLevel - 1];
            }
        }
    }
    DEBUGONLY(
    else if (nLevel == 1) {
        if (pCluster->Count() <= LOWER_BOUND)
                g_nClustersTooSmall++;
        else
            g_nClustersTooBig++;
                //DEBUGONLY(cout << "Cluster not mapped to word: " << pCluster->Count() << " descriptors\n");
        });
}

int CBoW::CBoWDictionary::CGetClusterNum::getClusterNum(const int nLevel, const int nDescriptors, const CBOWParams::CBOWClusteringParams & BOWCLUSTERPARAMS) const {
    int nClusters = -1;

    if (BOWCLUSTERPARAMS.BRANCH_METHOD == CBOWParams::CBOWClusteringParams::eFixedBF) {
        nClusters = nBranchFactor;
    } else if (BOWCLUSTERPARAMS.BRANCH_METHOD == CBOWParams::CBOWClusteringParams::eFixedClusterSizeTarget) {
        int nTotalWordsAtThisLevel = nBranchFactor;

        for (int nLevelsDown = BOWCLUSTERPARAMS.LEVELS; nLevelsDown > nLevel; nLevelsDown--)
            nTotalWordsAtThisLevel *= nBranchFactor;

        nClusters = 1 + doubleToInt(((double) nTotalWordsAtThisLevel * (double) nDescriptors) / ((double) nTotalDescriptors)); //prevent overflow
    } else
        THROW("Unhandled clustering strategy");

    const int MIN_DESCRIPTORS_PER_WORD = min<int>(nImages, 8); //otherwise they're just stopwords (although at first we need lots of words with not so many descriptors)
    //Kind-of corresponds to seeing runs of 10 scenes where a distinctive feature is visible.
    int nMaxClusters = max<int>(nDescriptors / MIN_DESCRIPTORS_PER_WORD, 1);

    if (nMaxClusters < nClusters) {
        /*char pc[100];
        sprintf4(pc, 100, "LIMITING CLUSTERS Level=%d Descriptors=%d Clusters=%d new Clusters=%d\n", nLevel, nDescriptors, nClusters, nMaxClusters );
        cout << pc;*/
        nClusters = nMaxClusters;
    }

    //char pc[100];
    //sprintf3(pc, 100, "Level=%d TotalWords=?? Descriptors=%d Clusters=%d\n", nLevel, nTotalDescriptors, nClusters );
    //cout << pc;

    return nClusters;
}

//Build new dictionary from a load of descriptors
//May be dispatched many times in seperate threads

CBoW::CBoWDictionary::CBoWDictionary(const CCluster * pParentCluster, CDescriptorSet * pDescriptors, const CGetClusterNum & getClusterCount, int nLevel_in, const int * anLevelWordCounts_in, const CBOWParams::CDescriptorBinningParams & BINNING_PARAMS, const CBOWParams::CBOWClusteringParams & BOWCLUSTERPARAMS, CClusterDrawer ** ppDrawCS)
: nLevel(nLevel_in), anWordLevelCounts(0), pClusters(0) {
    if(IS_DEBUG) CHECK(nLevel <= 0 || nLevel_in >= MAX_BOW_LEVELS || !anLevelWordCounts_in || (*anLevelWordCounts_in < 0), "CBoW::CBoW::CBoWDictionary: Bad parameters");

    if(IS_DEBUG) CHECK(Words.size(), "CBoW::CBoW::CBoWDictionary: Dictionary already contains words");

    int nClusters = getClusterCount.getClusterNum(nLevel, pDescriptors->Count(), BOWCLUSTERPARAMS);

    anWordLevelCounts = new int[nLevel]; //count of words at this level
    //memcpy(anWordLevelCounts, anLevelWordCounts_in, sizeof(int)*nLevel);
    //memset(anWordLevelCounts, 0, sizeof(int)*nLevel);
    setZero(anWordLevelCounts, nLevel);

    /*boost::scoped_ptr<CClusterSet>*/ pClusters = pDescriptors->Cluster(nClusters, pParentCluster);
    CHECK(!pClusters->Assigned(), "CBoW::CBoW::CBoWDictionary: Descriptors not assigned to centres");

    if (nClusters != pClusters->size()) {
        cout << nClusters << " != " << pClusters->size() << ": Less clusters than requested found\n";
        if (nClusters + 5 < pClusters->size() / 2) {
            cout << "Warning: Big reduction in number of clusters found. Consider increasing number of levels (BOW.BOWClustering.LEVELS). Currently ";
            BOWCLUSTERPARAMS.LEVELS.printParam();
            cout << "\n";
        }
    }

    nClusters = pClusters->size();

    anWordLevelCounts[0] += nClusters;
    nFirstWordIdx = anLevelWordCounts_in[0];

    //pClusters->CheckDisimilarity();

    Words.reserve(nClusters);

    //Todo: do some research into tree balancing -- look at num of words/cluster here
    /*int nWordsBelowEachWord = -1;
    if(nLevel>1)
        nWordsBelowEachWord=doubleToInt(nBranchFactor/dWordsInThisDictionary);*/ //=nBranchFactor ^ ((level-1)/level)

    if (BINNING_PARAMS.CLUSTER_THREADS > 1 && nLevel == BOWCLUSTERPARAMS.LEVELS) //Spawn threads here. Probably 1 or 2 levels under here so worth clustering
    {
        const int NUM_WORDS_PER_THREAD = nClusters / BINNING_PARAMS.CLUSTER_THREADS;
        CDynArrayOwner<boost::thread> apClusterThreads(BINNING_PARAMS.CLUSTER_THREADS);

        for (int nThread = 0; nThread < BINNING_PARAMS.CLUSTER_THREADS; nThread++) {
            int nClusterStart = nThread*NUM_WORDS_PER_THREAD;
            int nClusterEnd = nThread * NUM_WORDS_PER_THREAD + NUM_WORDS_PER_THREAD;
            if (nThread == BINNING_PARAMS.CLUSTER_THREADS - 1) nClusterEnd = nClusters;

            //char pc[100];
            //sprintf4(pc, 100, "Spawning thread %d nLevel=%d start=%d end=%d\n", nThread, nLevel, nClusterStart, nClusterEnd);
            //	cout << pc;

            apClusterThreads[nThread] = new boost::thread(boost::bind(&CBoW::CBoWDictionary::ClusterRangeToWords_recurse,
                    this, boost::ref(BINNING_PARAMS), pClusters, boost::ref(getClusterCount), nClusterStart, nClusterEnd, boost::ref(BOWCLUSTERPARAMS), ppDrawCS));
        }
        for (int nThread = 0; nThread < BINNING_PARAMS.CLUSTER_THREADS; nThread++) {
            apClusterThreads[nThread]->join();
        }
    } else
        ClusterRangeToWords_recurse(BINNING_PARAMS, pClusters, getClusterCount, 0, nClusters, BOWCLUSTERPARAMS, ppDrawCS);

    //Now all words are created in this dictionary. Sort them by cluster * to speed-up assignment
    std::sort(Words.begin(), Words.end(), CBoW::CBoWWord::sortWordsByDescPtr);

    checkDictionaryException();
}

//May be dispatched in a seperate thread
//Everyone else should be locked out now

templateCompMethod
void CBoW::RecreateWordBags(constImIt ppImStart, constImIt ppImEnd, boost::barrier & wordWeightBarrier) {
    try {
        RecreateWordBags_int<eCompMethod > (ppImStart, ppImEnd, wordWeightBarrier);
    } catch (CException pEx) {
        pClusteringException = pEx;
    } catch (...) {
        pClusteringException = CException("Unknown exception type caught in RecreateWordBags thread");
    }
}
//May be dispatched in a seperate thread
//Everyone else should be locked out now
//Delete existing words and Re-map descriptors to words

templateCompMethod
void CBoW::RecreateWordBags_int(constImIt ppImStart, constImIt ppImEnd, boost::barrier & wordWeightBarrier) {
    for (constImIt ppIm = ppImStart; ppIm < ppImEnd; ppIm++) {
        //Delete existing words and Re-map descriptors to words
        CBoWWordBag * pWB = *ppIm;
        if (pWB)
            pWB->RecreateWordBag();
        //else
        //cout << "One word bag doesn't exist\n"; Now ok, its been erased
    }

    if (eCompMethod == CBOWParams::eNisterDist || eCompMethod == CBOWParams::eVectorDistFast) {
        //Now calc weights for WB vectors. Reset when we recreated dictionary. Might be more descriptors now than when we started but will count occurances of those too
        //**New images do not affect occurance counts until we recluster**
        try {
            if (PARAMS.WB_WEIGHT_METHOD == CBOWParams::eTF_IDF || PARAMS.WB_WEIGHT_METHOD == CBOWParams::eTF_IDF_Wikipedia) {
                //Now count words. All threads wait here and take turns
                LOCK_WHILE_COUNTING_OCCURANCES;
                cout << "Counting word occurances...";
                for (constImIt ppIm = ppImStart; ppIm < ppImEnd; ppIm++) {
                    CBoWWordBag * pWB = *ppIm;
                    if (pWB)
                        pWB->CountWordOccurances < false > ();
                }
                //cout << "done-unlocking\n";
            }//unlock
        } catch (...) {
            if (PARAMS.RWB_THREADS > 1)
                wordWeightBarrier.wait();

            throw;
        }
        //Concurrent again
        if (PARAMS.RWB_THREADS > 1)
            wordWeightBarrier.wait();

        //Now re-score
        //cout << "Scoring words in images...\n";
        for (constImIt ppIm = ppImStart; ppIm < ppImEnd; ppIm++) {
            CBoWWordBag * pWB = *ppIm;
            if (pWB)
                (pWB->*pfn_WeightWordBag)();
        }
    }
}

//May be dispatched in one separate clustering thread
//Everyone else is (should be) locked out when this is called, so can reset + count occurances safely

void CBoW::RecreateWB() {
    resetBeforeRecreateWB();

    //Re-create every image's wordBag AND rescore if necessary
    cout << "Rebuilding word bags...\n";

    boost::barrier wordWeightBarrier(PARAMS.RWB_THREADS); //won't be locked if PARAMS.RWB_THREADS==0

    if (PARAMS.RWB_THREADS > 1) {
        const int NUM_IMAGES_PER_THREAD = (int) (vImages.size_th()) / PARAMS.RWB_THREADS;
        CDynArrayOwner<boost::thread> apRWBThreads(PARAMS.RWB_THREADS);

        for (int nThread = 0; nThread < PARAMS.RWB_THREADS; nThread++) {
            constImIt ppImStart = vImages.begin() + nThread*NUM_IMAGES_PER_THREAD;
            constImIt ppImEnd = ppImStart + NUM_IMAGES_PER_THREAD;
            if (nThread == PARAMS.RWB_THREADS - 1) ppImEnd = vImages.end();

            apRWBThreads[nThread] = new boost::thread(boost::bind(pfn_RecreateWordBags,
                    this, ppImStart, ppImEnd, boost::ref(wordWeightBarrier)));
        }
        for (int nThread = 0; nThread < PARAMS.RWB_THREADS; nThread++) {
            apRWBThreads[nThread]->join();
        }

        checkClusteringException();
    } else
        (this->*pfn_RecreateWordBags)(vImages.begin(), vImages.end(), wordWeightBarrier);

    testCorrespondences();
    //DEBUGONLY(testCorrespondences());

    doOR();
}

void CBoW::testCorrespondences() const {
    //Get first 2 WBs. Check
    if (vImages.Count() < 2) return;

    const CBoWWordBag * pWB0 = 0;
    const CBoWWordBag * pWB1 = 0;
    for (int i = 0; !pWB1; i++) {
        if (vImages.exists(i)) {
            if (!pWB0)
                pWB0 = vImages[i];
            else
                pWB1 = vImages[i];
        }
    }

    CStopWatch s;

    s.startTimer();
    CBOWMatchingParams MATCHING_PARAMS(0, 0);
    boost::scoped_ptr<const CBoWCorrespondences> pCorrBoW(pWB1->getCorrespondences(pWB0, MATCHING_PARAMS));
    s.stopTimer();
    const double dBoWTime = s.getElapsedTime();

    s.startTimer();
    MATCHING_PARAMS.BF_CORRESPONDENCES = CBOWMatchingParams::eOldBoW_BF_Correspondences;

    boost::scoped_ptr<const CBoWCorrespondences> pOldCorrBF(pWB0->getCorrespondences(pWB1, MATCHING_PARAMS));

    s.stopTimer();
    const double dOldBoWBFTime = s.getElapsedTime();

    s.startTimer();

    MATCHING_PARAMS.BF_CORRESPONDENCES = CBOWMatchingParams::eBoW_BF_Correspondences;

    boost::scoped_ptr<const CBoWCorrespondences> pCorrBF(pWB0->getCorrespondences(pWB1, MATCHING_PARAMS));
    s.stopTimer();
    const double dBFTime = s.getElapsedTime();

    s.startTimer();

    const int BML = PARAMS.BOWClustering.BF_LEVEL;
    const_cast<CNumParam<int> &> (PARAMS.BOWClustering.BF_LEVEL) = 0;
    boost::scoped_ptr<const CBoWCorrespondences> pCorrBest(pWB0->getBruteForceCorrespondences(pWB1, MATCHING_PARAMS.BF_CORNER_CONDITION, MATCHING_PARAMS.MATCH_NN));
    const_cast<CNumParam<int> &> (PARAMS.BOWClustering.BF_LEVEL) = BML;
    s.stopTimer();
    const double dOldBFTime = s.getElapsedTime();

    s.startTimer();

    const CMatchableDescriptors::CMatchSettings MS(MATCHING_PARAMS.BF_CORNER_CONDITION, 0.7, MATCHING_PARAMS.MATCH_NN, -1, 1, MATCHING_PARAMS.OI_NEARBY_ANGLE);
    boost::scoped_ptr<const CBoWCorrespondences> pCorrNewBF(pWB1->DescriptorSet()->getBruteForceCorrespondenceSet(pWB0->DescriptorSet(), MS, 0));

    s.stopTimer();
    const double dNewBFTime = s.getElapsedTime();

    int nBF = 0, nBoW = 0, nOldBF = 0, nOldBoWBF = 0;
    for (CBoWCorrespondences::const_iterator pCorr = pCorrNewBF->begin(); pCorr != pCorrNewBF->end(); pCorr++) {
        if (pCorrBF->contains(pCorr))
            nBF++;

        if (pCorrBest->contains(pCorr))
            nOldBF++;

        if (pCorrBoW->contains(pCorr))
            nBoW++;

        if (pOldCorrBF->contains(pCorr))
            nOldBoWBF++;
    }
    cout << pCorrNewBF->size() << " new BF correspondences, time=" << dNewBFTime << "\n";
    cout << pCorrBF->size() << " BF correspondences, prop that are prob good: " << (double) nBF / pCorrBF->size() << " (" << nBF << ") time=" << dBFTime << "\n";
    cout << pCorrBoW->size() << " BoW correspondences, prop that are prob good: " << (double) nBoW / pCorrBoW->size() << " (" << nBoW << ") time=" << dBoWTime << "\n";
    cout << pCorrBest->size() << " old full BF correspondences, prop that are prob good: " << (double) nOldBF / pCorrBest->size() << " (" << nOldBF << ") time=" << dOldBFTime << "\n";
    cout << pOldCorrBF->size() << " old BoW BF correspondences, prop that are prob good: " << (double) nOldBoWBF / pOldCorrBF->size() << " (" << nOldBoWBF << ") time=" << dOldBoWBFTime << "\n";
}

void CBoW::copyAllDescriptors(CDescriptorSet & myDescriptorSetCopy) const {
    for (constImIt ppIm = vImages.begin(); ppIm != vImages.end(); ppIm++) {
        const CBoWWordBag * pIm = *ppIm;
        if (pIm)
            myDescriptorSetCopy.Push(pIm->DescriptorSet());
    }
    if (myDescriptorSetCopy.Count() != vImages.totalDescriptorCount())
        cout << "POSSIBLE ERROR: CBoW::copyAllDescriptors: Descriptor count mismatch (ok if we have just erased an image)";
}

void CBoW::copyDescriptorsTS(CDescriptorSet & myDescriptorSetCopy, unsigned int & nTotalWordsTarget) {
    if (PARAMS.BOWClustering.CLUSTER_IN_SEPERATE_THREAD) {
        WRITE_LOCK;
        copyAllDescriptors(myDescriptorSetCopy);
    } else
        copyAllDescriptors(myDescriptorSetCopy);

    nTotalWordsTarget = myDescriptorSetCopy.Count() / PARAMS.BOWClustering.DESCRIPTORS_PER_WORD;
    
    if(nTotalWordsTarget == 0)
    {
        nTotalWordsTarget = 1;
        if(myDescriptorSetCopy.Count() < 50)
        {
            cout << "WARNING: Target number of image words is 1. Very few image features have been detected ("
                    << myDescriptorSetCopy.Count() 
                    << "). Possible causes: featureless images, parameter Im.SCALE_DOWN set too high, BOW.BOWClustering.DESCRIPTORS_PER_WORD too high, Corner.MAX_FEATURES too low." << endl;
        }
    }

    //Actually don't want too many clusters at first, we're brute-force matching... unsigned int nMinWords =  myDescriptorSetCopy.Count() / (2*vImages.Count()); //For well cond. correspondences in first few images (where most words are found in each image) we need about this many words
    //	unsigned int nMinWords = (*vImages.begin())->TotalWordCount()/2 + 2*vImages.Count(); //For well cond. correspondences in first few images (where most words are found in each image) we need about this many words

    //if (nTotalWordsTarget<nMinWords) nTotalWordsTarget=nMinWords; //We need about 0.45*descriptorsperimage for the first few frames to get good correspondences

    RELEASE_WRITE_LOCK;
}

enum ePriority {
    eRTPriority, eHighPriority, eAboveNormalPriority, eNormalPriority, eBelowNormalPriority, eIdlePriority
};

#ifdef __GNUC__

int getPthreadPriority(ePriority priority) //1-100
{
    int nMax = sched_get_priority_max(SCHED_OTHER);
    int nMin = sched_get_priority_min(SCHED_OTHER);
    switch (priority) {
        case eRTPriority:
            return nMax;
        case eHighPriority:
            return nMin + 3 * (nMax - nMin) / 4;
        case eAboveNormalPriority:
            return nMin + 2 * (nMax - nMin) / 3;
        case eNormalPriority:
            return nMin + (nMax - nMin) / 2;
        case eBelowNormalPriority:
            return nMin + (nMax - nMin) / 4;
        case eIdlePriority:
            return nMin;
        default:
            THROW("Not handled");
    }
}
#else

unsigned int getWinPriority(ePriority priority) {
    switch (priority) {
        case eRTPriority:
            return REALTIME_PRIORITY_CLASS;
        case eHighPriority:
            return HIGH_PRIORITY_CLASS;
        case eAboveNormalPriority:
            return ABOVE_NORMAL_PRIORITY_CLASS;
        case eNormalPriority:
            return NORMAL_PRIORITY_CLASS;
        case eBelowNormalPriority:
            return BELOW_NORMAL_PRIORITY_CLASS;
        case eIdlePriority:
            return IDLE_PRIORITY_CLASS;
        default:
            THROW("Not handled");
    }
}
#endif

void applyPriority(boost::thread & thread, ePriority priority) {
#ifndef __GNUC__
    HANDLE th = thread.native_handle();

    BOOL res = SetPriorityClass(th, getWinPriority(priority));

    if (res == FALSE) {
        cout << "Error setting thread priority (Windows)\n";
        //int err = GetLastError();
    }
#else
    /*
    struct sched_param param;
    param.sched_priority = getPthreadPriority(priority);
    int policy = SCHED_OTHER; //normal
    int retcode = pthread_setschedparam(thread.native_handle(), policy, &param))
if (retcode != 0)
{
            cout << "Error setting thread priority (pthread)\n";
}
     */
#endif
}

void setMyPriority(ePriority priority) {
#ifdef __GNUC__
    struct sched_param param;
    param.sched_priority = getPthreadPriority(priority);
    int policy = SCHED_OTHER; //normal
    int retcode = pthread_setschedparam(0, policy, &param);
    if (retcode != 0) {
        cout << "Error setting my priority (pthread)\n";
    }
#endif
}

void CBoW::doClustering() {
    cout << "Started clustering thread\n";

    if (PARAMS.BOWClustering.CLUSTER_IN_SEPERATE_THREAD && false) //not working...
        setMyPriority(eBelowNormalPriority);

    const int LEVELS = PARAMS.BOWClustering.LEVELS;
    ARRAY(int, anTopLevelWordCount, LEVELS); //static NOT ok because only *1* separate thread here
    CBoWDictionary * pNewDictionary = 0;
    try {
        if(IS_DEBUG) CHECK(LEVELS >= MAX_BOW_LEVELS || LEVELS <= 0 || !vImages.totalDescriptorCount(), "CBoW::doClustering: Bad parameters/state");

        memset(PTR(anTopLevelWordCount), 0, LEVELS * sizeof (int)); //Set current running totals of words at each level to 0

        boost::scoped_ptr<CDescriptorSet> myDescriptorSetCopy(vImages.makeNewDS());
        unsigned int nTotalWordsTarget=0;

        //First lock here, then unlock
        copyDescriptorsTS(*myDescriptorSetCopy.get(), nTotalWordsTarget);

        int nBranchFactor = doubleToInt(pow((double) nTotalWordsTarget, 1. / (double) LEVELS));
        cout << "Clustering: Total descriptors=" << myDescriptorSetCopy->Count() << " branch factor=" << nBranchFactor << '\n';

        DEBUGONLY(g_nClustersTooBig = g_nClustersTooSmall = g_nClustersJustRight = g_nWordsOutsideRad = g_nWordsInsideRad = 0);
        pDrawCS_temp = pDrawCS;
        CBoWDictionary::CGetClusterNum getClusterNum(nBranchFactor, myDescriptorSetCopy->Count(), (int) vImages.Count());

        pNewDictionary = new CBoWDictionary(0, myDescriptorSetCopy.get(), getClusterNum, LEVELS, PTR(anTopLevelWordCount), PARAMS.DescriptorBinning, PARAMS.BOWClustering, pDrawCS_temp ? &pDrawCS_temp : 0); //This becomes the top level dictionary. All other dictionaries are members of CBoWNodeWord

        //Now we need to lock everything again while we replace the old dictionary: Further down we may be multi-threaded again, so be careful locking again
        ReplaceDictionary(&pNewDictionary);

        if(IS_DEBUG) CHECK(pNewDictionary, "New dictionary not given away--will leak");
        cout << "Clustering completed successfully\n";

#ifdef _DEBUG
        if (g_nClustersTooBig > g_nClustersJustRight / 20) {
            cout << "WARNING: Many clusters are too big: " << g_nClustersTooBig << " of " << g_nClustersJustRight << endl
                    << "Consider increasing ";
            PARAMS.DescriptorBinning.UPPER_BOUND.printParam();
        } else if (g_nClustersTooBig < g_nClustersJustRight / 100) {
            cout << "WARNING: Few clusters are too big: " << g_nClustersTooBig << " of " << g_nClustersJustRight << endl
                    << "Consider decreasing ";
            PARAMS.DescriptorBinning.UPPER_BOUND.printParam();
        }

        if (g_nClustersTooSmall > g_nClustersJustRight / 10) {
            cout << "WARNING: Many clusters are too small: " << g_nClustersTooSmall << " of " << g_nClustersJustRight << endl
                    << "Consider decreasing ";
            PARAMS.DescriptorBinning.LOWER_BOUND.printParam();
        } else if (g_nClustersTooSmall < g_nClustersJustRight / 25) {
            cout << "WARNING: Few clusters are too small: " << g_nClustersTooSmall << " of " << g_nClustersJustRight << endl
                    << "Consider increasing ";
            PARAMS.DescriptorBinning.LOWER_BOUND.printParam();
        }

        if (g_nWordsOutsideRad > g_nWordsInsideRad / 10) {
            cout << "WARNING: Many descriptors are outside radius of closest cluster: " << g_nWordsOutsideRad << " of " << g_nWordsInsideRad << endl
                    << "Consider increasing ";
            PARAMS.DescriptorBinning.RADIUS.printParam();
        } else if (g_nWordsOutsideRad < g_nWordsInsideRad / 50) {
            cout << "WARNING: Almost all descriptors are inside radius of closest cluster: " << g_nWordsOutsideRad << " of " << g_nWordsInsideRad << endl
                    << "Consider decreasing ";
            PARAMS.DescriptorBinning.RADIUS.printParam();
        }
#endif
    } catch (CException pEx) {
        delete pNewDictionary;
        pClusteringException = pEx;
        //..and unlock
    } catch (...) {
        delete pNewDictionary;
        pClusteringException = CException("Unknown exception type caught in seperate clustering thread");
        //..and unlock
    }

    mxClusteringCantDelete.unlock();
    bClustering = false; //done--can start again now
}

void CBoW::ReplaceDictionary_int(CBoW::CBoWDictionary ** ppNewDictionary) {
    cout << "Now erasing images while have a write lock before RWB, no more can be erased until we are done...\n";
    //vImages.cleanUpErasedImages(false); //before clustering starts again AND while mxClusteringCantDelete is set

    const int LEVELS = PARAMS.BOWClustering.LEVELS;

    if (pDictionary) {
        delete pDictionary;
        pDictionary = 0;
    }
    pDictionary = *ppNewDictionary;
    *ppNewDictionary = 0;

    const int * anWordCounts = pDictionary->WordCountArray();
    cout << anWordCounts[LEVELS - 1] << " words in total. ";
    if (anWordCounts[LEVELS - 1])
        cout << vImages.totalDescriptorCount() / anWordCounts[LEVELS - 1] << " descriptors per word, target " << (int) PARAMS.BOWClustering.DESCRIPTORS_PER_WORD << "\n";
    else
        cout << vImages.totalDescriptorCount() << " descriptors and no clusters; target " << (int) PARAMS.BOWClustering.DESCRIPTORS_PER_WORD << "\n";

    RecreateWB();
}

void CBoW::ReplaceDictionary(CBoW::CBoWDictionary ** ppNewDictionary) {
    if(IS_DEBUG) CHECK(!ppNewDictionary || !*ppNewDictionary, "ReplaceDictionary: No dictionary supplied")
    if (PARAMS.BOWClustering.CLUSTER_IN_SEPERATE_THREAD) {
        WRITE_LOCK;

        if (!bDeleting)
            ReplaceDictionary_int(ppNewDictionary);
        else {
            cout << "Not replacing dictionary because deleting everything" << endl;
            delete *ppNewDictionary;
            *ppNewDictionary = 0;
        }
        RELEASE_WRITE_LOCK;
    } else {
        ReplaceDictionary_int(ppNewDictionary);
    }
}

void CBoW::ClusterDescriptorsIntoWords() {
    if (PARAMS.BOWClustering.CLUSTER_IN_SEPERATE_THREAD) {
        cout << "Launching clustering thread...\n";
        if (!bDeleting) {
            mxClusteringCantDelete.lock(); //launch here rather than after the thread starts...

            if (pClusterThread) //We shouldn't be here if the thread is still running! bClustering flag should stop us
            {
                delete pClusterThread;
            }
            pClusterThread = new boost::thread(boost::bind(&CBoW::doClustering, this));
            cout << "Forked clustering thread\n";
            applyPriority(*pClusterThread, eBelowNormalPriority);
        }
    } else {
        mxClusteringCantDelete.lock();
        doClustering();
    }
}

CBoW::CBoWWord::CBoWWord(CCluster * pCluster) {
    pDescriptor = pCluster->Centre();
    bIOwnMyDescriptor = pCluster->IWantCentreMemory(); //Cluster may retain memory briefly...
    nTotalFrequency = pCluster->Count();
    if(IS_DEBUG) CHECK(nTotalFrequency <= 0, "CBoW::CBoWWord::CBoWWord: Cluster has no descriptors assigned");
    if(IS_DEBUG) CHECK(!pDescriptor, "CBoW::CBoWWord::CBoWWord: No cluster centre");

    nTotalOccurances = 0; //how many images does it occur in
}

CBoW::CBoWNodeWord::CBoWNodeWord(CCluster * pCluster, const CBoW::CBoWDictionary::CGetClusterNum & getClusterCount, int nLevel, const int * anLevelWordCounts_in, const CBOWParams::CDescriptorBinningParams & BINNING_PARAMS, const CBOWParams::CBOWClusteringParams & BOWCLUSTERPARAMS, CClusterDrawer ** ppDrawCS)
: CBoWWord(pCluster) {
    pSubDictionary = new CBoW::CBoWDictionary(pCluster, pCluster->Members(), getClusterCount, nLevel, anLevelWordCounts_in, BINNING_PARAMS, BOWCLUSTERPARAMS, ppDrawCS);
}

CBoW::CBoWWord::~CBoWWord() {
    if (bIOwnMyDescriptor)
        delete pDescriptor;
}

CBoW::CBoWWord * CBoW::findBinaryChop(const CDescriptor * pDesc, CBoW::CBoWWord * const* apWords, const int nLength) {
    if (nLength < 4) {
        CBoW::CBoWWord * const* ppWords = apWords;
        for (int i = 0; i < nLength; i++, ppWords++)
            if ((*ppWords)->Descriptor() == pDesc)
                return *ppWords;

        return 0;
    } else {
        const int nMiddle = nLength / 2;
        if (apWords[nMiddle]->Descriptor() > pDesc)
            return findBinaryChop(pDesc, apWords, nMiddle);
        else
            return findBinaryChop(pDesc, apWords + nMiddle, nLength - nMiddle);
    }
}

//Pass in an array apWords[nLevelsremaining] that will be filled with pointers to words.

void CBoW::CBoWDictionary::LookupWordAllLevels(const CDescriptor * pDescriptor, int nLevel, CBoWWord ** apWords, const CDescriptor ** apClosestDescriptors, const CBOWParams::CDescriptorBinningParams & BINNING_PARAMS) const {
    //if(IS_DEBUG) CHECK( !BINNING_PARAMS, "CBoW::CBoWDictionary::LookupWordAllLevels: Empty dictionary/No binning params");

    if (Words.size() == 0) {
        setZero(apWords, nLevel + 1);
        return;
    }
    CBoWWord * pClosestWord = 0;
    const CDescriptor * pClosestDesc = *apClosestDescriptors;
    if (pClosestDesc) {
        /*cout << "Level " << nLevel << ": Searching for " << pClosestDesc << " in ";
        for(CDynArray<CBoWWord *>::const_iterator ppWord = Words.begin();  ppWord != Words.end(); ppWord++)
                cout << (*ppWord)->Descriptor() << ',';
        cout << endl;*/

        pClosestWord = findBinaryChop(pClosestDesc, Words.begin(), Words.size());

        DEBUGONLY(
                static int s_nWordsNotFound = 0, s_nWordsFound = 0;
        if (IS_DEBUG && !pClosestWord)
                s_nWordsNotFound++;
                //cout << "Warning: word not found. Maybe its cluster was too small\n";
        else
            s_nWordsFound++;

            if (s_nWordsNotFound > s_nWordsFound / 5) {
            cout << "Warning: " << s_nWordsNotFound << " words not found, " << s_nWordsFound << " words total. Prob from clusters being too small\n";
        }
        )
    }

    if (nLevel > 0) {
        if (!pClosestWord)
            pClosestWord = ClosestWord(pDescriptor, MAX_ALLOWED_DIST);
        *apWords = pClosestWord;

        CBoWNodeWord * pNodeWord = static_cast<CBoWNodeWord *> (pClosestWord);
        return pNodeWord->SubDictionary()->LookupWordAllLevels(pDescriptor, nLevel - 1, apWords + 1, apClosestDescriptors + 1, BINNING_PARAMS);
    } else {
        //Todo Look at radius, etc if we have pClosestWord
        if (!pClosestWord)
            pClosestWord = ClosestWord(pDescriptor, BINNING_PARAMS.RADIUS);
        else {
            CDescriptor::TDist dist = pClosestDesc->assignmentDist();
            if (dist == MAX_ALLOWED_DIST) //This means the actual distance wasn't cached
                dist = pClosestDesc->distance(pDescriptor);

            if (dist > BINNING_PARAMS.RADIUS)
                pClosestWord = 0;
        }
        *apWords = pClosestWord;
    }
}

void CBoW::CreateTemplateFns() {
    setZero(apfn_getBoWMatches, MAX_NUM_COMP_METHODS);
    setZero(apfn_WeightWordBag, MAX_NUM_COMP_METHODS);
    setZero(apfn_RecreateWordBags, MAX_NUM_COMP_METHODS);

    apfn_getBoWMatches[(int) CBOWParams::eBruteForce] = &CBoWWordBag::getBoWMatches<CBOWParams::eBruteForce>;
    apfn_getBoWMatches[(int) CBOWParams::eNisterDist] = &CBoWWordBag::getBoWMatches<CBOWParams::eNisterDist>;
    apfn_getBoWMatches[(int) CBOWParams::eVectorDistFast] = &CBoWWordBag::getBoWMatches<CBOWParams::eVectorDistFast>;
    /*apfn_getBoWMatches[eLogMBayes] = &CBoWWordBag::getBoWMatches<eLogMBayes>;
    apfn_getBoWMatches[eLogBayes] = &CBoWWordBag::getBoWMatches<eLogBayes>;
    apfn_getBoWMatches[eMBayes] = &CBoWWordBag::getBoWMatches<eMBayes>;
    apfn_getBoWMatches[eLambda] = &CBoWWordBag::getBoWMatches<eLambda>;*/

    apfn_WeightWordBag[(int) CBOWParams::eBruteForce] = &CBoWWordBag::WeightWordBag<CBOWParams::eBruteForce>;
    apfn_WeightWordBag[(int) CBOWParams::eNisterDist] = &CBoWWordBag::WeightWordBag<CBOWParams::eNisterDist>;
    apfn_WeightWordBag[(int) CBOWParams::eVectorDistFast] = &CBoWWordBag::WeightWordBag<CBOWParams::eVectorDistFast>;
    //apfn_RecreateWordBag[eLogMBayes] = &CBoWWordBag::RecreateWordBag<eLogMBayes>;
    //apfn_RecreateWordBag[eLogBayes] = &CBoWWordBag::RecreateWordBag<eLogBayes>;
    //apfn_RecreateWordBag[eMBayes] = &CBoWWordBag::RecreateWordBag<eMBayes>;
    //apfn_RecreateWordBag[eLambda] = &CBoWWordBag::RecreateWordBag<eLambda>;
    apfn_RecreateWordBags[(int) CBOWParams::eNisterDist] = &CBoW::RecreateWordBags<CBOWParams::eNisterDist>;
    apfn_RecreateWordBags[(int) CBOWParams::eVectorDistFast] = &CBoW::RecreateWordBags<CBOWParams::eVectorDistFast>;
    apfn_RecreateWordBags[(int) CBOWParams::eBruteForce] = &CBoW::RecreateWordBags<CBOWParams::eBruteForce>;

}

void CBoW::SetupMinMatchThreshholds(CBOWParams::eCOMP_METHOD eCompMethod) {
    const int LEVELS = PARAMS.BOWClustering.LEVELS;
    anMinMatchStrength = new int[LEVELS];

    switch (eCompMethod) {
            DEPRECATED(
        case eLogMBayes:
            anMinMatchStrength[LEVELS - 1] = doubleToInt(dMinMatchStrength * eCompMethod);
            for (int i = LEVELS - 2; i >= 0; i--)
                    //anMinMatchStrength[i]=dMinMatchStrength_in*anMinMatchStrength[i+1];
                    anMinMatchStrength[i] = doubleToInt(0.75 * anMinMatchStrength[i + 1]);
                break;)
                case CBOWParams::eNisterDist:
                case CBOWParams::eVectorDistFast:
                for (int i = 0; i < LEVELS - 1; i++)
                    anMinMatchStrength[i] = -3000000 * 100;
            anMinMatchStrength[LEVELS - 1] = MIN_INT;

            break;
        default:
            memset(anMinMatchStrength, 0, LEVELS * sizeof (int));
    }
}

//Internal

CBoW::CBoWWord * CBoW::CBoWDictionary::ClosestWord(const CDescriptor * pDescriptor, CDescriptor::TDist dShortestDistance) const {
    wordNum nClosestWord = INVALID_WORD;

    //Iterate through words, choose one with closest centre
    for (wordNum nWord = 0; nWord < (wordNum) Words.size(); nWord++) {
        CDescriptor::TDist dThisDist = Words[nWord]->Descriptor()->distance(pDescriptor);
        if (dThisDist < dShortestDistance) {
            dShortestDistance = dThisDist;
            nClosestWord = nWord;
        }

        static int nShow100 = 5;
        if (nShow100 > 0) {
            if (nShow100 == 5)
                cout << "Distances: ";
            cout << dThisDist << ' ';
            if (nShow100 == 1) cout << '\n';

            nShow100--;
        }
    }
    if (nClosestWord == INVALID_WORD) return 0;
    //if(IS_DEBUG) CHECK(nClosestWord==INVALID_WORD, "CBoW::CBoWDictionary::ClosestWord: No word found");

    return Words[nClosestWord]; //Used to Need to add an offset so have unique indices for all words at this level
}

CBoW::CBoWDictionary::~CBoWDictionary() {
    if (Words.size()) {
        CDynArray<CBoWWord *>::iterator ppEnd = Words.end();
        for (CDynArray<CBoWWord *>::iterator ppWord = Words.begin(); ppWord < ppEnd; ppWord++) {
            if (nLevel == 1)
                delete *ppWord;
            else {
                CBoWNodeWord * pNodeWord = static_cast<CBoWNodeWord *> (*ppWord);
                delete pNodeWord;
            }
        }
    }

    delete [] anWordLevelCounts;
    delete pClusters;
}

void CBoW::recreateDictionary() {
    ClusterDescriptorsIntoWords();
    //Clustering HASN'T NECESSARILY HAPPENED YET--may be happening in a seperate thread.
}

/*CBoW::CBoWWordBag::TWordSet * CBoW::CBoWWordBag::GetWordSet(int LEVELS)
{
    CBoW::CBoWWordBag::TWordSet * aWordSet = new CBoW::CBoWWordBag::TWordSet[LEVELS];
    return aWordSet;
}

void CBoW::CBoWWordBag::DeleteWordSet(CBoW::CBoWWordBag::TWordSet * aWordSet, int LEVELS)
{
    delete [] aWordSet;
}*/

void CBoW::CBoWWordBag::RecreateWordBag() {
    typedef set2<CBoW::CBoWWordBag::CBoWImageWord, CBoW::CBoWWordBag::WordPointerSort> TWordSet;

    if(IS_DEBUG) CHECK(!pParent->pDictionary, "CBoW::CBoWWordBag::RecreateWordBag: Dictionary not setup yet");

    const CBoW::CBoWDictionary * pDic = pParent->pDictionary;
    if(IS_DEBUG) CHECK(!pDic, "CBoW::CBoWWordBag::RecreateWordBag: Cannot build word bag--no words created yet");

    const int LEVELS = pParent->PARAMS.BOWClustering.LEVELS;
    const int BF_LEVEL = pParent->PARAMS.BOWClustering.BRUTEFORCE_MATCHING_LEVEL();

    //   //Re-map descriptors to words and insert into set
    ARRAY(CBoWWord *, apWords, LEVELS);
    ARRAY(TWordSet, aWordSet, LEVELS);

    ARRAY(const CDescriptor *, apClosestDescriptors, LEVELS);

    aDescriptorWordsBF.clear();
    const CBOWParams::CDescriptorBinningParams & pBinParams = pParent->PARAMS.DescriptorBinning;
    for (int i = 0; i < pDescriptors->Count(); i++) {
        const CDescriptor * pDescriptor = pDescriptors->get_const(i);

        setZero(PTR(apClosestDescriptors), LEVELS);

        CCluster const * pCluster = pDescriptor->assignment();
        if (pCluster) {
            int nChainLength = 1;
            while (pCluster->parentCluster()) {
                pCluster = pCluster->parentCluster();
                nChainLength++;
            }

            pCluster = pDescriptor->assignment();
            for (int i = nChainLength - 1; pCluster; i--) {
                if(IS_DEBUG) CHECK(!pCluster->Centre(), "Cluster has no centre");
                apClosestDescriptors[i] = pCluster->Centre();
                pCluster = pCluster->parentCluster();
            }
        }

        pDic->LookupWordAllLevels(pDescriptor, LEVELS - 1, PTR(apWords), PTR(apClosestDescriptors), pBinParams);

        if (apWords[LEVELS - 1]) //We might not find a place for this word if it's too far from any existing words, or in too big/small a cluster
        {
            for (int nLevel = 0; nLevel < LEVELS; nLevel++) {
                CLocation loc;
                if (nLevel == LEVELS - 1 /*|| nLevel == LEVELS-2*/) //### If we add locations at other levels, also need to delete them! (see below...)
                    loc = pDescriptor->location(); //otherwise will be 0

                CBoW::CBoWWordBag::CBoWImageWord myImWord(apWords[nLevel], loc); //has 1 location exactly
                pair < TWordSet::iterator, bool> insertResult = aWordSet[nLevel].insert(myImWord);
                if (!insertResult.second) {
                    const CBoW::CBoWWordBag::CBoWImageWord * imWord = &(*insertResult.first);
                    CBoW::CBoWWordBag::CBoWImageWord * imWord2 = const_cast<CBoW::CBoWWordBag::CBoWImageWord *> (imWord);
                    imWord2->IncrementFrequency(loc);
                }
            }
            if (BF_LEVEL > 0) {
                CBoWWordDescMatch match(apWords[BF_LEVEL - 1], i);
                aDescriptorWordsBF.push_back(match);
            }
            DEBUGONLY(
                    g_nWordsInsideRad++;
        } else {
            g_nWordsOutsideRad++);
        }
    }

    if (aWordSet[LEVELS - 1].size() == 0) {
        cout << "No words in word bag! Meaningless to represent this image so will now be ignored (TODO--check this)\nCheck radius is sensible size\n";
    }

    if (BF_LEVEL > 0) // Every word bag has this sorted lookup list
        std::sort(aDescriptorWordsBF.begin(), aDescriptorWordsBF.end(), CBoWWordDescMatch::CBoWWordDescMatchSort());

    //Now size a safeVector appropriately and copy in the sorted words and
    for (int nLevel = 0; nLevel < LEVELS; nLevel++) {
        TWordSet & wordset = aWordSet[nLevel];

        //Exactly resize the safeVector
        aWords[nLevel].clear(); //todo check max_size isn't excessively big

        aWords[nLevel].reserve(wordset.size());

        for (TWordSet::iterator pImWord = wordset.begin(); pImWord != wordset.end(); pImWord++)
            aWords[nLevel].push_back(*pImWord); //Location ptr (or single location) should now be moved here.

        //Triggers 'unsafe' warning with bare ptr: copy(pWordSet->begin(), pWordSet->end(), aWords[nLevel].begin()); //Locations might be duplicate here
        wordset.clear(); //Locations might be duplicate here

        //if(IS_DEBUG) CHECK((*(pWordSet->begin()))->Word() > (*(pWordSet->end()-1))->Word(), "CBoW::CBoWWordBag::RecreateWordBag: Sort error");
        if(IS_DEBUG) CHECK(aWords[nLevel].size() > 1 && aWords[nLevel][0].Word() > aWords[nLevel][1].Word(), "CBoW::CBoWWordBag::RecreateWordBag: Sort error");
    }

    cout << "WB count: " << aWords[LEVELS - 1].size() << endl;
}

templateCompMethod
void CBoW::CBoWWordBag::WeightWordBag()//TODO duplication
{
    const int LEVELS = pParent->PARAMS.BOWClustering.LEVELS;
    if ((eCompMethod == CBOWParams::eNisterDist || eCompMethod == CBOWParams::eVectorDistFast)) {
        //Occurances should already be counted globally, possibly not for this WB tho

        //cout << "Normalising word vectors...\n";
        int nLogTotal = intintLog(pParent->vImages.Count());
        //cout << "pParent->vImages.Count()=" << pParent->vImages.Count() << endl;
        int nLogAllWordCount = intintLog(pParent->vImages.totalDescriptorCount()); //Occasionallly wrong (? ? when erasing) but shouldn't be used anyway...
        //cout << pParent->nAllDescriptorsCount << " descriptors, log= " << nLogAllWordCount << endl;
        for (int nLevel = LEVELS - 1; nLevel < LEVELS; nLevel++) {
            int nNormalisingConst = 0; //can be 0 if all words are in all images
            TWordBag::iterator pWEnd = aWords[nLevel].end();
            for (TWordBag::iterator pWord = aWords[nLevel].begin(); pWord != pWEnd; pWord++) {
                int nWordWeight = -1, nTotalOccurances = -1;
                switch (pParent->PARAMS.WB_WEIGHT_METHOD) {
                    case CBOWParams::eTF_IDF_Wikipedia:
                        //Word frequency/total words in this im   *   log(total # images / images containing this word)
                        //Difference between this and Nister's is the normalising term.
                    case CBOWParams::eTF_IDF:
                        nTotalOccurances = pWord->Word()->TotalOccurances();
                        //cout << nTotalOccurances << "...";
                        if (nTotalOccurances <= 0) {
                            nTotalOccurances = 1;
                            cout << "WARNING: Word occurances not counted--probably exception thrown or threading issue or images deleted. Assuming 1\n";
                        }
                        nWordWeight = nLogTotal - intintLog(nTotalOccurances); //TF-IDF weight, ala Nister-Stewenius-2006
                        if (nWordWeight < 0) nWordWeight = 0; //This happens when images have been deleted--word weights aren't updated until later
                        //cout << nWordWeight << endl;
                        break;
                    case CBOWParams::eDF_ITDF:
                        nWordWeight = nLogAllWordCount - intintLog(pWord->Word()->TotalFrequency());
                        break;
                    default:
                        THROW("Unhandled word weighting method");
                }

                int nWordScore = pWord->FrequencyInImage() * nWordWeight;
                nNormalisingConst += (pParent->PARAMS.WB_WEIGHT_METHOD == CBOWParams::eTF_IDF_Wikipedia) ? pWord->FrequencyInImage() : nWordScore;
                pWord->SetWeight(nWordScore);
            }
            /*For TF-IDF nWordScore = 100 log( # images / # occurances ) * Freq
             *                     = about 100 * 10 * 20 = 20 000 max
             *                     so NISTERVECTOR_LENGTH should be < 2 000 000 000/20 000 = 100 000
             */

            CHECK(nNormalisingConst < 0, "WeightWordBag: Overflow");

            if (nNormalisingConst == 0) {
                //Every word is uninformative
                nNormalisingConst = 1;
                if (pParent->vImages.Count() > 1) {
                    CTSOut cout;
                    cout << "WARNING: Every word appears in every image! Word bag cannot be normalised. This is ok as long as we only have 1 image.\n";
                    cout << aWords[nLevel].size() << " words in this word bag\n";
                    cout << pParent->vImages.Count() << " total images\n";
                }
            } else {
                DEBUGONLY(int nZeroWords = 0);
                //nNormalisingConst is sum-of-weights so NISTERVECTOR_LENGTH>>num descriptors will keep numbers here working well.
                for (TWordBag::iterator pWord = aWords[nLevel].begin(); pWord != pWEnd; pWord++) {
                    int nNormalisedWeight = ((pWord->Weight() * NISTERVECTOR_LENGTH) / nNormalisingConst);

                    if(IS_DEBUG) CHECK(nNormalisedWeight < 0, "CBoW::CBoWWordBag::RecreateWordBag: Overflow");

                    if (pWord->Weight() > 0) {
                        if (nNormalisedWeight <= 0) {
                            cout << "Warning: underflow calculating word weight\n";
                            cout << pWord->Weight() << "= word weight, " << nNormalisingConst << " = normalising const\n";
                        } else if (nNormalisedWeight <= 3)
                            cout << "Warning: poor precision calculating word weight\n";
                    }
                    DEBUGONLY(
                    else {
                        nZeroWords++;
                                if(IS_DEBUG) CHECK((int) pWord->Word()->TotalOccurances() > (int) pParent->vImages.size_th(), "Error counting occurances");
                    });

                    nNormalisedWeight = max<int>(1, nNormalisedWeight);
                    pWord->SetWeight(nNormalisedWeight);
                }

                if (pParent->vImages.Count() > 2) {
                    DEBUGONLY(if (nZeroWords > 1) cout << nZeroWords << " of " << aWords[nLevel].size() << " have zero weight (occur as many times than there are images)\n");
                }
            }
        }

    }
}

//Is it a valid probability in [0,1]?

bool isProb(double p) {
    return p >= 0.0 && p <= 1.0;
}

//int CBoW::s_nLevels = -1;

CBoW::CBoW(const CBOWParams & PARAMS) : PARAMS(PARAMS), anMinMatchStrength(0), pDictionary(0), pClusterThread(0), bDeleting(false), pDrawCS(0), pDrawCS_temp(0) {
    WRITE_LOCK;
    try {
        //TODO: Make library functions
        CBoW::setupIntLogCache();
        CBoW::setupIntIntLogCache();

        SetupMinMatchThreshholds(PARAMS.COMP_METHOD);

        CreateTemplateFns();
        pfn_WeightWordBag = apfn_WeightWordBag[PARAMS.COMP_METHOD];
        pfn_getBoWMatches = apfn_getBoWMatches[PARAMS.COMP_METHOD];
        pfn_RecreateWordBags = apfn_RecreateWordBags[PARAMS.COMP_METHOD];

        nNextClusterCount = 1; //Cluster as soon as we have some descriptors

        bClustering = false;
    } catch (...) {
        RELEASE_WRITE_LOCK;
        throw;
    }
    RELEASE_WRITE_LOCK;
}

CBoW::~CBoW() {
    cout << "Deleting CBoW...";
    bDeleting = true;
    MX_NO_DELETE_WHILE_CLUSTER_IN_SEPERATE_THREAD;

    if (pClusterThread) {
        //sleep(1);
        if (bClustering) //NOT TS
        {
            cout << "rejoining cluster thread...";
            pClusterThread->join();
            cout << "Done";
        }
        delete pClusterThread;
        pClusterThread = 0;
        cout << "Deleted cluster thread\n";
        cout.flush();
    }

    WRITE_LOCK; //either order could deadlock, not quite safe now but should be ok...
    cout << "...locks obtained\n";

    cout << "Deleting dictionary...";
    delete pDictionary;
    pDictionary = 0;
    cout << "done deleting anMinMatchStrength...";
    delete [] anMinMatchStrength;
    anMinMatchStrength = 0;
    cout << "done" << endl;

    RELEASE_WRITE_LOCK;
}

templateCompMethod
void CBoW::CBoWWordBag::getBoWMatches_Loop(TBoWMatchVector *pvMatches, int nReturnMax, const int * anScoreAgainstBackground, constImIt imageIterBegin, constImIt imageIterEnd) {
    try {
        getBoWMatches_Loop_int<eCompMethod > (pvMatches, nReturnMax, anScoreAgainstBackground, imageIterBegin, imageIterEnd);
    } catch (CException pEx) {
        pWBException = pEx;
    } catch (...) {
        pWBException = CException("Unknown exception type caught in query thread");
    }
}

templateCompMethod
void CBoW::CBoWWordBag::getBoWMatches_Loop_int(TBoWMatchVector *pvMatches, int nReturnMax, const int * anScoreAgainstBackground, constImIt imageIterBegin, constImIt imageIterEnd) {
    const int LEVELS = pParent->PARAMS.BOWClustering.LEVELS;
    int nNoMatch = pParent->anMinMatchStrength[LEVELS - 1];

    TSortedIntSet * pSubvecScores = (eCompMethod == CBOWParams::eVectorDistFast) ? new TSortedIntSet() : 0;
    int nVecReturnMax = (pParent->PARAMS.FastVectorComparison.SUBVEC_CACHE_SIZE * nReturnMax) / pParent->PARAMS.QUERY_THREADS;

    for (constImIt imageIter = imageIterBegin; imageIter != imageIterEnd; imageIter++) {
        CBoWWordBag * pCompareAgainst = *imageIter;
        if (!pCompareAgainst) continue; //has been erased

        TWordBag & pBottomLevelWords = pCompareAgainst->aWords[LEVELS - 1];
        if (!pBottomLevelWords.size()) continue; //Don't match against images with empty word bags

        //CBoWMatch * pMatch;
        int nMatch;
        if (eCompMethod == CBOWParams::eBruteForce)
            nMatch = BruteForceMatch(pCompareAgainst);
        else
            nMatch = getBoWMatch<eCompMethod > (pCompareAgainst, anScoreAgainstBackground, pSubvecScores, nVecReturnMax);

        if (nMatch > nNoMatch) {
            boost::mutex::scoped_lock scoped_lock(mxQuery);
            TBoWMatchVector::iterator last = pvMatches->end();
            if ((int) pvMatches->size() < nReturnMax || (nMatch > (--last)->MatchStrength())) {
                if ((int) pvMatches->size() >= nReturnMax)
                    pvMatches->erase(last);

                CBoWMatch<int> Match(pCompareAgainst->id(), nMatch);
                pvMatches->insert(Match);

                if ((int) pvMatches->size() >= nReturnMax)
                    nNoMatch = (--(pvMatches->end()))->MatchStrength();
            }
        }
    }
    delete pSubvecScores;
}

template<bool DECREMENT>
void CBoW::CBoWWordBag::CountWordOccurances() const {
    const int LEVELS = pParent->PARAMS.BOWClustering.LEVELS;
    //Increment occurance count of each word this WB contains (could build Inverse File here too)
    for (int nLevel = LEVELS - 1; nLevel < LEVELS; nLevel++) {
        TWordBag::iterator ppEnd = aWords[nLevel].end();
        for (TWordBag::iterator ppWord = aWords[nLevel].begin(); ppWord < ppEnd; ppWord++) {
            if (DECREMENT)
                ppWord->Word()->DecrementOccurances();
            else
                ppWord->Word()->IncrementOccurances();
        }
    }
}

templateCompMethod
TBoWMatchVector * CBoW::CBoWWordBag::getBoWMatches(int nReturnMax) {
    TBoWMatchVector * pvMatches = new TBoWMatchVector();

    const int LEVELS = pParent->PARAMS.BOWClustering.LEVELS;

    if (!pParent->pDictionary || !this->aWords[LEVELS - 1].size() || nReturnMax <= 0) //If the wordbag is empty return no matches
        return pvMatches;

    ARRAY(int, anScoreAgainstBackground, LEVELS);

    //First calculate the score against the 'Background' category
    for (int nLevel = 0; nLevel < LEVELS; nLevel++) {
        if ((eCompMethod == CBOWParams::eNisterDist || eCompMethod == CBOWParams::eVectorDistFast) || eCompMethod == CBOWParams::eBruteForce)
            anScoreAgainstBackground[nLevel] = 0;
        else
            DEPRECATED(anScoreAgainstBackground[nLevel] = CompareWithAll<eCompMethod > (aWords[nLevel], nTotalDescriptorCount););

        //cout << anScoreAgainstBackground[nLevel] << "=background score, level=" << nLevel << '\n' ;
    }

    //This loop is the outer loop around the critical section
    if (pParent->PARAMS.QUERY_THREADS > 1) {
        const int NUM_IMAGES_PER_THREAD = (int) (pParent->vImages.size_th()) / pParent->PARAMS.QUERY_THREADS;
        //		ARRAY(boost::thread *, apQueryThreads, (const int)pParent->PARAMS.QUERY_THREADS); //TODO: array of smart ptrs
        CDynArrayOwner<boost::thread> apQueryThreads(pParent->PARAMS.QUERY_THREADS);

        for (int nThread = 0; nThread < pParent->PARAMS.QUERY_THREADS; nThread++) {
            constImIt ImStart = pParent->vImages.begin() + (nThread * NUM_IMAGES_PER_THREAD);
            constImIt ImEnd = ImStart + NUM_IMAGES_PER_THREAD;
            if (nThread == pParent->PARAMS.QUERY_THREADS - 1) ImEnd = pParent->vImages.end();

            apQueryThreads[nThread] = new boost::thread(boost::bind(&CBoW::CBoWWordBag::getBoWMatches_Loop<eCompMethod>,
                    this, pvMatches, nReturnMax, PTR(anScoreAgainstBackground), ImStart, ImEnd));
        }
        for (int nThread = 0; nThread < pParent->PARAMS.QUERY_THREADS; nThread++) {
            apQueryThreads[nThread]->join();
            //delete apQueryThreads[nThread];
        }
    } else {
        constImIt imageIterBegin = pParent->vImages.begin();
        constImIt imageIterEnd = pParent->vImages.end();
        getBoWMatches_Loop<eCompMethod > (pvMatches, nReturnMax, PTR(anScoreAgainstBackground), imageIterBegin, imageIterEnd);
    }

    checkWBException();

    //DEBUGONLY(
    //    if(pvMatches->size()>1)
    //    {
    //        CHECK((*pvMatches)[0]->MatchStrength() <= (*pvMatches)[pvMatches->size()-1]->MatchStrength(), "CBoW::CBoWWordBag::getBoWMatches: Match strengths equal or not sorted");
    //    }
    //);

    return pvMatches;
}

inline double CBoW::intLog(int i) {
    if(IS_DEBUG) CHECK(i <= 0, "CBoW::intLog: Can't log negative number");
    return (i < INT_LOG_CACHE_LIM) ? adIntLogCache[i] : log((double) i);
}

//same as intlog but returns int

inline int CBoW::intintLog(int i) {
    if(IS_DEBUG) CHECK(i <= 0, "CBoW::intintLog: Can't log negative number"); //This is tripped when a word isn't counted
    if (i < INTINT_LOG_CACHE_LIM) {
        COUNTHITS(nHit++);
        return anIntLogCache[i];
    } else if (i == nCache_i) {
        COUNTHITS(nCache++);
        //cout << 'h';
        return nLog_i;
    } else {
        COUNTHITS(nMiss++);
        int nLog = doubleToInt(INTINTLOG_SCALE * log((double) i));
        //cout << 'm';
        nCache_i = i;
        nLog_i = nLog; //Should really be atomic/threadsafe
        return nLog;
    }
}

double CBoW::intPow(double d, int i) {
    double x = 1;
    unsigned int bit = 1;
    if (i) {
        for (;;) {
            if (i & bit) {
                x *= d;
                i ^= bit;
                if (!i) return x;
            }
            d *= d;
            bit <<= 1;
        }
    }
    return x;
}
//Should we adjust frequency of words that occur once only, to avoid huge update factors (factor=number of vImages)?
//We need a Laplace m-estimate apparently. This sounds about the same as what I've done.

double getBayesUpdateFast(int nTotalWords, int nImageWords, int nWordFrequencyTotal, int nWordFrequencyImage) {
    //Update factor = P(this word|image)/P(this word)
    //Use P(this word|image)=(1-lambda)P(this word|image)+lambda*P(this word) as we expect lambda of our words to not occur in the original image
    const double lambda = 0.4;

    //Faster version than below (1 division)
    return lambda + ((1 - lambda) * nWordFrequencyImage * nTotalWords) / (nWordFrequencyTotal * nImageWords);
}

double getBayesUpdate(int nTotalWords, int nImageWords, int nWordFrequencyTotal, int nWordFrequencyImage) {
    //Update factor = P(this word|image)/P(this word)
    //Use P(this word|image)=(1-lambda)P(this word|image)+lambda*P(this word) as we expect lambda of our words to not occur in the original image

    const double lambda = 0.4;

    double dPw_givenFromImage = nWordFrequencyImage / (double) nImageWords;
    double dPw = nWordFrequencyTotal / (double) nTotalWords;

    return ((1 - lambda) * dPw_givenFromImage + lambda * dPw) / dPw;
}

double CBoW::getBayesPosteriorFromPriorWithLambda(double dPrior, int nImageWords, int nWordFrequencyTotal, int nWordFrequencyImage, int nWordFrequencyNewImage) {
    //See Maple: Naive Bayes time-independence justification.mws

    //Update factor = P(this word|image)/P(this word|prior assumption)
    //Use P(this word|image)=(1-lambda)P(this word|image)+lambda*P(this word) as we expect lambda of our words to not occur in the original image
    DEPRECATED(
            const double lambda = dProportionBackground; //proportion of words from background that we expect to see
            const unsigned int nTotalWords = pAllDescriptors->Count();

            double dPw_givenFromImage = nWordFrequencyImage / (double) nImageWords;
            double dPw = nWordFrequencyTotal / (double) nTotalWords;
            double dPw_givenFromImageAndBackground = dPw_givenFromImage * (1 - lambda) + dPw*lambda;
            double dNumerator = dPrior * intPow(dPw_givenFromImageAndBackground, nWordFrequencyNewImage);
            double dPw_givenPrior = (1 - dPrior) * intPow(dPw, nWordFrequencyNewImage) + dNumerator;

    return dNumerator / dPw_givenPrior;);
}

double CBoW::getBayesPosteriorFromPriorMEstimate(double dPrior, int nImageWords, int nWordFrequencyTotal, int nWordFrequencyImage, int nWordFrequencyNewImage) const {
    //See Maple: Naive Bayes time-independence justification.mws

    //Update factor = P(this word|image)/P(this word|prior assumption)
    //Use P(this word|image)=(1-lambda)P(this word|image)+lambda*P(this word) as we expect lambda of our words to not occur in the original image
    DEPRECATED(
            const double lambda = dProportionBackground; //proportion of words from background that we expect to see
            const unsigned int nTotalWords = pAllDescriptors->Count();

            double dPw_givenFromImage = nWordFrequencyImage / (double) nImageWords;
            double dPw = nWordFrequencyTotal / (double) nTotalWords;
            double dPw_givenFromImageAndBackground = dPw_givenFromImage * (1 - lambda) + dPw*lambda;
            double dNumerator = dPrior * intPow(dPw_givenFromImageAndBackground, nWordFrequencyNewImage);
            double dPw_givenPrior = (1 - dPrior) * intPow(dPw, nWordFrequencyNewImage) + dNumerator;

    return dNumerator / dPw_givenPrior;);
}

/*double CBoW::getLogBayesPosteriorFromPriorMEstimate(double dLogPrior, int nImageWords, int nWordFrequencyTotal, int nWordFrequencyImage, int nWordFrequencyNewImage) const
{
    //See Maple: Naive Bayes time-independence justification.mws

    //Update factor = P(this word|image)/P(this word|prior assumption)
    //Use P(this word|image)=(1-lambda)P(this word|image)+lambda*P(this word) as we expect lambda of our words to not occur in the original image

    const double lambda=dProportionBackground; //proportion of words from background that we expect to see
    const unsigned int nTotalWords = pAllDescriptors->Count();

    double dLogPw_givenFromImage = intLog(nWordFrequencyImage) - intLog(nImageWords);
    double dLogPw = intLog(nWordFrequencyTotal) - intLog(nTotalWords);
    double dPw_givenFromImageAndBackground = dLogPw_givenFromImage*(1-lambda)+dPw*lambda;
    double dNumerator = dPrior*intPow(dPw_givenFromImageAndBackground, nWordFrequencyNewImage);
    double dPw_givenPrior = (1-dPrior)*intPow(dPw, nWordFrequencyNewImage) + dNumerator;

    return dNumerator/dPw_givenPrior;
}*/

/*double CBoW::getBayesPFromPMEstimate(double dPrior, int nImageWords, int nWordFrequencyTotal, char nWordFrequencyImage, char nWordFrequencyNewImage) const
{
    //Add 1 ('dummy') to cat word counts. Add n to total cat word count
    const double m=1.0; //dummy for m estimate
    const unsigned int nTotalWords = pAllDescriptors->Count();

    double dPw_givenFromImage = (nWordFrequencyImage+1)/(double)(nImageWords+nTotalWords);
    double dPw = nWordFrequencyTotal/(double)nTotalWords;
    double dNumerator = dPrior*intPow(dPw_givenFromImage, nWordFrequencyNewImage);
    double dPw_givenPrior = (1-dPrior)*intPow(dPw, nWordFrequencyNewImage) + dNumerator;

    return dNumerator/dPw_givenPrior;
}*/

double getBayesUpdateFromPriorWithLambda(double dPrior, int nTotalWords, int nImageWords, int nWordFrequencyTotal, int nWordFrequencyImage) {
    //Update factor = P(this word|image)/P(this word|prior assumption)
    //Use P(this word|image)=(1-lambda)P(this word|image)+lambda*P(this word) as we expect lambda of our words to not occur in the original image

    const double lambda = 0.4; //proportion of words from background that we expect to see

    double dPw_givenFromImage = nWordFrequencyImage / (double) nImageWords;
    double dPw = nWordFrequencyTotal / (double) nTotalWords;
    double dPw_givenFromImageAndBackground = dPw_givenFromImage * (1 - lambda) + dPw*lambda;
    double dPw_givenPrior = (1 - dPrior) * dPw + dPrior*dPw_givenFromImageAndBackground;

    return dPw_givenFromImageAndBackground / dPw_givenPrior;
}

double getBayesUpdateFromPrior(double dPrior, int nTotalWords, int nImageWords, int nWordFrequencyTotal, int nWordFrequencyImage) {
    //Update factor = P(this word|image)/P(this word|prior assumption)

    double dPw_givenFromImage = nWordFrequencyImage / (double) nImageWords;
    double dPw = nWordFrequencyTotal / (double) nTotalWords;
    double dPw_givenPrior = (1 - dPrior) * dPw + dPrior*dPw_givenFromImage;

    return dPw_givenFromImage / dPw_givenPrior;
}

//L0 is the top of the tree

templateCompMethod
int CBoW::CBoWWordBag::getBoWMatch(const CBoWWordBag * pwbExisting, const int * anScoreAgainstBackground, TSortedIntSet * pSubvecScores, int nVecReturnMax) const {
    if(IS_DEBUG) CHECK(!pwbExisting, "CBoW::getBoWMatch: Bad parameter");

    int p = 0;
    int * anMinMatchStrength = pParent->anMinMatchStrength;
    int LEVELS = pParent->PARAMS.BOWClustering.LEVELS;

    for (int nLevel = LEVELS - 1; nLevel < LEVELS; nLevel++) //TODO go back to higher levels (??)
    {
        if (eCompMethod == CBOWParams::eNisterDist)
            p = VectorCompare(pwbExisting, nLevel);
        else if (eCompMethod == CBOWParams::eVectorDistFast)
            p = VectorCompareFast(pwbExisting, nLevel, pSubvecScores, nVecReturnMax);
        else {
            DEPRECATED(
                    bool bUseLogs = (eCompMethod == eLogBayes || eCompMethod == eLogMBayes);

                    p = Compare<eCompMethod > (aWords[nLevel], pwbExisting->aWords[nLevel], pwbExisting->TotalWordCount());

            if (bUseLogs)
                    p -= anScoreAgainstBackground[nLevel];
            else
                p /= anScoreAgainstBackground[nLevel];)
            }

#ifndef __GNUC__
        static int nPrintScoreCount = 10;

        if (nPrintScoreCount > 0 && nLevel < 2) {
            nPrintScoreCount--;
            char pc[50];
            sprintf_s(pc, 50, "%d: %d\n", nLevel, p);
            cout << pc;
        }
#endif

        //if(nLevel > -1) cout << p << "=p level=" << nLevel << '\n' ;
        if (anMinMatchStrength[nLevel] > p) return anMinMatchStrength[LEVELS - 1];
    }

    return p;
}

/*inline double CBoW::getBayesScoreM(int nTimesWordAppearsInCat, double dTotalWordsInCat_inv)
{
    //Return prob of seeing this word in this category
    return (nTimesWordAppearsInCat+1)*dTotalWordsInCat;
}*/

inline double CBoW::getLogBayesScoreM(int nTimesWordAppearsInCat, int nTotalWordsInCat) {
    //Return prob of seeing this word in this category
    return CBoW::intLog(nTimesWordAppearsInCat + 1) - CBoW::intLog(nTotalWordsInCat);
}

inline double CBoW::getBayesScoreM(int nTimesWordAppearsInCat, int nTotalWordsInCat) {
    //Return prob of seeing this word in this category
    return (nTimesWordAppearsInCat + 1) / (double) (nTotalWordsInCat);
}

/*templateCompMethod
inline int CBoW::CBoWWordBag::CompareWithAll(TWordBag * pWordsInImage, int nWords) const
{
        DEPRECATED(
        double p=0;
    if(IS_DEBUG) CHECK(!pWordsInImage, "CBoW::CBoWWordBag::CompareWithAll: Bad parameters");

    bool bUseLogs = (eCompMethod==eLogBayes || eCompMethod==eLogMBayes);

    if(bUseLogs)
        p = log(1 - pParent->dPriorProbOfAnyMatch) ; //=prior probability we're going to match this image; //TODO no log
    else
        p = 1 - pParent->dPriorProbOfAnyMatch ; //=prior probability we're going to match this image; //TODO no log

    TWordBag::iterator pImWordInImEnd = pWordsInImage->end();

    for(TWordBag::iterator pImWordInIm = pWordsInImage->begin();
        pImWordInIm != pImWordInImEnd;
        pImWordInIm++)
    {
        const CBoWWord * pWordInIm = pImWordInIm->Word();

        switch(eCompMethod)
        {
        case eLambda://NB Repeated params TODO
            p = pParent->getBayesPosteriorFromPriorMEstimate(p, nWords, pWordInIm->TotalFrequency(), pWordInIm->TotalFrequency(), pImWordInIm->FrequencyInImage());
            break;
        case eLogMBayes:
        case eLogBayes:
            p += pImWordInIm->FrequencyInImage() * getLogBayesScoreM(pWordInIm->TotalFrequency(), nWords);
            break;
        case eMBayes:
            p *= CBoW::intPow(getBayesScoreM(pWordInIm->TotalFrequency(), nWords), pImWordInIm->FrequencyInImage());
            break;
        case CBOWParams::eNisterDist:
        case CBOWParams::eVectorDistFast:
        case CBOWParams::eBruteForce:
        default:
            throw new CBoWException("Shouldn't call CompareWithAll when safeVector comparison used");
        }
    }
    return doubleToInt(p*1000);
        )
}*/

inline int CBoW::CBoWWordBag::VectorCompare_int(TWordBag::iterator catWordList, TWordBag::iterator catWordListEnd, TWordBag::iterator imageWordList, TWordBag::iterator imageWordListEnd) const {
    const CBoWWord * pWordInIm = imageWordList->Word();
    const CBoWWord * pWordInCat = catWordList->Word();

    //    const bool bTF_IDF = (pParent->eVecWeightMethod==CBOWParams::eTF_IDF);

    //    int nLogCount = intintLog(bTF_IDF ? pParent->vImages.Count() : pParent->pAllDescriptors->Count());
    int p = 0;
    bool bImageListEnd = false;
    bool bCatListEnd = false;
    //This is the inner critical section
    do {
        if ((pWordInIm > pWordInCat || bImageListEnd) && !bCatListEnd) {
            //the word in the cat but not the im:
            //int nWordWeight = nLogCount - intintLog(bTF_IDF ? pWordInCat->TotalOccurances() : pWordInCat->TotalFrequency());

            p += /*nWordWeight*/catWordList->Weight();

            catWordList++;
            if (catWordList != catWordListEnd) {
                pWordInCat = catWordList->Word();
            } else
                bCatListEnd = true;
        } else if ((pWordInIm < pWordInCat || bCatListEnd) && !bImageListEnd) {
            //the word in the im but not the cat:
            //int nWordWeight = nLogCount - intintLog(bTF_IDF ? pWordInIm->TotalOccurances() : pWordInIm->TotalFrequency());

            p += /*nWordWeight*/imageWordList->Weight();

            imageWordList++;
            if (imageWordList != imageWordListEnd) {
                pWordInIm = imageWordList->Word();
            } else
                bImageListEnd = true;
        } else {
            //int nWordWeight = nLogCount - intintLog(bTF_IDF ? pWordInIm->TotalOccurances() : pWordInIm->TotalFrequency());

            p += /*nWordWeight*/ abs(imageWordList->Weight()
                    - catWordList->Weight());

            imageWordList++;
            if (imageWordList != imageWordListEnd) {
                pWordInIm = imageWordList->Word();
            } else
                bImageListEnd = true;

            catWordList++;
            if (catWordList != catWordListEnd) {
                pWordInCat = catWordList->Word();
            } else
                bCatListEnd = true;
        }

    } while (!(bCatListEnd && bImageListEnd));
    return -p;
}

//Faster version using fact that SUM(xi)=SUM(yi)=const
//Returns 2const-VectorCompare_int

template<class TImWordIterator>
inline int CBoW::CBoWWordBag::NormalisedVectorCompare_int(TImWordIterator catWordList, const TImWordIterator catWordListEnd, TImWordIterator imageWordList, const TImWordIterator imageWordListEnd) const {
    //TImWordIterator catWordList = *pcatWordList;
    //TImWordIterator imageWordList = *pimageWordList;

    const CBoWWord * pWordInIm = imageWordList->Word();
    const CBoWWord * pWordInCat = catWordList->Word();

    //    const bool bTF_IDF = (pParent->eVecWeightMethod==CBOWParams::eTF_IDF);

    //    int nLogCount = intintLog(bTF_IDF ? pParent->vImages.Count() : pParent->pAllDescriptors->Count());
    int p = 0;
    bool bImageListEnd = false;
    bool bCatListEnd = false;
    //This is the inner critical section
    do {
        if ((pWordInIm > pWordInCat || bImageListEnd) && !bCatListEnd) {
            //the word in the cat but not the im:
            catWordList++;
            if (catWordList != catWordListEnd) {
                pWordInCat = catWordList->Word();
            } else
                bCatListEnd = true;
        } else if ((pWordInIm < pWordInCat || bCatListEnd) && !bImageListEnd) {
            //the word in the im but not the cat:
            imageWordList++;
            if (imageWordList != imageWordListEnd) {
                pWordInIm = imageWordList->Word();
            } else
                bImageListEnd = true;
        } else {
            //int nWordWeight = nLogCount - intintLog(bTF_IDF ? pWordInIm->TotalOccurances() : pWordInIm->TotalFrequency());

            //14-12-09: Weight includes everything
            int nImWordComponent = imageWordList->Weight();
            int nCatWordComponent = catWordList->Weight();

            p += /*nWordWeight* */ min(nImWordComponent, nCatWordComponent);

            imageWordList++;
            if (imageWordList != imageWordListEnd) {
                pWordInIm = imageWordList->Word();
            } else
                bImageListEnd = true;

            catWordList++;
            if (catWordList != catWordListEnd) {
                pWordInCat = catWordList->Word();
            } else
                bCatListEnd = true;
        }
    } while (!(bCatListEnd && bImageListEnd));

    /*/Return where iterators have got to, in case we call again for fast vector comp...
     *pcatWordList = catWordList;
     *pimageWordList = imageWordList;*/

    return p;
}

//For Nister's safeVector dist comparison

inline int CBoW::CBoWWordBag::VectorCompare(const CBoWWordBag * pWB, int nLevel) const {
    TWordBag & pWordsInImage = aWords[nLevel]; //todo rename
    TWordBag & pWordsMatchAgainst = pWB->aWords[nLevel];

    const TWordBag::iterator catWordList = pWordsMatchAgainst.begin();
    const TWordBag::iterator catWordListEnd = pWordsMatchAgainst.end();

    const TWordBag::iterator imageWordList = pWordsInImage.begin();
    const TWordBag::iterator imageWordListEnd = pWordsInImage.end();

    return fnVectorCompare(catWordList, catWordListEnd, imageWordList, imageWordListEnd);
}

inline int CBoW::CBoWWordBag::VectorCompareFast(const CBoWWordBag * pWB, int nLevel, TSortedIntSet * pSubvecScores, int nVecReturnMax) const {
    const CBOWParams::CFastVectorComparisonParams & FV_PARAMS = pParent->PARAMS.FastVectorComparison;

    TWordBag & pWordsInImage = aWords[nLevel];
    TWordBag & pWordsMatchAgainst = pWB->aWords[nLevel];

    TWordBag::iterator catWordListIt = pWordsMatchAgainst.begin();
    TWordBag::iterator imageWordListIt = pWordsInImage.begin();

    const TImWordIt catWordList = *((TImWordIt *) (void *) (&catWordListIt)); //Iterator to ptr conversion now not needed
    const TImWordIt imageWordList = *((TImWordIt *) (void *) (&imageWordListIt));

    int nImWords = (int) pWordsInImage.size(), nCatWords = (int) pWordsMatchAgainst.size();

    //cout << "Testing iterator conv: " << catWordListIt << " =? " << catWordList << endl;
    //only if enough points to make worthwhile

    if (nImWords > FV_PARAMS.VEC_MIN_LENGTH && nCatWords > FV_PARAMS.VEC_MIN_LENGTH) {
        const TImWordIt catWordListStop = catWordList + FV_PARAMS.SUBVEC_LENGTH;
        const TImWordIt imageWordListStop = imageWordList + FV_PARAMS.SUBVEC_LENGTH;

        int p_subset = fnVectorCompare(catWordList, catWordListStop, imageWordList, imageWordListStop);

        if ((int) pSubvecScores->size() < nVecReturnMax || p_subset > *(pSubvecScores->begin())) {
            pSubvecScores->insert(p_subset);
        } else return MIN_INT; //not big enough
    }

    const TImWordIt catWordListEnd = catWordList + nCatWords;
    const TImWordIt imageWordListEnd = imageWordList + nImWords;

    //Could re-use earlier score--an iterator may have been incremented past a match in the other category but is unlikely... Actually probably can if we track where the iterators have got to...
    return fnVectorCompare(catWordList, catWordListEnd, imageWordList, imageWordListEnd);

    /* Test approximation
     *
    catWordList = *((TImWordIt *)(void *)(&catWordListIt));
    imageWordList = *((TImWordIt *)(void *)(&imageWordListIt));
    int nTest = fnVectorCompare(&catWordList, catWordListEnd, &imageWordList, imageWordListEnd);
    if(nTest != p_subset)
    {
        cout << nTest << "=correct score, " << p_subset << " = approx score\n";
    } Actually not that good an approximation */
}

templateCompMethod
inline int CBoW::CBoWWordBag::Compare(TWordBag * pWordsInImage, TWordBag * pWordsMatchAgainst, int nWordsInCat) const {
    DEPRECATED(
            bool bUseLogs = (eCompMethod == eLogBayes || eCompMethod == eLogMBayes);

            double p = 0;

    if (bUseLogs)
            p = log(pParent->dPriorProbOfAnyMatch) - intLog((int) pParent->vImages.Count()); //=prior probability we're going to match this image;
    else
        p = pParent->dPriorProbOfAnyMatch / pParent->vImages.Count(); //=prior probability we're going to match this image;

            TWordBag::iterator pImWordInCat = pWordsMatchAgainst->begin();
            TWordBag::iterator pImWordInCatEnd = pWordsMatchAgainst->end();

            TWordBag::iterator pImWordInImEnd = pWordsInImage->end();

            //const CBoWImageWord * pImWordInCat = (*catWordList);
            const CBoWWord * pWordInCat = pImWordInCat->Word();

            int nFreqInCatTotal = 0;

            //This is the inner critical section
        for (TWordBag::iterator pImWordInIm = pWordsInImage->begin();
                pImWordInIm != pImWordInImEnd;
                pImWordInIm++) {

            //const CBoWImageWord * pImWordInIm = (*imageWordList);
            const CBoWWord * pWordInIm = pImWordInIm->Word();

                    //increment cat word iterator until greater than or = (or end).
            if (pImWordInCat != pImWordInCatEnd) {
                while (pWordInIm > pWordInCat) {
                    pImWordInCat++;
                    if (pImWordInCat != pImWordInCatEnd) {
                        //pImWordInCat = (*catWordList);
                        pWordInCat = pImWordInCat->Word();
                    } else
                        break;
                    }
            }

            int nFreqInCat = (pWordInIm == pWordInCat) ? pImWordInCat->FrequencyInImage() : 0;
                    nFreqInCatTotal += nFreqInCat;
                    /*cout << pWordInIm ;
                            if(nFreqInCat)
                                    cout << " Image and cat-cat freq=" << nFreqInCat << ", im freq = " << pImWordInIm->FrequencyInImage() << '\n';
                            else
                                cout << " Image but not cat, im freq = " << pImWordInIm->FrequencyInImage() << '\n';*/

            switch (eCompMethod) {
                case eLambda:
                    p = pParent->getBayesPosteriorFromPriorMEstimate(p, nWordsInCat, pWordInIm->TotalFrequency(), nFreqInCat, pImWordInIm->FrequencyInImage());
                    break;
                case eLogBayes:
                case eLogMBayes:
                    if (eCompMethod != eLogBayes || nFreqInCat)
                            p += CBoW::getLogBayesScoreM(nFreqInCat, nWordsInCat) * pImWordInIm->FrequencyInImage();
                        break;
                        case eMBayes:
                        p *= CBoW::intPow(CBoW::getBayesScoreM(nFreqInCat, nWordsInCat), pImWordInIm->FrequencyInImage());
                        break;
                        default: //CBOWParams::eNisterDist, CBOWParams::eBruteForce
                        THROW("Compare shouldn't be called--we're using safeVector distance");
                    }
        }
    if (eCompMethod == eLogBayes)
            p *= (-nFreqInCatTotal / (double) nWordsInCat);

            //cout << "Score:" << p << " Proportion matched = " << nFreqInCatTotal << '/' << nWordsInCat << '\n';
        return doubleToInt(p * 1000););
}

void CBoW::addImage(CBoWWordBag * pWB) {
    if(IS_DEBUG) CHECK(!pWB || pWB->id() == CBoW::DONT_ADD, "CBoW::Compare: Bad parameter");

    cout << "About to add " << pWB->id() << "...\n";

    vImages.Push(pWB);

    if (pDictionary)
        pWB->CountWordOccurances < false > ();

    //Todo: when do we recreate dictionary? Every frame at first (doesn't matter how slow at first anyway)
    if (//vImages.Count() >=2 && //otherwise dont get good corresp. between 1 and 2 --doesn't seem to matter
            vImages.totalDescriptorCount() > nNextClusterCount && !bClustering
            && PARAMS.BOWClustering.RECLUSTER_FREQUENCY > 0 //dReclusterFrequency<0 == no auto. clustering
            ) {
        bClustering = true; //This is safe as we currently have a lock, and it will be unset when the clustring code has a lock.
        nNextClusterCount = doubleToInt(vImages.totalDescriptorCount() * PARAMS.BOWClustering.RECLUSTER_FREQUENCY);
        recreateDictionary();
    }
    cout << "Added " << pWB->id() << "\n";
}

//Always takes memory

void CBoW::addImage(CDescriptorSet ** ppDescriptors, int nNewImageId) {
    WRITE_LOCK;

    CHECK(!ppDescriptors, "CBoW::addImage: Null parameter");
    CDescriptorSet * pDescriptors = *ppDescriptors;
    CHECK(!pDescriptors, "CBoW::addImage: Null parameter");

    CHECK(!pDescriptors->Count(), "CBoW::addImage: No descriptors");
    CHECK(nNewImageId == DONT_ADD, "CBoW::addImage: Not adding image");
    CHECK(vImages.exists(nNewImageId), "CBoW::addImage: Id already in use");

    /*	if(!pAllDescriptors)
            {
                    pAllDescriptors = pDescriptors->makeNewDS();
            }*/

    addImage(new CBoWWordBag(pDescriptors, this, nNewImageId));
    *ppDescriptors = 0; //I now own the memory;
}

TBoWMatchVector * CBoW::getMatches(CDescriptorSet * pDescriptors, int nReturnMax) //takes pDescriptors memory from caller. Is it freed?
{
    //UPGRADEABLE_LOCK; //Might wait here for clustering to finish
    WRITE_LOCK; //Use a write lock cos we need one for recreatewordbag.

    if(IS_DEBUG) CHECK(!pDescriptors || pDescriptors->Count() == 0, "CBoW::getMatches: Bad/empty descriptor set");
    //if(IS_DEBUG) CHECK(!pAllDescriptors, "CBoW::getMatches: Error initialising all descriptors' set");
    if(IS_DEBUG) CHECK(!pDescriptors || !pDescriptors->Count(), "CBoW::getMatches: Bad parameter");

    CBoWWordBag wb(pDescriptors, this, DONT_ADD);

    return getMatches_int(&wb, nReturnMax);
}

TBoWMatchVector * CBoW::getMatches_int(CBoWWordBag * pwb, int nReturnMax) {
    if (!nReturnMax) return 0;

    if (nReturnMax == RETURN_ALL) nReturnMax = (int) vImages.Count();
    return (pwb->*pfn_getBoWMatches)(nReturnMax); //will be empty for first image(s)
};

TBoWMatchVector * CBoW::getMatches(int nId_in, int nReturnMax) {
    WRITE_LOCK; //Might wait here for clustering to finish
    cout << "Getting matches with id " << nId_in << endl;
    TBoWMatchVector * pRet = 0;
    try {
        CBoWWordBag * pwb = vImages.Find(nId_in);
        pRet = getMatches_int(pwb, nReturnMax);
    }    catch (...) {
        RELEASE_WRITE_LOCK;
        throw;
    }
    RELEASE_WRITE_LOCK;

    return pRet;
};

/*extern double statA,statB,statT1,statT2,statScale;
extern int statN, statMinCorr, bScale;

double CBoW::geomAdd(double stat)
{
    return statScale * geomScale(stat);
};

double CBoW::geomScale(double stat)
{
    if(stat < statT1) return statA;
    if(stat > statT2) return statB;
    return statA + (stat-statT1)*(statB-statA)/(statT2-statT1);
};

//Apply a geometric compat. update to top nProcess matches
TBoWStatMatchVector * CBoW::getGeomMatches(TBoWMatchVector * pMatches, int id, int nProcess)
{
    //READ_LOCK; //Might wait here for clustering to finish

    TBoWStatMatchVector * pNewMatches = new TBoWStatMatchVector();
    TBoWMatchVector::iterator pMatchEnd = pMatches->end();

    for(TBoWMatchVector::iterator pMatch = pMatches->begin(); pMatch != pMatchEnd && nProcess>0; pMatch++, nProcess--)
    {
        CBoWCorrespondences * pCorr = getCorrespondences(pMatch->id(), id, statN, statN); //this will lock

        double dAdd = 0, dScale = 1;

        if((int)pCorr->size() >= statMinCorr / *CBoWCorrespondences::MIN_CORRESPONDENCES* /)
        {
            CBoWCorrespondences::TNumberType stat=pCorr->CompStat();
            if(!bScale)
                dAdd = geomAdd(stat);
            else
                dScale = geomScale(stat);
        }

        delete pCorr;

        double dMatchStrength = dAdd + dScale * pMatch->MatchStrength();
        CBoWMatch<double> newMatch(pMatch->id(), dMatchStrength);
        pNewMatches->insert(newMatch);
    }
    return pNewMatches;
};*/

int CBoW::CBoWWordBag::BruteForceMatch(CBoW::CBoWWordBag * pWB) const {
    const CDescriptorSet * pDS = pWB->DescriptorSet();
    double dDist = 0;
    for (int i = 0; i < pDS->Count(); i++)
        dDist += pDescriptors->Closest(pDS->get_const(i));

    dDist /= pDS->Count();

    return -doubleToInt(1000 * dDist);
}

CBoW::CBoWWordBag::CBoWWordBag(CDescriptorSet * pDescriptors_in, const CBoW * pParent_in, int nId_in) : pParent(pParent_in), pDescriptors(pDescriptors_in), nId(nId_in) {
    if(IS_DEBUG) CHECK(!pDescriptors_in || !pDescriptors_in->Count() || !pParent_in, "CBoW::CBoWWordBag::CBoWWordBag: Bad parameter");

    if (pParent->PARAMS.BOWClustering.BRUTEFORCE_MATCHING_LEVEL() > 0)
        aDescriptorWordsBF.reserve(pDescriptors_in->Count());

    if(IS_DEBUG) CHECK(pParent->vImages.exists(nId), "Id already in use");

    const int LEVELS = pParent->PARAMS.BOWClustering.LEVELS;
    if(IS_DEBUG) CHECK(LEVELS <= 0 || LEVELS > MAX_BOW_LEVELS, "CBoW::CBoWWordBag::CBoWWordBag: Uninitialised parent");

    aWords = new TWordBag[LEVELS];

    if (pParent->pDictionary) {
        RecreateWordBag();
        (this->*(pParent->pfn_WeightWordBag))();
    }
}

CBoW::CBoWWordBag::~CBoWWordBag() {
    if (pParent->pDictionary) {
        const int LEVELS = pParent->PARAMS.BOWClustering.LEVELS;

        if(IS_DEBUG) CHECK(LEVELS <= 0 || LEVELS >= MAX_BOW_LEVELS, "Access OOB/overwrite");

        for (int nLevel = 0; nLevel < LEVELS; nLevel++) {
            //Delete word locations SHOULD HAPPEN IN D'TOR
            TWordBag::iterator pEnd = aWords[nLevel].end();
            for (TWordBag::iterator pImWord = aWords[nLevel].begin(); pImWord < pEnd; pImWord++)
                pImWord->DestroyLoc();
        }

    }
    delete [] aWords;

    if (id() != DONT_ADD) {
        CDescriptorSet::deleteDS(&pDescriptors);
    }
}

void CBoW::IDToWBIterator(int nId1, int nId2, constImIt &ppWb1, constImIt &ppWb2) const {
    constImIt end = vImages.end();
    ppWb1 = end;
    ppWb2 = end;
    for (constImIt ppIm = vImages.begin(); ppIm != end; ppIm++) {
        int nId = (*ppIm)->id();
        if (nId == nId1)
            ppWb1 = ppIm;
        /*else*/if (nId == nId2)
            ppWb2 = ppIm;
    }
    if (ppWb1 == end || ppWb2 == end) THROW("CBoW::IDToWBIterator: Id not found");
}

void CBoW::IDToWBIterator(int nId1, int nId2, constImIt &ppWb1, constImIt &ppWb2) {
    constImIt end = vImages.end();
    ppWb1 = end;
    ppWb2 = end;
    //typedef pair<constImIt, constImIt> TWBPair;
    //constImIt ppWb = find<>(begin, end, nId1);
    for (constImIt ppIm = vImages.begin(); ppIm != end; ppIm++) {
        int nId = (*ppIm)->id();
        if (nId == nId1)
            ppWb1 = ppIm;
        else if (nId == nId2)
            ppWb2 = ppIm;
    }
    if (ppWb1 == end || ppWb2 == end) THROW("CBoW::IDToWBIterator: Id not found");
}

//Merge 2 word bags into one category:
/*void CBoW::merge(int nId1, int nId2)
{
        if(IS_DEBUG) CHECK(nId1 == nId2 || nId1 == DONT_ADD || nId2 == DONT_ADD, "CBoW::merge: Bad parameters");

        WRITE_LOCK;
        try {

                constImIt ppWb1 = vImages.Find(nId1);
                constImIt ppWb2 = vImages.Find(nId2);
                //IDToWBIterator(nId1, nId2, ppWb1, ppWb2);
                CBoWWordBag * pWb1 = *ppWb1;
                CBoWWordBag * pWb2 = *ppWb2;

                pWb1->merge(pWb2);

                //Erase Wb2 from vImages
                vImages.Erase(nId2);
        } catch(...)
        {
                RELEASE_WRITE_LOCK;
                throw;
        }
        RELEASE_WRITE_LOCK;
}*/

//Merge 2 word bags into one category:

void CBoW::remove(int nId) {
    if(IS_DEBUG) CHECK(nId == DONT_ADD, "CBoW::merge: Bad parameters");
    if(IS_DEBUG) CHECK(!vImages.exists(nId), "Erasing image that doesn't exist")

    WRITE_LOCK;
    const bool bLocked_notClustering = mxClusteringCantDelete.try_lock();

    try {
        vImages.Erase(nId, bLocked_notClustering);
    } catch (...) {
        if (bLocked_notClustering)
            mxClusteringCantDelete.unlock();
        throw;
    }

    if (bLocked_notClustering)
        mxClusteringCantDelete.unlock();

    RELEASE_WRITE_LOCK;

}

CBoW::CBoWWordBag::CBoWImageWord::CBoWImageWord(CBoWWord * pWord_in, const CLocation & Loc) : pWord(pWord_in), nFrequency(1), nWeight(0) {
    if(IS_DEBUG) CHECK(!pWord_in, "CBoW::CBoWWordBag::CBoWImageWord: No word in c'tor");
    aLocations = Loc;
};

void CBoW::CBoWWordBag::CBoWImageWord::IncrementFrequency(const CLocation & Loc) {
    nFrequency++;

    if (!aLocations.loc().zero()) {
        if(IS_DEBUG) CHECK(Loc.xAsIntFast() <= 0 || Loc.x() > 4000 || Loc.yAsIntFast() <= 0 || Loc.y() > 3000, "CBoW::CBoWWordBag::CBoWImageWord::IncrementFrequency: trying to increment location with 0 or multiple locations");
        if (nFrequency == 2) {
            //New array
            CLocation * aNewLocationArray = new CLocation[2];
            aNewLocationArray[0] = aLocations.loc();
            aNewLocationArray[1] = Loc;
            aLocations.pLoc = aNewLocationArray; /*.CLocationPtr;*/
        }
#if LOCATION_STORE_LIM>2
        else if (nFrequency == 3) {
            //New array (resize)
            CLocation * aNewLocationArray = new CLocation[LOCATION_STORE_LIM];
            aNewLocationArray[0] = aLocations.pLoc[0];
            aNewLocationArray[1] = aLocations.pLoc[1];
            aNewLocationArray[2] = Loc;
            delete [] aLocations.pLoc;
            aLocations.pLoc = aNewLocationArray;
        }
#if LOCATION_STORE_LIM>3
        else if (nFrequency <= LOCATION_STORE_LIM) {
            //Add to array
            aLocations.pLoc[nFrequency - 1] = Loc;
        }
#endif
#endif
        else //(nFrequency > LOCATION_STORE_LIM)
        {
            //Don't bother storing anything
            delete [] aLocations.pLoc;
            aLocations.pLoc = 0;
        }
    }
};

/*inline const CLocation & CBoW::CBoWWordBag::CBoWImageWord::DupLocation() const
{
    return aLocations;
    //if(aLocations.CLocationPtr == 0)
    //    return aLocations; //equiv return 0

    //if(nFrequency == 1)
    //    return aLocations;

    ////else
    //if(IS_DEBUG) CHECK(nFrequency>LOCATION_STORE_LIM, "CBoW::CBoWWordBag::CBoWLocation: We have more locations that we should do");

    //CLocation * aDupLocations = new CLocation[nFrequency == 2 ? 2 : 4];
    //
    //for(int i=0; i<nFrequency; i++)
    //    aDupLocations[i] = Location(i);

    ////CLocation dupedLoc(aDupLocations);

    //return CLocation(aDupLocations);
};*/

//inline const CLocation & CBoW::CBoWWordBag::CBoWImageWord::Location(int n) const

CBoW::CBoWWordBag::CBoWImageWord::~CBoWImageWord() {
    //DestroyLoc(); TB 14-5-2011 Remove this line to stop crash on uninitialised variable. Todo: Debug properly, is prob leaking now.
    //Location should already be destroyed or duplicated
};

void CBoW::CBoWWordBag::CBoWImageWord::DestroyLoc() {
    if (nFrequency > 1 && aLocations.pLoc) {
        delete [] aLocations.pLoc;
        aLocations.pLoc = 0;
    }
};

//Assimilate descriptors and words from pOtherWB
/*void CBoW::CBoWWordBag::merge(CBoWWordBag * pOtherWB)
{
    if(pParent->pDictionary)
    {
        pOtherWB->CountWordOccurances<true>();
        CountWordOccurances<true>();
    }

    pDescriptors->Push(pOtherWB->DescriptorSet());

    if(pParent->pDictionary)
    {
        RecreateWordBag();
        CountWordOccurances<false>(); //if needed??
        (this->*(pParent->pfn_WeightWordBag))();
    }
    delete pOtherWB;
}*/

//CBoWCorrespondences * CBoW::getCorrespondences(int nId1, int nId2, int nN, int nM)
//{
//	if(!pDictionary)
//		return 0; //Need to wait for dictionary to be created...
//
//	WRITE_LOCK;
//
//    CBoWCorrespondences * pCorr = 0;
//
//    try
//    {
//		if(IS_DEBUG) CHECK(nN > nM || nN<=0 || nM <= 0 || nM>LOCATION_STORE_LIM || nId1 == DONT_ADD || nId2 == DONT_ADD, "CBoW::getCorrespondences: Bad params");
//
//		constImIt ppWb1 = vImages.Find(nId1);
//		constImIt ppWb2 = vImages.Find(nId2);
//		//IDToWBIterator(nId1, nId2, ppWb1, ppWb2);
//		CBoWWordBag * pWb1 = *ppWb1;
//		CBoWWordBag * pWb2 = *ppWb2;
//
//		pCorr = pWb1->getCorrespondences(pWb2, nN, nM);
//	} catch(...)
//	{
//		RELEASE_WRITE_LOCK;
//		throw;
//	} 
//	RELEASE_WRITE_LOCK;
//
//	return pCorr;
//}

void CBoW::ensureClusteredOnce() //should NOT have a write lock here
{
    if (pDictionary == 0) {
        MX_NO_DELETE_WHILE_CLUSTER_IN_SEPERATE_THREAD;
    }

    if (pDictionary == 0) {
        cout << PARAMS.BOWClustering.CLUSTER_IN_SEPERATE_THREAD;
        cout << PARAMS.BOWClustering.RECLUSTER_FREQUENCY;
        THROW("Before correspondences or matches are available you must create dictionary by calling recreateDictionary(), or by setting BOW.BOWClustering.CLUSTER_IN_SEPERATE_THREAD=true and BOW.BOWClustering.RECLUSTER_FREQUENCY > 1")
    }
}

const CBoWCorrespondences * CBoW::getCorrespondences(CDescriptorSet * pDS1, int nId2, const CBOWMatchingParams & MATCHING_PARAMS) {
    ensureClusteredOnce();

    SCOPED_WRITE_LOCK; //Use STATICS while matching

    if(IS_DEBUG) CHECK(nId2 == DONT_ADD || !pDS1, "CBoW::getCorrespondences: Bad params");

    CBoWWordBag * pWb2 = vImages.Find(nId2);
    CBoWWordBag Wb1(pDS1, this, DONT_ADD);

    return Wb1.getCorrespondences(pWb2, MATCHING_PARAMS);
}

const CBoWCorrespondences * CBoW::getCorrespondences(int nId1, CDescriptorSet * pDS2, const CBOWMatchingParams & MATCHING_PARAMS) {
    ensureClusteredOnce();

    SCOPED_WRITE_LOCK; //Use STATICS while matching

    if(IS_DEBUG) CHECK(nId1 == DONT_ADD || !pDS2, "CBoW::getCorrespondences: Bad params");

    CBoWWordBag * pWb1 = vImages.Find(nId1);
    CBoWWordBag Wb2(pDS2, this, DONT_ADD);

    return pWb1->getCorrespondences(&Wb2, MATCHING_PARAMS);
}

const CBoWCorrespondences * CBoW::getCorrespondences(CDescriptorSet * pDS1, CDescriptorSet * pDS2, const CBOWMatchingParams & MATCHING_PARAMS) {
    ensureClusteredOnce();

    SCOPED_WRITE_LOCK; //Use STATICS while matching
    if(IS_DEBUG) CHECK(!pDS1 || !pDS2, "CBoW::getCorrespondences: Bad params");

    CBoWWordBag Wb1(pDS1, this, DONT_ADD);
    CBoWWordBag Wb2(pDS2, this, DONT_ADD);

    return Wb1.getCorrespondences(&Wb2, MATCHING_PARAMS);
}

const CBoWCorrespondences * CBoW::getCorrespondences(int nId1, int nId2, const CBOWMatchingParams & MATCHING_PARAMS) {
    ensureClusteredOnce();

    SCOPED_WRITE_LOCK; //Use STATICS while matching
    if(IS_DEBUG) CHECK(nId1 == DONT_ADD || nId2 == DONT_ADD, "CBoW::getCorrespondences: Bad params");
    CBoWWordBag * pWb1 = vImages.Find(nId1);
    CBoWWordBag * pWb2 = vImages.Find(nId2);

    return pWb1->getCorrespondences(pWb2, MATCHING_PARAMS);
}

#define QUIET(x)

void CBoW::CBoWWordBag::BFCorrFromRange(const CDescriptorSet * pDS1, const CDescriptorSet * pDS2,
        const CDynArray<CBoWWordDescMatch>::const_iterator WB1WordList,
        const CDynArray<CBoWWordDescMatch>::const_iterator WB1WordListRangeEnd,
        const CDynArray<CBoWWordDescMatch>::const_iterator WB2WordList,
        const CDynArray<CBoWWordDescMatch>::const_iterator WB2WordListRangeEnd,
        double dCondition, CBoWCorrespondences * pCorrespondences) {
    typedef std::map<CLocation, pair<CLocation, CDescriptor::TDist> > TMatchMap;
    TMatchMap aStrongestMatch;

    for (CDynArray<CBoW::CBoWWordBag::CBoWWordDescMatch>::const_iterator pW1 = WB1WordList; pW1 != WB1WordListRangeEnd; pW1++) {
        CDescriptor::TDist nClosest = MAX_INT, n2ndClosest = MAX_INT;
        const CDescriptor * pClosestDesc = 0;

        int nDesc1 = pW1->DescriptorSetIdx();
        const CDescriptor * pDesc1 = pDS1->get_const(nDesc1);

        for (CDynArray<CBoW::CBoWWordBag::CBoWWordDescMatch>::const_iterator pW2 = WB2WordList; pW2 != WB2WordListRangeEnd; pW2++) {
            int nDesc2 = pW2->DescriptorSetIdx();
            const CDescriptor * pDesc2 = pDS2->get_const(nDesc2);
            CDescriptor::TDist dist = pDesc1->distance(pDesc2);
            if (dist < nClosest) {
                pClosestDesc = pDesc2;
                n2ndClosest = nClosest;
                nClosest = dist;
            } else if (dist < n2ndClosest) {
                n2ndClosest = dist;
            }
        }

        if (dCondition * n2ndClosest > (double) nClosest) {
            TMatchMap::iterator pStrongestMatch = aStrongestMatch.find(pClosestDesc->location());
            if (pStrongestMatch != aStrongestMatch.end()) {
                //We've already matched a different one to this descriptor
                if (pStrongestMatch->second.second > nClosest) //but this one is closer
                    pStrongestMatch->second = make_pair<CLocation, CDescriptor::TDist > (pDesc1->location(), nClosest);
            } else
                aStrongestMatch[pClosestDesc->location()] = make_pair<CLocation, CDescriptor::TDist > (pDesc1->location(), nClosest);
        }
    }
    //Now turn map into correspondences:
    const double dPriorProb = 0.7;
    for (TMatchMap::const_iterator pBestMatch = aStrongestMatch.begin(); pBestMatch != aStrongestMatch.end(); pBestMatch++)
        pCorrespondences->push_back(CCorrespondence(pBestMatch->first, pBestMatch->second.first, dPriorProb));
}

void CBoW::CBoWWordBag::ensureBiDiMatches_remove(TSimpleMatchMap & aStrongestMatchesA, TSimpleMatchMap & aStrongestMatchesB) {
    for (TSimpleMatchMap::iterator pMatchesWithA = aStrongestMatchesA.begin(); pMatchesWithA != aStrongestMatchesA.end(); pMatchesWithA++) {
        TSimpleLocSet & matchesInB = pMatchesWithA->second;
        TSimpleLocSet matchesWithAToRemove;
        CLocation locInA = pMatchesWithA->first;
        QUIET(cout << "Size before = " << matchesInB.size() << endl;)
        for (TSimpleLocSet::const_iterator pMatchInB = matchesInB.begin(); pMatchInB != matchesInB.end(); ++pMatchInB) {
            if (!aStrongestMatchesB[*pMatchInB].exists(locInA))
                //Remove match from A
                matchesWithAToRemove.insert(*pMatchInB);
        }
        for (TSimpleLocSet::const_iterator pMatchInBToRemove = matchesWithAToRemove.begin(); pMatchInBToRemove != matchesWithAToRemove.end(); ++pMatchInBToRemove)
            matchesInB.erase(*pMatchInBToRemove);

        QUIET(cout << "Size after = " << matchesInB.size() << endl;)
    }
}

//void CBoW::CBoWWordBag::ensureBiDiMatches(TSimpleMatchMap & aStrongestMatchesA, TSimpleMatchMap & aStrongestMatchesB)
//{
//	for(TSimpleMatchMap::const_iterator pMatchesWithA = aStrongestMatchesA.begin(); pMatchesWithA != aStrongestMatchesA.end(); pMatchesWithA++)
//	{
//		const TSimpleLocSet & matchesInB = pMatchesWithA->second;
//		CLocation locInA = pMatchesWithA->first;
//		QUIET(cout << "Size before = " << matchesInB.size() << endl;)
//		for(TSimpleLocSet::const_iterator pMatchInB = matchesInB.begin(); pMatchInB != matchesInB.end(); pMatchInB++)
//		{
//			aStrongestMatchesB[*pMatchInB].insert(locInA);
//		}
//		QUIET( cout << "Size after = " << matchesInB.size() << endl; )
//	}
//	//B should contain a forall ab.
//	/*for(TMatchMap::iterator pMatchesWithB = aStrongestMatchesB.begin(); pMatchesWithB != aStrongestMatchesB.end(); pMatchesWithB++)
//	{
//		cout << "Match with " << pMatchesWithB->first.x() << ',' << pMatchesWithB->first.y() << ":";
//		for(TMatchSet::iterator pMatchWithB = pMatchesWithB->second.begin(); pMatchWithB != pMatchesWithB->second.end(); pMatchWithB++)
//		{
//			cout << pMatchWithB->first.x() << ',' << pMatchWithB->first.y() << "-" << pMatchWithB->second << ',';
//		}
//		cout << endl;
//	}*/
//}

void CBoW::CBoWWordBag::cullMatches(const TMatchMap & aStrongestMatchesLeft, TSimpleMatchMap & aStrongestMatchesLeftGood, double dCondition) {
    for (TMatchMap::const_iterator pMatches = aStrongestMatchesLeft.begin(); pMatches != aStrongestMatchesLeft.end(); pMatches++) {
        const TMatchSet & aMatchesRight = pMatches->second;

        QUIET(for (TMatchSet::iterator pMatch = aMatchesRight.begin(); pMatch != aMatchesRight.end(); pMatch++) {
            cout << pMatch->second << ',';
        }
        cout << endl;)

        //bool bFirst = true;

        CDescriptor::TDist nLastScore = MAX_INT;
        for (TMatchSet::const_iterator pMatch = aMatchesRight.begin(); pMatch != aMatchesRight.end(); pMatch++) {
            CDescriptor::TDist nThisScore = pMatch->second;

            //if(bFirst && nThisScore > 2000) break;

            if (nLastScore < nThisScore * dCondition) {
                //aMatchesRight.erase(pMatch, aMatchesRight.end());
                break;
            }
            aStrongestMatchesLeftGood[pMatches->first].insert(pMatch->first);
            nLastScore = nThisScore;
            //bFirst = false;
        }
    }
}

void CBoW::CBoWWordBag::addMatchesToLeftPoint(const TSimpleMatchMap & aStrongestMatchesLeft, const TSimpleMatchMap & aStrongestMatchesRight, const CLocation locLeft, TSimpleLocSet & matchesLeft, TSimpleLocSet & matchesRight) {
    const TSimpleLocSet & aMatchesRight = aStrongestMatchesLeft.find(locLeft)->second;
    for (TSimpleLocSet::const_iterator pMatch = aMatchesRight.begin(); pMatch != aMatchesRight.end(); ++pMatch) {
        CLocation locRight = *pMatch;
        if (!matchesRight.exists(locRight)) {
            //if((int)matchesRight.size() >= BF_NN) return false;
            matchesRight.insert(locRight);
            addMatchesToRightPoint(aStrongestMatchesRight, aStrongestMatchesLeft, locRight, matchesRight, matchesLeft);
            //if(!bSuccess) return false;
        }
    }
}

void CBoW::CBoWWordBag::addMatchesToRightPoint(const TSimpleMatchMap & aStrongestMatchesRight, const TSimpleMatchMap & aStrongestMatchesLeft, const CLocation locRight, TSimpleLocSet & matchesRight, TSimpleLocSet & matchesLeft) {
    const TSimpleLocSet & aMatchesLeft = aStrongestMatchesRight.find(locRight)->second;
    for (TSimpleLocSet::const_iterator pMatch = aMatchesLeft.begin(); pMatch != aMatchesLeft.end(); ++pMatch) {
        CLocation locLeft = *pMatch;
        if (!matchesLeft.exists(locLeft)) {
            //if((int)matchesLeft.size() >= BF_NN) return false;
            matchesLeft.insert(locLeft);
            addMatchesToLeftPoint(aStrongestMatchesLeft, aStrongestMatchesRight, locLeft, matchesLeft, matchesRight);
            //if(!bSuccess) return false;
        }
    }
}

/*void CBoW::CBoWWordBag::printMap(TSimpleMatchMap & aStrongestMatches)
{
        for(TSimpleMatchMap::iterator pMatchSet = aStrongestMatches.begin(); pMatchSet != aStrongestMatches.end(); pMatchSet++)
        {
                cout << pMatchSet->first.x() << ',' << pMatchSet->first.y() << " matches: ";
                TSimpleLocSet & aMatchSet = pMatchSet->second;
                for(TSimpleLocSet::const_iterator pMatch = aMatchSet.begin(); pMatch != aMatchSet.end(); pMatch++)
                        cout << pMatch->x() << ',' << pMatch->y() << ", ";
                cout << endl;
        }
}*/


void CBoW::CBoWWordBag::BFCorrFromRangeNN(const CDescriptorSet * pDS1, const CDescriptorSet * pDS2,
        const CDynArray<CBoWWordDescMatch>::const_iterator WB1WordList,
        const CDynArray<CBoWWordDescMatch>::const_iterator WB1WordListRangeEnd,
        const CDynArray<CBoWWordDescMatch>::const_iterator WB2WordList,
        const CDynArray<CBoWWordDescMatch>::const_iterator WB2WordListRangeEnd,
        double dCondition, const int BF_NN, CBoWCorrespondences * pCorrespondences, const CBOWParams::CBOWCorrespondenceProbParams & CORR_PROB_PARAMS) {
    /*
    TMatchMap aStrongestMatchesLeft, aStrongestMatchesRight;

    for(CDynArray<CBoW::CBoWWordBag::CBoWWordDescMatch>::const_iterator pW1 = WB1WordList; pW1 != WB1WordListRangeEnd; pW1++)
    {
            int nDesc1 = pW1->DescriptorSetIdx();
            const CDescriptor * pDesc1 = (*pDS1)[nDesc1];

            for(CDynArray<CBoW::CBoWWordBag::CBoWWordDescMatch>::const_iterator pW2 = WB2WordList; pW2 != WB2WordListRangeEnd; pW2++)
            {
                    int nDesc2 = pW2->DescriptorSetIdx();
                    const CDescriptor * pDesc2 = (*pDS2)[nDesc2];
                    CDescriptor::TDist dist = pDesc1->distance(pDesc2);

                    //if(dist<1500)
                    {
                            aStrongestMatchesLeft[pDesc1->location()].insert(make_pair<CLocation, CDescriptor::TDist>(pDesc2->location(), dist));
                            aStrongestMatchesRight[pDesc2->location()].insert(make_pair<CLocation, CDescriptor::TDist>(pDesc1->location(), dist));
                    }
            }
    }

    //Now cull weaker matches
    TSimpleMatchMap aStrongestMatchesLeftGood, aStrongestMatchesRightGood;
    cullMatches(aStrongestMatchesLeft, aStrongestMatchesLeftGood, dCondition);
    cullMatches(aStrongestMatchesRight, aStrongestMatchesRightGood, dCondition);
     */
    static CMatchMap aStrongestMatches(CORR_PROB_PARAMS); //Not TS, that's ok cos we're locked out

    aStrongestMatches.init(WB1WordListRangeEnd - WB1WordList, WB2WordListRangeEnd - WB2WordList);
    int nL = 0;
    for (CDynArray<CBoW::CBoWWordBag::CBoWWordDescMatch>::const_iterator pW1 = WB1WordList; pW1 != WB1WordListRangeEnd; pW1++) {
        int nDesc1 = pW1->DescriptorSetIdx();
        const CDescriptor * pDesc1 = pDS1->get_const(nDesc1);
        int nR = 0;

        for (CDynArray<CBoW::CBoWWordBag::CBoWWordDescMatch>::const_iterator pW2 = WB2WordList; pW2 != WB2WordListRangeEnd; pW2++) {
            int nDesc2 = pW2->DescriptorSetIdx();
            const CDescriptor * pDesc2 = pDS2->get_const(nDesc2);
            CDescriptor::TDist dist = pDesc1->distance(pDesc2);

            aStrongestMatches.add(nL, nR, pDesc1->location(), pDesc2->location(), dist);

            nR++;
        }
        nL++;
    }
    aStrongestMatches.sort();

    //Now cull weaker matches
    TSimpleMatchMap aStrongestMatchesLeftGood, aStrongestMatchesRightGood;
    //cullMatches(aStrongestMatchesLeft, aStrongestMatchesLeftGood, dCondition);
    //cullMatches(aStrongestMatchesRight, aStrongestMatchesRightGood, dCondition);
    nL = 0;
    for (CDynArray<CBoW::CBoWWordBag::CBoWWordDescMatch>::const_iterator pW1 = WB1WordList; pW1 != WB1WordListRangeEnd; pW1++) {
        int nDesc1 = pW1->DescriptorSetIdx();
        const CDescriptor * pDesc1 = pDS1->get_const(nDesc1);
        aStrongestMatches.getLRMatchesWC(aStrongestMatchesLeftGood[pDesc1->location()], nL, dCondition);
        nL++;
    }
    int nR = 0;
    for (CDynArray<CBoW::CBoWWordBag::CBoWWordDescMatch>::const_iterator pW2 = WB2WordList; pW2 != WB2WordListRangeEnd; pW2++) {
        int nDesc2 = pW2->DescriptorSetIdx();
        const CDescriptor * pDesc2 = pDS2->get_const(nDesc2);
        aStrongestMatches.getRLMatchesWC(aStrongestMatchesRightGood[pDesc2->location()], nR, dCondition);
        nR++;
    }

    QUIET(cout << "LEFT before: ";
            printMap(aStrongestMatchesLeftGood);
            cout << "RIGHT before: ";
            printMap(aStrongestMatchesRightGood);)

            ensureBiDiMatches_remove(aStrongestMatchesLeftGood, aStrongestMatchesRightGood);
    ensureBiDiMatches_remove(aStrongestMatchesRightGood, aStrongestMatchesLeftGood);

    QUIET(cout << "LEFT: ";
            printMap(aStrongestMatchesLeftGood);
            cout << "RIGHT: ";
            printMap(aStrongestMatchesRightGood);)

            //For each point recurse and find a L-R set. If N-N or less add corresps.
            nL = 0;
    TSimpleLocSet leftPointsAlreadyMatched;
    for (CDynArray<CBoW::CBoWWordBag::CBoWWordDescMatch>::const_iterator pW1 = WB1WordList; pW1 != WB1WordListRangeEnd; pW1++, nL++) {
        int nDesc1 = pW1->DescriptorSetIdx();
        const CDescriptor * pDesc1 = pDS1->get_const(nDesc1);

        if (leftPointsAlreadyMatched.exists(pDesc1->location())) {
            QUIET(cout << pDesc1->location().x() << ',' << pDesc1->location().y() << " already matched\n";)
            continue;
        }
        QUIET(cout << pDesc1->location().x() << ',' << pDesc1->location().y() << " not yet excluded\n";)

        TSimpleLocSet matchesLeft, matchesRight;
        matchesLeft.insert(pDesc1->location());
        addMatchesToLeftPoint(aStrongestMatchesLeftGood, aStrongestMatchesRightGood, pDesc1->location(), matchesLeft, matchesRight);

        //Why doesn't this work?! leftPointsAlreadyMatched.insert(matchesLeft.begin(), matchesLeft.end());
        for (TSimpleLocSet::const_iterator pLocLeft = matchesLeft.begin(); pLocLeft != matchesLeft.end(); ++pLocLeft) {
            QUIET(cout << pLocLeft->x() << ',' << pLocLeft->y() << " added to xclude set\n";)
            leftPointsAlreadyMatched.insert(*pLocLeft);
        }

        int nLikelihoodInv = max((int) matchesRight.size(), (int) matchesLeft.size());
        if (nLikelihoodInv <= BF_NN) {
            const double dPriorProb = aStrongestMatches.getPP(nL) / nLikelihoodInv;

            QUIET(cout << "New set\n";)
            for (TSimpleLocSet::const_iterator pLocLeft = matchesLeft.begin(); pLocLeft != matchesLeft.end(); ++pLocLeft) {
                for (TSimpleLocSet::const_iterator pLocRight = matchesRight.begin(); pLocRight != matchesRight.end(); ++pLocRight) {
                    pCorrespondences->push_back(CCorrespondence(*pLocRight, *pLocLeft, dPriorProb));
                    QUIET(cout << "Adding corr. from " << pLocLeft->x() << ',' << pLocLeft->y() << " to " << pLocRight->x() << ',' << pLocRight->y() << '-' << aStrongestMatchesLeft[*pLocLeft].size() << '-' << aStrongestMatchesRight[*pLocRight].size() << '-' << nLikelihoodInv << endl;)
                }
            }
        }
    }
}

const CBoWCorrespondences * CBoW::CBoWWordBag::getNewBruteForceCorrespondences(const CBoW::CBoWWordBag * pWB, const CBOWMatchingParams & MATCHING_PARAMS) const {
    CBoWCorrespondences * pCorrespondences = new CBoWCorrespondences();
    const CDescriptorSet * pDS1 = DescriptorSet();
    const CDescriptorSet * pDS2 = pWB->DescriptorSet();
    const double dPP = pParent->PARAMS.BOWCorrespondenceProb.MAX_PRIOR;

    int NEARBY = abs(pWB->id() - id()) > 20 ? -1 : 1;
    const CMatchableDescriptors::CMatchSettings MS(MATCHING_PARAMS.BF_CORNER_CONDITION, dPP, MATCHING_PARAMS.MATCH_NN, -1 /*TODO*/, NEARBY, MATCHING_PARAMS.OI_NEARBY_ANGLE);

    if (pParent->PARAMS.BOWClustering.BRUTEFORCE_MATCHING_LEVEL() == 0) {
        return pDS2->getBruteForceCorrespondenceSet(pDS1, MS, pCorrespondences);
    } else {
        //const int EST_SIZE = 32;
        CSimpleMatchableDescriptors DS1Small; // ( pDS1->makeNewDS(EST_SIZE) );
        CSimpleMatchableDescriptors DS2Small; // ( pDS1->makeNewDS(EST_SIZE) );

        CDynArray<CBoWWordDescMatch>::const_iterator WB1WordList = aDescriptorWordsBF.begin();
        CDynArray<CBoWWordDescMatch>::const_iterator WB1WordListEnd = aDescriptorWordsBF.end();

        CDynArray<CBoWWordDescMatch>::const_iterator WB2WordList = pWB->aDescriptorWordsBF.begin();
        CDynArray<CBoWWordDescMatch>::const_iterator WB2WordListEnd = pWB->aDescriptorWordsBF.end();

        if (WB1WordList != WB1WordListEnd && WB2WordList != WB2WordListEnd) {
            const CBoWWord * pWordInWB2 = WB2WordList->Word();
            const CBoWWord * pWordInWB1 = WB1WordList->Word();

            bool bWB2ListEnd = false;
            bool bWB1ListEnd = false;

            do {
                if (pWordInWB2 > pWordInWB1) {
                    //the word in WB1 but not WB2:
                    do {
                        WB1WordList++;
                        bWB1ListEnd = WB1WordList == WB1WordListEnd;

                    } while (!bWB1ListEnd && WB1WordList->Word() == pWordInWB1);

                    if (!bWB1ListEnd) {
                        pWordInWB1 = WB1WordList->Word();
                    }
                } else if (pWordInWB2 < pWordInWB1) {
                    //the word in WB2 but not the WB1:
                    do {
                        WB2WordList++;
                        bWB2ListEnd = WB2WordList == WB2WordListEnd;

                    } while (!bWB2ListEnd && WB2WordList->Word() == pWordInWB2);

                    if (!bWB2ListEnd) {
                        pWordInWB2 = WB2WordList->Word();
                    }
                } else {
                    //Possible correspondences. Prob a range of each.
                    DS1Small.Clear();
                    DS2Small.Clear();

                    do {
                        DS1Small.Push_const(pDS1->get_const(WB1WordList->DescriptorSetIdx()));

                        WB1WordList++;
                        bWB1ListEnd = WB1WordList == WB1WordListEnd;

                    } while (!bWB1ListEnd && WB1WordList->Word() == pWordInWB1);

                    do {
                        DS2Small.Push_const(pDS2->get_const(WB2WordList->DescriptorSetIdx()));

                        WB2WordList++;
                        bWB2ListEnd = WB2WordList == WB2WordListEnd;

                    } while (!bWB2ListEnd && WB2WordList->Word() == pWordInWB2);

                    //Now brute-force match these ranges' descriptors.
                    DS2Small.getBruteForceCorrespondenceSet(&DS1Small, MS, pCorrespondences);

                    if (!bWB1ListEnd)
                        pWordInWB1 = WB1WordList->Word();

                    if (!bWB2ListEnd)
                        pWordInWB2 = WB2WordList->Word();
                }

            } while (!(bWB1ListEnd || bWB2ListEnd));
        }

    }

    return pCorrespondences;
}

const CBoWCorrespondences * CBoW::CBoWWordBag::getBruteForceCorrespondences(const CBoW::CBoWWordBag * pWB, const double dCondition, const int BF_NN) const {
    CBoWCorrespondences * pCorrespondences = new CBoWCorrespondences();

    if (pParent->PARAMS.BOWClustering.BRUTEFORCE_MATCHING_LEVEL() == 0) {
        const CDescriptorSet * pDS = pWB->DescriptorSet();
        //int nID = 0;
        for (int i = 0; i < pDS->Count(); i++) {
            const CDescriptor * pThisDesc = pDS->get_const(i);
            const int nRadius = pParent->PARAMS.DescriptorBinning.RADIUS;
            const CDescriptor * pClosestDesc = pDescriptors->ClosestDescriptor(pThisDesc, dCondition, nRadius);
            if (pClosestDesc) {
                //DO NOT NEED TO Check symmetry cos used condition
                //const CDescriptor * pClosestDescSym = pWB->pDescriptors->ClosestDescriptor(pClosestDesc, dCondition, nRadius);
                //if(pClosestDescSym == pThisDesc)
                pCorrespondences->push_back(CCorrespondence(pThisDesc->location(), pClosestDesc->location(), 0.25)); //Keep PPs low...
            }
        }
    } else {
        const CDescriptorSet * pDS1 = DescriptorSet();
        const CDescriptorSet * pDS2 = pWB->DescriptorSet();

        CDynArray<CBoWWordDescMatch>::const_iterator WB1WordList = aDescriptorWordsBF.begin();
        CDynArray<CBoWWordDescMatch>::const_iterator WB1WordListEnd = aDescriptorWordsBF.end();

        CDynArray<CBoWWordDescMatch>::const_iterator WB2WordList = pWB->aDescriptorWordsBF.begin();
        CDynArray<CBoWWordDescMatch>::const_iterator WB2WordListEnd = pWB->aDescriptorWordsBF.end();

        if (WB1WordList != WB1WordListEnd && WB2WordList != WB2WordListEnd) {

            const CBoWWord * pWordInWB2 = WB2WordList->Word();
            const CBoWWord * pWordInWB1 = WB1WordList->Word();

            bool bWB2ListEnd = false;
            bool bWB1ListEnd = false;

            do {
                if (pWordInWB2 > pWordInWB1) {
                    //the word in WB1 but not WB2:
                    do {
                        WB1WordList++;
                        bWB1ListEnd = WB1WordList == WB1WordListEnd;

                    } while (!bWB1ListEnd && WB1WordList->Word() == pWordInWB1);

                    if (!bWB1ListEnd) {
                        pWordInWB1 = WB1WordList->Word();
                    }
                } else if (pWordInWB2 < pWordInWB1) {
                    //the word in WB2 but not the WB1:
                    do {
                        WB2WordList++;
                        bWB2ListEnd = WB2WordList == WB2WordListEnd;

                    } while (!bWB2ListEnd && WB2WordList->Word() == pWordInWB2);

                    if (!bWB2ListEnd) {
                        pWordInWB2 = WB2WordList->Word();
                    }
                } else {
                    //Possible correspondences. Prob a range of each.
                    CDynArray<CBoWWordDescMatch>::const_iterator WB1WordListRangeEnd = WB1WordList;
                    do {
                        WB1WordListRangeEnd++;

                    } while (WB1WordListRangeEnd != WB1WordListEnd && WB1WordListRangeEnd->Word() == pWordInWB1);

                    CDynArray<CBoWWordDescMatch>::const_iterator WB2WordListRangeEnd = WB2WordList;
                    do {
                        WB2WordListRangeEnd++;

                    } while (WB2WordListRangeEnd != WB2WordListEnd && WB2WordListRangeEnd->Word() == pWordInWB2);

                    //Now brute-force match these ranges' descriptors.
                    if (BF_NN == 1) {
                        BFCorrFromRange(pDS1, pDS2, WB1WordList,
                                WB1WordListRangeEnd,
                                WB2WordList,
                                WB2WordListRangeEnd,
                                dCondition, pCorrespondences);
                    } else {
                        BFCorrFromRangeNN(pDS1, pDS2, WB1WordList,
                                WB1WordListRangeEnd,
                                WB2WordList,
                                WB2WordListRangeEnd,
                                dCondition, BF_NN, pCorrespondences, pParent->PARAMS.BOWCorrespondenceProb);
                    }

                    WB1WordList = WB1WordListRangeEnd;
                    bWB1ListEnd = WB1WordList == WB1WordListEnd;
                    if (!bWB1ListEnd)
                        pWordInWB1 = WB1WordList->Word();

                    WB2WordList = WB2WordListRangeEnd;
                    bWB2ListEnd = WB2WordList == WB2WordListEnd;
                    if (!bWB2ListEnd)
                        pWordInWB2 = WB2WordList->Word();
                }

            } while (!(bWB1ListEnd || bWB2ListEnd));
        }
    }

    return pCorrespondences;
}

const CBoWCorrespondences * CBoW::CBoWWordBag::getCorrespondences(const CBoWWordBag * pWB2, const CBOWMatchingParams & MATCHING_PARAMS) const {
    if (pParent->pDictionary == 0) {
        cout << "Error: no dictionary\n";
        return 0;
    }

    if (MATCHING_PARAMS.BF_CORRESPONDENCES == CBOWMatchingParams::eBF_Correspondences) {
        int NEARBY = abs(pWB2->id() - id()) > 20 ? -1 : 1;
        const CMatchableDescriptors::CMatchSettings MS(MATCHING_PARAMS.BF_CORNER_CONDITION, pParent->PARAMS.BOWCorrespondenceProb.MAX_PRIOR, MATCHING_PARAMS.MATCH_NN, -1, NEARBY, MATCHING_PARAMS.OI_NEARBY_ANGLE);
        return pWB2->DescriptorSet()->getBruteForceCorrespondenceSet(DescriptorSet(), MS, 0);
    } else if (MATCHING_PARAMS.BF_CORRESPONDENCES == CBOWMatchingParams::eBoW_BF_Correspondences) {
        return getNewBruteForceCorrespondences(pWB2, MATCHING_PARAMS);
    } else if (MATCHING_PARAMS.BF_CORRESPONDENCES == CBOWMatchingParams::eOldBoW_BF_Correspondences) {
        return getBruteForceCorrespondences(pWB2, MATCHING_PARAMS.BF_CORNER_CONDITION, MATCHING_PARAMS.MATCH_NN);
    } else if (MATCHING_PARAMS.BF_CORRESPONDENCES == CBOWMatchingParams::eBoWCorrespondences) {
        const int nBottomLevelIdx = pParent->PARAMS.BOWClustering.LEVELS - 1;
        return getCorrespondences(aWords[nBottomLevelIdx], pWB2->aWords[nBottomLevelIdx], MATCHING_PARAMS.MATCH_NN, MATCHING_PARAMS.MATCH_NN);
    } else
        THROW("Unhandled correspondence-finding method")
    }

/* Find all correspondences where there are N possiblilties in one image and M in the other
 * N<=M. If we assume N are correct then the proportion that are correct is N/(N*M)
 * => prob of one being correct = 1/M
 */
const CBoWCorrespondences * CBoW::CBoWWordBag::getCorrespondences(const TWordBag & WB1, const TWordBag & WB2, int nN, int nM) {
    if(IS_DEBUG) CHECK(nN > nM || nN <= 0 || nM > LOCATION_STORE_LIM || nM <= 0, "CBoW::CBoWWordBag::getCorrespondences: Bad params");

    CBoWCorrespondences * pCorrespondences = new CBoWCorrespondences();

    TWordBag::const_iterator WB1WordList = WB1.begin();
    TWordBag::const_iterator WB1WordListEnd = WB1.end();

    TWordBag::const_iterator WB2WordList = WB2.begin();
    TWordBag::const_iterator WB2WordListEnd = WB2.end();

    //const CBoWImageWord * WB2WordList = (*WB2WordList);
    const CBoWWord * pWordInWB2 = WB2WordList->Word();
    //const CBoWImageWord * WB1WordList = (*WB1WordList);
    const CBoWWord * pWordInWB1 = WB1WordList->Word();

    bool bWB2ListEnd = false;
    bool bWB1ListEnd = false;
    //This is the inner critical section
    do {
        if ((pWordInWB2 > pWordInWB1 || bWB2ListEnd) && !bWB1ListEnd) {
            //the word in WB1 but not WB2:
            WB1WordList++;
            if (WB1WordList != WB1WordListEnd) {
                pWordInWB1 = WB1WordList->Word();
            } else
                bWB1ListEnd = true;
        } else if ((pWordInWB2 < pWordInWB1 || bWB1ListEnd) && !bWB2ListEnd) {
            //the word in WB2 but not the WB1:
            WB2WordList++;
            if (WB2WordList != WB2WordListEnd) {
                pWordInWB2 = WB2WordList->Word();
            } else
                bWB2ListEnd = true;
        } else {
            //Possible correspondences
            int nF1 = WB1WordList->FrequencyInImage();
            int nF2 = WB2WordList->FrequencyInImage();
            int nMinF, nMaxF;
            if (nF1 < nF2) {
                nMinF = nF1;
                nMaxF = nF2;
            } else {
                nMinF = nF2;
                nMaxF = nF1;
            }

            if (nMaxF <= nM && nMinF <= nN) {
                //it's an N:M or less correspondence
                //iterate thru possibilities
                //if(IS_DEBUG) CHECK(!aLoc1 || !aLoc2, "CBoW::CBoWWordBag::getCorrespondences: No point locations");
                //LikelihoodInv()static const int MAX_STRENGTH = 3*4;//*5*7;
                //int nStrength= MAX_STRENGTH/nMaxF;
                for (int i1 = 0; i1 < nF1; i1++)
                    for (int i2 = 0; i2 < nF2; i2++) {
                        //Todo--use: int nWordWeight = nLogCount - intintLog(bTF_IDF ? WB1WordList->Word()->TotalOccurances() : WB1WordList->Word()->TotalFrequency());
                        CLocation Loc1 = WB1WordList->Location(i1);
                        CLocation Loc2 = WB2WordList->Location(i2);
                        int nLikelihoodInv = nMaxF;
                        const double dPriorProb = 0.6;
                        pCorrespondences->push_back(CCorrespondence(Loc1, Loc2, dPriorProb / nLikelihoodInv));
                        /*if(i1*i2)
                            if(Loc1.x()>Loc2.y())
                                nMinF=Loc2.y(); //for VG to find uninitialised locations
                        else
                                                        if(Loc1.x()>Loc2.y())
                                                                nMinF=Loc1.y(); //for VG to find uninitialised locations*/
                    }
            }

            WB2WordList++;
            if (WB2WordList != WB2WordListEnd) {
                pWordInWB2 = WB2WordList->Word();
            } else
                bWB2ListEnd = true;

            WB1WordList++;
            if (WB1WordList != WB1WordListEnd) {
                pWordInWB1 = WB1WordList->Word();
            } else
                bWB1ListEnd = true;
        }
    } while (!(bWB1ListEnd && bWB2ListEnd));

    return pCorrespondences;
}

/*CBoWCorrGroup * CBoW::CBoWWordBag::new_cBoWCorrGroup(TWordBag::iterator word1, TWordBag::iterator word2)
{
        CBoWCorrGroup *  pCorrGrp = new CBoWCorrGroup();
        pCorrGrp->nNumCorr1 = word1->FrequencyInImage();
        pCorrGrp->nNumCorr2 = word2->FrequencyInImage();
        for (int i=0; i<pCorrGrp->nNumCorr1; i++)
        {
                pCorrGrp->Loc1[i] = word1->Location(i);
        }
        for (int i=0; i<pCorrGrp->nNumCorr2; i++)
        {
                pCorrGrp->Loc2[i] = word2->Location(i);
        }
        return pCorrGrp;
};*/

inline void CBoW::CBoWImageDB::Push(CBoWWordBag * pWB) {
    if(IS_DEBUG) CHECK(exists(pWB->id()), "Push: Word bag id already in use");
    if(IS_DEBUG) CHECK(!pWB, "Push: Null Word bag");
    push_back(pWB);
    //cout << "Adding id " << pWB->id() << "," << size()-1 << " to lookup map\n";
    imIdLookupMap[pWB->id()] = size() - 1; //size-1 is the index into the safeVector of what we've just inserted

    const CDescriptorSet * pDS = pWB->DescriptorSet();
    nAllDescriptorsCount += pDS->Count();

    if(IS_DEBUG) CHECK(!exists(pWB->id()), "Push: Word bag push failed");
};

void CBoW::CBoWImageDB::cleanUpErasedImages() {
    //Don't bother until the end--tends to delete cluster centres. Shouldn't waste too much mem
    for (iterator ppDeleteMe = vImagesToErase.begin(); ppDeleteMe != vImagesToErase.end(); ppDeleteMe++) {
        cout << "Erasing " << (*ppDeleteMe)->id() << endl;
        /*if(bUpdateWordCounts)
        {
                const CBoWWordBag::TWordBag & words = (*ppDeleteMe)->bottomLevelWords();
                for(int i=0; i < words.size(); i++) //May be concurrent with clustering but shouldn't be a probelm
                {
                        words[i].Word()->DecrementOccurances();
                }
        }*/

        delete *ppDeleteMe;
    }
    vImagesToErase.clear();
}

inline void CBoW::CBoWImageDB::Erase(imageNum nID, const bool bLocked_notClustering) {
    //Otherwise the clustering thread will also have a ptr to this image, so will remap words

    iterator ppWBPos = Find_int(nID);

    CBoWWordBag * pWB = *ppWBPos;

    const CDescriptorSet * pDS = pWB->DescriptorSet();
    nAllDescriptorsCount -= pDS->Count(); //this is not used while actually clustering

    imIdLookupMap.erase(nID);
    vImagesToErase.push_back(pWB);
    nErased++;
    *ppWBPos = 0;

    if (bLocked_notClustering) {
        cout << "Erasing images (" << nID << ") and updateing word freq. counts when not clustering\n";

        const CBoWWordBag::TWordBag & words = pWB->bottomLevelWords();
        for (int i = 0; i < words.size(); i++) //May be concurrent with clustering but shouldn't be a probelm
        {
            words[i].Word()->DecrementOccurances(/*words[i].FrequencyInImage()*/);
        }

        //BUT DO NOT DELETE YET OR MAY DELETE CLUSTER CENTRES
        //cleanUpErasedImages(true);
    } else
        cout << "Will erase image " << nID << " after clustering, don't need to update occurance counts\n";
}

/*inline CBoW::constImIt CBoW::CBoWImageDB::Find(imageNum nID)
{
    return Find_int(nID);
}

inline CBoW::CBoWImageDB::iterator CBoW::CBoWImageDB::Find_int(imageNum nID)
{
        return TImageDB::begin() + imIdLookupMap[nID];
    / *std::map<imageNum, int>::iterator pFound = imIdLookupMap.find(nID);
    if(IS_DEBUG) CHECK(pFound == imIdLookupMap.end(), "ID not found");
    return TImageDB::begin() + pFound->second;* /
}*/

CBoW::CBoWImageDB::~CBoWImageDB() {
    constImIt ppEnd = end();
    for (constImIt ppIm = begin(); ppIm != ppEnd; ppIm++) {
        CBoWWordBag * pWB = *ppIm;
        if (pWB)
            delete pWB;
    }

    cleanUpErasedImages();
}

//CBoWException::CBoWException(const char * szError_in) : CException(szError_in) {};

pragma_warning(pop)
