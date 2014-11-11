/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

/*
 * findBestMatchingDisjoint.h
 *
 *  Created on: 8/10/2009
 *      Author: tom
 */

#ifndef FINDBESTMATCHINGDISJOINT_H_
#define FINDBESTMATCHINGDISJOINT_H_

#include "findBestMatching.h"

class CFindBestMatchingDisjoint : public CFindBestMatching {

    class CPossibleMatches {
        int numMatches, anMatches[NN_MAX];
    public:

        CPossibleMatches() :
        numMatches(0) {
        }

        void reset() {
            numMatches = 0;
        }

        void add(int matchIdx) {
            if(IS_DEBUG) CHECK(numMatches >= NN_MAX, "CPossibleMatches::add: too many collisions!! (synthetic data or bad val for NN_MAX?)");
            anMatches[numMatches++] = matchIdx;
        }

        int operator[](int idx) const {
            if(IS_DEBUG) CHECK(idx >= numMatches, "CPossibleMatches::add: idx OOB");
            return anMatches[idx];
        }

        int getNumMatches() const {
            return numMatches;
        }

    };

    const int count;
    //const CvPoint2D64f * m0, * m1;
    const CPointIdentifiers & corr;

    std::vector<int> path;
    typedef std::vector<int>::iterator iter;

    static const int NO_MATCH = -1;

    /*struct cvPointCmp
     {
     bool operator()(const CvPoint2D64f & p0, const CvPoint2D64f & p1) const { return tooClose(p0,p1) ? 0 : p0<p1; };
     };*/
    typedef int TLocId;
    typedef map2_NF<TLocId, int> TPointIdxMap;
    TPointIdxMap leftIndices, rightIndices; //These store the index of the correspondence using each unique location
    typedef TPointIdxMap::iterator pointIter;
    map2_NF<std::pair<int, int>, int> pointIndices;
    int * anLeft;
    int * anRight;

    CPossibleMatches * rightIndicesMatchingLeft, *leftIndicesMatchingRight;

    bool extendRight(std::vector<int> &path) {
        int nLastLeftIdx = path[path.size() - 1]; //*(path.end() - 1);

        int nMatchesRight = rightIndicesMatchingLeft[nLastLeftIdx].getNumMatches();

        if(IS_DEBUG) CHECK(nMatchesRight == 0, "We've followed an edge the wrong way, or our map is corrupted");

        for (int nMatch = 0; nMatch < nMatchesRight; nMatch++) {
            int rightIdx = rightIndicesMatchingLeft[nLastLeftIdx][nMatch];

            if (anLeft[rightIdx] == NO_MATCH) {
                path.push_back(rightIdx);
                return true;
            }
        }
        for (int nMatch = 0; nMatch < nMatchesRight; nMatch++) {
            int rightIdx = rightIndicesMatchingLeft[nLastLeftIdx][nMatch];

            for (int oldRightIdx = path.size() - 2; oldRightIdx > 0; oldRightIdx -= 2)
                if (rightIdx == path[oldRightIdx])
                    //this right vertex is already in the path
                    goto next;

            //This right vertex is connected to something, but isn't in the path
            path.push_back(rightIdx);
            //if(IS_DEBUG) CHECK(anRight[anLeft[rightIdx]]);
            path.push_back(anLeft[rightIdx]);
            //The left vertex isn't in the path because we could only get there by following this right vertex back there

            if (extendRight(path)) return true; //We've extended it right again and eventually found an unmatched node on the right

            path.pop_back();
            path.pop_back();

next:
            ;
        }
        return false; //This vertex is a dead end
    }

    bool extendLeft(std::vector<int> &path) {
        int nLastRightIdx = path[path.size() - 1]; //*(path.end() - 1);

        int nMatchesLeft = leftIndicesMatchingRight[nLastRightIdx].getNumMatches();

        if(IS_DEBUG) CHECK(nMatchesLeft == 0, "We've followed an edge the wrong way, or our map is corrupted");

        for (int nMatch = 0; nMatch < nMatchesLeft; nMatch++) {
            int LeftIdx = leftIndicesMatchingRight[nLastRightIdx][nMatch];

            if (anRight[LeftIdx] == NO_MATCH) {
                path.push_back(LeftIdx);
                return true;
            }
        }
        for (int nMatch = 0; nMatch < nMatchesLeft; nMatch++) {
            int LeftIdx = leftIndicesMatchingRight[nLastRightIdx][nMatch];

            for (int oldLeftIdx = path.size() - 2; oldLeftIdx > 0; oldLeftIdx -= 2)
                if (LeftIdx == path[oldLeftIdx])
                    //this Left vertex is already in the path
                    goto next;

            //This Left vertex is connected to something, but isn't in the path
            path.push_back(LeftIdx);
            path.push_back(anRight[LeftIdx]);
            //The Right vertex isn't in the path because we could only get there by following this Left vertex back there

            if (extendLeft(path)) return true; //We've extended it Left again and eventually found an unmatched node on the Left

            path.pop_back();
            path.pop_back();

next:
            ;
        }
        return false; //This vertex is a dead end
    }

    double corrErr(int nLeftIdx, int nRightIdx, const double * adErrs) {
        int nCorrIdx = pointIndices[std::pair<int, int> (nLeftIdx, nRightIdx)];
        return adErrs[nCorrIdx];
    }

    bool findCheaperCycle_Left(std::vector<int> &path, double & dCostImprovementSoFar, const double * adErrs) {
        int nLastRightIdx = path[path.size() - 1]; //*(path.end() - 1);

        int nMatchesLeft = leftIndicesMatchingRight[nLastRightIdx].getNumMatches();

        if(IS_DEBUG) CHECK(nMatchesLeft == 0, "We've followed an edge the wrong way, or our map is corrupted");

        for (int nMatch = 0; nMatch < nMatchesLeft; nMatch++) {
            int LeftIdx = leftIndicesMatchingRight[nLastRightIdx][nMatch];

            if (anRight[LeftIdx] == NO_MATCH) {
                return false; //we've failed this time, maybe starting from here would work
            }
        }
        for (int nMatch = 0; nMatch < nMatchesLeft; nMatch++) {
            int LeftIdx = leftIndicesMatchingRight[nLastRightIdx][nMatch];

            double dErrRL = corrErr(LeftIdx, nLastRightIdx, adErrs);

            if (LeftIdx == path[0]) {
                //It's a cycle--Maybe it's an augmenting cycle?
                if (path.size() > 3) {
                    dCostImprovementSoFar -= dErrRL;
                    if (dCostImprovementSoFar > 0) {
                        path.push_back(LeftIdx);
                        return true;
                    }
                } //otherwise continue with other links
            } else {
                //extend
                for (int oldLeftIdx = path.size() - 2; oldLeftIdx > 0; oldLeftIdx -= 2)
                    if (LeftIdx == path[oldLeftIdx])
                        //this Left vertex is already in the path
                        goto next;

                int nNewRightIdx = anRight[LeftIdx];

                for (int oldRightIdx = path.size() - 1; oldRightIdx > 0; oldRightIdx -= 2)
                    if (nNewRightIdx == path[oldRightIdx])
                        //this Right vertex is already in the path
                        goto next;

                //This Left vertex is connected to something, but isn't in the path
                path.push_back(LeftIdx);
                dCostImprovementSoFar -= dErrRL;
                path.push_back(nNewRightIdx);
                double dErrLR = corrErr(LeftIdx, nNewRightIdx, adErrs);
                dCostImprovementSoFar += dErrLR;
                //The Right vertex isn't in the path because we could only get there by following this Left vertex back there

                if (findCheaperCycle_Left(path, dCostImprovementSoFar, adErrs)) return true; //We've extended it Left again and eventually matched back to the start

                dCostImprovementSoFar -= dErrLR;
                path.pop_back();
                dCostImprovementSoFar += dErrRL;
                path.pop_back();
            }

next:
            ;
        }
        return false; //This vertex is a dead end
    }

public:

    ~CFindBestMatchingDisjoint() {
        delete[] anLeft;
        delete[] anRight;
        delete[] rightIndicesMatchingLeft;
        delete[] leftIndicesMatchingRight;
    }

    CFindBestMatchingDisjoint(const CPointIdentifiers & pCorr) :
    count(pCorr.size()), corr(pCorr) {
        anLeft = new int[count];
        anRight = new int[count];
        rightIndicesMatchingLeft = new CPossibleMatches[count];
        leftIndicesMatchingRight = new CPossibleMatches[count];

        //		int numCollisionsRight=0, numCollisionsLeft=0;

        //Need to map CvPoint2D64f's to indices, and pairs of indices to indices into m0,m1
        for (int i = 0; i < count; i++) {
            const int leftIdx = leftIndices.initOrGet(corr[i].id1(), i);
            const int rightIdx = rightIndices.initOrGet(corr[i].id2(), i);

            //BLOBS cause this to fail!! Duplicate features in the same place
            if (IS_DEBUG && pointIndices.find(std::pair<int, int> (leftIdx, rightIdx)) != pointIndices.end())
                std::cout << leftIdx/*pLeft->first*/ << ',' << rightIdx /* pRight->first*/ << std::endl;
            if(IS_DEBUG) CHECK(pointIndices.find(std::pair<int, int> (leftIdx, rightIdx)) != pointIndices.end(), "Error indexing positions. This usually means 2 points in the same location have different point IDs");
            pointIndices.init(std::pair<int, int> (leftIdx, rightIdx), i);
            if(IS_DEBUG) CHECK((pointIndices[std::pair<int, int> (leftIdx, rightIdx)] != i), "Error indexing positions");
        }
        if(IS_DEBUG) CHECK((int) pointIndices.size() != count, "CFindBestMatching: Failed to add indices to map (duplicates?)");

        //cout << leftIndices.size() << " points indexed left\n";
        //cout << rightIndices.size() << " points indexed right\n";

        /*for(TPointIdxMap::const_iterator pIndex = leftIndices.begin(); pIndex != leftIndices.end(); pIndex++ )
         {
         cout << pIndex->first << "=>" << pIndex->second << endl;
         }*/
    }

    virtual bool supplyResiduals() const {
        return true;
    }

    virtual void refine(CMask & mask, int & nInliers, const double * adResiduals) {
        if(IS_DEBUG) CHECK((supplyResiduals() && !adResiduals) || (!supplyResiduals() && adResiduals), "refine: adResiduals should match whether they are needed");
        refine_int(mask, adResiduals);
        nInliers = mask.countInliers();
    }

    void refine_int(CMask & mask, const double * adResiduals) {
        for (int i = 0; i < count; i++) {
            rightIndicesMatchingLeft[i].reset();
            leftIndicesMatchingRight[i].reset(); //because i always >= idx fir this point
            anLeft[i] = NO_MATCH;
            anRight[i] = NO_MATCH;

            if (mask[i]) {
                if(IS_DEBUG) CHECK(leftIndices.find(corr[i].id1()) == leftIndices.end(), "Refine: point was never added to map");
                int leftIdx = leftIndices[corr[i].id1()];
                int rightIdx = rightIndices[corr[i].id2()];
                if(IS_DEBUG) CHECK(leftIdx > i || rightIdx > i, "refine: Index assumptions invalid");
                //int n=rightIndicesMatchingLeft[leftIdx].getNumMatches();
                rightIndicesMatchingLeft[leftIdx].add(rightIdx);
                //if(IS_DEBUG) CHECK( n+1!=rightIndicesMatchingLeft[leftIdx].getNumMatches(), "Matches aren't being counted properly");
                leftIndicesMatchingRight[rightIdx].add(leftIdx);

                mask[i] = 0;
            }
        }

        int nPathsFound = 0;
        do {
            nPathsFound = 0;

            for (int i = 0; i < count; i++) {
                if (rightIndicesMatchingLeft[i].getNumMatches() && anRight[i] == NO_MATCH) {
                    path.clear();
                    path.push_back(i);

                    if (extendRight(path)) //odd path
                    {
                        if(IS_DEBUG) CHECK(path.size() % 2 == 1, "Augmenting path isn't odd");
                        //iterate over path, update links
                        //if( path.size()>2 ) cout << "found path of length " << path.size() << endl;
                        for (int idx = 0; idx < (int) path.size() - 1; idx += 2) {
                            anLeft[path[idx + 1]] = path[idx];
                            anRight[path[idx]] = path[idx + 1];
                        }
                        nPathsFound++;
                    }
                }
                if (leftIndicesMatchingRight[i].getNumMatches() && anLeft[i] == NO_MATCH) {
                    path.clear();
                    path.push_back(i);

                    if (extendLeft(path)) //odd path
                    {
                        if(IS_DEBUG) CHECK(path.size() % 2 == 1, "Augmenting path isn't odd");
                        //iterate over path, update links
                        for (int idx = 0; idx < (int) path.size() - 1; idx += 2) {
                            anRight[path[idx + 1]] = path[idx];
                            anLeft[path[idx]] = path[idx + 1];
                        }
                        nPathsFound++;
                    }
                }
            }
            //cout << "Found " << nPathsFound << " paths\n";

        } while (nPathsFound > 0);
        int nInliersBefore = 0;
        double dErrBefore = 0;

        for (int i = 0; i < count; i++) {
            if (anLeft[i] != NO_MATCH) {
                if(IS_DEBUG) CHECK(anRight[anLeft[i]] != i, "Error finding consistent link set");
                        nInliersBefore++;
                if (adResiduals)
                        dErrBefore += corrErr(anLeft[i], i, adResiduals);
                }
            if(IS_DEBUG) CHECK(anRight[i] != NO_MATCH && anLeft[anRight[i]] != i, "Error finding consistent link set");
        }

        //http://valis.cs.uiuc.edu/~sariel/teach/courses/473/notes/27_matchings_notes.pdf
        //Now improve this matching by finding maximum weight augmenting paths by Dijkstra's alg
        if (adResiduals) {
            //For each unselected match:
            for (int i = 0; i < count; i++) {
                int nNumMatchesRight = rightIndicesMatchingLeft[i].getNumMatches();
                if (nNumMatchesRight > 1 && anRight[i] != NO_MATCH) {
                    int nRightMatch = anRight[i];
                    path.clear();
                    path.push_back(i);
                    path.push_back(nRightMatch);
                    double dTotalErrChange = corrErr(i, nRightMatch, adResiduals);

                    if (findCheaperCycle_Left(path, dTotalErrChange, adResiduals)) //maximum weight augmenting path
                    {
                        if(IS_DEBUG) CHECK(dTotalErrChange <= 0, "No improvement found");
                        if(IS_DEBUG) CHECK(path.size() % 2 == 0, "Augmenting cycle is odd");
                        if(IS_DEBUG) CHECK(path[path.size() - 1] != i, "Augmenting cycle isn't a cycle");
                        /*for (int idx = 0; idx < (int) path.size(); idx ++)
                                cout << path[idx] << ' ';
                        cout << "=path" << endl;
                        cout << "BEFORE:\n";

                        for (int idx = 0; idx < (int) path.size(); idx += 2)
                                cout << path[idx] << "L matches " << anRight[path[idx]] << endl;
                        for (int idx = 1; idx < (int) path.size(); idx += 2)
                                cout << path[idx] << "R matches " << anLeft[path[idx]] << endl;*/

                        for (int idx = 1; idx < (int) path.size() - 1; idx += 2) {
                            anLeft[path[idx]] = path[idx + 1];
                            anRight[path[idx + 1]] = path[idx];
                        }

                        /*cout << "AFTER:\n";
                        for (int idx = 0; idx < (int) path.size(); idx += 2)
                                cout << path[idx] << "L matches " << anRight[path[idx]] << endl;
                        for (int idx = 1; idx < (int) path.size(); idx += 2)
                                cout << path[idx] << "R matches " << anLeft[path[idx]] << endl;*/
                    }
                }
            }
        }
        int nInliersAfter = 0;
        double dErrAfter = 0;

        for (int i = 0; i < count; i++) {
            if (anLeft[i] != NO_MATCH) {
                if(IS_DEBUG) CHECK(anRight[anLeft[i]] != i, "Error finding consistent link set");
                int nCorrIdx = pointIndices[std::pair<int, int> (anLeft[i], i)];
                mask[nCorrIdx] = 1;
                nInliersAfter++;
                if (adResiduals)
                   dErrAfter += corrErr(anLeft[i], i, adResiduals);
            }
            if(IS_DEBUG) CHECK(anRight[i] != NO_MATCH && anLeft[anRight[i]] != i, "Error finding consistent link set");
        }
        if(IS_DEBUG) CHECK(nInliersAfter != nInliersBefore, "Inlier count has changed");
        if(IS_DEBUG) CHECK(dErrAfter > dErrBefore, "Total err has increased");

        //cout << nInliersAfter << "=after, err=" << dErrAfter << endl;

        /*
         #ifdef _DEBUG
         TPointIdxMap left;
         for (int i = 0; i < count; i++) {
         if(mask[i])
         {
         if(IS_DEBUG) CHECK(left.find(corr[i]->Location1()) != left.end(), "Error--have selected same L-point twice");
         left[corr[i]->Location1()] = i;
         }
         }

         #endif*/

    }
};

/*void testFindBestMatching()
 {
 CvPoint2D64f m0[5];
 CvPoint2D64f m1[5];

 m0[0] = cvPoint2D64f(0,0);
 m0[1] = cvPoint2D64f(0,0);
 m0[2] = cvPoint2D64f(1,1);
 m0[3] = cvPoint2D64f(1,1);
 m0[4] = cvPoint2D64f(23,1);

 m1[0] = cvPoint2D64f(0,0);
 m1[1] = cvPoint2D64f(1,1);
 m1[2] = cvPoint2D64f(0,0);
 m1[3] = cvPoint2D64f(1,1);
 m1[4] = cvPoint2D64f(0,0);

 uchar abMask[5]={1,1,1,1,1};

 CFindBestMatching match(5,m0,m1);
 match.refine(abMask);

 int inliers = 0;
 for(int i=0;i<5;i++)
 if(abMask[i]) inliers++;
 if(IS_DEBUG) CHECK(inliers != 2, "Failed to find good matching");
 }

 //A list of indices with which this point is incompatible
 class CPointCollisions {
 const int count;
 const CvPoint2D64f * m0, * m1;

 public:
 CPointCollisions(int count, const CvPoint2D64f * m0,
 const CvPoint2D64f * m1) :
 count(count), m0(m0), m1(m1) {

 }
 ;
 };*/

class CFindFirstMatchingDisjoint : public CFindBestMatchingDisjoint {
public:

    CFindFirstMatchingDisjoint(const CPointIdentifiers & pointIds) : CFindBestMatchingDisjoint(pointIds) {
    }

    virtual bool supplyResiduals() const {
        return false;
    }
};

#endif /* FINDBESTMATCHINGDISJOINT_H_ */
