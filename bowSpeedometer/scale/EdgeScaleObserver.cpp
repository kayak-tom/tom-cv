/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * EdgeScaleObserver.cpp
 *
 *  Created on: 24/08/2009
 *      Author: tom
 */
#include <boost/filesystem/operations.hpp>   // includes path.hpp

#include "EdgeScaleObserver.h"

//Todo: clean up and library-fy these
#include "geom/pointAlignment.h"
#include "../../BoWSLAM/BoWSLAM_Main/threeDPointCollection.h"

#include "imageSource/imageSourceFromDir.h"
#include "image/imageAccess.h"
#include "stdio.h"

#include "util/opencv.h"
#include "util/opencv_highgui.h"

#include "ScaleMLE.h"

#ifndef USE_OLD_OPENCV
#include "opencv2/opencv.hpp"
#endif

using namespace std;

//extern int IM_WIDTH, IM_HEIGHT, IM_CHANNELS;
CEdgeScaleObserver::~CEdgeScaleObserver() {

}
void getFile(int nId, char * szFile) {
    sprintf(szFile, "objects/%d.png", nId);
}
int val(int & nShift) {
    nShift++;
    if (nShift == 5)
        nShift = 1;

    return 50 + (205 * nShift) / 5;
}
CvScalar getCol(int nId) {
    int nShift = nId % 5;
    nId %= 6;
    nId++; //3 bits, 1 or 2 set
    int R = (nId & 1) ? val(nShift) : 0;
    int G = (nId & 2) ? val(nShift) : 0;
    int B = (nId & 4) ? val(nShift) : 0;
    int m = max(max(R, G), B);
    m = 255 - m;

    return CV_RGB(R + m, G + m, B + m);
}
void CEdgeScaleObserver::markOneObjectObservation(int id, IplImage * pFrameCol1, const CBoWSpeedo::CBoWObjectOccurance *pObservedObj, IplImage * pFrameCol2) const {
    static CvFont font_small;
    REPEAT(1, cvInitFont(&font_small, CV_FONT_HERSHEY_PLAIN, 1.5, 1.5, 0, 2));

    char szId[10];
    sprintf_s(szId, 10, "%d", id);
    CvScalar col = getCol(id);
    if (pFrameCol1) {
        cvCircle(pFrameCol1, locToCvPoint(pObservedObj->l1Im1()), 2, col, 2);
        cvCircle(pFrameCol1, locToCvPoint(pObservedObj->l2Im1()), 2, col, 2);
        cvLine(pFrameCol1, locToCvPoint(pObservedObj->l1Im1()), locToCvPoint(pObservedObj->l2Im1()), col, 2);
        cvPutText(pFrameCol1, szId, locToCvPoint(pObservedObj->l1Im1()), &font_small, col);
    }
    if (pFrameCol2) {
        cvCircle(pFrameCol2, locToCvPoint(pObservedObj->l1Im2()), 2, col, 2);
        cvCircle(pFrameCol2, locToCvPoint(pObservedObj->l2Im2()), 2, col, 2);
        cvLine(pFrameCol2, locToCvPoint(pObservedObj->l1Im2()), locToCvPoint(pObservedObj->l2Im2()), col, 2);
        cvPutText(pFrameCol2, szId, locToCvPoint(pObservedObj->l1Im2()), &font_small, col);
    }
}
void CEdgeScaleObserver::drawObjects(const TObsCounter & bestObjects, CBoWSpeedo::CObservedObjectsVec & observedObjects) const {
    IplImage * pFrame1 = 0;
    IplImage * pFrame2 = 0;
    IplImage * pFrameCol1 = 0;
    IplImage * pFrameCol2 = 0;
    char szFile1[100], szFile2[100];
    getFile(observedObjects.id1(), szFile1);
    getFile(observedObjects.id2(), szFile2);

    bool bDrawnObjects = false;

    for (CBoWSpeedo::CObservedObjectsVec::const_iterator ppObservedObj = observedObjects.begin(); ppObservedObj != observedObjects.end(); ppObservedObj++) {
        const CBoWSpeedo::CBoWObjectOccurance * pObservedObj = *ppObservedObj;
        if (bestObjects.find(pObservedObj->object()) != bestObjects.end()) {
            if (!pObservedObj->l1Im1().zero()) {
                if (!bDrawnObjects) {
                    pFrameCol1 = cvCreateImage(IM_PARAMS.SIZE(), IPL_DEPTH_8U, 3);
                    pFrameCol2 = cvCreateImage(IM_PARAMS.SIZE(), IPL_DEPTH_8U, 3);

                    if (boost::filesystem::exists(szFile1)) {
                        pFrame1 = cvLoadImage(szFile1);
                    } else {
                        pFrame1 = cvCreateImage(IM_PARAMS.SIZE(), IPL_DEPTH_8U, IM_PARAMS.IM_CHANNELS);
                        pImageLoader->loadImage(observedObjects.id1(), pFrame1);
                    }

                    if (boost::filesystem::exists(szFile2)) {
                        pFrame2 = cvLoadImage(szFile2);
                    } else {
                        pFrame2 = cvCreateImage(IM_PARAMS.SIZE(), IPL_DEPTH_8U, IM_PARAMS.IM_CHANNELS);
                        pImageLoader->loadImage(observedObjects.id2(), pFrame2);
                    }
                    //pImageLoader->loadImage(observedObjects.id2(), pFrame2);

                    if (pFrame1->nChannels == 1)
                        cvCvtColor(pFrame1, pFrameCol1, CV_GRAY2RGB);
                    else
                        //cvConvertImage(pFrame1, pFrameCol1);
                        cvCopy(pFrame1, pFrameCol1);

                    if (pFrame2->nChannels == 1)
                        cvCvtColor(pFrame2, pFrameCol2, CV_GRAY2RGB);
                    else
                        cvCopy(pFrame2, pFrameCol2);
                    //						cvConvertImage(pFrame1, pFrameCol1);

                    bDrawnObjects = true;

                }

                int id = bestObjects.find(pObservedObj->object())->second;

                if(IS_DEBUG) CHECK(!pFrame1 || !pFrame2, "Missing images");
                cout << "Drawing Obj " << id << " which has log-size " << pObservedObj->object()->estimatedObjLogSize() << " and var " << pObservedObj->object()->estimatedObjLogVar() << " from " << pObservedObj->object()->count() << endl;
                pObservedObj->object()->pp();
                markOneObjectObservation(id, pFrameCol1, pObservedObj, pFrameCol2);
            }
        }
    }

    if (bDrawnObjects) {
        cvSaveImage(szFile1, pFrameCol1);
        cvSaveImage(szFile2, pFrameCol2);
        cvReleaseImage(&pFrame1);
        cvReleaseImage(&pFrame2);
        cvReleaseImage(&pFrameCol1);
        cvReleaseImage(&pFrameCol2);
    }
}
void CEdgeScaleObserver::resetScales() {
    for (TEdgeVec::iterator ppEdge = vEdges.begin(); ppEdge != vEdges.end(); ppEdge++) {
        CEdge * pEdge = *ppEdge;
        if (pEdge)
            pEdge->reset();
    }
}
void CEdgeScaleObserver::markObjects(const int nIdLast, const int nIdCurrent, IplImage * pCurrentFrame) const {
    if (vEdges.size() == 0)
        return;

    const CEdge * pCorrectEdge = 0;
    int nEdgeIdx = -1;
    for (TEdgeVec::const_iterator ppEdge = std::max<TEdgeVec::const_iterator > (vEdges.begin(), vEdges.end() - 10); ppEdge != vEdges.end(); ppEdge++) {
        const CEdge * pEdge = *ppEdge;
        if (pEdge) {
            if ((int) pEdge->id1() == nIdLast && (int) pEdge->id2() == nIdCurrent) {
                pCorrectEdge = pEdge;
                nEdgeIdx = ppEdge - vEdges.begin();
                break;
            }
        }
    }

    if (!pCorrectEdge) {
        cout << "Failed to find correct edge...Error?\n";
    } else {
        CBoWSpeedo::CObservedObjectsVec * pObservedObjects = pCorrectEdge->getObservations();
        //observedObjects.init(nIdLast, nIdCurrent);
        //measureObjectOccurancesInEdge(observedObjects, const_cast<CEdge *>(pCorrectEdge), nEdgeIdx, 0);
        if (!pObservedObjects) {
            cout << "Error? No cached observations for edge\n";
        } else {
            for (CBoWSpeedo::CObservedObjectsVec::const_iterator ppObservedObj = pObservedObjects->begin(); ppObservedObj != pObservedObjects->end(); ppObservedObj++) {
                const CBoWSpeedo::CBoWObjectOccurance * pObservedObj = *ppObservedObj;
                //if(bestObjects.find(pObservedObj->object()) != bestObjects.end())
                {
                    if (pObservedObj->object()->isInit())
                        if (!pObservedObj->l1Im2().zero()) {
                            //int id=bestObjects.find(pObservedObj->object())->second;
                            int id = pObservedObj->object()->getObjectId();
                            markOneObjectObservation(id, 0, pObservedObj, pCurrentFrame);
                        }
                }
            }
        }
    }
}
void CEdgeScaleObserver::printScales() const {
    double dSum = 0;
    int nCount = 0;
    for (TEdgeVec::const_iterator ppEdge = vEdges.begin(); ppEdge != vEdges.end(); ppEdge++) {
        const CEdge * pEdge = *ppEdge;
        if (pEdge && pEdge->hasScale()) {
            double dScale = pEdge->SLAMscale().scaleEstimate();
            double dTrueScale = pEdge->id2() - pEdge->id1();
            double dDeviation = dScale / dTrueScale;
            if (dScale > 0) {
                dSum += dDeviation;
                nCount++;
            }
            cout << dDeviation << ',';
        }
    }
    if (nCount > 0) {
        double dMean = dSum / nCount, dVar = 0;
        ;

        for (TEdgeVec::const_iterator ppEdge = vEdges.begin(); ppEdge != vEdges.end(); ppEdge++) {
            const CEdge * pEdge = *ppEdge;
            if (pEdge && pEdge->hasScale()) {
                double dScale = pEdge->SLAMscale().scaleEstimate();
                double dTrueScale = pEdge->id2() - pEdge->id1();
                double dDeviation = dScale / dTrueScale;
                if (dScale > 0) {
                    dVar += sqr(dMean - dDeviation);
                }
            }
        }

        dVar /= nCount;
        cout << endl << "Scales mean = " << dMean << " var = " << dVar << endl;
    }
}
void CEdgeScaleObserver::findObjects(const CBoWSpeedo & bow, CBoWSpeedo::CObjectFinder & objectFinder) {
    //CBoWSpeedo::CBoWSpeedometer::CObjectFinder localObjectFinder;
    CBoWSpeedo::TObjectsAndLocations localObjectFinder;
    for (TEdgeVec::const_iterator ppEdge = vEdges.begin(); ppEdge != vEdges.end(); ppEdge++) {
        const CEdge * pEdge = *ppEdge;
        if (pEdge) {
            localObjectFinder.clear();

            bow.getObjectsAndLocationsInIntersection(localObjectFinder, pEdge->id1(), pEdge->id2());

            pEdge->findReconstructedPairs(localObjectFinder);

            for (CBoWSpeedo::TObjectsAndLocations::const_iterator ppLocalObj = localObjectFinder.begin(); ppLocalObj != localObjectFinder.end(); ppLocalObj++) {
                const CBoWSpeedo::CObjAndLocations * pLocalObj = *ppLocalObj;
                for (CBoWSpeedo::TObjectsAndLocations::const_iterator ppLocalObj2 = ppLocalObj + 1; ppLocalObj2 != localObjectFinder.end(); ppLocalObj2++) {
                    const CBoWSpeedo::CObjAndLocations * pLocalObj2 = *ppLocalObj2;
                    if (pLocalObj->isReconstructed() && pLocalObj2->isReconstructed())
                        objectFinder.addCoOccurance(pLocalObj->pW1, pLocalObj2->pW1, pLocalObj->nReconCount, pLocalObj2->nReconCount);
                }
            }
        }
    }
}
void CEdgeScaleObserver::measureObjectOccurancesInEdge(CBoWSpeedo::CObservedObjectsVec & observedObjects, CEdge *pEdge, const int nEdgeIdx, TObsCounter *pCountObservations) const {
    pSpeedo->getObjectsInBothFrames(observedObjects);
    for (CBoWSpeedo::CObservedObjectsVec::iterator ppObservedObj = observedObjects.begin(); ppObservedObj != observedObjects.end(); ppObservedObj++) {
        CBoWSpeedo::CBoWObjectOccurance * pObservedObj = *ppObservedObj;
        double dObjLogLength = pEdge->measureObject_LN(pObservedObj); //might well be 0 because none of these points were reconstructed

        if (dObjLogLength > -HUGE) {
            pObservedObj->observeOccurance(dObjLogLength, pEdge->SLAMscale(), nEdgeIdx);
            if (pCountObservations)
                (*pCountObservations)[pObservedObj->object()]++;
        }
    }
}
void CEdgeScaleObserver::reobserveScales(bool bDrawObjects) {
    resetScales();

    cout << "Scales before: ";
    printScales();

    //Find the new set of objects in frame-pairs
    CDynArray<CBoWSpeedo::CObservedObjectsVec> aaObservedObjects(vEdges.size());

    //Also count the most successfully observed objects
    TObsCounter countObservations, bestObservations, *pCountObservations = (bDrawObjects ? &countObservations : 0);

    int nEdgeIdx = 0;
    for (TEdgeVec::iterator ppEdge = vEdges.begin(); ppEdge != vEdges.end(); ppEdge++, nEdgeIdx++) {
        CEdge * pEdge = *ppEdge;

        if (pEdge && pEdge->hasScale()) {
            CBoWSpeedo::CObservedObjectsVec & observedObjects = aaObservedObjects[nEdgeIdx];
            observedObjects.init(pEdge->id1(), pEdge->id2());
            measureObjectOccurancesInEdge(observedObjects, pEdge, nEdgeIdx, pCountObservations);
        }
    }

    //For each object calculate mean+var
    pSpeedo->updateObjectSizeDistns();

    if (bDrawObjects) {
        boost::filesystem::remove_all("objects");
        boost::filesystem::create_directory("objects");

        //Copy to vector and sort to find best objects
        int id = 0, max = 0;
        for (TObsCounter::iterator pObs = countObservations.begin(); pObs != countObservations.end(); pObs++)
            if (pObs->first->isInit() && pObs->second > max)
                max = pObs->second;
        for (TObsCounter::iterator pObs = countObservations.begin(); pObs != countObservations.end(); pObs++)
            if (pObs->first->isInit() && pObs->second > (max * 3 / 4)) {
                bestObservations.insert(pair<const CBoWSpeedo::CBoWObject *, int>(pObs->first, id++));
            }


        if (bestObservations.size() == 0) {
            for (TObsCounter::iterator pObs = countObservations.begin(); pObs != countObservations.end(); pObs++)
                if (pObs->first->isInit() && pObs->second > 10)
                    bestObservations.insert(pair<const CBoWSpeedo::CBoWObject *, int>(pObs->first, id++));
        }

        /*if(bestObservations.size() == 0)
        {
                for(TObsCounter::iterator pObs = countObservations.begin(); pObs != countObservations.end(); pObs++)
                        if(pObs->second > 4)
                                bestObservations.insert(pair<const CBoWSpeedo::CBoWObject *, int>( pObs->first, id++ ));
        }*/
    }

    //Apply a scale update to every edge

    nEdgeIdx = 0;
    for (TEdgeVec::iterator ppEdge = vEdges.begin(); ppEdge != vEdges.end(); ppEdge++, nEdgeIdx++) {
        CEdge * pEdge = *ppEdge;
        if (pEdge) {
            pEdge->updateScale(aaObservedObjects[nEdgeIdx], nEdgeIdx);

            if (bDrawObjects)
                drawObjects(bestObservations, aaObservedObjects[nEdgeIdx]);
        }
    }

    cout << "Scales after: ";
    printScales();
}
void CEdgeScaleObserver::removeEdge(CEdge * pEdge) {
    //vEdges.erase(pEdge); //its probably near the back
    //Just zero it so we don't invalidate indices
    for (TEdgeVec::iterator ppEdge = vEdges.end() - 1; ppEdge >= vEdges.begin(); ppEdge--) {
        if (*ppEdge == pEdge) {
            *ppEdge = 0;
            break;
        }
    }
}
void CEdgeScaleObserver::addEdge(CEdge * pEdge) {
    CBoWSpeedo::CObservedObjectsVec * pObservedObjects = new CBoWSpeedo::CObservedObjectsVec(pEdge->id1(), pEdge->id2());
    //pSpeedo->getObjectsInBothFrames(aObservedObjects);
    const int nEdgeIdx = vEdges.size();
    vEdges.push_back(pEdge);

    //Will update parameters for this object...and the measurement of the object in this frame will act to damp any update in scale
    measureObjectOccurancesInEdge(*pObservedObjects, pEdge, nEdgeIdx, 0);

    if (pObservedObjects->size() > 5)
        cout << "Observed " << pObservedObjects->size() << " objects in edge " << pEdge->id1() << ',' << pEdge->id2() << endl;

    pEdge->setObservedObjVec(pObservedObjects);

    pEdge->updateScale(*pObservedObjects, nEdgeIdx);
}
void CEdge::updateScale(CBoWSpeedo::TObservedObjectsVec & aObservedObjects, int nEdgeIdx) {
    CScaleFromObservingAndMeasuringKnownObjs scaleObserver;

    for (CBoWSpeedo::TObservedObjectsVec::iterator ppObservedObj = aObservedObjects.begin(); ppObservedObj != aObservedObjects.end(); ppObservedObj++) {
        CBoWSpeedo::CBoWObjectOccurance * pObservedObj = *ppObservedObj;
        if (pObservedObj->object()->isInit()) {
            double dObservedScale = NO_SCALE();

            if (nEdgeIdx == -1) {
                dObservedScale = measureObject_LN(pObservedObj); //now should actually already have measured object
            } else {
                //This object remembers the log-measurement plus scale from edge, and var from this edge
                if (SLAMscale().hasScale()) {
                    //This is the baseline 1 scale
                    dObservedScale = pObservedObj->object()->rememberScaleObservation(nEdgeIdx);

                    /*if(dObservedScale != NO_SCALE())
                    {
                            dObservedScale -= SLAMscale().getD();
                    }*/
                }
            }

            if (dObservedScale != NO_SCALE()) {
                CBoWSpeedo::CBoWObjectOccurance * pObservedObj = *ppObservedObj;

                double dObjScale = pObservedObj->object()->estimatedObjLogSize();
                double dObjVar = pObservedObj->object()->estimatedObjLogVar();

                scaleObserver.addObjectMeasurement(dObservedScale, dObjScale, dObjVar);
            }
        }
    }

    double dScaleMLE = 0, dVar = -1;
    scaleObserver.computeMeanVar(dScaleMLE, dVar);

    if (dVar > 0)
        updateScaleWithOR(dScaleMLE, dVar);
}
void CBoWSpeedo::CBoWObject::calculateScale() {
    //if(IS_DEBUG) CHECK(bInit, "Already initialised"); Might be re-computing after making more observations
    const int MIN_OBSERVATIONS = 3;
    if (scaleObservations.size() < MIN_OBSERVATIONS) {
        //std::cout << "Object " << getObjectId() << " has insufficient observations (" << scaleObservations.size() << ")\n";
        bInit = false;
        return; //Don't want to fit 1 scale to essentially 1 other
    }

    //check we actually have 3 good measurements...
    int nGood = 0;
    for (TLengths::const_iterator pScale = scaleObservations.begin(); pScale != scaleObservations.end(); pScale++) {
        if (pScale->hasSLAMLength()) {
            nGood++;
            if (nGood == MIN_OBSERVATIONS)
                break;
        }
    }
    
    const bool bVerbose = false;

    if (nGood >= MIN_OBSERVATIONS && CScaleMLE::getBayesianScaleParams_NG(scaleObservations, dMean, dVariance, false)) {
        if (bVerbose) {
            std::cout << "Object " << getObjectId() << " created with " << scaleObservations.size() << " observations, log mean=" << dMean << " s.d.=" << sqrt(dVariance) << " mean=" << exp(dMean) << std::endl;

            for (TLengths::iterator pLength = scaleObservations.begin(); pLength != scaleObservations.end(); pLength++) {
                cout << pLength->getSLAMLength() << ",";
            }
            cout << endl;
        }

        bInit = true;
    } else {
        std::cout << "Object " << getObjectId() << " not init from " << scaleObservations.size() << " (" << nGood << ") observations" << std::endl;
        bInit = false;
        //CScaleMLE::getMLEScaleParams_2ndOrderGD(scaleObservations, dMean, dVariance, true);
    }

    /*dMean = 0; dVariance = 0;
    dSumWeights=0;
    for(TLengths::iterator pLength = scaleObservations.begin(); pLength != scaleObservations.end(); pLength++)
    {
            double dWeight = pLength->getWeight();
            dMean += pLength->getLength()*dWeight;
            dSumWeights += dWeight;
    }
    if(zero(dSumWeights))
    {
            std::cout << "Warning: object " << getObjectId() << " has insufficient total-weight (" << scaleObservations.size() << " objects)\n";
            return; //Don't want to fit 1 scale to essentially 1 other
    }

    dVariance = 1.0 / dSumWeights;
    dMean *= dVariance;
    / * Old weighted var: double dSumWeights_inv = 1.0/dSumWeights;
    dMean *= dSumWeights_inv;

    double dV2 = 0;
    for(TLengths::iterator pLength = scaleObservations.begin(); pLength != scaleObservations.end(); pLength++)
    {
            double dWeight = pLength->getWeight() * dSumWeights_inv;
            dVariance += dWeight * sqr(pLength->getLength() - dMean);
            dV2 += sqr(dWeight);
    }

    dVariance *= 1/(1-dV2);* /
    bInit = true;
    //dVariance += sqr(dMean); //Variance as big as mean is not entirely unexpected.
    std::cout << "Object " << getObjectId() << " created with " << scaleObservations.size() << " observations, log mean=" << dMean << " s.d.="<< sqrt(dVariance)<< " mean="<< exp(dMean) << std::endl;
    bInit = true;*/
}
