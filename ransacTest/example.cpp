/* ****** See ransacAndRefine/main.cpp for a simple example of use of the ransac library. This file is mostly for me to test things. *******

Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * example.cpp
 *
 *  Created on: 7/10/2009
 *      Author: tom
 */
#include <iostream>
#include "ransac/ransac.h"
#include "ransac/findBestMatchingDisjoint.h"
#include "ransac/essentialMat_Eigen.h"
#include "ransac/openCV8PtEssentialMatModified.h"
#include <boost/smart_ptr.hpp>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include "ransac/ransacParams.h"
#include "geom/geom.h"
#include "geom/geom_eigen.h"
#include "time/SpeedTest.h"
#include "util/stats.h"
#include <fstream>
#include <Eigen/SVD>
#include "ransac/houghForE.h"
#include "ransac/refineEOnRTManifold.h"

//USING_PART_OF_NAMESPACE_EIGEN
using namespace Eigen;
using namespace std;
using namespace boost;

/* Alternative way to do ransac, a few more options exposed...
void ransacForE(T2dPoints & points0, T2dPoints & points1, CInlierProbs & adPriorProbs, CPointIdentifiers & pointIds, CModel & E)
{
        const int nSampleSize = 5, nMaxIters = 1000;
        const int nCount = points0.size();
        const double pProbSuccess = 0.99, d8ptCutoff = 5, dInlierThresh= 0.01;

        CMask inliers(nCount);

    scoped_ptr<CSampler> pSampler(new CBaySACSampler_discretePP(nSampleSize, points0, points1, adPriorProbs, true));
    scoped_ptr<CRansacTerminator> pTerminator(new CSimpleTerminator(nCount));
    scoped_ptr<CIterTerminator> pIterTerminator(new CPPIterTerminator(nMaxIters, adPriorProbs, nSampleSize, pProbSuccess));
    scoped_ptr<CModelHypothesiser> pHypothesise(new C5ptEssentialMat);
    scoped_ptr<CModelRefiner> pRefiner(new COpenCV8PtEssentialMatModified(d8ptCutoff));
    scoped_ptr<CInlierCounter> pInlierCounter(new CInlierCounter(dInlierThresh, points0, points1));
    scoped_ptr<CFindBestMatching> pRefineMatches(new CFindBestMatchingDisjoint(pointIds));

    int nInliers = doRansac(points0, points1, pSampler.get(), pHypothesise.get(), pRefiner.get(), pTerminator.get(), pIterTerminator.get(), pInlierCounter.get(), pRefineMatches.get(), E, inliers, 1, false);
    cout << nInliers << " inliers\n";
}*/
double threeDPointSSD(const T3dPoints & points3d, const Matrix<double, 4, Dynamic> & points) {
    double dSS = 0;
    for (int i = 0; i < points3d.size(); i++) {
        C3dPoint original(points(0, i), points(1, i), points(2, i)),
                recon = points3d[i];
        original.normalise();
        recon.normalise();
        if (dotproduct(original, recon) < 0)
            recon = -recon;
        dSS += (recon - original).sum_square();
        //cout << original << '-' << recon << endl ;
    }
    return dSS;
}

extern int MODELS; //Model counter, used in RANSAC

/*class CRTErrorStats {
    CStats statsR, statsT;
    int nFail, nLargeErrors, nTotalRuns;
    ostream & results;
public:
    const string strName;

    CRTErrorStats(ostream & results, const char * szName) : nFail(0), nLargeErrors(0), nTotalRuns(0), results(results), strName(szName) {
    }

    void addFail() {
        nFail++;
        nTotalRuns++;
    }

    void addRT(double dErrInR, double dErrInT) {
        nTotalRuns++;
        double dSuccessThresh = 0.25;
        if (dErrInR > dSuccessThresh)
            nLargeErrors++;
        //else {
        statsR.add(dErrInR);
        statsT.add(dErrInT);
        //}
    }

    ~CRTErrorStats() {
        pp(results);
    }

    void pp(ostream & os) {
        if (nTotalRuns <= 1)
            return;

        double dSuccessRate = 1 - ((double) nFail / (double) nTotalRuns);
        double dSolutionsWithErrors = 1 - ((double) (nLargeErrors + nFail) / (double) nTotalRuns);
        std::ostringstream ss;

        ss << "\"" << strName << "\"\tRuns:" << nTotalRuns << '\t' << statsR.mean() << '\t' << statsR.median() << '\t' << statsT.mean() << '\t' << statsT.median() << "\tSuccesses" << dSuccessRate << "\tUsefulSuccesses" << dSolutionsWithErrors << endl;
        cout << ss.str();
        os << ss.str();
    }
};*/

void printErrors(CRTErrorStats & stats, const C3dRotation & R_exact, const C3dPoint & T_exact, const C3dRotation & R_fromRANSAC, const C3dPoint & T_fromRANSAC) {
    double dErrTheta = diff(R_fromRANSAC, R_exact);
    double dErrT = angle(T_fromRANSAC, T_exact);
    cout << dErrTheta << "=error in angle (" << stats.strName << ")" << endl;
    cout << dErrT << "=error in t (" << stats.strName << ")" << endl;
    stats.addRT(dErrTheta, dErrT);
}

class COutFile : public std::ofstream {
public:

    COutFile(const char * FN) : std::ofstream(FN) {
    }

    ~COutFile() {
        close();
    }
};

int compareMasks(const CMask & mask, const CMask & maskExact, bool bVerbose) {
    int nInliers = mask.countInliers();
    int nCorrectInliers = 0, nIncorrectInliers = 0;
    int nMissedInliers = 0;
    for (int i = 0; i < (int) mask.size(); i++) {
        if (mask[i] && maskExact[i])
            nCorrectInliers++;
        else if (mask[i] && !maskExact[i])
            nIncorrectInliers++;
        else if (!mask[i] && maskExact[i])
            nMissedInliers++;
    }
    const int nTotalInliers = nCorrectInliers + nMissedInliers;
    if(IS_DEBUG) CHECK(nTotalInliers != maskExact.countInliers(), "Inlier counting failed");
    if (bVerbose) {
        cout << "Correct Inliers" << nCorrectInliers << " (" << 100. * nCorrectInliers / (double) (nTotalInliers) << "% of true inliers)\n";
        cout << "Incorrect Inliers" << nIncorrectInliers << " (" << 100. * nIncorrectInliers / (double) (nInliers) << "% of inliers found)\n";
        cout << "Missed Inliers" << nMissedInliers << " (" << 100. * nMissedInliers / (double) (nTotalInliers) << "% of true inliers)\n";
    }
    return nIncorrectInliers;
}

int mainFindE(bool bUprightCam, const int nCount, const int nIters, double dNoiseSD /*e.g. = 0.000005*/, const double dInlierRate, bool bVerbose, ostream & results, ostream & resultsTrans) //if bUprightCam cam is assumed approximately upright
{
    MODELS = 0;

    CRANSACParams PARAMS(0, 0); //Create parameter object
    PARAMS.VERBOSE = false; //IS_DEBUG;
    PARAMS.PROB_SUCCESS = 0.9;
    PARAMS.E_INLIER_THRESH_PX = 0.01; //Low because not using uncalibrated coordinates -- camera has FOV of 1
    PARAMS.RANSACSampler = CRANSACParams::eRANSAC;
    PARAMS.TOPDOWN_ITERS = 0;
    PARAMS.RANSACTerminator = CRANSACParams::eSimpleTerminator; // CRANSACParams::eWaldSAC;// CRANSACParams::eBrownianBridge;
    PARAMS.MAX_ITERS = 1000;
    PARAMS.HypothesiseAlg = CRANSACParams::e5Pt_GradientDesc; //e5PtE, e7PtE, e7PtF, e2PtE, eMCE, e1PtE, e5Pt_GradientDesc, e4Pt_GradientDesc, e3Pt_GradientDesc e7PtEFast
    PARAMS.RANSACIterTerminator = CRANSACParams::eTerminateOnPropInliers;

    CHoughForE * pHough = CHoughForE::makeHough();

    T2dPoints points0(nCount), points1(nCount);
    CPointIdentifiers pointIds(nCount);
    CInlierProbs adPriorProbs(nCount);

    CStats stats8pt, statsRefineP, statsRefine5pt, statsRefineHZ;
    CStats stats8ptTrans, statsRefinePTrans, statsRefine5ptTrans, statsRefineHZTrans;

    COutFile newResults("refinerAlgs.tsv");
    CRTErrorStats RANSACstats(newResults, "RANSAC"), LSMaskStats(newResults, "LS on inliers"), LSStatsLinear(newResults, "Linear LS on inliers"), robustStats(newResults, "Robust optimisation"), robustStatsLinear(newResults, "Robust linear optimisation"),
            robustStatsLinearF(newResults, "Robust linear optimisation F");

    //ofstream testOut("testDepth.tsv");
    double dTotalTime = 0;

    CStopWatch s;
    s.startTimer();
    const int MAX_NUM_OUTLIERS = (int)(0.9 * nCount);

    int nSuccesses = 0;

    for (int iter = 0; iter < nIters; iter++) {
        //Choose a random essential mat:
        double dAngle = CRandom::Uniform(-.75, .75); //Will put inliers behind cam if rotate too far
        Vector3d axis;
        axis << CRandom::Normal(), CRandom::Normal(), CRandom::Normal();

        Vector3d temp, translation;
        translation << CRandom::Normal(), CRandom::Normal(), CRandom::Normal();

        enum eTransModes {
            eRotAndTrans, ePureRot, eForwardMotion, eSidewaysMotion
        };
        eTransModes eMode = eRotAndTrans;

        if (eMode == ePureRot)
            translation *= 0.001;
        else if (eMode == eForwardMotion) {
            translation << 0, 0, 1;
            dAngle *= 0.01;
        } else if (eMode == eSidewaysMotion) {
            translation(2) *= 0.01;
            dAngle *= 0.01;
        }

        if (bUprightCam) {
            axis(1) = 20; //Makes approx. upright
            translation(1) *= 0.05; //Constrains approx. to plane
        }
        
        //Really degenerate if plane rotates away too far
        const bool bPlanarPoints = false;
        if(bPlanarPoints)
        {
            axis(2) = 10;
        }
        
        translation.normalize();
        C3dPoint T_exact(translation);
        T_exact.normalise();
        
        

        C3dRotation R_exact(axis, dAngle);

        Matrix3d rotMat;
        R_exact.asMat(rotMat);

        if (bVerbose) cout << "#######################################\nRotation=" << R_exact << ", trans=" << translation.transpose() << endl;

        (R_exact * translation).asVector(temp);
        CHECK((temp - rotMat * translation).squaredNorm() > 0.0000001, "Rotation test failed");
        C3dRotation q_temp(rotMat);
        (q_temp * translation).asVector(temp);
        CHECK((temp - rotMat * translation).squaredNorm() > 0.0000001, "Convert mat back to quat failed");

        Matrix3d E_exact;
        makeE(R_exact, translation, E_exact);

        Matrix<double, 3, 4 > cam2;
        cam2.block < 3, 3 > (0, 0) = rotMat;
        cam2.block < 3, 1 > (0, 3) = translation;
        Matrix<double, 3, 4 > cam1;
        cam1.block < 3, 3 > (0, 0).setIdentity();
        cam1.block < 3, 1 > (0, 3) = 0 * translation;

        Matrix<double, 4, Dynamic> points = Matrix<double, 4, Dynamic>::Random(4, nCount);
        if(bPlanarPoints)
        {
            points.row(2).array() = points.row(2).array()*0.1 + 3; //Force PLANAR
            //points.row(2).array().setConstant(3);
        }
        else
            points.row(2).array() += 3; //Force in front

        points.row(3).setOnes();

        //points.Random();
        Matrix<double, 3, Dynamic> pointsLeft = cam1*points;
        Matrix<double, 3, Dynamic> pointsRight = cam2*points;

        double dAssumedInlierProb = 0.8;
        int nReusedPointCount = 1;

        for (int i = 0; i < nCount; i++) {
            pointsLeft.col(i) /= pointsLeft(2, i);
            pointsRight.col(i) /= pointsRight(2, i);
        }

        //what's the range of points?
        double dApproxRangeRads = pointsLeft.row(0).maxCoeff() - pointsLeft.row(0).minCoeff();

        //What's the variance in radians?
        double dMeanX = pointsLeft.row(0).sum() / nCount;
        double dMeanY = pointsLeft.row(1).sum() / nCount;
        double dVar = 0;

        for (int i = 0; i < nCount; i++) {
            dVar += sqr(pointsLeft(0, i) - dMeanX);
            dVar += sqr(pointsLeft(1, i) - dMeanY);
        }
        dVar /= nCount;
        if (bVerbose) cout << "Range (rads): " << dApproxRangeRads << ", sd = " << sqrt(dVar) << endl;
        //dNoiseSD /= sqrt(dVar) ;

        //Also check we have enough point movement to calculate a relative position
        Matrix<double, 3, Dynamic> pointsTotalMovement = pointsLeft - pointsRight;
        double dAverageMotion = 0;
        for (int i = 0; i < nCount; i++) {
            dAverageMotion += sqrt((pointsTotalMovement.col(i)).squaredNorm());
        }
        dAverageMotion /= nCount;

        if (bVerbose) cout << "Average motion: " << dAverageMotion << endl;
        //cout << "Average noise: " << 2.0 * dNoiseSD << endl; //2 because noise added to both points
        double dNoiseSDUse = dNoiseSD; //* dAverageMotion / 2.0; //Noise is a fraction of average point motion.

        double dMean2X = pointsRight.row(0).sum() / nCount;
        double dMean2Y = pointsRight.row(1).sum() / nCount;
        if (bVerbose) cout << "Means1: " << dMeanX << ", " << dMeanY << endl;
        if (bVerbose) cout << "Means2: " << dMean2X << ", " << dMean2Y << endl;

        double dTotalD = 0, dTotalID = 0;

        //Check these points are in front of both cameras!
        CCamera camP, camPp = R_exact | C3dPoint(translation);
        int i = 0;
        for (; i < nCount; i++) {
            C3dPoint X(points(0, i), points(1, i), points(2, i));
            if (!X.testInFront(camP, camPp)) {
                REPEAT(1, cout << "ERROR: Points not all in front of cameras\n");
                //break;
            }

            double depth = X.depth(camP);
            dTotalID += 1 / depth;
            dTotalD += depth;
        }

        CStats checkMotionMeasurement;

        int numOutliers = 0;
        CMask maskExact(nCount);
        //Introduce outliers
        for (int i = 0; i < nCount; i++) {
            maskExact[i] = 1;
            points0[i] = (CSimple2dPoint(pointsLeft(0, i) + CRandom::Normal(0, dNoiseSDUse), pointsLeft(1, i) + CRandom::Normal(0, dNoiseSDUse)));
            points1[i] = (CSimple2dPoint(pointsRight(0, i) + CRandom::Normal(0, dNoiseSDUse), pointsRight(1, i) + CRandom::Normal(0, dNoiseSDUse)));

            double dTotalMotion = sqrt(pointsTotalMovement.col(i).squaredNorm());
            double dTotalErr = sqrt(sqr(CRandom::Normal(0, dNoiseSDUse) + CRandom::Normal(0, dNoiseSDUse)) + sqr(CRandom::Normal(0, dNoiseSDUse) + CRandom::Normal(0, dNoiseSDUse)));
            checkMotionMeasurement.add(dTotalErr / dTotalMotion);

            pointIds[i] = CPointIds(i, i);
            adPriorProbs[i] = dAssumedInlierProb;

            const bool bMismatchOutliers = false;

            const int MAX_REUSED_POINT_COUNT = NN_MAX;
            if (i > 0 && CRandom::Uniform(1.0) > dInlierRate && numOutliers < MAX_NUM_OUTLIERS && nReusedPointCount < MAX_REUSED_POINT_COUNT) {
                numOutliers++;
                //Add an outlier:
                maskExact[i] = 0;

                if (bMismatchOutliers) //This helps RANSAC a lot
                {
                    nReusedPointCount++;
                    points1[i] = CSimple2dPoint(points1[i - 1].getX(), points1[i - 1].getY());
                    pointIds[i] = CPointIds(i, pointIds[i - 1].id2()); //TODO: Can make prior probs break law of total prob.

                    //This point is now used in nReusedPointCount correspondences
                    for (int j = i - nReusedPointCount + 1; j <= i; j++)
                        adPriorProbs[j] = dAssumedInlierProb / nReusedPointCount;
                } else
                    points1[i] = CSimple2dPoint(points1[i - 1].getX() + CRandom::Normal(0.01, dNoiseSDUse), points1[i - 1].getY() + CRandom::Normal(0, dNoiseSDUse)); //A random point *close* to previous points

            } else {
                nReusedPointCount = 1;
            }

            //cout << pointsLeft.col(i).transpose() * E * pointsRight.col(i) << ',';
            //cout << pointsRight.col(i).transpose() * E * pointsLeft.col(i) << endl;
        }

        if (bVerbose) {
            cout << numOutliers << " outliers\n";
            cout << "Err as prop of motion = " << checkMotionMeasurement.mean() << " sd=" << checkMotionMeasurement.sd() << endl;
        }

        /*if(bVerbose)
    {
        cout << "Testing refinement...\n";
                C3dRotation R_refined = R_exact;
                C3dPoint T_refined = C3dPoint( translation(0),  translation(1),  translation(2) );
        refinePp(0, (const CPointVec2d&)points0, (const CPointVec2d&)points1, R_refined, T_refined, bVerbose);
        //refineE((const CPointVec2d&)points0, (const CPointVec2d&)points1, R_refined, T_refined, true, bVerbose);
    }*/

        C3x3MatModel RANSACmodel;
        CMask mask(points0.size());

        const bool USE_ALL_POINTS = false;
        if (USE_ALL_POINTS) //Hack to make RANSAC accept all points as inliers
            PARAMS.E_INLIER_THRESH_PX = HUGE;

        int nInliers = 0;

        const bool bHough = false;
        if (bHough) {
            nInliers = pHough->findRTHough(points0, points1, E_exact, maskExact, RANSACmodel, mask);
            mask.pp();
        } else {
            //for(int i=0;i<9;i++)
            //getE(points0, points1, adPriorProbs, pointIds, PARAMS, model, mask, bUprightCam ? 0.2 : -1, 1, 1);
            CStopWatch s; s.startTimer();
            nInliers = getE(points0, points1, adPriorProbs, pointIds, PARAMS, RANSACmodel, mask, bUprightCam ? 0.2 : -1, 1, 1);
            s.stopTimer();
            cout << "RANSAC took " << s.getElapsedTime() << " seconds" << endl;
            mask.pp();
        }
        if (bVerbose)
            cout << nInliers << " inliers\n";

        compareMasks(mask, maskExact, bVerbose);

        //bool bSuccess = (nIncorrectInliers / (double) (nInliers)) < 0.1 && nCorrectInliers / (double) (nTotalInliers) > 0.5;
        //if (bSuccess)
        //  nSuccesses++;

        Eigen::Matrix3d E_fromRANSAC(RANSACmodel.asDouble9());
        E_fromRANSAC.transposeInPlace();

        const bool bTestStableAtSolution = false, bApplyRobustNormsToAll = true;

        C3dRotation R_fromRANSAC;
        C3dPoint T_fromRANSAC;
        if (!RTFromE(mask, E_fromRANSAC, points0, points1, pointIds, R_fromRANSAC, T_fromRANSAC, sqr(PARAMS.E_INLIER_THRESH_PX))) {
            cout << "Error recovering R,T after RANSAC" << endl;
            RANSACstats.addFail();
            continue;
        }
        
        //AFTER RANSAC masking, so can choose whether we want the exact mask or not
        if (bTestStableAtSolution)
        {
            E_fromRANSAC = E_exact;
            //maskExact.copyInto(mask);
        }
        //for(int i=0;i<mask.size(); i++)
            //if(!maskExact[i]) mask[i] = false;            
          //  if(maskExact[i]) mask[i] = true;            

        int nIncorrectInliers = compareMasks(mask, maskExact, bVerbose);


        //%%%%%%%%%%%%%%%%
        //if(diff(R_fromRANSAC, R_exact) > 0.25)
        //  continue;
        //if(nIncorrectInliers != 5)
          //continue;
        //%%%%%%%%%%%%%%%%

        printErrors(RANSACstats, R_exact, T_exact, R_fromRANSAC, T_fromRANSAC);

        //First-pass refinement, LS on inliers:
        C3dPoint T_fromLS = T_fromRANSAC;
        C3dRotation R_fromLS = R_fromRANSAC;
        CMask mask_LSInliers(mask);

        const bool bChainRefinementAlgs = false;
        
        bool bLSOnInliersSucceeded = true;

        {
            CRefineEOnRTManifold::refineLSOnMask(points0, points1, R_fromLS, T_fromLS, mask_LSInliers);
            printErrors(LSMaskStats, R_exact, T_exact, R_fromLS, T_fromLS);
            if (!RTFromRT(mask_LSInliers, points0, points1, pointIds, R_fromLS, T_fromLS, sqr(PARAMS.E_INLIER_THRESH_PX))) {
                cout << "Error recovering R,T after LS on mask" << endl;
                if(bChainRefinementAlgs)
                    bLSOnInliersSucceeded = false;
            }
        }
        
        //continue;
        
        if (bLSOnInliersSucceeded) {
            if (bTestStableAtSolution) {
                T_fromLS = T_exact;
                R_fromLS = R_exact; //NOT USED
            }

            C3dPoint T_robust = bChainRefinementAlgs ? T_fromLS : T_fromRANSAC;
            C3dRotation R_robust = bChainRefinementAlgs ? R_fromLS : R_fromRANSAC;
            CMask mask_robustInliers(mask);
            if(bChainRefinementAlgs)
                mask_LSInliers.copyInto(mask_robustInliers);

            if (bApplyRobustNormsToAll)
                CRefineEOnRTManifold::refineRobustOnAll(points0, points1, R_robust, T_robust);
            else
                CRefineEOnRTManifold::refineRobustOnMask(points0, points1, R_robust, T_robust, mask_robustInliers);

            printErrors(robustStats, R_exact, T_exact, R_robust, T_robust);
        }

        //Start again evaluating linear algs
        Eigen::Matrix3d E_fromLinearLS = E_fromRANSAC;
        CMask mask_LLSInliers(mask);

        { //Weighted linear LS
            refineWeightedLinearLSOnMask(points0, points1, E_fromLinearLS, mask_LLSInliers, false);
            C3dPoint T_fromLLS;
            C3dRotation R_fromLLS;
            if (!RTFromE(mask_LLSInliers, E_fromLinearLS, points0, points1, pointIds, R_fromLLS, T_fromLLS, sqr(PARAMS.E_INLIER_THRESH_PX))) {
                cout << "Error recovering R,T after linear weighted LS" << endl;
                LSStatsLinear.addFail();
                if (!bApplyRobustNormsToAll && bChainRefinementAlgs) continue; //Can't use linear input for robust linear if this fails
            } else
                printErrors(LSStatsLinear, R_exact, T_exact, R_fromLLS, T_fromLLS);
        }

        { //Weighted linear robust
            CMask mask_LLSInliers1(mask);
            if(bChainRefinementAlgs)
                mask_LLSInliers.copyInto(mask_LLSInliers1);
            Eigen::Matrix3d E_fromLinearRobust = bChainRefinementAlgs ? E_fromLinearLS : E_fromRANSAC;
            if (bApplyRobustNormsToAll)
                refineWeightedLinearRobustOnAll(points0, points1, E_fromLinearRobust, false);
            else
                refineWeightedLinearRobustOnMask(points0, points1, E_fromLinearRobust, mask_LLSInliers1, false);
            C3dPoint T_fromLinearRobust;
            C3dRotation R_fromLinearRobust;
            if (!RTFromE(mask_LLSInliers1, E_fromLinearRobust, points0, points1, pointIds, R_fromLinearRobust, T_fromLinearRobust, sqr(PARAMS.E_INLIER_THRESH_PX))) {
                cout << "Error recovering R,T after linear weighted LS" << endl;
                robustStatsLinear.addFail();
            } else
                printErrors(robustStatsLinear, R_exact, T_exact, R_fromLinearRobust, T_fromLinearRobust);
        }

        { //Weighted linear robust for F
            CMask mask_LLSInliers2(mask);
            if(bChainRefinementAlgs)
                mask_LLSInliers.copyInto(mask_LLSInliers2);
            
            Eigen::Matrix3d E_fromLinearRobustF = bChainRefinementAlgs ? E_fromLinearLS : E_fromRANSAC;
            if (bApplyRobustNormsToAll)
                refineWeightedLinearRobustOnAll(points0, points1, E_fromLinearRobustF, true);
            else
                refineWeightedLinearRobustOnMask(points0, points1, E_fromLinearRobustF, mask_LLSInliers2, true);

            C3dPoint T_fromLinearRobustF;
            C3dRotation R_fromLinearRobustF;
            if (!RTFromE(mask_LLSInliers2, E_fromLinearRobustF, points0, points1, pointIds, R_fromLinearRobustF, T_fromLinearRobustF, sqr(PARAMS.E_INLIER_THRESH_PX))) {
                cout << "Error recovering R,T after linear weighted LS" << endl;
                robustStatsLinearF.addFail();
            } else
                printErrors(robustStatsLinearF, R_exact, T_exact, R_fromLinearRobustF, T_fromLinearRobustF);
        }

        continue;

        /*double dTotalResidExactE = E_sumSquaredResiduals(E_exact, points0, points1, mask);
        if(bVerbose) cout << dTotalResidExactE << " = total residual in exact E\n";

        double dTotalResid = E_sumSquaredResiduals(Efound, points0, points1, mask);
        if(bVerbose) cout << dTotalResid << " = total residual\n";*/

        if (bVerbose) cout << "Now refining E:\n";

        CPointVec2d aTestPoints1;
        aTestPoints1.reserve(nInliers);
        CPointVec2d aTestPoints2;
        aTestPoints2.reserve(nInliers);
        mask2Vectors<C2dPoint > (mask, points0, points1, aTestPoints1, aTestPoints2);

        CCamera Pp;
        if (!chooseCamFromE(E_fromRANSAC, aTestPoints1, aTestPoints2, Pp)) {
            cout << "Error: Camera selection failed\n";
            continue;
        }

        if (bVerbose) cout << aTestPoints2.size() << " points remaining\n";
        if (aTestPoints2.size() < 8)
            continue;

        const C3dRotation R_8pt = Pp.rotation();
        const C3dPoint T_8pt = Pp.translation();

        double dCalcAngle = R_8pt.angle(); //could move earlier?
        if (bVerbose) cout << dCalcAngle << "=Chosen angle\n";
        double d8ptAngleError = (R_8pt * R_exact.t()).angle();

        //if (d8ptAngleError > 0.1)
        //  continue;

        stats8pt.add(d8ptAngleError);
        stats8ptTrans.add(angle(T_8pt, translation));

        continue;

        C3dRotation R_refined = R_8pt;
        C3dPoint T_refined = T_8pt;

        C3dRotation R_refined_old = R_8pt;
        C3dPoint T_refined_old = T_8pt;
        C3dRotation R_refined_HZ = R_8pt;
        C3dPoint T_refined_HZ = T_8pt;

        if (bVerbose) cout << "Refining camera position\n";
        //refineE(aTestPoints1, aTestPoints2, R_refined, T_refined, true, bVerbose);
        //cout << "Chosen cam before refinement:\n" << (R_refined | T_refined) << endl;

        refineRT_LM_SampsonError(aTestPoints1, aTestPoints2, R_refined, T_refined, true, bVerbose);

        dCalcAngle = R_refined.angle(); //could move earlier?
        if (bVerbose) cout << dCalcAngle << "=Refined angle\n";
        double dRefinePpAngleError = (R_refined * R_exact.t()).angle();
        double dRefinePpTransError = angle(T_refined, translation);
        statsRefineP.add(dRefinePpAngleError);
        statsRefinePTrans.add(dRefinePpTransError);

        if (d8ptAngleError < 0.9 * dRefinePpAngleError) {
            cout << "Error higher after refinement\n";
        } else if (d8ptAngleError > 1.1 * dRefinePpAngleError) {
            cout << "Error lower after refinement\n";
        }

        //continue;
        /*{
                C5ptEssentialMat get5pt(points0, points1);
                TSubSet HS(points0.size());
                CMask mask5(points0.size()); mask5.setZero();
                for(int i=0; i<5; i++)
                {
                        HS[i]=i;
                        mask5[i]=1;
                }
                const T3x3MatModels * m = dynamic_cast<const T3x3MatModels *>( get5pt.getModels(HS) );

                C3x3MatModel bestModel;
                double dErrBest=HUGE;

                double dRbest=100, dTbest=100;

                CPointVec2d aTestPoints15; aTestPoints15.reserve(nInliers);
                CPointVec2d aTestPoints25; aTestPoints25.reserve(nInliers);
                mask2Vectors<C2dPoint>(mask5, points0, points1, aTestPoints15, aTestPoints25);

                Matrix3d Efound5;
                for(int i=0;i<m->numModels();i++)
                {
                        C3x3MatModel & mm= reinterpret_cast<C3x3MatModel &>( const_cast<CModel &>( m->getData(i)));
                Matrix3d Efound5_temp(mm.asDouble9()); Efound5_temp.transposeInPlace();
                / *double dErrInE5 = EssentialMat_SSD(E_exact, Efound5);
                        if(dErrBest > dErrInE5)
                        {
                                dErrBest = dErrInE5;
                                bestModel = mm;
                                Efound5 = Efound5_temp;
                        }* /
                        CCamera Pp5;
                        if(!chooseCamFromE(Efound5_temp, aTestPoints15, aTestPoints25, Pp5) )
                        {
                                cout << "Error: 5pt Camera selection failed\n";
                                continue;
                        }

                        if(bVerbose) cout << aTestPoints25.size() << " 5points remaining\n";
                        if(aTestPoints25.size() < 5)
                                continue;

                        const C3dRotation R_5pt_temp =Pp5.rotation();
                        const C3dPoint T_5pt_temp = Pp5.translation();

                        double dRtemp = (R_exact * R_5pt_temp.t()).angle();
                        double dTtemp = angle(T_5pt_temp, translation);
                        if(dRbest > dRtemp)
                        {
                                dRbest=dRtemp;
                                dTbest=dTtemp;
                        }
                }
                delete m;

                if(dRbest < 100)
                {
                        statsRefine5pt.add(dRbest);
                        statsRefine5ptTrans.add(dTbest);
                }
        }*/

        /*testOut << dRefinePpAngleError << '\t';
        testOut << dRefinePpTransError << "\t\n";

        {
                refineE(aTestPoints1, aTestPoints2, R_refined_old, T_refined_old, bVerbose);
                statsRefine5pt.add((R_refined_old * R_exact.t()).angle());
                statsRefine5ptTrans.add(angle(T_refined_old, translation));
                T3dPoints points3d;
                CCamera P_old, Pp_old = R_refined_old | T_refined_old;
                for(int i=0;i<nCount; i++)
                        points3d.push_back(reconstruct(P_old, Pp_old, aTestPoints1[i], aTestPoints2[i]));
				
                cout << threeDPointSSD(points3d, points) << " = 3d point SSD from refine Sampson err" << endl;
        }*/
        {
            T3dPoints points3d;
            refine3d(aTestPoints1, aTestPoints2, points3d, R_refined_HZ, T_refined_HZ, bVerbose);
            statsRefineHZ.add((R_refined_HZ * R_exact.t()).angle());
            statsRefineHZTrans.add(angle(T_refined_HZ, translation));
            cout << angle(T_refined_HZ, translation) << "=HZ error" << endl;
            cout << threeDPointSSD(points3d, points) << " = 3d point SSD from HZ" << endl;
        }

        CCamera Pp_refined = R_refined | T_refined;

        //cout << "Refined cam back:\n" << Pp << endl;

        makeE(R_refined, T_refined, E_fromRANSAC);
        //double dTotalResidNew = E_sumSquaredResiduals(Efound, points0, points1, mask);
        //if(bVerbose) cout << dTotalResidNew << " = total residual after refinement\n";
        //CHECK(dTotalResidNew > dTotalResid, "refinement has increased residual\n");

        if (bVerbose) cout << "E from refinement SE " << E_fromRANSAC << endl << endl;
        //		cout << Pp << "= refined camera" << endl;
        //cout << "Conversion to/from E is failing...??\n";

        /*getCamsFromE(Efound, aPp);
        // chooseCamMat also removes points behind camera
        if(bVerbose) cout << aTestPoints2.size() << " points...";
        nCamMat = chooseCamMat(aPp, aTestPoints1, aTestPoints2);
        if(bVerbose) cout << aTestPoints2.size() << " points remaining\n";

        //cout << "Chosen cam again:\n" << aPp[nCamMat] << endl;

        CHECK (nCamMat == CAM_SELECT_ERROR, "Error selecting cam")

        //cout << "Cam mat: " << aPp[nCamMat] << endl;
        C3dRotation R_fromcamAgain = aPp[nCamMat].rotation();
        C3dPoint T_fromcamAgain = aPp[nCamMat].translation();
        if(bVerbose) cout << "R and T back again from P\n" << R_fromcamAgain << endl << T_fromcamAgain << endl;*/

        if ((R_refined * R_exact.t()).angle() > 5 * (dNoiseSD + 0.01)) {
            cout << "Error recovering rotation\n";
        }

        if ((T_refined - translation).sum_square() > 5 * (dNoiseSD + 0.0001)) {
            cout << "Error recovering translation\n";
        }

        /*double dErrInE2 = EssentialMat_SSD(E, Efound);
        double dChangeInRotation = (R_refined * R_8pt.t()).angle();
        if(dErrInE1+1e-5 < 0.8*dErrInE2)
        {
                if(bVerbose) cout << "Possible error refining E--solution is worse. probably just due to the inlier set found\n" << "Rotation change: " << dChangeInRotation << endl;
                CHECK(dChangeInRotation > 1e-2, "Significant change in rotation despite error being worse")
        }

CHECK(dErrInE1+2e-1 < dErrInE2, "Error refining E--solution is significantly worse")*/
    }
    s.stopTimer();
    dTotalTime = s.getElapsedTime();

    cout << dInlierRate << ',' << dTotalTime / (nIters) << endl;

    //cout << stats8pt.mean() << ',' << statsRefineP.mean() << ',' << statsRefineHZ.mean() << endl;

    //statsRefine5pt.writeTSVdata(results);
    stats8pt.writeTSVdata(results, true);
    statsRefineP.writeTSVdata(results);
    statsRefineHZ.writeTSVdata(results);

    //statsRefine5ptTrans.writeTSVdata(resultsTrans);
    stats8ptTrans.writeTSVdata(resultsTrans);
    statsRefinePTrans.writeTSVdata(resultsTrans);
    statsRefineHZTrans.writeTSVdata(resultsTrans);

    //testOut.close();
    /*for(int i=0; i<nCount; i++)
        {
                cout << pointsLeft.col(i).transpose() * Efound * pointsRight.col(i) << ',';
                cout << pointsRight.col(i).transpose() * Efound * pointsLeft.col(i) << endl;
        }
        char c;
        cin >> c;	*/

    cout << "Hypotheses per iteration: " << MODELS / (double) nIters << endl;

    ofstream hypPerIt("hypPerIt.tsv", ios_base::app);
    //hypPerIt << PARAMS.HypothesiseAlg << '\t';
    hypPerIt << dInlierRate << "\t" << MODELS / (double) nIters << "\t" << nSuccesses / (double) nIters << endl;
    hypPerIt.close();

    return 0;
}

/*Matrix3d makeH(const Matrix3d & rotMat, const Vector3d & planeNormal, const Vector3d & translation)
{
        return rotMat - translation * planeNormal.transpose();
}*/
bool sizeEqual(const C3dPoint & v1, const C3dPoint & v2) {
    return zero((v1 - v2).sum_square()) || zero((v1 + v2).sum_square());
}

void normaliseS(Eigen::Matrix3d & H) {
    Eigen::JacobiSVD< Eigen::Matrix3d > svdH(H);
    const double s = svdH.singularValues()(1);
    H /= s;
}

void checkS(const Eigen::Matrix3d & H) {
    Eigen::JacobiSVD< Eigen::Matrix3d > svdH(H);
    const double s = svdH.singularValues()(1);
    const Eigen::Matrix3d Hnorm = H / s;

    cout << Hnorm << endl;
    cout << "s = " << svdH.singularValues()(1) << endl;
}

void simulateRD(C2dPoint & p) {
    //Assume centred at 0
    const double K = -0.01;
    const double shift = 1 + K * p.sum_square();
    p *= shift;
}

int mainFindH(const int nCount) {
    CRANSACHomographyParams PARAMS(0, 0);

    T2dPoints points0(nCount), points1(nCount);
    cout << points0.size() << endl;
    CPointIdentifiers pointIds(nCount);
    CInlierProbs adPriorProbs(nCount);
    //	CStopWatch s;

    for (int iter = 0; iter < 1; iter++) {
        double dAngle = CRandom::Uniform(-0.1, 0.1);
        C3dPoint axis(CRandom::Uniform(-0.5, 0.5), 2, CRandom::Uniform(-0.5, 0.5));
        C3dPoint translation(3, 2, 1); //translation.normalise();
        Eigen::Vector3d translationVec;
        translationVec << translation.getX(), translation.getY(), translation.getZ();
        C3dRotation q(axis, dAngle);
        CHECK(!zero(q.angle() - fabs(dAngle)), "Quaternion axis-angle failed")

        Matrix3d rotMat;
        q.asMat(rotMat);

        /*	    Matrix3d transMat;
                    Xmat(translation, transMat);
                    Matrix3d E = rotMat * transMat;*/

        C3dPoint planeNormal(0, 0, 1);
        //cout << planeNormal * translation.transpose() << endl << "=normal x t";
        Matrix3d H;
        makeH(rotMat, planeNormal, translation, H); // see http://en.wikipedia.org/wiki/Homography#Computer_vision_applications
        checkS(H);

        C3dRotation rotation;
        C3dPoint planeNormal_est;
        C3dPoint camMotion;
        decomposeHomography(H, rotation, planeNormal_est, camMotion);
        //cout << rotMat << "=rotMat original" << endl << endl;

        if(IS_DEBUG) CHECK(!zero(diff(rotation, q)), "rotation recovery from original H failed");
        if(IS_DEBUG) CHECK(!sizeEqual(planeNormal_est, planeNormal), "planeNormal_est recovery from original H failed");
        if(IS_DEBUG) CHECK(!sizeEqual(camMotion, translation), "camMotion recovery from original H failed");

        if (false) //all seems to work fine
        {
            cout << H << "=H original" << endl << endl;

            rotation.asMat(rotMat);
            cout << rotMat << "=rotMat original" << endl << endl;
            Matrix3d H_recovered;
            makeH(rotMat, planeNormal_est, camMotion, H_recovered);
            cout << H_recovered << "=H from fitted decomp." << endl << endl;
            if(IS_DEBUG) CHECK(!zero((H_recovered - H).squaredNorm()), "Decomposing and recovering H failed");
        }

        CHECK(!zero((rotation.t() * q).angle()), "Rotation estimation failed");
        CHECK(!zero(fabs(dotproduct(planeNormal_est, planeNormal)) - 1), "Normal estimation failed");
        CHECK(!zero((camMotion - translation).sum_square()) && !zero((camMotion + translation).sum_square()), "Translation estimation failed");

        Matrix<double, 3, 4 > cam2;
        cam2.block < 3, 3 > (0, 0) = rotMat;
        cam2.block < 3, 1 > (0, 3) = -translationVec;
        Matrix<double, 3, 4 > cam1;
        cam1.block < 3, 3 > (0, 0).setIdentity();
        cam1.block < 3, 1 > (0, 3) = 0 * translationVec;

        //cout << cam1 << endl << endl << cam2 << endl << endl << H << endl << endl << rotMat <<endl << endl << transMat << endl << endl << " = 1, 2, E\n";
        //cout << cam2.transpose()*E*cam1 << " should be skew-sym\n";

        Matrix<double, 4, Dynamic> points = Matrix<double, 4, Dynamic>::Random(4, nCount);
        const double dDepth = 3; // CRandom::Uniform(1,20);
        points.row(2).setConstant(dDepth); //Force planar. Todo: more sophisticated
        points.row(3).setOnes();

        //points.Random();
        Matrix<double, 3, Dynamic> pointsLeft = cam1*points;
        Matrix<double, 3, Dynamic> pointsRight = cam2*points;

        Matrix<double, 3, Dynamic> pointsRightNoise = Matrix<double, 3, Dynamic>::Random(3, nCount);
        pointsRight += pointsRightNoise * 0.01;

        for (int i = 0; i < nCount; i++) {
            pointsLeft.col(i) /= pointsLeft(2, i);
            pointsRight.col(i) /= pointsRight(2, i);

            points0[i] = (CSimple2dPoint(pointsLeft(0, i), pointsLeft(1, i)));
            points1[i] = (CSimple2dPoint(pointsRight(0, i), pointsRight(1, i)));

            //simulateRD((C2dPoint &)points0[i]);
            //simulateRD((C2dPoint &)points1[i]);

            pointIds[i] = CPointIds(i, i);
            adPriorProbs[i] = 0.3;

            if (i % 2 == 1) {
                //Outlier:
                points1[i] = points1[i - 1];
                pointIds[i] = CPointIds(i, i - 1);
            }

            //cout << pointsLeft.col(i).transpose() * E * pointsRight.col(i) << ',';
            //cout << pointsRight.col(i).transpose() * E * pointsLeft.col(i) << endl;
        }

        C3x3MatModel model;
        CMask inliers(nCount);
        PARAMS.VERBOSE = true;
        int nInliers = getH(points0, points1, adPriorProbs, pointIds,
                PARAMS, model, inliers, 1, 1);
        if (nInliers < 500)
            cout << nInliers << " inliers" << endl;

        Matrix3d Hfound(model.asDouble9());
        Hfound.transposeInPlace();
        //normaliseS(Hfound);

        cout << "H:\n" << H << endl << endl;
        cout << "\nH_ransac:\n" << Hfound << endl << endl;
        //Should decompose correctly with no scale ambiguity:
        decomposeHomography(Hfound, rotation, planeNormal_est, camMotion);
        makeH(rotation, planeNormal_est, camMotion*dDepth, Hfound);
        cout << "\nH_ransac normalised:\n" << Hfound << endl << endl;

        double dDepthFromH = translation.length() / camMotion.length();
        cout << dDepthFromH << "=depth from H\n" << endl;

        CHECK((H - Hfound).squaredNorm() > 0.01, "Error recovering Homography")
    }

    /*for(int i=0; i<nCount; i++)
        {
                cout << pointsLeft.col(i).transpose() * Efound * pointsRight.col(i) << ',';
                cout << pointsRight.col(i).transpose() * Efound * pointsLeft.col(i) << endl;
        }*/

    return 0;
}

template<int N, typename T>
void testInversion_int() {
    Matrix<T, N, N> I;
    I.setIdentity();

    CStopWatch s;
    s.startTimer();
    for (int i = 0; i < 1000000; i++) {
        Matrix<T, N, N> M;
        do {
            M.setRandom();
        } while (fabs(M.determinant()) < 0.01);

        const Matrix<T, N, N> & M_inv = M.inverse();

        Matrix<T, N, N> I_test = M_inv * M;
        I_test -= I;
        if (I_test.squaredNorm() > 1e-1) {
            cout << I_test << endl;
            CHECK(1, "Inversion failed");
        }
    }
    s.stopTimer();
    cout << N << "=N, Time=" << s.getElapsedTime() << endl;
}

template<int N>
void testInversion() {
    testInversion_int<N, double>();
    testInversion_int<N, float>();
    testInversion_int<N, complex<double> >();
    testInversion_int<N, complex<float> >();
}
void findK();
void testELM();

void test5ptNumerical(const double dErrRads, const int nIters, ofstream & results) {
    //CModelHypothesiser * getHypothesiser(const CRANSACParams::eHypothesiseAlg alg, const T2dPoints & points0, const T2dPoints & points1, const double dUprightThresh)

    const int NUM_POINTS = 7;
    T2dPoints points0(NUM_POINTS), points1(NUM_POINTS);
    TSubSet anHypSet(NUM_POINTS);

    for (int nPoint = 0; nPoint < NUM_POINTS; nPoint++)
        anHypSet[nPoint] = nPoint;

    CRANSACParams::eHypothesiseAlg alg = CRANSACParams::e5PtRt;
    //CRANSACParams::e7PtEFast;
    //CRANSACParams::e5Pt_GradientDesc;
    //CRANSACParams::e5PtE;
    const double LENGTH = sqrt(2.0);

    boost::scoped_ptr<CModelHypothesiser> pHypothesiseAlg(getHypothesiser(alg, points0, points1, -1));
    CStats stats, statsAngle;

    double dAvModels = 0, nRejected = 0, closestThresh = 0.5, //1e-5, 
            targetPropRejected = 0.2;

    for (int nIter = 0; nIter < nIters; nIter++) {
        double dAngle = 0.3;
        C3dRotation q;
        q.setRandom(dAngle);
        C3dPoint t;
        t.setRandomNormal();
        t *= 0.25;

        Matrix3d E_exact, E_closest;
        makeE(q, t, E_exact);
        E_exact *= LENGTH / E_exact.norm();

        //Make 2 cameras
        CCamera P, Pp = q | t;

        C3dPoint aPoints[NUM_POINTS];
        for (int nPoint = 0; nPoint < NUM_POINTS; nPoint++) {
            aPoints[nPoint].setRandomPlanar(1);
            aPoints[nPoint].addNoise(0.0001);
            points0[nPoint] = aPoints[nPoint].photo(P);
            points1[nPoint] = aPoints[nPoint].photo(Pp);

            points0[nPoint] = CSimple2dPoint(points0[nPoint].getX() + CRandom::Normal(0, dErrRads), points0[nPoint].getY() + CRandom::Normal(0, dErrRads));
        }

        const T3x3MatModels * pModels = (const T3x3MatModels *) pHypothesiseAlg->getModels(anHypSet);
        double dClosest = HUGE;

        //cout << endl << E_exact << endl; 
        dAvModels += pModels->numModels();
        for (int i = 0; i < pModels->numModels(); i++) {
            Matrix3d E = Matrix3d((dynamic_cast<const C3x3MatModel &> (pModels->getData(i))).asDouble9()).transpose();

            double dDist = (E - E_exact).squaredNorm();
            double dDist2 = (E + E_exact).squaredNorm();
            if (dDist < dClosest || dDist2 < dClosest) {
                dClosest = min<double>(dDist, dDist2);
                E_closest = E;
            }

            //cout << endl << i << E << endl;
        }

        //cout << dClosest << endl;
        //cout << dAvModels/nIter << endl;
        double dClosestAngle = HUGE, dClosestTrans = HUGE;

        if (dClosest < closestThresh) {
            stats.add(sqrt(dClosest));

            const int N_CAM_MATS = 4;
            CCamera aPp[N_CAM_MATS];
            getCamsFromE(E_closest, aPp);

            for (int nCam = 0; nCam < N_CAM_MATS; nCam++) {
                const C3dRotation R_8pt = aPp[nCam].rotation();
                const C3dPoint T_8pt = aPp[nCam].translation();
                double dAngle = diff(R_8pt, q);
                if (dAngle < dClosestAngle)
                    dClosestAngle = dAngle;
            }

            statsAngle.add(dClosestAngle);
        }
    }
    //stats.printData();
    ofstream numErrors("allNumericalErrors");
    stats.pp("\n", numErrors);
    numErrors.close();

    //stats.cropTop(0.2);
    REPEAT(1, stats.writeTSVheader("", cout); cout << endl);
    statsAngle.writeTSVdata(cout, true);
    cout << endl;
    statsAngle.writeTSVdata(results, true);
    results << endl;
    //stats.writeTSVdata(cout, true); cout << endl;
    //stats.writeTSVdata(results, true); results << endl;


    //cout << "Av models" << dAvModels/nIters << endl;
    //    cout << "Prop rejected" << nRejected/nIters << endl;
}

void test5ptRoots() {

    const int NUM_POINTS = 5;
    T2dPoints points0(NUM_POINTS), points1(NUM_POINTS);
    TSubSet anHypSet(NUM_POINTS);

    for (int nPoint = 0; nPoint < NUM_POINTS; nPoint++)
        anHypSet[nPoint] = nPoint;

    CRANSACParams::eHypothesiseAlg alg = //CRANSACParams::e5PtRt;
    //CRANSACParams::e7PtEFast;
    //CRANSACParams::e5Pt_GradientDesc;
    CRANSACParams::e5PtE;
    const double LENGTH = sqrt(2.0);

    boost::scoped_ptr<CModelHypothesiser> pHypothesiseAlg(getHypothesiser(alg, points0, points1, -1));

    double dAngle = 0.3;
    C3dRotation q;
    q.setRandom(dAngle);
    C3dPoint t;
    t.setRandomNormal();
    t *= 0.25;

    Matrix3d E_exact, E_closest;
    makeE(q, t, E_exact);
    E_exact *= LENGTH / E_exact.norm();

    //Make 2 cameras
    CCamera P, Pp = q | t;

    C3dPoint aPoints[NUM_POINTS];
    for (int nPoint = 0; nPoint < NUM_POINTS; nPoint++) {
        aPoints[nPoint].setRandom(1);
        aPoints[nPoint].addNoise(0.0001);
        points0[nPoint] = aPoints[nPoint].photo(P);
        points1[nPoint] = aPoints[nPoint].photo(Pp);
    }

    for (double dx = -1; dx < 1; dx+=0.02) {
        for (double dy = -1; dy < 1; dy+=0.02) {
            points0[0] = CSimple2dPoint(dx, dy);
            
            const T3x3MatModels * pModels = (const T3x3MatModels *) pHypothesiseAlg->getModels(anHypSet);
            double dClosest = HUGE;

            //cout << endl << E_exact << endl; 
            cout << dx << '\t' << dy << '\t' << pModels->numModels() << endl;
            for (int i = 0; i < pModels->numModels(); i++) {
                Matrix3d E = Matrix3d((dynamic_cast<const C3x3MatModel &> (pModels->getData(i))).asDouble9()).transpose();

            }
        }
    }
}

template<typename TFloat>
void testMatMulSpeed() {
    const int ROWS = 3, COLS = 3;
    typedef Eigen::Matrix<TFloat, ROWS, 1 > TTestVec;
    typedef Eigen::Matrix<TFloat, ROWS, COLS> TMat;

    TTestVec v;
    v.setOnes();
    TMat aM[1000];
    for (int i = 0; i < 1000; i++) {
        aM[i].setRandom();
        aM[i] *= 0.9;
    }

    CStopWatch s;
    s.startTimer();
    for (int i = 0; i < 100000000; i++) {
        v = aM[i % 1000] * v;
    }
    s.stopTimer();
    cout << v.transpose() << " time " << s.getElapsedTime() << endl;

}

int main() {

    //test5ptRoots();
    //return 0;

    //testMatMulSpeed<float>();
    //testMatMulSpeed<double>();
    //return 0;

    ofstream resultsNumerical("numericalResults.tsv");

    //for(double dErrRads=0; dErrRads <= 0.025; dErrRads += 0.0005)
    //    test5ptNumerical(dErrRads, 20000, resultsNumerical);
    //return 0;

    //testELM(); return 0;

    //Can normalise E either by setting SVs to 1 or by setting SSDs equal. Anyway, should be minimising the *reprojection error*
    ofstream results("camPosRefinementResults.tsv");
    ofstream resultsTrans("camPosRefinementResultsTrans.tsv");


    results << "Points\tNoise\t";
    //CStats::writeTSVheader("5pt", results);
    CStats::writeTSVheader("Refine P", results);
    CStats::writeTSVheader("Refine SE", results);
    CStats::writeTSVheader("Refine HZ", results);
    results << endl;

    resultsTrans << "Points\tNoise\t";
    //CStats::writeTSVheader("5pt", resultsTrans);
    CStats::writeTSVheader("8pt", resultsTrans);
    CStats::writeTSVheader("Refine SE", resultsTrans);
    CStats::writeTSVheader("Refine HZ", resultsTrans);
    resultsTrans << endl;

    for (double IR = 0.25; IR >= 0.0; IR -= 0.05)
        for (int nPoints = 500;; nPoints++)
            for (double dNoise = 0.0025; dNoise <= 0.0025; dNoise += 0.02) {
                results << nPoints << '\t' << dNoise << '\t';
                resultsTrans << nPoints << '\t' << dNoise << '\t';
                mainFindE(false, nPoints, 50, dNoise, IR, IS_DEBUG, results, resultsTrans);
                results << endl;
                resultsTrans << endl;

                return 0;
            }

    results.close();
    resultsTrans.close();
    return 0;

    //Test faster inversion

    /*testInversion<5>();
    testInversion<6>();
    testInversion<7>();
    testInversion<8>();
    testInversion<9>();
    testInversion<10>();
    testInversion<11>();
    testInversion<12>();
    testInversion<13>();
    testInversion<14>();
    testInversion<15>();
    testInversion<16>();*/
    return 0;

    /*double dAngle = 0.72;//CRandom::Uniform(-1.0,1.0);
    //    Vector3d axis; axis << 1, 2, 6;
C3dPoint axis(1, 2, 6);
Vector3d translation; translation << 5.5,-.5, 4;

//AngleAxis<double> q(dAngle, axis);
//Matrix3d rotMat = q.toRotationMatrix();

    C3dRotation q(axis, dAngle);
    Matrix<double, 3,3,2,3,3> rotMat;
    q.asMat(rotMat);

    cout << rotMat << endl;

    Vector3d planeNormal; planeNormal << 1, 2, 3;
    planeNormal /= sqrt(planeNormal.squaredNorm());
    //cout << planeNormal * translation.transpose() << endl << "=normal x t";
    Matrix3d H = makeH(rotMat, planeNormal, translation);
    cout << H << endl;

    C3dRotation rotation; C3dPoint planeNormal_est; C3dPoint camMotion;
    decomposeHomography( H, rotation, planeNormal_est, camMotion );

    cout << H << " before" << endl;
    Matrix<double, 3,3,2,3,3> RotAfter;
    rotation.asMat(RotAfter);

    planeNormal << planeNormal_est.getX(), planeNormal_est.getY(), planeNormal_est.getX();
    translation << camMotion.getX(), camMotion.getY(), camMotion.getX();

    cout << RotAfter << " RotAfter" << endl;
    Matrix3d Hafter = makeH(RotAfter, planeNormal, translation);
    cout << H << " after" << endl;
    return 0; */
}
