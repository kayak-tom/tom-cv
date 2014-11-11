/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * ransacWithParams.cpp
 *
 *  Created on: 19/10/2009
 *      Author: tom
 */
#include "ransac.h"
#include "ransacParams.h"

#include <iostream>
#include "findBestMatchingDisjoint.h"
#include "essentialMat_Eigen.h"
#include "openCV8PtEssentialMatModified.h"
#include "openCV4PtHomography.h"
#include <boost/smart_ptr.hpp>
#include "topdownRLRefiner.h"
#include "openCV7PtFundamentalMat.h"
#include "essentialMatLM.h"


#include "geom/geom.h"
#include "geom/geom_eigen.h"


using namespace boost;

CModelHypothesiser * getHypothesiser(const CRANSACParams::eHypothesiseAlg alg, const T2dPoints & points0, const T2dPoints & points1, const double dUprightThresh)
{
    switch (alg) {
        case CRANSACParams::eMCE:
            return new CMCEssentialMat(points0, points1, 10);
        case CRANSACParams::e2PtE:
            return new C2ptMCEssentialMat(points0, points1, 100);
        case CRANSACParams::e5PtE:
            return new C5ptEssentialMat(points0, points1, dUprightThresh);
        case CRANSACParams::e7PtF:
        case CRANSACParams::e7PtE:
            return new COpenCV7PtFundamentalMat(points0, points1);
        case CRANSACParams::e7PtEFast:
            return new C7ptEssentialMat_GS(points0, points1);
        case CRANSACParams::e5Pt_GradientDesc:
            return new CEssentialMatLMbasis(points0, points1, 5);
        case CRANSACParams::e4Pt_GradientDesc:
            return new CEssentialMatLM(points0, points1, 4, 10);
        case CRANSACParams::e3Pt_GradientDesc:
            return new CEssentialMatLM(points0, points1, 3, 10);
        case CRANSACParams::e5PtRt:
            return new CEssentialMatLM(points0, points1, 5, 10);
        default:
            THROW("HypothesiseAlg not implemented yet");
    }
}

int getE(const T2dPoints & points0, const T2dPoints & points1, CInlierProbs & adPriorProbs, const CPointIdentifiers & pointIds,
        const CRANSACParams & PARAMS,
        C3x3MatModel & E, CMask & inliers, const double dUprightThresh, const double dFocalLength, const int nThreads) {
    const double E_INLIER_THRESH = PARAMS.E_INLIER_THRESH_PX / dFocalLength;

    const int nCount = points0.size();

    scoped_ptr<CModelHypothesiser> pHypothesiseAlg ( getHypothesiser(PARAMS.HypothesiseAlg, points0, points1, dUprightThresh) );
    bool bFindE = (PARAMS.HypothesiseAlg != CRANSACParams::e7PtF);

    const int nSampleSize = pHypothesiseAlg->numPoints();

    if (nCount < nSampleSize)
        return 0;

    scoped_ptr<CSampler> pSampler;

    switch (PARAMS.RANSACSampler) {
        case CRANSACParams::eBaySAC:
            pSampler.reset(new CBaySACSampler_discretePP(nSampleSize, pointIds, adPriorProbs, false /* NOT disjoint */));
            break;
        case CRANSACParams::eSimSAC:
            pSampler.reset(new CSimSACSampler(nSampleSize, pointIds, adPriorProbs, true, true /*not the old approx one*/, 25 /* T. Todo: make a param */));
            break;
        case CRANSACParams::eRANSAC:
            pSampler.reset(new CDisjointRANSACSampler(nSampleSize, pointIds));
            break;
        default:
            THROW("RANSAC sampler parameter val not handled")
    }

    const double dMinPropInliers = 0.1;

    scoped_ptr<CRansacTerminator> pTerminator;
    switch (PARAMS.RANSACTerminator) {
        case CRANSACParams::eSimpleTerminator:
            pTerminator.reset(new CSimpleTerminator(nCount));
            break;
        case CRANSACParams::eWaldSAC:
            //pTerminator.reset(new CWaldSACTerminator(nCount,nSampleSize,4,1000,0.1));
            pTerminator.reset(new CLogWaldSACTerminator(nCount, nSampleSize, pHypothesiseAlg->modelsPerIteration(), pHypothesiseAlg->timePerIteration(), 0.1));
            break;
        case CRANSACParams::eBrownianBridge://TODO: Move to homography too.
            pTerminator.reset(new CBBTerminator(nCount, dMinPropInliers, 0.01));
            break;
        case CRANSACParams::eBrownianBridgeLinear:
            pTerminator.reset(new CBBLinearTerminator(nCount, dMinPropInliers, 0.01));
            break;
        case CRANSACParams::eNoTestTerminator:
            pTerminator.reset(new CRansacTerminator());
            break;
        default:
            THROW("Terminator not recognised");
    }

    scoped_ptr<CRansacTerminator> pTerminatorLogger;
    if (PARAMS.TERMINATOR_LOG)
        pTerminatorLogger.reset(new CTerminatorLogger(pTerminator.get(), nCount));

    scoped_ptr<CIterTerminator> pIterTerminator;
    switch (PARAMS.RANSACIterTerminator) {
        case CRANSACParams::eMaxIters:
            pIterTerminator.reset(new CIterTerminator(PARAMS.MAX_ITERS));
            break;
        case CRANSACParams::eClassicRANSAC:
            pIterTerminator.reset(new CPPIterTerminator(PARAMS.MAX_ITERS, adPriorProbs, nSampleSize, PARAMS.PROB_SUCCESS));
            break;
        case CRANSACParams::eTerminateOnPropInliers:
            pIterTerminator.reset(new CPropIterTerminator(PARAMS.MAX_ITERS, nCount /2));
            break;
    }

    scoped_ptr<CModelRefiner> pSimpleRefiner(new COpenCV8PtEssentialMatModified(points0, points1, PARAMS.E_8PT_CUTOFF, !bFindE));
    scoped_ptr<CInlierCounter> pInlierCounter(new CEssentialMatInlierCounter(E_INLIER_THRESH, points0, points1));
    scoped_ptr<CFindBestMatching> pRefineMatches(new CFindBestMatchingDisjoint(pointIds));
    scoped_ptr<CModelRefiner> pRefiner;
    if (PARAMS.TOPDOWN_ITERS == 0)
        pRefiner.reset(new CNoModelRefiner(nSampleSize));
    else
        pRefiner.reset(new CTopdownRLRefiner(pSimpleRefiner.get(), pRefineMatches.get(), pInlierCounter.get(), PARAMS.TOPDOWN_ITERS, PARAMS.TOPDOWN_EXPAND, PARAMS.TOPDOWN_SCALEDOWN));

    return doRansac(pSampler.get(), pHypothesiseAlg.get(), pRefiner.get(), pTerminatorLogger ? pTerminatorLogger.get() : pTerminator.get(), pIterTerminator.get(), pInlierCounter.get(), pRefineMatches.get(), E, inliers, nThreads, PARAMS.VERBOSE);
}

int getF(const T2dPoints & points0, const T2dPoints & points1, CInlierProbs & adPriorProbs, const CPointIdentifiers & pointIds,
        const CRANSACParams & PARAMS,
        C3x3MatModel & F, CMask & inliers, const int nThreads) {
    if (PARAMS.HypothesiseAlg != CRANSACParams::e7PtF) {
        cout << PARAMS.HypothesiseAlg;
        THROW("To find F hypothesis generator must be set to e7PtF (RANSAC.HypothesiseAlg=7PtF)");
    }

    return getE(points0, points1, adPriorProbs, pointIds,
            PARAMS, F, inliers, -1, 1, nThreads);
}

/*int getF(const T2dPoints & points0, const T2dPoints & points1, CInlierProbs & adPriorProbs, const CPointIdentifiers & pointIds,
                const CRANSACParams & PARAMS,
                CModel & F, CMask & inliers, const int nThreads)
{
        const int nCount = points0.size();

        //Could declare on stack but maybe will parameterise more:
        COpenCV7PtFundamentalMat hypothesise(points0, points1);

        if(nCount < hypothesise.numPoints())
                return 0;

        const int nSampleSize = hypothesise.numPoints();
    //scoped_ptr<CSampler> pSampler(new CRANSACSampler(nSampleSize, points0, points1));
        scoped_ptr<CSampler> pSampler;

        switch(PARAMS.RANSACSampler)
        {
        case CRANSACParams::eBaySAC:
                pSampler.reset(new CBaySACSampler_discretePP(nSampleSize, pointIds, adPriorProbs, true));
                break;
        case CRANSACParams::eSimSAC:
                pSampler.reset(new CSimSACSampler(nSampleSize, pointIds, adPriorProbs, true, true / *not the old approx one* /, 25 / * T. Todo: make a param * /));
                break;
        case CRANSACParams::eRANSAC:
                pSampler.reset(new CDisjointRANSACSampler(nSampleSize, pointIds));
                break;
        default:
                THROW("RANSAC sampler parameter val not handled")
        }

        scoped_ptr<CRansacTerminator> pTerminator(new CSimpleTerminator(nCount));
        scoped_ptr<CIterTerminator> pIterTerminator(new CPPIterTerminator(PARAMS.MAX_ITERS, adPriorProbs, nSampleSize, PARAMS.PROB_SUCCESS));
        scoped_ptr<CModelRefiner> pSimpleRefiner(new COpenCV8PtEssentialMatModified(points0, points1, -1, true));
        scoped_ptr<CInlierCounter> pInlierCounter(new CEssentialMatInlierCounter(PARAMS.E_INLIER_THRESH_PX, points0, points1));
        scoped_ptr<CFindBestMatching> pRefineMatches(new CFindBestMatchingDisjoint(pointIds));
        scoped_ptr<CModelRefiner> pRefiner(new CTopdownRLRefiner(pSimpleRefiner.get(), pRefineMatches.get(), pInlierCounter.get(), PARAMS.TOPDOWN_ITERS, PARAMS.TOPDOWN_EXPAND, PARAMS.TOPDOWN_SCALEDOWN));

        return doRansac(pSampler.get(), &hypothesise, pRefiner.get(), pTerminator.get(), pIterTerminator.get(), pInlierCounter.get(), pRefineMatches.get(), F, inliers, nThreads, PARAMS.VERBOSE);
}*/

int getH(const T2dPoints & points0, const T2dPoints & points1, CInlierProbs & adPriorProbs, const CPointIdentifiers & pointIds,
        const CRANSACHomographyParams & PARAMS,
        C3x3MatModel & H, CMask & inliers, const double dFocalLength, const int nThreads) {
    const double H_INLIER_THRESH = PARAMS.H_INLIER_THRESH_PX / dFocalLength;
    const int nCount = points0.size();

    COpenCV4ptHomography hypothesise(points0, points1);

    if (nCount < hypothesise.numPoints())
        return 0;

    const int nSampleSize = hypothesise.numPoints();
    scoped_ptr<CSampler> pSampler;
    switch (PARAMS.RANSACSampler) {
        case CRANSACHomographyParams::eBaySAC:
            pSampler.reset(new CBaySACSampler_discretePP(nSampleSize, pointIds, adPriorProbs, false));

            /*if(padMMLikelihoods)
            {
                    dynamic_cast<CBaySACSampler_discretePP *>(pSampler.get())->updateProbs(*padMMLikelihoods);
            }*/

            break;
        case CRANSACHomographyParams::eSimSAC:
            pSampler.reset(new CSimSACSampler(nSampleSize, pointIds, adPriorProbs, true, true /*not the old approx one*/, 25 /* T. Todo: make a param */));
            break;
        case CRANSACHomographyParams::eRANSAC:
            pSampler.reset(new CDisjointRANSACSampler(nSampleSize, pointIds));
            break;
        default:
            THROW("RANSAC sampler parameter val not handled")
    }//	(new CBaySACSampler_discretePP(nSampleSize, points0, points1, adPriorProbs, true));

    scoped_ptr<CRansacTerminator> pTerminator;
    switch (PARAMS.RANSACTerminator) {
        case CRANSACHomographyParams::eSimpleTerminator:
            pTerminator.reset(new CSimpleTerminator(nCount));
            break;
        case CRANSACHomographyParams::eWaldSAC:
            pTerminator.reset(new CWaldSACTerminator(nCount, nSampleSize, 4, 1000, 0.1));
            break;
        case CRANSACHomographyParams::eBrownianBridge:
            //pTerminator.reset(new CBBTerminator(nCount));
            break;
        default:
            THROW("Terminator not recognised");
    }

    //scoped_ptr<CRansacTerminator> pTerminator(new CSimpleTerminator(nCount));
    scoped_ptr<CIterTerminator> pIterTerminator(new CPPIterTerminator(PARAMS.MAX_ITERS, adPriorProbs, nSampleSize, PARAMS.PROB_SUCCESS));
    scoped_ptr<CModelRefiner> pSimpleRefiner(new COpenCVHomography(points0, points1));
    scoped_ptr<CInlierCounter> pInlierCounter(new CHomographyInlierCounter(H_INLIER_THRESH, points0, points1));
    scoped_ptr<CFindBestMatching> pRefineMatches(new CFindBestMatchingDisjoint(pointIds));
    scoped_ptr<CModelRefiner> pRefiner(new CTopdownRLRefiner(pSimpleRefiner.get(), pRefineMatches.get(), pInlierCounter.get(), PARAMS.TOPDOWN_ITERS, PARAMS.TOPDOWN_EXPAND, PARAMS.TOPDOWN_SCALEDOWN));

    return doRansac(pSampler.get(), &hypothesise, pRefiner.get(), pTerminator.get(), pIterTerminator.get(), pInlierCounter.get(), pRefineMatches.get(), H, inliers, nThreads, PARAMS.VERBOSE);
}

//Heuristic estimate of what the minimum inlier count from a good model would be.
//dMinPropInliersGood == min overlap between images?
//Maybe dMinPropInliersGood = 0.25 is realistic if there;s only 25% overlap

int getMinInliers(CInlierProbs & adPriorProbs, const double dMinPropInliersGood, const int nNumOfPointsForModel) {
    double dExpectedInliers = adPriorProbs.sum();
    double dPotentialInliers = dExpectedInliers - nNumOfPointsForModel;
    double dMinInliers = nNumOfPointsForModel + (dPotentialInliers > 0 ? (dPotentialInliers * dMinPropInliersGood) : 0); //Heuristic
    return (int) dMinInliers;
}

inline double sampsonsErr(const C2dPoint & px, const C2dPoint & pxp, Eigen::Matrix3d & E) {
    Eigen::Vector3d x, xp;
    px.asVector(x);
    pxp.asVector(xp);

    Eigen::Vector3d xpE = xp.transpose() * E;
    Eigen::Vector3d Ex = E * x;

    double numerator = xpE.dot(x);
    double denom = xpE.segment(0, 2).squaredNorm() + Ex.segment(0, 2).squaredNorm();

    return sqr(numerator) / denom;
}

//Integrate analysis of what is in front/behind camera into RANSAC:
bool RTFromE(CMask & mask, const C3x3MatModel & model, const T2dPoints & points0, const T2dPoints & points1, const CPointIdentifiers & ids, C3dRotation & R, C3dPoint & T, const double THRESH_SQ) {
    int nInliers = mask.countInliers();
    CPointVec2d a1stPassInlierPoints1;
    a1stPassInlierPoints1.reserve(nInliers);
    CPointVec2d a1stPassInlierPoints2;
    a1stPassInlierPoints2.reserve(nInliers);
    mask2Vectors<C2dPoint > (mask, points0, points1, a1stPassInlierPoints1, a1stPassInlierPoints2);

    CCamera Pp_new;
    if (!chooseCamFromE(model.asDouble9(), a1stPassInlierPoints1, a1stPassInlierPoints2, Pp_new)) {
        cout << "Error: Camera selection failed\n";
        return false;
    }

    R = Pp_new.rotation();
    T = Pp_new.translation();
    Eigen::Matrix3d E;
    if (T.sum_square() == 0)
    {
        cout << "Pure rotation detected" << endl;
        //T = C3dPoint(-1,-1,-1); //Should essentially be arbitrary
        return false;
    }
    //makeE(R, C3dPoint(0,0,1), E);
    //else
    makeE(R, T, E);

    CCamera P;

    CDynArray<double> adSE(mask.size());
    for (int i = 0; i < (int) mask.size(); i++) {

        double dErr = sampsonsErr(points0[i], points1[i], E);
        CHECKNAN(dErr);
        adSE[i] = dErr;
        //cout << dErr << endl;
        mask[i] = false;
        if (dErr < THRESH_SQ) {
            if (testPair(P, Pp_new, points0[i], points1[i]))
                mask[i] = true;
        }
    }
    cout << "Inliers in front: " << mask.countInliers() << endl;

    CFindBestMatchingDisjoint disjointMatch(ids);
    disjointMatch.refine(mask, nInliers, adSE.begin());

    cout << "Inliers in front-disjoint: " << nInliers << endl;

    return nInliers >= 8;
}

bool RTFromE(CMask & mask, const Eigen::Matrix3d & E, const T2dPoints & points0, const T2dPoints & points1, const CPointIdentifiers & ids, C3dRotation & R, C3dPoint & T, const double THRESH_SQ) {
    C3x3MatModel model;
    for (int r = 0; r < 3; r++)
        for (int c = 0; c < 3; c++)
            model(r, c) = E(r, c);
    return RTFromE(mask, model, points0, points1, ids, R, T, THRESH_SQ);
}

//inefficient. Constrict E, reassess inliers/outliers and in front/behind (R,t)

bool RTFromRT(CMask & mask, const T2dPoints & points0, const T2dPoints & points1, const CPointIdentifiers & ids, C3dRotation & R, C3dPoint & T, const double THRESH_SQ) {
    Eigen::Matrix3d E;
    makeE(R, T, E);
    return RTFromE(mask, E, points0, points1, ids, R, T, THRESH_SQ);
}
