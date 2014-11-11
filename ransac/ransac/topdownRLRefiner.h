/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * topdownRLRefiner.h
 *
 *  Created on: 29/10/2009
 *      Author: tom
 */

#ifndef TOPDOWNRLREFINER_H_
#define TOPDOWNRLREFINER_H_

#include "util/exception.h"
#include "refiner.h"
#include "cRansacTerminator.h"

//Topdown refinement by fitting a model, discarding least compatible inliers and adding most compatible outliers, then iterating
//See Robust regression and outlier detection, Rousseeuw, P. J. and Leroy, A. M., 1987

class CTopdownRLRefiner : public CModelRefiner {
    CModelRefiner * pSimpleRefiner;
    CFindBestMatching * pMatchRefiner;
    CInlierCounter * pGetResiduals;
    const int nNumIters;
    const double dExpandThresh, dScaleThresh;
    //const T2dPoints & p1, & p2;

    virtual bool fitModel_int(CMask & mask, CModel & pModel, bool bVerbose) {
        CRansacTerminator doNotTerminate;
        ARRAY(double, adResiduals, (const int) mask.size());
        DEBUGONLY(setConstant(PTR(adResiduals), -1.0, mask.size()));

        double *pdResiduals = 0;
        if (pMatchRefiner->supplyResiduals())
            pdResiduals = PTR(adResiduals);

        int nInliers = mask.countInliers();

        if (bVerbose)
            mask.pp();

        int nOldChecksum = mask.checksum();
        double dThresholdScale = dExpandThresh;

        for (int i = 0; i < nNumIters; i++) //Todo: terminate on convergence or when mask is unchanged
        {
            if (nInliers < pSimpleRefiner->minNumPoints()) return false;
            if (!pSimpleRefiner->fitModel(mask, pModel, bVerbose))
                return (i > 0); //Success if we refined one previously (should still be in pModel, actually probably isn't), otherwise fail
            
            //cout << "Model r1 during refinement:\n" << endl;

            mask.setZero();

            pGetResiduals->countInliers(pModel, &doNotTerminate, 0, 0, mask, pdResiduals, nInliers, dThresholdScale);

            if (bVerbose)
            {
                cout << "Iter " << i << " ";
                cout << "(" << mask.countInliers() << ")\n";
                printModel(pModel);
            }

            pMatchRefiner->refine(mask, nInliers, pdResiduals);

            int nNewChecksum = mask.checksum();

            if (nNewChecksum == nOldChecksum) {
                //Mask hasn't changed, have fitted model to mask
                break;
            } else
                nOldChecksum = nNewChecksum;

            //if (bVerbose)
            //    mask.pp();

            dThresholdScale *= dScaleThresh;
        }

        DEBUGONLY(
        if (pdResiduals) {
            for (int i = 0; i < (int) mask.size(); i++)
                if (pdResiduals[i] == -1)
                        std::cout << "pdResiduals[" << i << "] uninitialised\n";
                })

        if (bVerbose && pdResiduals) {
            //std::cout << "Inlier residuals (should be approx. Gaussian): [";
            int nPoints = mask.size();

            double dInlierMean = 0, dOutlierMean = 0; //, dOutlierMin = HUGE;

            for (int nPoint = 0; nPoint < nPoints; nPoint++) {
                if (mask[nPoint]) {
                    //std::cout << pdResiduals[nPoint] << ", ";
                    dInlierMean += pdResiduals[nPoint];
                } else {
                    dOutlierMean += pdResiduals[nPoint];
                    //if(dOutlierMin > pdResiduals[nPoint]) dOutlierMin = pdResiduals[nPoint];
                }
            }
            //std::cout << "]\n";

            int nInliers = mask.countInliers();
            if (nInliers != nPoints && nInliers > 0) {
                dInlierMean /= nInliers;
                dOutlierMean /= nPoints - nInliers;

                double dRatio = dOutlierMean / dInlierMean;

                std::cout << "Inlier mean residual = " << dInlierMean << std::endl;
                std::cout << "Outlier mean residual (approx) = " << dOutlierMean << std::endl;
                std::cout << "Ratio = " << dRatio << std::endl;

                if (dRatio < 5)
                    std::cout << "Warning: inliers have high residual compared with outliers. Consider dropping RANSAC threshhold??\n";

                double dAreInliersBellshaped = dInlierMean / (dThresholdScale * pGetResiduals->threshold_sq());
                if (dAreInliersBellshaped < 0.1)
                    std::cout << "Outliers probably included in soln; ratio of inlier mean to thresh=" << dAreInliersBellshaped << std::endl;
                if (dAreInliersBellshaped > 0.4)
                    std::cout << "Inliers probably missed from soln; ratio of inlier mean to thresh=" << dAreInliersBellshaped << std::endl;
            }
        }

        return true;
    }
public:
    //dExpandThresh increases the inlier threshold by this scale factor on the first pass. It is progressively shrunk by dScaleThresh afterwards.

    CTopdownRLRefiner(CModelRefiner * pSimpleRefiner, CFindBestMatching * pMatchRefiner, CInlierCounter * pGetResiduals, const int nNumIters, const double dExpandThresh = 1.0, const double dScaleThresh = 1.0) : CModelRefiner(pSimpleRefiner->minNumPoints()), pSimpleRefiner(pSimpleRefiner), pMatchRefiner(pMatchRefiner), pGetResiduals(pGetResiduals), nNumIters(nNumIters), dExpandThresh(dExpandThresh), dScaleThresh(dScaleThresh) {

    }

    virtual bool fitModel(CMask & mask, CModel & pModel, bool bVerbose) {
        return fitModel_int(mask, pModel, bVerbose && IS_DEBUG);
    }

};

#endif /* TOPDOWNRLREFINER_H_ */
