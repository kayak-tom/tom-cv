/*
 * CFullRelPose: A relative position including a scale estimate and covariance estimate. Can now write it to file or make a map
 *
 *  Created on: 14 Jan 2010
 *      Author: tom
 */

#ifndef FULLRELPOSE_H_
#define FULLRELPOSE_H_

#include "bow/scales.h"
#include "geom/geom.h"
#include "scale/EdgeScaleObserver.h"

class CFullRelPose {
    const C3dNormalisedPoseWithSD normalisedPose;
    CScale SLAMscale, //Calculated from previous combined scale + rel scale
    SLAMonlyScale, //Calculated from previous SLAM scale + rel scale
    SCOREscale, //Also used for hacky scale from constraining the scale of the world
    combinedScale; //Combining SLAMscale and SCOREscale

    void updateCombinedScale() {
        bool bNewCompOrigin = SLAMscale.isNewComponentOrigin();
        if (!SCOREscale.hasScale() || (SLAMscale.isOrigin() && !bNewCompOrigin))
            combinedScale = SLAMscale;
        else if (!SLAMscale.hasScale()
                || bNewCompOrigin)// Use SCORE scale if we're starting a new component
            combinedScale = SCOREscale;
        else {
            double dNewMean = CEdge::NO_SCALE(), dNewVar = CEdge::NO_SCALE();

            CVariablesFromMultiDistnsObserved::combineTwoObservations(SLAMscale.getD(), SLAMscale.getG_sq(), SCOREscale.getD(), SCOREscale.getG_sq(), dNewMean, dNewVar);
            if(IS_DEBUG) CHECK(!(dNewMean > CEdge::NO_SCALE() && dNewMean < HUGE && dNewVar >= 0 && dNewVar < HUGE), "Combining distns failed");

            //if(!within(SLAMscale.getD(), dNewMean, 0.2))
            //cout << "SCORE updated log-scale from " << SLAMscale.getD() << " to " << dNewMean << "; "<< SLAMscale.getG_sq() << " to " << dNewVar << " variance change\n";

            combinedScale = CScale(dNewMean, dNewVar);
        }
    }
public:

    CFullRelPose(const C3dNormalisedPoseWithSD & normalisedPose) : normalisedPose(normalisedPose) {
        if (normalisedPose.SD.havePointDepth()) {
            const double dD_baseline = normalisedPose.SD.baselineMean();
            const double dG_sq_baseline = normalisedPose.SD.baselineVar();
            if(IS_DEBUG) CHECK(dG_sq_baseline <= 0, "MM needs positive var")

            /*static const double dD_baseLineToDepthRatio = 4.4;
            static const double dG_sq_baseLineToDepthRatio = 1.5; //Can actually set this pretty high as certainty drops rght off away from origin...

            const double dD_baseline = dD_pointDepth - dD_baseLineToDepthRatio,
                                     dG_sq_baseline = dG_sq_pointDepth + dG_sq_baseLineToDepthRatio;*/

            SCOREscale = CScale(dD_baseline, dG_sq_baseline);
        }
    }

    virtual ~CFullRelPose() {
    }

    bool hasPosition() const {
        return combinedScale.hasScale();
    }

    //Unset SLAM scale, before Dijkstra update

    void setUnused() {
        SLAMscale.setUninit();
        SLAMonlyScale.setUninit();
        //SCOREscale.setUninit();
        //combinedScale.setUninit();
        updateCombinedScale();
    }

    C3dPose position(bool bReverse) const {
        if (!bReverse)
            return normalisedPose.scale(combinedScale.scaleEstimate());
        else
            return normalisedPose.reverse().scale(combinedScale.scaleEstimate()); //Todo: not entirely sure this shouldn't be 1/scale...
    }

    void setSCOREScale(const CScale & s) {
        if(IS_DEBUG) CHECK(normalisedPose.SD.usePointDepths(), "Point depth constraint on motion is incompatible with SCORE")
        SCOREscale = s;
        updateCombinedScale();
    }

    const CScale & getScale() const //unused?
    {
        return combinedScale;
    }

    const CScale & getSLAMscale() const //Includes earlier combined scales
    {
        return SLAMonlyScale;
    }

    //Combine SCORE and accumulated scales

    double length() const {
        return combinedScale.scaleEstimate();
    }

    double variance() const {
        //Use SCORE scale, should also use loop-optimised var (not available atm)
        return combinedScale.scaleVarGeometric();
    }

    void setRoot(bool bIsFirstComponent) {
        SLAMscale.setOrigin(bIsFirstComponent); //Should set one true origin for 1st component, then uncertain origins later
        SLAMonlyScale.setOrigin(true);
        //SCOREscale.setUninit(); //should probably already be uninit... In a new component should really get scale from SCORE,
        if (SCOREscale.hasScale()) {
            cout << "Setting origin where we have a SCORE scale\n";
        }
        updateCombinedScale();
    }

    void updateScale(const CScale & prevScale, const CScale & prevSLAMOnlyScale, const CRelScale & relScale) {
        SLAMscale = prevScale /*includes SCORE*/ + relScale;
        SLAMonlyScale = prevSLAMOnlyScale + relScale;
        updateCombinedScale();
    }

    void writeEdge(std::ostream & toroGraphFile, const bool b2d) const;
};

std::ostream& operator<<(std::ostream& s, const CFullRelPose & X);

#endif /* FULLRELPOSE_H_ */
