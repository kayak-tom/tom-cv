#pragma once
#ifndef _TRANSFORM_ERR_FN
#define _TRANSFORM_ERR_FN

#include "util/set2.h"
#include "TransformSet.h"
#include "util/opencv.h"

namespace grc {

    //! Represents a function that measures the inconsistency of a transform set. Refine using LM to make more mutually consistent.
    class TransformErrFunction
    {
        //! Evaluate function
        double evaluateReprojErrInt();

        Transform::TParamLocations indexedParamLocations_; //! < Locations of parameters to refine.
        
        TransformSet * transformSet_; //!<This transform set is actually adjusted. Adjust by small amounts to get derivatives, then update vals here.
        double reprojErrVal_; //! < Cache reprojection error value.
        bool RECalculated_; //!< Reproj. error has been calculated since params were last changed (consistency check really0
        enum eTEWeightFunction { eEqualWeight, eWeightByInlierCountSquared, eWeightByInlierCount};
        static const eTEWeightFunction eWeightByInlierCountSquared_ = eWeightByInlierCountSquared; 

        typedef std::vector<CvPoint2D64f> TPointVec;
        TPointVec pointVec_; //!< Set of points to use to measure alignment error.
        int alignTo_; //!< Frame that all frames have a transformation to.

        double scaleCost_;//!< Set scaleCost_ so that costs are around 1 and we get sensible well-conditioned derivatives
        const int minIdToRefine_; //!< Sometimes (incremental rendering) we don't want to refine all parameters--just refine the last few

        //! Return inconsistency between two transforms
        double getErr(const Transform * T1, const Transform * T2, int inlierCount);
        
        //! Return pointer to a parameter
        double * getParamRef(int idx) const;
        
        //! Return scaled value of a parameter
        double getParamVal(int idx) const;
        
        //! Set parameter val
        void setParamVal(int idx, double val);

    public:
        //! Represents a set of parameter values
        class cParamVals : public std::vector<double> 
        {
        public:
            cParamVals(size_t len);

            //! Make new parameter set by adding delta to an old parameter set
            cParamVals(const cParamVals * oldParams, const CvMat * delta);
            
            //!Return sum of square of parameter vals
            double sumSquare() const;
            
            //!Return Frobenius norm of parameter vals
            double frobeniusNorm() const;

            //! Output how much two parameter sets have changed
            void printDiff(const cParamVals * params2) const;
        };

        //! Output how much two parameter sets have changed
        size_t numParams() const { return indexedParamLocations_.size(); };

        ~TransformErrFunction()
        {
            //Does not own the transform set, although will have modified values in it during minimisation

        };

        //! Assign particular vals to the parameters
        void assignAllVals(cParamVals * valVector);

        //! Get current parameter vals
        cParamVals * getVals() const;

        //! Calculate partial derivatives by adjust..evaluate..adjust back
        double pd(int idx);
        
        //! Setup function from a transform set
        /* \param latestNOnly Only refine last N transforms (for speed) */
        TransformErrFunction(TransformSet * transformSet, int alignTo, int latestImageId, int latestNOnly);

        //! Evaluate error function (using a cached value if possible)
        double evaluateReprojErr();
    };

}


#endif // _TRANSFORM_ERR_FN
