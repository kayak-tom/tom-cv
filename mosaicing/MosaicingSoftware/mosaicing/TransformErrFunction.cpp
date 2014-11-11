#include "TransformErrFunction.h"
#include "GRCException.h"
#include "cvGeom.h"
#include <iostream>

namespace grc {

TransformErrFunction::TransformErrFunction(TransformSet * transformSet, int alignTo, int latestImageId, int latestNOnly)
    : transformSet_(transformSet), reprojErrVal_(-HUGE), alignTo_(alignTo), RECalculated_(false), scaleCost_(1.0), minIdToRefine_(latestImageId - latestNOnly)
{
    //Choose some points to project (4 is best I think--necessary+sufficient)
    double pointSquareSideLen = 400; //roughly 4 corners of an image
    pointVec_.push_back(cvPoint2D64f(0, 0));
    pointVec_.push_back(cvPoint2D64f(pointSquareSideLen, 0));
    pointVec_.push_back(cvPoint2D64f(0, pointSquareSideLen));
    pointVec_.push_back(cvPoint2D64f(pointSquareSideLen, pointSquareSideLen));

    //Add all the parameters to the param vector
    for(TransformSet::iterator ppTrans = transformSet_->begin(); ppTrans != transformSet_->end(); ppTrans++)
    {
        int j = ppTrans->first.im1Id();
        int k = ppTrans->first.im2Id();

        if(k == alignTo_ && j != alignTo_) //Refine all transform T_j0 : j!=0
        {
            if(minIdToRefine_ < j)
            {
                Transform * Tj0 = ppTrans->second->transform();
                Tj0->addParams(&indexedParamLocations_);
            }
        }
    }

    //Set scaleCost_ so that costs are around 1 and we get sensible well-conditioned derivatives
    evaluateReprojErr();
    if(reprojErrVal_ > 0.0001) scaleCost_ = 1.0 / reprojErrVal_;
    RECalculated_ = false;
}


double TransformErrFunction::getErr(const Transform * T1, const Transform * T2, int inlierCount)
{
    double projErr = 0;
    for(TPointVec::iterator pPoint = pointVec_.begin(); pPoint != pointVec_.end(); pPoint++)
    {
        CvPoint2D64f point = *pPoint;
        CvPoint2D64f transformedPoint1 = T1->applyToPoint(point);
        CvPoint2D64f transformedPoint2 = T2->applyToPoint(point);
        projErr += SSD(transformedPoint1, transformedPoint2);
    }

    if(eWeightByInlierCountSquared_ == eWeightByInlierCount)
    {
        projErr *= (double)(inlierCount);
    }
    else if(eWeightByInlierCountSquared_ == eWeightByInlierCountSquared)
    {
        projErr *= (double)(inlierCount*inlierCount);
    }

    return projErr;
}

double TransformErrFunction::evaluateReprojErr()
{
    if(RECalculated_)
        return reprojErrVal_;

    reprojErrVal_ = evaluateReprojErrInt();
    RECalculated_=true;
    return reprojErrVal_;
};

double TransformErrFunction::evaluateReprojErrInt()
{
    double reprojErr = 0;
    for(TransformSet::const_iterator ppTrans = transformSet_->begin(); ppTrans != transformSet_->end(); ppTrans++)
    {
        int j = ppTrans->first.im1Id();
        int k = ppTrans->first.im2Id();

        if(j != alignTo_ && k != alignTo_ && (j > minIdToRefine_ || k > minIdToRefine_))
        {
            IdPair j0(j, alignTo_);
            IdPair k0(k, alignTo_);
            const Transform * Tj0 = (*transformSet_)[j0]->transform();
            const Transform * Tk0 = (*transformSet_)[k0]->transform();

            const Transform * Tjk = ppTrans->second->transform();
            int inlierCount = ppTrans->second->correspondences().size();

            Transform * Tj0_from_jk = Tk0->accumulate(Tjk);
            reprojErr += getErr(Tj0, Tj0_from_jk, inlierCount);
            delete Tj0_from_jk;
        }
        /*else if(j != alignTo_ ) //Here ppTrans->first.second == alignTo_. Don't calc reproj err for identity!
        {
            const Transform * Tj0 = ppTrans->second->transform();
            reprojErr += getErr(Tj0);

        }*/
    }
    return reprojErr * scaleCost_; //scaleCost_ makes errors initially about 1
}
double * TransformErrFunction::getParamRef(int idx) const
{
    return indexedParamLocations_[idx].paramRef_;
}

double TransformErrFunction::getParamVal(int idx) const
{
    return *(indexedParamLocations_[idx].paramRef_) * indexedParamLocations_[idx].scale_inv_;
}

void TransformErrFunction::setParamVal(int idx, double val)
{
    *(indexedParamLocations_[idx].paramRef_) = val*indexedParamLocations_[idx].scale_;
}

void TransformErrFunction::assignAllVals(cParamVals * valVector)
{
    RECalculated_ = false;
    size_t paramCount = indexedParamLocations_.size();
    if(valVector->size() != paramCount)
        throw new GRCException("TransformErrFunction::assignAllVals: Wrong number of params");

    for(size_t i=0; i < paramCount; i++)
    {
        //double * param = indexedParamLocations_[i];
        //*param =  * getParamFactor(i);
        setParamVal(i, (*valVector)[i]);
    }
}

TransformErrFunction::cParamVals * TransformErrFunction::getVals() const
{
    cParamVals * vals = new cParamVals(indexedParamLocations_.size());

    size_t paramCount = indexedParamLocations_.size();
    for(size_t i=0; i < paramCount; i++)
    {
        (*vals)[i] = getParamVal(i);
    }
    return vals;
}

/*double TransformErrFunction::getParamFactor(int idx) const//Todo: properly!!
{
    idx %= 8;
    if(idx==2 || idx == 5) return 1000.0;
    if(idx == 6 || idx == 7) return 0.01;
    return 1.0;
}
double TransformErrFunction::getParamFactorInv(int idx) const//Todo: properly!!
{
    idx %= 8;
    if(idx == 2 || idx == 5) return 0.001;
    if(idx == 6 || idx == 7) return 100.0;
    return 1.0;
}*/

double TransformErrFunction::pd(int idx)
{
    if(!RECalculated_)
        throw new GRCException("TransformErrFunction::pd: Reproj. Err. not calculated");
    if(idx < 0 || idx >= (int)indexedParamLocations_.size())
        throw new GRCException("TransformErrFunction::pd: Idx OOB");

    /*double * param = indexedParamLocations_[idx];

    double paramFactor = getParamFactor(idx); //Scale

    double paramVal = *param;
	double delta_inv, delta = 1e-4*fabs(paramVal); //Todo: constants...
	if (delta < 1e-6) {
		delta = 1e-6;
        delta_inv = 1e+6;
	}
    else
        delta_inv = 1.0/delta;

    *param = paramVal + delta;
    double f_x_plus_delta = evaluateReprojErrInt(); //Todo: only need to measure change in reproj error--look only at the few transforms this val affects (all transforms including the image containing this point)
    *param = paramVal;

    return (f_x_plus_delta - reprojErrVal_) * delta_inv * paramFactor;*/

    double paramVal = getParamVal(idx);

	double delta_inv, delta = 1e-4*fabs(paramVal); //Todo: constants...
	if (delta < 1e-6) {
		delta = 1e-6;
        delta_inv = 1e+6;
	}
    else
        delta_inv = 1.0/delta;

    setParamVal(idx, paramVal + delta);
    double f_x_plus_delta = evaluateReprojErrInt(); //Todo: only need to measure change in reproj error--look only at the few transforms this val affects (all transforms including the image containing this point)
    setParamVal(idx, paramVal);

    return (f_x_plus_delta - reprojErrVal_) * delta_inv;
};

/*void TransformErrFunction::incVal(int idx, double delta)
{
    double * param = indexedParamLocations_[idx];
    *param += delta;
};

void TransformErrFunction::decVal(int idx, double delta)
{
    double * param = indexedParamLocations_[idx];
    *param -= delta;
};*/

//! Equivalent to frobenius norm for Lev. Mar.
double TransformErrFunction::cParamVals::sumSquare() const
{
    double sum_square=0;
    size_t paramCount = size();
    for(size_t i=0; i < paramCount; i++)
    {
        double param = (*this)[i];
        sum_square += param*param;
    }
    return sum_square;
}

double TransformErrFunction::cParamVals::frobeniusNorm() const
{
    return sqrt(sumSquare());
}

void TransformErrFunction::cParamVals::printDiff(const cParamVals * params2) const
{
    size_t paramCount = size();
    for(size_t i=0; i < paramCount; i++)
    {
        if(i % 8 == 2 || i % 8 == 5) //translation
        {
            double param1 = fabs((*this)[i]);
            double param2 = fabs((*params2)[i]);
            //if((param1 < param2*0.9 || param1 > param2*1.1) && param1 > .1)
            if(fabs(param1-param2) > 0.001)
                std::cout << i << ": param1=" << param1 << ", param2=" << param2 << std::endl;
        }
    }
}

TransformErrFunction::cParamVals::cParamVals(size_t len) :  std::vector<double>(len) {};

TransformErrFunction::cParamVals::cParamVals(const cParamVals * oldParams, const CvMat * delta) : std::vector<double>(oldParams->size())
{
    const_iterator pOldParam = oldParams->begin();
    double * pUpdate = delta->data.db;
    for(iterator pParam = begin(); pParam < end(); pParam++, pOldParam++, pUpdate++)
    {
        *pParam = *pOldParam + *pUpdate;
    };
}


}


