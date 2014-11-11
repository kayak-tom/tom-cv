//! Implements transform engine: Most logic is done here. Requests and organises sets of image+transformation (from 0) pairs that the rendering engine will draw the latest of.
#include "TransformSet.h"
#include "GRCException.h"
#include <iostream>

namespace grc {

//! 
TransformSet::TransformSet() :
onlyRenderLatestTransform_(false), mosaicTransform_(0), latestTransform_(0,0) {};

TransformSet::~TransformSet()
{
    for(TTransformSet::iterator ppTI = begin(); ppTI != end(); ppTI++)
        delete ppTI->second;
    delete mosaicTransform_;
};

void TransformSet::erase(const IdPair & id)
{
	iterator pToErase = find(id);
    if(pToErase != end())
    {
	    delete pToErase->second;
	    TTransformSet::erase(pToErase);
    } 
};

bool TransformSet::exists(IdPair ids) const 
{ 
    return find(ids) != end(); 
};

bool TransformSet::onlyRenderLatestTransform() const 
{
    return onlyRenderLatestTransform_; 
};
const TransformInfo2 * TransformSet::latestTransform() const
{ 
    //if(!onlyRenderLatestTransform_)
        //throw new GRCException("Transform::latestTransform: We should be rendering all transforms");

	if(onlyRenderLatestTransform_)
		return find(latestTransform_)->second;
	else
	{
		size_t nLatestId = 0;
		TTransformSet::const_iterator pLatest = end();
	    for(TTransformSet::const_iterator ppTI = begin(); ppTI != end(); ppTI++)
	    {
	        size_t thisId = ppTI->first.im2Id();
	        if(thisId > nLatestId)
	        {
	        	pLatest = ppTI;
	        	nLatestId = thisId;
	        }
	    }
	    if(nLatestId>0)
	    	return pLatest->second;

    	return 0;
	}
};

const Transform * TransformSet::mosaicTransform() const 
{
    if(!onlyRenderLatestTransform_)
        throw new GRCException("Transform::mosaicTransform: We should be rendering all transforms");

    if(!mosaicTransform_)
        throw new GRCException("Transform::mosaicTransform: mosaicTransform_ should be initialised");

    return mosaicTransform_;
};

void TransformSet::setLatestTransform(IdPair latestTrans)
{
    latestTransform_ = latestTrans;
};

void TransformSet::setFastRender(Transform * mosaicTrans, IdPair latestTrans)
{
    //todo check
    onlyRenderLatestTransform_ = true;
    mosaicTransform_ = mosaicTrans;
    latestTransform_ = latestTrans;
};
	
}
