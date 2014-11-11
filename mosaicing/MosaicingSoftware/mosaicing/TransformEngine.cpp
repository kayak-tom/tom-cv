#include "util/exception.h"
//! Implements transform engine: Most logic is done here. Requests and organises sets of image+transformation (from 0) pairs that the rendering engine will draw the latest of.
pragma_warning(push)
pragma_warning(disable:4996)  // warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)

#include "TransformEngine.h"
#include "TransformSet.h"
#include "TransformErrFunction.h"
#include "LevenbergMarquardt.h"
#include "GRCException.h"
#include "boost/bind.hpp"
#include "cvGeom.h"

namespace grc {

//! Will spawn a thread that generates a new transformation set every time we are notified that a new image is available
TransformEngine::TransformEngine(TransformEstimator2 & transformEstimator, int maxFrames, bool useIncrementalRendering, int LMIterations, CvSize mosaicSize, Enums::eFrameChoiceMethod frameChoiceMethod, int fullFrameUpdateFreq, int maxSearchForTransform)
    : transformEstimator_(transformEstimator), latestFrameId_(NO_FRAME), semWaitForFrame_(0), semWaitForTransformSet_(0), semWaitForRenderer_(1), semRendererFinished_(0), finished_(false), latestTS_(0),
    maxFrames_(maxFrames), doFastRender_(useIncrementalRendering), maxLMIterations_(LMIterations), mosaicSize_(mosaicSize), imageSize_(cvSize(0,0)),
    positionMethod_(useIncrementalRendering ? Enums::ePosCentre : Enums::ePosLatest),
    frameChoiceMethod_(frameChoiceMethod),
    fullFrameUpdateFreq_(fullFrameUpdateFreq),
    maxSearchForTransform_(maxSearchForTransform)
{
    //Start thread
    mainEngineThread_ = new boost::thread(boost::bind(&grc::TransformEngine::mainTransformLoop, this));
};

TransformEngine::~TransformEngine()
{

    if(!finished_)
        throw new GRCException("TransformEngine::~TransformEngine: Should 1) be finished, and 2) my thread should have returned by now");
};

//!We don't know image sizes until capture has started
bool TransformEngine::notifySize(CvSize s)
{
    imageSize_ = s; //todo checks
    return true;
}

void TransformEngine::mainTransformLoop()
{
    try
    {
        for(;;)
        {
            size_t imageId = NO_FRAME;

            semWaitForFrame_.wait();

            {
                boost::mutex::scoped_lock renderDisplayLock(notifyLock_);
                if(finished_)
                {
                    setCurrentTransform(0); //Tell renderer
                    break; //and return
                }

                //We have a new image
                imageId = latestFrameId_;
            } //stop locking--new images can arrive while we compute transforms

            TransformSet * TS = computeTransforms(imageId);
            setCurrentTransform(TS);
        }
    }
    catch(grc::GRCException * pEx)
    {
        std::cout << "ERROR: Unhandled exception in transform engine main thread: " << *(pEx->GetErrorMessage()) << std::endl ;
        delete pEx;
    }

}

//!Called when new images are available to render
void TransformEngine::notifyImageAvailable(size_t id)
{
    boost::mutex::scoped_lock renderDisplayLock(notifyLock_);
    if(id == NO_FRAME || (latestFrameId_ != NO_FRAME && id != (latestFrameId_ + 1)))
        throw new GRCException("TransformEngine::notifyImageAvailable: Id not sequential (as was expected)");

    latestFrameId_ = id;
    semWaitForFrame_.post();
}

void TransformEngine::notifyRendererFinished()
{
    semRendererFinished_.post();
}

void TransformEngine::notifyFinished()
{
    {
        boost::mutex::scoped_lock renderDisplayLock(notifyLock_);
        finished_ = true;
        semWaitForFrame_.post();
    }
    //Now block until the renderer is finished
    semRendererFinished_.wait();

    //and my thread is finished
    mainEngineThread_->join();
    delete mainEngineThread_; mainEngineThread_=0;
}

//! Apply transform to 4 corners of source image. If any are further away than current max bounds then adjust max bounds
void TransformEngine::findDistantCorners(const Transform * pTrans, double & globalMaxX, double & globalMaxY, double & globalMinX, double & globalMinY) const
{
    static const CvMat * sourceCorners = getSourceCorners(imageSize_);
    static CvMat * projectedCorners = getSourceCorners(imageSize_);

    pTrans->applyToPoints(sourceCorners, projectedCorners);

    for(int i=0; i<4; i++)
    {
        double x = cvmGet(projectedCorners, 0, i);
        if(x<globalMinX) globalMinX=x;
        else if(x>globalMaxX) globalMaxX=x;

        double y = cvmGet(projectedCorners, 1, i);
        if(y<globalMinY) globalMinY=y;
        else if(y>globalMaxY) globalMaxY=y;
    }
}

//!Calculate shift so that shift+localMin and shift+localMax fall within mosaic, and hopefully global ones do too
double TransformEngine::position1d(double localMin, double localMax, double globalMin, double globalMax, double width)
{
    double shiftMin = -localMin;
    double shiftMax = width-localMax;
    if(shiftMin>=shiftMax)
        return shiftMin; //Image doesn't fit in mosaic, can't do much

    if(globalMax-globalMin < width) //Plenty of room--centre:
        return -globalMin + (width - (globalMax-globalMin))*0.5;
    else if(globalMax-localMin > width)
        return shiftMin; //Mosaic is far too wide--stick latest frame on left
    else if(localMax-globalMin > width)
        return shiftMax; //Mosaic is far too wide--stick latest frame on right

    //Mosaic is too wide, but we have room for the latest im + some past frames on both sides: fit all of the smaller side on
    if(localMin-globalMin < globalMax-localMax)//More on the RH side
        return -globalMin;

    //More on the LH side:
    return width-globalMax;
};

//! Translate all transforms so latest image is visible and as many previous images as possible are visible
/*! TODO: Cull images that may be out of view */
void TransformEngine::positionMosaic(TransformSet * TS, int alignToId, int latestImageId)
{
    if(positionMethod_ == Enums::ePosCentre)
    {
        double top = mosaicSize_.height/2 - imageSize_.height/2;
        double left = mosaicSize_.width/2 - imageSize_.width/2;
        for(TransformSet::iterator ppTrans = TS->begin(); ppTrans != TS->end(); ppTrans++)
        {
            Transform * pTrans = ppTrans->second->transform();

            pTrans->shift(left, top);
        }
    }
    else if(positionMethod_ == Enums::ePosLatest)
    {
        //Get positions of corners of most recent and most distant frames. Move most recent to be at the best corner
        IdPair latestTransId(latestImageId, alignToId);
        const Transform * transLatest = (*TS)[latestTransId]->transform();

        double latestMaxX=-HUGE, latestMaxY=-HUGE, latestMinX=HUGE, latestMinY=HUGE;
        findDistantCorners(transLatest, latestMaxX, latestMaxY, latestMinX, latestMinY);

        double globalMaxX=-HUGE, globalMaxY=-HUGE, globalMinX=HUGE, globalMinY=HUGE;

        //Iterate thro transforms and compute all corners:
        for(TransformSet::iterator ppTrans = TS->begin(); ppTrans != TS->end(); ppTrans++)
        {
            const Transform * pTrans = ppTrans->second->transform();
            findDistantCorners(pTrans, globalMaxX, globalMaxY, globalMinX, globalMinY);
        }

        //Now choose shift:
        double shiftX = position1d(latestMinX, latestMaxX, globalMinX, globalMaxX, mosaicSize_.width);
        double shiftY = position1d(latestMinY, latestMaxY, globalMinY, globalMaxY, mosaicSize_.height);

        //and apply
        for(TransformSet::iterator ppTrans = TS->begin(); ppTrans != TS->end(); ppTrans++)
        {
            Transform * pTrans = ppTrans->second->transform();

            pTrans->shift(shiftX, shiftY);
        }
    }

    //And cull frames out-of-view (don't bother as they're almost free to render anyway)
}

//! Set the latest set of image+transformations. Set 0 to kill renderer.
void TransformEngine::setCurrentTransform(TransformSet * TS)
{
    if(renderAllFrames_)
        semWaitForRenderer_.wait();

    boost::mutex::scoped_lock giveAwayTransformLock(giveAwayTransformLock_);

    if(!finished_ && TS==0) throw new GRCException("TransformEngine::setCurrentTransform: Null param received before we have finished");

    if(latestTS_)
    {
        delete latestTS_;
        //std::cout << "Warning: TransformEngine::setCurrentTransform: A new transform set has been computed before the old one was rendered. Deleting old transform set.\n";
    }
    latestTS_ = TS;

    semWaitForTransformSet_.post();
};

//! Returns the latest set of image+transformations.
/*! The rendering engine calls this. It will block if we haven't
 *  finished computing a new one. Will always return the latest transform: it
 *  will discard a transformation if a new one is calculated before the previous one has been geto by the renderer.
 *  this class has generated a new transform before the last one has been rendered.
 */
TransformSet * TransformEngine::getLatestTransforms()
{
    semWaitForTransformSet_.wait();

    boost::mutex::scoped_lock giveAwayTransformLock(giveAwayTransformLock_);

    TransformSet * TStoReturn = latestTS_;
    latestTS_ = 0;

    if(TStoReturn && renderAllFrames_)
        semWaitForRenderer_.post();

    return TStoReturn;
}

//! Delete each transformation x->y with no x->0 (or y->0)
void TransformEngine::removeDisconnectedTransforms(TransformSet * TS, int alignToImageId)
{
    std::vector<IdPair> idPairsToDelete;
    for(TransformSet::const_iterator ppTrans = TS->begin(); ppTrans != TS->end(); ppTrans++)
    {
        IdPair idPair = ppTrans->first;

        //Should always have added both if connected, or neither if disconnected
        IdPair idConnectedPair1(idPair.im1Id(), alignToImageId);
        IdPair idConnectedPair2(idPair.im1Id(), alignToImageId);

        if(!TS->exists(idConnectedPair1))
        {
            idPairsToDelete.push_back(idPair); //don't actually delete while iterating
        }

        if(TS->exists(idConnectedPair1) != TS->exists(idConnectedPair2))
            throw new GRCException("TransformEngine::removeDisconnectedTransforms: Transform from connected image to 0 missing");
    }
    for(std::vector<IdPair>::iterator pId = idPairsToDelete.begin(); pId != idPairsToDelete.end(); pId++)
    {
        //cout << "erasing " << pId->im1Id() << ',' << pId->im2Id() << endl;
        TS->erase(*pId);
    }
}

//! Delete each transformation x->y with y != latestImageId
void TransformEngine::cleanupTransformSet(TransformSet * TS, int alignToImageId)
{

    std::vector<IdPair> idPairsToDelete;
    for(TransformSet::const_iterator ppTrans = TS->begin(); ppTrans != TS->end(); ppTrans++)
    {
        IdPair idPair = ppTrans->first;
        if((int)idPair.im2Id() != alignToImageId)
            idPairsToDelete.push_back(idPair); //don't actually delete while iterating
    }
    for(std::vector<IdPair>::iterator pId = idPairsToDelete.begin(); pId != idPairsToDelete.end(); pId++)
    {
    	//delete pId->second;
        TS->erase(*pId);
    }
}

//! For each transformation x->pivotId add x->latestImageId and recurse with x as pivot
void TransformEngine::alignTransformToLastImage(TransformSet * TS, int pivotId, int alignToImageId, bool bFirstTime)
{
    IdPair idPivToLatest(pivotId, alignToImageId);
    //std::pair<IdPair, TransformInfo2 *>
    TransformSet::const_iterator transPivotToLatest = TS->find(idPivToLatest);
    if(transPivotToLatest == TS->end())
        throw new GRCException("TransformEngine::alignTransformToLastImage: Pivot->latest transform does not exist");


    for(TransformSet::const_iterator ppTrans = TS->begin(); ppTrans != TS->end(); ppTrans++)
    {
        IdPair idPair = ppTrans->first;
        if((int)idPair.im2Id() == pivotId)
        {
            IdPair idPairNew(idPair.im1Id(), alignToImageId); //This is the transform we're after

            bool recurse = bFirstTime;
            if(!TS->exists(idPairNew))
            {
                TS->init(idPairNew, ppTrans->second->accumulate(transPivotToLatest->second));
                recurse = true;
            }
            if(recurse)
                //and recurse (DO NOT need to recurse if this transform already existed, or will get infinite loops)
                alignTransformToLastImage(TS, idPair.im1Id(), alignToImageId, false); //would be complex, and would miss some in loops if we tried deleting here
        }
        //Do the same for pivotId->x transforms (take inverse transformations)
        //(alternatively could make sure transforms are always directed lower-higher id. This wouldn't work if we ever want to align all images to one a few frames back from the latest)
        if((int)idPair.im1Id() == pivotId)
        {
            IdPair idPairNew(idPair.im2Id(), alignToImageId); //This is the transform we're after
            bool recurse = bFirstTime;
            if(!TS->exists(idPairNew))
            {
                TS->init(idPairNew, transPivotToLatest->second->accumulateInverse(ppTrans->second));
                recurse = true;
            }
            if(recurse)
                //and recurse (DO NOT need to recurse if this transform already existed, or will get infinite loops)
                alignTransformToLastImage(TS, idPair.im2Id(), alignToImageId, false); //would be complex, and would miss some in loops if we tried deleting here
        }
    }
}

/*! First add all transformations x->alignToImageId for all images
 *  Then delete all transformations not pointing to alignToImageId
 *  alignToImageId is normally the latest image
 */
void TransformEngine::alignTransformsToOneImage(TransformSet * TS, int alignToImageId)
{
    //This will recurse and add a transform x->alignToImageId for every x where this path exists
    alignTransformToLastImage(TS, alignToImageId, alignToImageId, true);

    //Delete all transforms not pointing to latest
    cleanupTransformSet(TS, alignToImageId);
}

/*//! Use bundle adjustment to refine transform set
void TransformEngine::refineTranformSetAndPoints(TransformSet * TS)
{

    //Todo: check correspondences copied ok (not likely atm)
    //Link points
    LinkedTransformSet linkedTS(TS); //Creates a wrapper with correspondences as links

    //setup error function
    ReprojErrorFunction reprojErr(linkedTS);

    // Now we have a function structure and params: reprojErr(p) = reprojection err
    // reprojErr(linkedTS.currentParamVal()) is the reprojection err at the current param values
    BAParamSet * newParams = minimiseLevenbergMarquardt(reprojErr, linkedTS.currentParamVal());

    TS->updateTransforms(newParams);
    delete newParams;

}*/

void TransformEngine::refineTranformSet(TransformSet * TS, int alignTo, int latestId)
{
    static int alignToIdLastTime = NO_FRAME;

    const int MAX_NUM_TRANS_TO_ADJUST = 25;
    int numTransToAdjust = MAX_NUM_TRANS_TO_ADJUST;
    //if(doFastRender_ && alignToIdLastTime == alignTo) //We're just doing an incremental render so only refine the latest transform
      //  numTransToAdjust = 1; //only adjust previous transform

    TransformErrFunction errFn(TS, alignTo, latestId, numTransToAdjust);

    if(errFn.numParams() > 0)
        LevenbergMarquardt LMminimise(errFn, maxLMIterations_);

    alignToIdLastTime = alignTo;
}

int TransformEngine::getAlignToId(int latestImageId, const std::set<size_t> & permittedAlignToIds)
{
    int alignToId = latestImageId > 4 ? latestImageId-4 : latestImageId;
    while(permittedAlignToIds.find(alignToId) == permittedAlignToIds.end() && alignToId != latestImageId)
        alignToId++;

    return alignToId;
}
int TransformEngine::getFastAlignToId(int latestImageId, const std::set<size_t> & permittedAlignToIds)
{
    //const int FAST_FULL_UPDATE_FREQ = 7; //We will generally align to images with ids dividing 7
    //int alignToId = FAST_FULL_UPDATE_FREQ * (latestImageId / FAST_FULL_UPDATE_FREQ);

	int alignToId = fullFrameUpdateFreq_ * (latestImageId / fullFrameUpdateFreq_);

    //while(permittedAlignToIds.find(alignToId) == permittedAlignToIds.end() && alignToId != latestImageId)
    //    alignToId++;
    if(permittedAlignToIds.find(alignToId) == permittedAlignToIds.end())
        alignToId = getAlignToId(latestImageId, permittedAlignToIds);

    return alignToId;
}

TransformInfo2 * TransformEngine::tryGetTransformCopy(int id, int lastId)
{
    //If we don't BA we're doing some unnecessary copying of correspondences here, doesn't really matter though
    TransformInfo2 * newTransInfo = 0;
    const TransformInfo2 * transBackToHere = transformEstimator_.getTransform(id, lastId);
    if(transBackToHere)
    {
        newTransInfo = new TransformInfo2(*transBackToHere);

        if(!newTransInfo)
            throw new GRCException("TransformEngine::tryGetTransformCopy: TransformInfo2 copy failed");
    }
    else
    {
    	cout << "No transform found between " << id << " and " << lastId << endl;
    }
    return newTransInfo;
}

//Every time an id is added link to one already included
void TransformEngine::addAllIds(TransformSet * TS, std::set<size_t> & permittedAlignToIds, const int latestImageId, TIntSetDesc & candidateAlignToIds, const int DEPTH)
{
	TCloseIdMap::const_iterator pIds = aBestBoWMatches.find(latestImageId);
	if(pIds != aBestBoWMatches.end())
	{
		const TCloseIdSet & aBestMatchingIds = pIds->second;
		for(TCloseIdSet::const_iterator pCloseMatch = aBestMatchingIds.begin(); pCloseMatch != aBestMatchingIds.end(); pCloseMatch++)
		{
			int nNewId = *pCloseMatch;
			if(permittedAlignToIds.find(nNewId) == permittedAlignToIds.end())
			{
	            TransformInfo2 * transBackToHere = tryGetTransformCopy(nNewId, latestImageId);
	            if(transBackToHere)
	            {
	                IdPair idPair(nNewId, latestImageId);
	                TS->init(idPair, transBackToHere);

					permittedAlignToIds.insert(nNewId);
					if(DEPTH)
						addAllIds(TS, permittedAlignToIds, *pCloseMatch, candidateAlignToIds, DEPTH-1);
	            }
			}
		}
	}
}

//! Choose pairs of images that we're likely to get pairs of transforms between, and are likely to fall within the rendered window
void TransformEngine::chooseTransforms(TransformSet * TS, int latestImageId, std::set<size_t> & permittedAlignToIds)
{
    const int NUM_SEQUENTIAL_TRANSFORMS = (maxLMIterations_ > 0) ? 3 : 1; // Try to get transforms with 3 previous frames if doing LM refinement
    const int MAX_SEQ_DIST = 3;

    //Choose twice as many if we're going to cull some:
    int tryMaxFrames = maxFrames_;
    if(frameChoiceMethod_==Enums::eChooseSequentialSkip)
        tryMaxFrames *= 3;

    TIntSetDesc candidateAlignToIds; //Actually unused now essentially

    permittedAlignToIds.insert(latestImageId);//We can align to this one!
     //Very simple: look at last n images...
    int lastId = latestImageId;//This ensures we get a continuous link, even if there are breaks

    if(frameChoiceMethod_==Enums::eChooseSequentialSkip || frameChoiceMethod_==Enums::eChooseSequential)
    {
		//For now just select last 30. Select backwards from here rather than forwards from 0 (so scale/skew doesn't accumulate)
		int selectStart = latestImageId > maxSearchForTransform_ ? latestImageId - maxSearchForTransform_ : 0;

		for(int id = latestImageId - 1; id >= selectStart; id--)
			candidateAlignToIds.insert(id);

	    for(TIntSetDesc::const_iterator pid = candidateAlignToIds.begin(); pid != candidateAlignToIds.end(); pid++)
	    {
	    	const int id = *pid;
	        bool foundOne = false;
	        for(int linkToId = lastId; linkToId <= min<int>(lastId+NUM_SEQUENTIAL_TRANSFORMS-1, latestImageId)
	            && ((frameChoiceMethod_==Enums::eChooseBoWRecurse) || (id>linkToId-MAX_SEQ_DIST)); //Don't try linking to a frame this far away
	            linkToId++)
	        {
	            TransformInfo2 * transBackToHere = tryGetTransformCopy(id, linkToId);
	            if(transBackToHere)
	            {
	                IdPair idPair(id, linkToId);
	                TS->init(idPair, transBackToHere);
	                foundOne = true;

	                if(permittedAlignToIds.find(linkToId) != permittedAlignToIds.end())
	                    permittedAlignToIds.insert(id); //Not necessarily a sequence of transforms up to the latest one now
	            }
	        }
	        if(foundOne)
	            lastId = id;
	    }
    }
    else if(frameChoiceMethod_==Enums::eChooseBoWRecurse)
    {
    	const int BF = 4, DEPTH = 6;

    	int nLastId = -1;
    	if(aBestBoWMatches.size() > 0) nLastId = aBestBoWMatches.backKey();

    	TCloseIdSet & aBestMatchingIds = aBestBoWMatches.init(latestImageId);

    	if(nLastId >= 0)
    		aBestMatchingIds.push_back(nLastId);

    	transformEstimator_.getSimilarIds(latestImageId, BF, aBestMatchingIds);
    	//Make matches bidirectional, for extra frames on LC
    	for(TCloseIdSet::const_iterator pCloseMatch = aBestMatchingIds.begin(); pCloseMatch != aBestMatchingIds.end(); pCloseMatch++)
    		aBestBoWMatches.initOrGet(*pCloseMatch).push_back(latestImageId);

    	addAllIds(TS, permittedAlignToIds, latestImageId, candidateAlignToIds, DEPTH);
    }

}

//! Discard transforms that are heavily overlapping neighbours
//These should be transforms to alignToId for all images, but we haven't refined yet so we want to maintain links between transforms
void TransformEngine::thinTransforms(TransformSet * TS, int alignToId, int latestId) const
{

    typedef std::map<size_t, std::vector<IdPair> > TTransformsById;

    TTransformsById transMap;

    for(TransformSet::const_iterator ppTrans = TS->begin(); ppTrans != TS->end(); ppTrans++)
    {
        int id1=ppTrans->first.im1Id();
        int id2=ppTrans->first.im2Id();

        /*if(transMap.find(id1) == transMap.end())
            transMap[id1];
        if(transMap.find(id2) == transMap.end())
            transMap[id2];
            */

        transMap[id1].push_back(ppTrans->first);

        if(id1 != id2)
            transMap[id2].push_back(ppTrans->first);
    }
    //Now for each image have a vector of transforms. It's a map so is sorted
    typedef std::set<size_t> TIds;
    TIds removeIds;

    TTransformsById::iterator pPreviousTrans = transMap.begin();
    if(pPreviousTrans == transMap.end()) return;

    TTransformsById::iterator pThisTrans = pPreviousTrans; pThisTrans++;
    if(pThisTrans == transMap.end()) return;

    TTransformsById::iterator pNextTrans = pThisTrans; pNextTrans++;
    for(; pNextTrans != transMap.end(); /*previousTrans++ not always*/ pThisTrans++, pNextTrans++)
    {
        //First use the transforms to alignTo to work out if it should be removed
        size_t idPrev = pPreviousTrans->first;
        size_t id = pThisTrans->first;
        size_t idNext = pNextTrans->first;

        if((int)id != alignToId || !doFastRender_) //Do NOT remove this frame during fast render cos we use the translation to work out how to move the entire mosaic next time
        {
            IdPair prevAlignPair(idPrev, alignToId), thisAlignPair(id, alignToId), nextAlignPair(idNext, alignToId);
            const Transform * prevTrans = (*TS)[prevAlignPair]->transform();
            const Transform * thisTrans = (*TS)[thisAlignPair]->transform();
            const Transform * nextTrans = (*TS)[nextAlignPair]->transform();

            double TPrevNext = SSD(prevTrans->translation(), nextTrans->translation());

            const double MIN_TRANSLATION_PX = 0.25 * imageSize_.height;

            if(TPrevNext < SQR(MIN_TRANSLATION_PX)) //1) These frames are closer together in the mosaic than is strictly necessary
            {

                double TPrevThis = SSD(prevTrans->translation(), thisTrans->translation());
                double TThisNext = SSD(thisTrans->translation(), nextTrans->translation());

                if(TPrevThis < TPrevNext && TThisNext < TPrevNext) // 2) The middle frame is actually between the previous and next frame
                {
                    removeIds.insert(id);//affects increments now

                    //Look thru the transforms for this image--make sure there's a prev-next transform for refinement
                    //We're throwing away some info here that could be used for BA

                    //cout << "Erasing " << id << endl;

                    if((int)id == alignToId) //we're deleting the frame we're aligning to--this is ok but don't want to delete transforms taking other images here
                    {
                        IdPair alignIds(alignToId, alignToId);
                        TS->erase(alignIds);
                    }
                    else
                    {
                        //If these exist then join them and add to transform set (make use of transformation info when we LM refine)
                        IdPair prevThisPair(idPrev, id), prevNextPair(idPrev, idNext), thisNextPair(id, idNext);

                        if(!TS->exists(prevNextPair)) //Commenting out this line improves performance slightly, but is kind of hacky. Should allow both prevNextPair transforms to be used.
                        {
                            if(TS->exists(prevThisPair) && TS->exists(thisNextPair))
                            {
                                const TransformInfo2 * prevThisTrans = (*TS)[prevThisPair];
                                const TransformInfo2 * thisNextTrans = (*TS)[thisNextPair];
                                Transform * pNewTrans = prevThisTrans->transform()->accumulate(thisNextTrans->transform());

                                //Choose the correspondence set from the transform with the least--the number of correspondences here is a good indicator of the accumulated transform's quality. THIS IS A HACK though--for full BA we must accumulate the correspondence sets and find tracks.
                                const CBoWCorrespondences & corrToUseHack = thisNextTrans->correspondences().size() <  prevThisTrans->correspondences().size() ? thisNextTrans->correspondences() : prevThisTrans->correspondences();
                                TransformInfo2 * prevNextTrans = new TransformInfo2( pNewTrans, &corrToUseHack, idPrev, idNext);

                                (*TS)[prevNextPair] = prevNextTrans;
                                transMap[idPrev].push_back(prevNextPair);
                                transMap[idNext].push_back(prevNextPair);
                                //cout << "Adding new transform: " << idPrev << ',' << idNext << endl;
                            }
                        } //TODO We're still throwing away some info here (sometimes) that could be used for BA

                        std::vector<IdPair> & transformToThisIm(transMap[id]);
                        for(std::vector<IdPair>::iterator pPair = transformToThisIm.begin(); pPair != transformToThisIm.end(); pPair++)
                        {
                            TS->erase(*pPair);
                        }
                    }
                }
            }
        }

        //May have deleted several ids--make sure previousTrans points to the last undeleted one
        pPreviousTrans = pThisTrans;
        while(removeIds.find(pPreviousTrans->first) != removeIds.end())
            pPreviousTrans--;

    }
}

//! Calculate transform set (image+transformations) to be rendered, including images up to latestImageId.
TransformSet * TransformEngine::computeTransforms(int latestImageId)
{
    TransformSet * TS = new TransformSet;

    Transform * idTrans = Transform::newTransform();

    //Choose and calculate some transformations
    std::set<size_t> permittedAlignToIds;
    chooseTransforms(TS, latestImageId, permittedAlignToIds);

    int alignToId = NO_FRAME;
	if(doFastRender_){
        alignToId = getFastAlignToId(latestImageId, permittedAlignToIds);
	} else {
        alignToId = getAlignToId(latestImageId, permittedAlignToIds);
	}
    IdPair latestPair(alignToId, alignToId);
    TS->init(latestPair, new TransformInfo2(idTrans, alignToId, alignToId));

    bool doLMRefinement = (maxLMIterations_ > 0); // Todo: decide (maybe only when there's a loop?)

    if(doLMRefinement)
    {
        //This will recurse and add a transform x->alignToId for every x where this path exists
        alignTransformToLastImage(TS, alignToId, alignToId, true);

        // Remove Tij where Ti0 not in TS
        removeDisconnectedTransforms(TS, alignToId);

        //if(frameChoiceMethod_==Enums::eChooseSequentialSkip)
        //    thinTransforms(TS, alignToId, latestImageId); //Sparsify transform set

        //Refine transforms by LM

        refineTranformSet(TS, alignToId, latestImageId);

        //Delete all transforms not pointing to latest
        cleanupTransformSet(TS, alignToId);
    }
    /*else if(doBundleAdjust)
    {

        //Refine with BA
        //refineTranformSetAndPoints(TS);

        //Change each transform to one taking the intermediate image to the last image
        alignTransformsToOneImage(TS, alignToId);
    }*/
    else
    {
        alignTransformsToOneImage(TS, alignToId);

    }

    if(frameChoiceMethod_==Enums::eChooseSequentialSkip)
        thinTransforms(TS, alignToId, latestImageId); //Sparsify transform set

    positionMosaic(TS, alignToId, latestImageId);

    if(doFastRender_ )
        setFastRender(TS, alignToId, latestImageId);
    /*else
    {
        IdPair idLatest(latestImageId, alignToId);
    	TS->setLatestTransform(idLatest);
    }*/

    return TS;
}

void TransformEngine::setFastRender(TransformSet * TS, int alignToId, int latestId) const
{
    static int lastAlignToId = NO_FRAME;
    static Transform * transToAlignToIdInLastRenderedFrame = 0;

    IdPair idAC(alignToId, alignToId);
    const Transform * transAC = (*TS)[idAC]->transform();//The one that should (will!) always exist
    if(alignToId == lastAlignToId)
    {
        //Last time mosaic was in frame b, transToAlignToIdInLastRenderedFrame takes a to b
        //Now it is in frame c, we've calculated a to c--it is (*TS)[ids(a,c)]. Want to shift mosaic b to c so want

        Transform * pGlobalMosaicTrans = 0;
        if(transToAlignToIdInLastRenderedFrame)
        {
            static Transform * IdentityTrans =Transform::newTransform();

            Transform * transBA = IdentityTrans->accumulateInverse(transToAlignToIdInLastRenderedFrame);

            pGlobalMosaicTrans = transBA->accumulate(transAC); //error if doesn't exist

            delete transBA;
        }

        IdPair idLatest(latestId, alignToId);
        TS->setFastRender(pGlobalMosaicTrans, idLatest);
    }

    delete transToAlignToIdInLastRenderedFrame;
    transToAlignToIdInLastRenderedFrame = Transform::copyTransform(transAC);//Will be transAB next time

    lastAlignToId = alignToId;
}


}

pragma_warning(pop)  // warning C4996: 'std::_Copy_backward_opt' was declared deprecated (used in boost 1_35_0)
