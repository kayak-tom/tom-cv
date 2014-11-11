#pragma once
#ifndef _CORRESPONDENCE_SET
#define _CORRESPONDENCE_SET

#include <map>
#include <vector>
#include "util/opencv.h"
#include "GRCException.h"

namespace grc {

//! A set of correspondences (x1,y1)-(x2,y2) between 2 images
class CorrespondenceSet
{
    //! For sorting CvPoint2D64f's.
	class CvPoint2D64fCompare : std::binary_function<const CvPoint2D64f &, const CvPoint2D64f &, bool> {
	public:
		bool operator()(const CvPoint2D64f & p1, const CvPoint2D64f & p2) const {
			if (p1.x != p2.x)
				return (p1.x < p2.x);
			else
				return (p1.y < p2.y);
		}
	};

    typedef std::pair<CvPoint2D64f, CvPoint2D64f> TCorrespondencePair ;
    class CorrespondencePair : private TCorrespondencePair //!< Wrapper around std::pair, with nicer names
    {
        public:
            CorrespondencePair(const CvPoint2D64f& p1, const CvPoint2D64f& p2) : TCorrespondencePair(p1, p2) {};
            CvPoint2D64f pointInIm1() const { return first; };
            CvPoint2D64f pointInIm2() const { return second; };
    };

    typedef std::map<CvPoint2D64f, CvPoint2D64f, CvPoint2D64fCompare> TCorrespondenceMap; //!<Map mapping points in one image to matching point in the other.
    typedef std::vector<CvPoint2D64f> TCorrespondenceVector ; //!< For random access--memory cost is at most 1/4 the size of the maps

    TCorrespondenceMap corr1to2_, corr2to1_; //!< Maps mapping points in each image to the corresponding point in the other. 2 for fast searching either way. 
    TCorrespondenceVector randAccessVec_;//!< Holds list of points in image 1 (for random access for RANSAC)
public:
    typedef TCorrespondenceMap::const_iterator const_iterator;

    const_iterator begin() const { return corr1to2_.begin(); };
    const_iterator end() const { return corr1to2_.end(); };

    //! Random access to correspondences
    CorrespondencePair operator[](int i) const
    {
        CvPoint2D64f p1 = randAccessVec_[i];
        CorrespondencePair pair(p1, im2Point(p1));
        return pair;
    };

    //! Number of correspondences
    size_t size() const { return randAccessVec_.size(); };

    //! Insert correspondence
    void insertCorresp(CvPoint2D64f p1, CvPoint2D64f p2)
    {
        corr1to2_[p1] = p2;
        corr2to1_[p2] = p1;
        randAccessVec_.push_back(p1);

        if(corr1to2_.size() != corr2to1_.size() || corr1to2_.size() != randAccessVec_.size())
            throw new GRCException("CorrespondenceSet::insertCorresp: Sizes out-of-sync--duplicate inserted?");
    };

    //! Is this point already matched to one in image 1?
    bool point2Exists(CvPoint2D64f p2) const
    {
    	return corr2to1_.find(p2) != corr2to1_.end();
    };

    //! Is this point already matched to one in image 2?
    bool point1Exists(CvPoint2D64f p1) const
    {
    	return corr1to2_.find(p1) != corr1to2_.end();
    };

    //! Return point in image 1 matching point in image 2 (or NO_CORRESPONDENCE if there isn't a match)
    CvPoint2D64f im1Point(CvPoint2D64f p2) const
    {
    	TCorrespondenceMap::const_iterator loc = corr2to1_.find(p2);

    	if(loc == corr2to1_.end())
        	throw new GRCException("CorrespondenceSet::im1Point: Correspondence doesn't exist");

    	return loc->second;
    };

    //! Return point in image 2 matching point in image 1 (or NO_CORRESPONDENCE if there isn't a match)
    CvPoint2D64f im2Point(CvPoint2D64f p1) const
    {
    	TCorrespondenceMap::const_iterator loc = corr1to2_.find(p1);

    	if(loc == corr1to2_.end())
        	throw new GRCException("CorrespondenceSet::im2Point: Correspondence doesn't exist");

    	return loc->second;
    };
};

}


#endif // _CORRESPONDENCE_SET
