/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

/*
 * findBestMatching.h
 *
 *  Created on: 8/10/2009
 *      Author: tom
 */

#ifndef FINDBESTMATCHING_H_
#define FINDBESTMATCHING_H_

class CMask;

class CFindBestMatching
{
public:
	virtual bool supplyResiduals() const = 0;
	virtual void refine(CMask & bInliers, int & nInliers, const double * adResiduals) = 0;
	virtual ~CFindBestMatching() {};
};

//Use when datapoints are independent (i.e. 1-1 correspondences rather than N--M correspondences)
class CNoMatchRefinement : public CFindBestMatching
{
public:
	virtual bool supplyResiduals() const { return false; }
	virtual void refine(CMask & , int & , const double * ) {};
};


#endif /* FINDBESTMATCHING_H_ */
