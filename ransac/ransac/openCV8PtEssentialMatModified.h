/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

/*
 * openCV8PtEssentialMatModified.h
 *
 *  Created on: 8/10/2009
 *      Author: tom
 */

#ifndef OPENCV8PTESSENTIALMATMODIFIED_H_
#define OPENCV8PTESSENTIALMATMODIFIED_H_

#include "refiner.h"

class COpenCV8PtEssentialMatModified: public CImCorrModelRefiner
{
	const double d8ptCutoff;
	const bool bFindF;
	virtual bool fitModel_int(CMask & mask, CModel & pModel, bool bVerbose);
public:
	COpenCV8PtEssentialMatModified(const T2dPoints & p1, const T2dPoints & p2, double d8ptCutoff, bool bFindF);
	virtual ~COpenCV8PtEssentialMatModified();
};

#endif /* OPENCV8PTESSENTIALMATMODIFIED_H_ */
