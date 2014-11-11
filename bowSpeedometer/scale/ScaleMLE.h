/*
 * ScaleMLE.h
 *
 *  Created on: 23 Apr 2010
 *      Author: tom
 */

#ifndef SCALEMLE_H_
#define SCALEMLE_H_

#include "BoWSpeedo.h"

class CScaleMLE {
	static bool getMLEScaleParams_2ndOrderGD_int(const CBoWSpeedo::CBoWObject::TLengths & aLengths, double & dMean, double & dVar, const bool bVerbose);
public:
	CScaleMLE();
	//virtual ~CScaleMLE();

	static bool getBayesianScaleParams_NG(const CBoWSpeedo::CBoWObject::TLengths & aLengths, double & dMean, double & dVar, const bool bVerbose);
	static bool getMLEScaleParams(const CBoWSpeedo::CBoWObject::TLengths & aLengths, double & dMean, double & dVar, const bool bVerbose);
	static bool getMLEScaleParams_GD(const CBoWSpeedo::CBoWObject::TLengths & aLengths, double & dMean, double & dVar, const bool bVerbose);
	static bool getMLEScaleParams_2ndOrderGD(const CBoWSpeedo::CBoWObject::TLengths & aLengths, double & dMean, double & dVar, const bool bVerbose);

	static void testScaleMLE();
};

#endif /* SCALEMLE_H_ */
