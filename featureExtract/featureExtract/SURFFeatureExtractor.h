/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * CSURFFeatureExtractor.h
 *
 *  Created on: 29/06/2009
 *      Author: tom
 */

#ifndef CSURFFEATUREEXTRACTOR_H_
#define CSURFFEATUREEXTRACTOR_H_

#include "featureExtractor.h"
#include "params/param.h"

class CSURFFeatureExtractor: public CFeatureExtractor
{
	const CCornerParams::CSURFParams & SURFParams;
protected:
	virtual CDescriptorSet * getDescriptors_int(const IplImage * pImage);
	CvMemStorage* pCvMemStorage;
public:
	CSURFFeatureExtractor(const int MAX_FEATURES, const CImParams & IM_PARAMS, const CDescriptorSetClusteringParams & DSC_PARAMS, const CCornerParams::CSURFParams & SURFParams);
	virtual ~CSURFFeatureExtractor();
};

#endif /* CSURFFEATUREEXTRACTOR_H_ */
