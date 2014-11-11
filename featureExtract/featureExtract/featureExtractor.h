/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * featureExtractor.h
 *
 *  Created on: 19/05/2009
 *      Author: tom
 */

#ifndef FEATUREEXTRACTOR_H_
#define FEATUREEXTRACTOR_H_

#include "cornerDetector.h"
#include "description/descriptor.h"
#include "image/imageAccess.h"
#include "image/convert_OpenCV.h"
#include "description/patchParams.h"

class CPatchDescriptorParams;
class CDescriptorSetClusteringParams;

class CFeatureExtractor
{
	const CGreyscaler greyScaler;
	const CDescriptorSetClusteringParams & DSC_PARAMS;
protected:
	virtual CDescriptorSet * newDescriptorSet(int nDescriptorEstimate);
	const int MAX_FEATURES;
	const CImParams & IM_PARAMS;
	IplImage * pGreyImg;
	virtual CDescriptorSet * getDescriptors_int(const IplImage * pImage) = 0;

public:
	CFeatureExtractor(const int MAX_FEATURES, const CImParams & IM_PARAMS, const CDescriptorSetClusteringParams & DSC_PARAMS);
	virtual ~CFeatureExtractor();
	CDescriptorSet * getDescriptors(const IplImage * pImage);

	static CFeatureExtractor * makeFeatureExtractor(const CImParams & IMPARAMS, const CCornerParams & CORNERPARAMS,
			const CPatchDescriptorParams & PATCHDESCRIPTORPARAMS, const CDescriptorSetClusteringParams & DSCPARAMS);
};

void markDescriptors(IplImage * pImage, const CDescriptorSet * pDesc);

#endif /* FEATUREEXTRACTOR_H_ */
