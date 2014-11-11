/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * PatchFeatureExtractor.h
 *
 *  Created on: 29/06/2009
 *      Author: tom
 */

#ifndef PATCHFEATUREEXTRACTOR_H_
#define PATCHFEATUREEXTRACTOR_H_

#include "featureExtractor.h"

class CPatchFeatureExtractor: public CFeatureExtractor
{
protected:
	CCornerDetector * pCornerDetector;
	CLocation * aCorners;
	const CPatchDescriptorParams & PATCH_PARAMS;

	virtual CDescriptorSet * getDescriptors_int(const IplImage * pImage);
public:
	CPatchFeatureExtractor(const int MAX_FEATURES, const CImParams & IM_PARAMS, CCornerDetector ** ppCornerDetector, const CPatchDescriptorParams & PATCH_PARAMS, const CDescriptorSetClusteringParams & DSC_PARAMS);
	virtual ~CPatchFeatureExtractor();
};

class CDSPatchFeatureExtractor: public CPatchFeatureExtractor
{
	IplImage * pDownsampledImage;
	virtual CDescriptorSet * getDescriptors_int(const IplImage * pImage);
	const int nScale;
public:
	CDSPatchFeatureExtractor(const int MAX_FEATURES, const CImParams & IM_PARAMS, CCornerDetector ** ppCornerDetector, CPatchDescriptorParams & PATCH_PARAMS, const CDescriptorSetClusteringParams & DSC_PARAMS);
	virtual ~CDSPatchFeatureExtractor();
};

#endif /* PATCHFEATUREEXTRACTOR_H_ */
