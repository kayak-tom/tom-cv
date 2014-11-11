/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * CImageSimulator.h
 *
 *  Created on: 10/05/2009
 *      Author: tom
 */

#ifndef CIMAGESIMULATOR_H_
#define CIMAGESIMULATOR_H_

//#include "../StereoNav/odometry_lib.h"
#include "description/descriptor.h"
#include <vector>
#include "bowslamParams.h"

class CCamCalibMatrix;

class CImageSimulator
{
	std::vector<CDescriptorSet *> aSimPhotos;
public:
	CImageSimulator(const CCamCalibMatrix & K, const CDescriptorSetClusteringParams & DSC_PARAMS);
	virtual ~CImageSimulator();

	CDescriptorSet * getSimDescriptors(int id);
};

#endif /* CIMAGESIMULATOR_H_ */
