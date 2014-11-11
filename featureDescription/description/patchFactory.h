/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * patchFactory.h
 *
 *  Created on: 13/10/2009
 *      Author: tom
 */

#ifndef PATCHFACTORY_H_
#define PATCHFACTORY_H_

#include "util/exception.h"
#include "patchDescriptor.h"

#define TDescriptorFactory CDescriptorFactory<unsigned char, unsigned char>
template<typename dataType, typename imDataType>
class CDescriptorFactory
{
	friend void LoadSettingsFromCfg(const char * szFN);
public:
	template<CPatchParams::eColourTypeSize ImType>
	static CDescriptor * makeDescriptor_int(const IplImage * pImage, const CLocation loc, double dScale, const CPatchDescriptorParams & PATCH_PARAMS)
	{
		switch(PATCH_PARAMS.PATCH_COMP_METHOD)
		{
		case CPatchDescriptorParams::ePatchEuclidFast:
			return makeDescriptor_int2<ImType, CPatchDescriptorParams::ePatchEuclidFast>(pImage, loc, dScale, PATCH_PARAMS);
		case CPatchDescriptorParams::ePatchL1Fast:
			return makeDescriptor_int2<ImType, CPatchDescriptorParams::ePatchL1Fast>(pImage, loc, dScale, PATCH_PARAMS);
		case CPatchDescriptorParams::ePatchMaxDist:
			return makeDescriptor_int2<ImType, CPatchDescriptorParams::ePatchMaxDist>(pImage, loc, dScale, PATCH_PARAMS);
		case CPatchDescriptorParams::ePatchCorrel:
			return makeDescriptor_int2<ImType, CPatchDescriptorParams::ePatchCorrel>(pImage, loc, dScale, PATCH_PARAMS);
		case CPatchDescriptorParams::ePatchL1:
			return makeDescriptor_int2<ImType, CPatchDescriptorParams::ePatchL1>(pImage, loc, dScale, PATCH_PARAMS);
		case CPatchDescriptorParams::ePatchEuclid:
			return makeDescriptor_int2<ImType, CPatchDescriptorParams::ePatchEuclid>(pImage, loc, dScale, PATCH_PARAMS);
		case CPatchDescriptorParams::ePatchEuclidParallel:
			return makeDescriptor_int2<ImType, CPatchDescriptorParams::ePatchEuclidParallel>(pImage, loc, dScale, PATCH_PARAMS);
		case CPatchDescriptorParams::ePatchL1Parallel:
			return makeDescriptor_int2<ImType, CPatchDescriptorParams::ePatchL1Parallel>(pImage, loc, dScale, PATCH_PARAMS);
		}
		THROW( "makeDescriptor_int: Norm not handled")
	}

	template<CPatchParams::eColourTypeSize ImType, CPatchDescriptorParams::ePATCH_COMP_METHOD PATCH_COMP_METHOD>
	static CDescriptor * makeDescriptor_int2(const IplImage * pImage, const CLocation loc, double dScale, const CPatchDescriptorParams & PATCH_PARAMS)
	{
		CDescriptor * pDescriptor = 0;
		switch(PATCH_PARAMS.PATCH_RAD)
		{
		case 1:
			pDescriptor = new CPatchAndLocationDescriptor<dataType, 1, PATCH_COMP_METHOD, ImType, imDataType>(dScale, pImage, loc, PATCH_PARAMS.Patch);
			break;
		case 2:
			pDescriptor = new CPatchAndLocationDescriptor<dataType, 2, PATCH_COMP_METHOD, ImType, imDataType>(dScale, pImage, loc, PATCH_PARAMS.Patch);
			break;
		case 3:
			pDescriptor = new CPatchAndLocationDescriptor<dataType, 3, PATCH_COMP_METHOD, ImType, imDataType>(dScale, pImage, loc, PATCH_PARAMS.Patch);
			break;
		case 4:
			pDescriptor = new CPatchAndLocationDescriptor<dataType, 4, PATCH_COMP_METHOD, ImType, imDataType>(dScale, pImage, loc, PATCH_PARAMS.Patch);
			break;
		case 5:
			pDescriptor = new CPatchAndLocationDescriptor<dataType, 5, PATCH_COMP_METHOD, ImType, imDataType>(dScale, pImage, loc, PATCH_PARAMS.Patch);
			break;
		case 6:
			pDescriptor = new CPatchAndLocationDescriptor<dataType, 6, PATCH_COMP_METHOD, ImType, imDataType>(dScale, pImage, loc, PATCH_PARAMS.Patch);
			break;
		default:
			THROW( "Bad PATCH_RAD")
		}
		return pDescriptor;
	}

	static CDescriptor * makeDescriptor(const IplImage * pImage, const CLocation loc, double dScale, const CPatchDescriptorParams & PATCH_PARAMS)
	{
		if(pImage->nChannels == 3)
		{
			return makeDescriptor_int<CPatchParams::RG>(pImage, loc, dScale, PATCH_PARAMS);
		}
		else
		{
			return makeDescriptor_int<CPatchParams::GREY>(pImage, loc, dScale, PATCH_PARAMS);
		}
	}
};

#endif /* PATCHFACTORY_H_ */
