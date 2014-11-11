/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * vectorDescriptorFactory.h
 *
 *  Created on: 14/10/2009
 *      Author: tom
 */

#ifndef VECTORDESCRIPTORFACTORY_H_
#define VECTORDESCRIPTORFACTORY_H_

#include "vectorDescriptor.h"
//typedef CVectorSpaceDescriptor<char, 64, eCosine> TSURFDescriptor;
//typedef CVectorSpaceDescriptor<unsigned char, 128, eCosine> TSIFTDescriptor;

typedef CVectorLocationDescriptor<char, 64, eEuclidSquared> TVectorDescriptor64;
typedef CVectorLocationDescriptor<char, 128, eEuclidSquared> TVectorDescriptor128;

typedef CVectorLocationOrientationDescriptor<char, 64, eEuclidSquared> TVectorOrientationDescriptor64;
typedef CVectorLocationOrientationDescriptor<char, 128, eEuclidSquared> TVectorOrientationDescriptor128;

class CVectorDescriptorFactory
{
public:
	static CDescriptor * makeDescriptor(const float * pfDesc, int size, const CLocation loc)
	{
		if(size==64)
			return new TVectorDescriptor64(pfDesc, 127.99f, loc);
		else if(size==128)
			return new TVectorDescriptor128(pfDesc, 127.99f, loc);
		else
			THROW( "Size not supported")
	}

	static CDescriptor * makeDescriptor(const float * pfDesc, int size, double orientation, const CLocation loc)
	{
		if(size==64)
			return new TVectorOrientationDescriptor64(pfDesc, 127.99f, orientation, loc);
		else if(size==128)
			return new TVectorOrientationDescriptor128(pfDesc, 127.99f, orientation, loc);
		else
			THROW( "Size not supported")
	}
};

#endif /* VECTORDESCRIPTORFACTORY_H_ */
