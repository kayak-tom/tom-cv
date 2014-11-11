/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * surfFeature.h
 *
 *  Created on: 14/10/2009
 *      Author: tom
 */

#ifndef SURFFEATURE_H_
#define SURFFEATURE_H_

class sSURFFeature
{
public:
	sSURFFeature() : x(0), y(0)
	{
		for(int i=0;i<128;i++)
			adDescriptor[i] = 0;
	}
    double x,y;
    double adDescriptor[128];
    typedef double xytype;
};

/*class sSIFTFeature
{
public:
	sSIFTFeature() : x(0), y(0)
	{
		for(int i=0;i<128;i++)
			afDescriptor[i] = 0;
	}
    double x,y;
    float afDescriptor[128];
    typedef float xytype;
};*/
#endif /* SURFFEATURE_H_ */
