/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

#include "descriptor.h"

#include <vector>
#include <cstdlib>			// C standard includes
#include <iostream>			// C++ I/O
#include <string>			// C++ strings
#include "util/exception.h"
#include "util/opencv.h"

#define CTVectorAndHistDescriptor CVectorAndHistDescriptor<elementType, nDescriptorLength, eNormToUse>

//////// Metric Spaces /////////////////////
templateVS
class CVectorAndHistDescriptor : public CTVectorSpaceDescriptor
{
public:
	CVectorAndHistDescriptor(CvHistogram * pHist, const elementType * aDescriptor);

	~CVectorAndHistDescriptor();

    double distance(const CDescriptor * pd) const; //still virtual

	inline const CvHistogram * Hist() const { return pHist; };

	virtual int size() const { return sizeof(this); }
};

templateVS
CTVectorAndHistDescriptor::CVectorAndHistDescriptor(CvHistogram * pHist_in, const elementType * aDescriptor) : CTVectorSpaceDescriptor(const elementType * aDescriptor)
{
	if(IS_DEBUG) CHECK(!pHist, "CTVectorAndHistDescriptor::CVectorAndHistDescriptor: Bad params");
	pHist = pHist_in;
}

templateVS
CTVectorAndHistDescriptor::~CVectorAndHistDescriptor()
{
	CvReleaseHist(pHist);
}

templateVS
double CTVectorAndHistDescriptor::distance(const CDescriptor * pd) const
{
	CTVectorAndHistDescriptor * pDesc = CAST<const CTVectorAndHistDescriptor *>(pd);
	double dHistDist = CvCompareHist(pDesc->Hist, pHist, COMP_METHOD);
	return dHistDist + CTVectorSpaceDescriptor::distance(pd);
}
