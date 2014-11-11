/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * FASTCornerDetector.h
 *
 *  Created on: 17/06/2009
 *      Author: tom
 */

#ifndef FASTCORNERDETECTOR_H_
#define FASTCORNERDETECTOR_H_

#include "../cornerDetector.h"

class CFASTCornerDetector: public CCornerDetector
{
	const CCornerParams::CCornerDetectorParams & CORNER_PARAMS;
	const int MARGIN;
public:
	CFASTCornerDetector(const CImParams& IM_PARAMS, int MARGIN, const CCornerParams::CCornerDetectorParams & CORNER_PARAMS);
	virtual ~CFASTCornerDetector();
	virtual void getCorners(IplImage * pImage, int & nCorners, CLocation * aCorners);
};

#endif /* FASTCORNERDETECTOR_H_ */
