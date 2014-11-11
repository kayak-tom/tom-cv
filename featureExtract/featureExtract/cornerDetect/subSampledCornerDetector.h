/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * openCVCornerDetector.h
 *
 *  Created on: 19/05/2009
 *      Author: tom
 */

#ifndef SSCORNERDETECTOR_H_
#define SSCORNERDETECTOR_H_

#include "goodFeaturesFast.h"
#include "../cornerDetector.h"

class CSubSampledCornerDetector: public CCornerDetector
{
	const double MARGIN;
	const int MARGIN_INT;
    IplImage * eigImg;
    IplImage * tempImg;
    CvPoint2D64f * aCvPointCorners;
    IplImage * pGreySSImage;
    CPointBin<32, 32, 16> pointBin;
    const CCornerParams::CCornerDetectorParams & CORNER_PARAMS;
public:
	CSubSampledCornerDetector(const CImParams & IM_PARAMS, int MARGIN, const CCornerParams &);
	virtual ~CSubSampledCornerDetector();

	virtual void getCorners(IplImage * pImage, int & nCorners, CLocation * aCorners);
};

#endif /* SSCORNERDETECTOR_H_ */
