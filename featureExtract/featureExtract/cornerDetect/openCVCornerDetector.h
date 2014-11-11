/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * openCVCornerDetector.h
 *
 *  Created on: 19/05/2009
 *      Author: tom
 */

#ifndef OPENCVCORNERDETECTOR_H_
#define OPENCVCORNERDETECTOR_H_

#include "goodFeaturesFast.h"
#include "../cornerDetector.h"

class COpenCVCornerDetector: public CCornerDetector
{
	const double MARGIN;
	IplImage * eigImg;
	IplImage * tempImg;
    CvPoint2D64f * aCvPointCorners;
    CPointBin<32, 32, 16> pointBin;
    const CCornerParams::CCornerDetectorParams & CORNER_PARAMS;
public:
	COpenCVCornerDetector(const CImParams & IM_PARAMS, int MARGIN, const CCornerParams &);
	virtual ~COpenCVCornerDetector();

	virtual void getCorners(IplImage * pImage, int & nCorners, CLocation * aCorners);
};

#endif /* OPENCVCORNERDETECTOR_H_ */
