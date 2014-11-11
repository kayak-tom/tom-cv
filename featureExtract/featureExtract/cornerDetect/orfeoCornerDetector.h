/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * openCVCornerDetector.h
 *
 *  Created on: 19/05/2009
 *      Author: tom
 */

#ifndef ORFEOCORNERDETECTOR_H_
#define ORFEOCORNERDETECTOR_H_

//#include "goodFeaturesFast.h"
#include "../cornerDetector.h"

class COpenCVGoodFeaturesCornerDetector: public CCornerDetector
{
	const double MARGIN;
    IplImage * eigImg;
    IplImage * tempImg;
    CvPoint2D32f * aCvPointCorners;
    const CCornerParams::CCornerDetectorParams & CORNER_PARAMS;
public:
	COpenCVGoodFeaturesCornerDetector(const CImParams & IM_PARAMS, int MARGIN, const CCornerParams &);
	virtual ~COpenCVGoodFeaturesCornerDetector();

	virtual void getCorners(IplImage * pImage, int & nCorners, CLocation * aCorners);
};

#endif /* ORFEOCORNERDETECTOR_H_ */
