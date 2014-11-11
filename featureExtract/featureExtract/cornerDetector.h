/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * cornerDetector.h
 *
 *  Created on: 19/05/2009
 *      Author: tom
 */

#ifndef CORNERDETECTOR_H_
#define CORNERDETECTOR_H_

#include "image/imageAccess.h"
#include "util/opencv.h"
#include "util/location.h"
#include "params/param.h"

class CCornerDetector
{
protected:
	const CImParams & IM_PARAMS;
public:
	CCornerDetector(const CImParams & IM_PARAMS) : IM_PARAMS(IM_PARAMS) {}
	virtual void getCorners(IplImage * pImage, int & nCorners, CLocation * aCorners) = 0;
	virtual ~CCornerDetector() {};
};

struct CvSURFParams;

PARAMCLASS(Corner)
	PARAME(SALIENT_FEATURE_TYPE, FastCorners, "Corner/blob type. FAST appear to be both faster and most repeatable")
	PARAM(MAX_FEATURES, 10, 100000, 400, "Max number of features detected per image. Usually reached, with top scoring corners (with min seperation) selected")
	CHILDCLASS(CornerDetector, "Parameters for FAST and Harris-like corners")
	CHILDCLASS(SURF, "SURF blob detector (from OpenCV). Not extensively tested.")
	{}

	PARAMCLASS(CornerDetector)
		PARAME(CORNER_MODE, OpenCVGoodFeatures, "Only OpenCVGoodFeatures working properly") //Harris set seperately
		PARAM(CORNER_QUAL, 0.004, 0.2, 0.006, "Min Shi-Tomasi corner score (not usually used--top MAX_FEATURES usually selected")
		PARAM(CORNER_MIN_DIST, 2, 25, 4, "Min seperation, in pixels")
		PARAM(CORNER_BLOCK_RAD, 1, 4, 2, "Gaussian blur size")
		PARAM(CORNER_HARRIS_K, 0.003, 0.1, 0.01, "k param for Harris corners")
		PARAMB(USE_SUBPIX, true, "Subpixel refinement. Improves localisation accuracy slightly but not to accuracy of FAST")
		PARAMB(USE_HARRIS, false, "Harris corners rather than OpenCV. Hardly actually faster.")
		PARAM(FASTCORNER_T, 1,255,21, "Min strength of Harris corners (min grey level difference between centre and circle)")
		PARAM(FASTCORNER_MIN_SEPERATION, 1,64,4, "5 works well, repeatability drops if worse corners let through")
		PARAMB(FASTCORNER_TWO_PASS_BINNING, false, "Slightly more repeatable false, true improves geometry")
		PARAM(FASTCORNER_LOCALISATION_SD, 1e-6, 5, 0.6 , "Expected localisation error in pixels. 0.6 derived from test imges")
		{}

		MAKEENUMPARAM3(CORNER_MODE, FasterOpenCVGoodFeatures, SubSampledCorners, OpenCVGoodFeatures);
		CNumParam<double> CORNER_QUAL;
		CNumParam<int> CORNER_MIN_DIST;
		CNumParam<int> CORNER_BLOCK_RAD;
		CNumParam<double> CORNER_HARRIS_K;
		CNumParam<bool> USE_SUBPIX, USE_HARRIS;
		CNumParam<int> FASTCORNER_T, FASTCORNER_MIN_SEPERATION;
		CNumParam<bool> FASTCORNER_TWO_PASS_BINNING;
		CNumParam<double> FASTCORNER_LOCALISATION_SD;
	};


	MAKEENUMPARAM2(SALIENT_FEATURE_TYPE, SURFBlobs /*deprecated , ShiTomasiCorners*/, FastCorners);
	CNumParam<int> MAX_FEATURES;
	MAKECHILDCLASS(CornerDetector);

	PARAMCLASS(SURF)
		PARAMB(EXTENDED, false, "")
		PARAM(HESSIAN_THRESH, 300, 500, 400, "")
		PARAM(OCTAVES, 1, 10, 3, "")
		PARAM(OCTAVE_LAYERS, 1, 10, 4, "")
		PARAM(BLOB_LOCALISATION_SD, 1e-6, 5, 0.64, "SD of error in pixels. 0.64 derived from test images")
		PARAM(CORNER_MIN_DIST, 1, 100, 5, "minimum seperation between blobs. Multiple blobs at same location never allowed")
		{}

		CNumParam<bool> EXTENDED;
		CNumParam<double> HESSIAN_THRESH;
		CNumParam<int> OCTAVES, OCTAVE_LAYERS;
		CNumParam<double> BLOB_LOCALISATION_SD;
		CNumParam<int> CORNER_MIN_DIST;

		/*Possibly should comment out... */void getSurfParams(CvSURFParams & surfParams) const;

		int SURF_LENGTH() const
		{
			int length = 64;
			return EXTENDED ? length*2 : length;
		}
	};
	MAKECHILDCLASS(SURF);

	double CORNER_LOCALISATION_SD() const
	{
		switch(SALIENT_FEATURE_TYPE)
		{
		case eFastCorners:
			return CornerDetector.FASTCORNER_LOCALISATION_SD;
		case eSURFBlobs:
			return SURF.BLOB_LOCALISATION_SD;

		default:
			std::cout << "Warning: CORNER_LOCALISATION_SD not known for this corner type";
			return 0.6;
		}
	}
};

#endif /* CORNERDETECTOR_H_ */
