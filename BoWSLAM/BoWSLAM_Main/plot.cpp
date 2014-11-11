/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#include "plot.h"
#include "util/convert.h"
#include "util/opencv.h"
#include "stdio.h"

using namespace std;

typedef std::map<int, double>::const_iterator it;

class CPlotConvert
{
	const double dMinX, dMaxX, dMinY, dMaxY, dWidth, dHeight, dMargin;
	double dScaleX, dScaleY;
public:
	CPlotConvert(double dMinX, double dMaxX, double dMinY, double dMaxY, IplImage * image, double dMargin) : dMinX(dMinX), dMaxX(dMaxX), dMinY(dMinY), dMaxY(dMaxY), dWidth(image->width), dHeight(image->height), dMargin(dMargin)
	{
		dScaleX = (dWidth-2*dMargin)/(dMaxX-dMinX);
		dScaleY = (dWidth-2*dMargin)/(dMaxY-dMinY);

		pair<int, double> pTop(0, dMaxY);
		pair<int, double> pBottom(0, dMinY);
		pair<int, double> pLeft(dMinX, 0);
		pair<int, double> pRight(dMaxX, 0);
		CvScalar col = CV_RGB(0,0,0);
		cvLine(image, convert(pTop), convert(pBottom), col);
		cvLine(image, convert(pLeft), convert(pRight), col);
		CvFont font;
		cvInitFont(&font, CV_FONT_HERSHEY_PLAIN, 0.9, 0.9);
		char num[20];
		sprintf(num, "%f", dMaxY);
		cvPutText(image, num, convert(pTop), &font, col);
		sprintf(num, "%f", dMinY);
		cvPutText(image, num, convert(pBottom), &font, col);
		sprintf(num, "%f", dMinX);
		cvPutText(image, num, convert(pLeft), &font, col);
		sprintf(num, "%f", dMaxX);
		cvPutText(image, num, convert(pRight), &font, col);
	}

	CvPoint convert(pair<int, double> datum)
	{
		return cvPoint(dMargin + (datum.first-dMinX) * dScaleX,
					   /*dHeight -*/ (dMargin + (datum.second-dMinY) * dScaleY));
	}
};

void CPlot::setWhite(IplImage * pIm)
{
	int * pnData = (int *)(void *)pIm->imageData;
	for(int nIntSize = (pIm->widthStep * pIm->height)/4; nIntSize>0; nIntSize--)
	{
		*pnData = -1;
		pnData++;
	}
}

void CPlot::plot(IplImage * image, const std::map<int, double> & aData, bool bPlot0)
{
	double dMax = MIN_INT, dMin = MAX_INT;
	int nMax = MIN_INT, nMin = MAX_INT;
	if(bPlot0)
	{
		dMax = 0, dMin = 0;
		nMax = 0, nMin = 0;
	}

	for(it pDatum = aData.begin(); pDatum != aData.end(); pDatum++)
	{
		if(nMin > pDatum->first) nMin = pDatum->first;
		if(dMin > pDatum->second) dMin = pDatum->second;
		if(dMax < pDatum->second) dMax = pDatum->second;

		nMax = pDatum->first;
	}

	setWhite(image);

	CPlotConvert PC(nMin, nMax, dMin, dMax, image, 30);

	for(it pDatum = aData.begin(); pDatum != aData.end(); pDatum++)
	{
		cvCircle(image, PC.convert(*pDatum), 2, CV_RGB(255,0,0));
	}
}
