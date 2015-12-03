/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*
 * convert_OpenCV.cpp
 *
 *  Created on: 12/10/2009
 *      Author: tom
 */

#include "convert_OpenCV.h"
#include "imageAccess.h"

template<>
CvPtr<IplImage>::~CvPtr() { cvReleaseImage(&ptr); }

template<>
CvPtr<CvMat>::~CvPtr() { cvReleaseMat(&ptr); }


template<int nScale, int nChannels>
void downSample(const IplImage * pSrc, IplImage * pDest)
{
    if(IS_DEBUG) CHECK(!(pSrc->nChannels ==nChannels && pDest->nChannels == nChannels), "downSample: Wring num channels");

    setZero(pDest);
    typedef uchar TDest;
    //Iterate over src rows incrementing dest px.
    TDest * pnDestRow = (TDest *)(void *)pDest->imageData;
    uchar * pcSrcRow = (uchar *)(void *)pSrc->imageData;
    for(int nDestRow=0; nDestRow < pDest->height; nDestRow++)
    {
        for(int nSrcSubRow=0; nSrcSubRow<nScale; nSrcSubRow++)
        {
            uchar * pcSrc = pcSrcRow;
            TDest * pnDest = pnDestRow;

            for(int nDestCol=0; nDestCol < pDest->width; nDestCol++)
            {
                volatile int anTotalVal[nChannels];
                for(int i=0; i<nChannels; i++)
                    anTotalVal[i] = 0;
                for(int nSrcSubCol=0; nSrcSubCol<nScale; nSrcSubCol++)
                {
                    for(int i=0; i<nChannels; i++)
                    {
                        anTotalVal[i] += (int)*pcSrc;
                        pcSrc++;
                    }
                }

                for(int i=0; i<nChannels; i++)
                {
                    *pnDest += (TDest)( anTotalVal[i]/(nScale*nScale) );
                    //cout << "Setting " << nDestRow << ',' << nDestCol << " to " << *pnDest << '\n';
                    pnDest++;
                }
            }

            pcSrcRow += pSrc->widthStep;
        }
        pnDestRow += pDest->widthStep/sizeof(TDest);
    }
    //cvSaveImage("Dest.bmp", pDest);
}

template<int nChannels>
void downSample(const int nScale, const IplImage * pSrc, IplImage * pDest)
{
    switch(nScale)
    {
    case 1:
        return downSample<1, nChannels>(pSrc, pDest);
    case 2:
        return downSample<2, nChannels>(pSrc, pDest);
    case 3:
        return downSample<3, nChannels>(pSrc, pDest);
    case 4:
        return downSample<4, nChannels>(pSrc, pDest);
    case 5:
        return downSample<5, nChannels>(pSrc, pDest);
    case 6:
        return downSample<6, nChannels>(pSrc, pDest);
    case 7:
        return downSample<7, nChannels>(pSrc, pDest);
    case 8:
        return downSample<8, nChannels>(pSrc, pDest);
    case 9:
        return downSample<9, nChannels>(pSrc, pDest);
    case 10:
        return downSample<10, nChannels>(pSrc, pDest);

    }
    THROW( "downSample: Bad scale");
}

void doDownSample(const int nScale, const IplImage * pSrc, IplImage * pDest)
{
    if(IS_DEBUG) CHECK(!(pSrc->nChannels == 1 || pSrc->nChannels == 3), "Unsupported num of channels in source image");
    if(IS_DEBUG) CHECK(!(pDest->nChannels == 1 || pDest->nChannels == 3) || pDest->nChannels != pSrc->nChannels, "Unsupported num of channels in dest image");

    if(pSrc->nChannels == 1)
        downSample<1>(nScale, pSrc, pDest);
    else if(pSrc->nChannels == 3)
        downSample<3>(nScale, pSrc, pDest);
}

void setZero(IplImage * pIm)
{
    int * pnData = (int *)(void *)pIm->imageData;
    for(int nIntSize = ((pIm->widthStep / sizeof(int)) * pIm->height); nIntSize>0; nIntSize--)
    {
        *pnData = 0;
        pnData++;
    }
}

void markCorrespondences(const CBoWCorrespondences * pCorr, const CMask & mask, const IplImage * pIm1, const IplImage * pIm2, IplImage ** ppImOut)
{
    CvSize size=cvSize(pIm1->width*2, pIm1->height);

    if(!*ppImOut)
        *ppImOut = cvCreateImage(size, pIm1->depth, pIm1->nChannels);

    IplImage * pImOut = *ppImOut;

    CvRect r1=cvRect(0, 0, pIm1->width, pIm1->height);
    CvRect r2=cvRect(pIm1->width, 0, pIm1->width, pIm1->height);

    IplImage subIm1 = *pImOut;
    IplImage subIm2 = *pImOut;

    CIplPx<uchar>::cropImageToRect(subIm1, r1);
    CIplPx<uchar>::cropImageToRect(subIm2, r2);
    
    cvCopy(pIm1, &subIm1);
    cvCopy(pIm2, &subIm2);

    int i=0;
    for(CBoWCorrespondences::const_iterator pC = pCorr->begin(); pC != pCorr->end(); pC++, i++)
    {
        CvScalar col = (mask[i]) ? CV_RGB(0,255,0) : CV_RGB(255,0,0);
        CvPoint p1 = locToCvPoint(pC->Location1());
        CvPoint p2 = locToCvPoint(pC->Location2());
        p2.x += pIm1->width;

        cvLine(pImOut, p1, p2, col, 2);
    }
}

void fastBlur(cv::Mat & im, const double dSigma/*, const double dSigma_box*/)
{
	//We need 3 integer filter radii n1,n2,n3 that sum to dSigma_box

	/*
	Best approximations: Sigma, Sigma_box, mean difference between blurred and fastBlurred
	5	12	0.0019651
10	26	0.122756
15	39	0.122226
20	55	0.364318 [I don't know why this is a worse approximation]
25	67	0.0178849
30	79	0.00146146
35	88	0.00202969
y = 2.6086x

  Code for testing:

  void testFilters(const cv::Mat & M)
  {
  cv::Mat gaussian, approx1, approx2;
  for (double dSigma = 5; dSigma <= 35; dSigma += 5)
  {
  ALWAYS_VERBOSE;
  SHOW(M)
  M.copyTo(gaussian);

  blur(gaussian, dSigma);
  SHOW(gaussian);

  //    fastBlur(approx1, dSigma, false);
  //  SHOW(approx1);

  for (double dSigma_box = std::floor(1.5*dSigma); dSigma_box <= 3*dSigma; dSigma_box++)
  {
  M.copyTo(approx2);
  fastBlur(approx2, dSigma, dSigma_box);

  const double dDiff = cv::mean(approx2 - gaussian)[0];

  SHOW2(approx2, toString(dDiff));

  cout << dSigma << '\t' << dSigma_box << '\t' << dDiff << endl;
  }
  }
  //SHOW(approx1 - approx2);
  //SHOW(approx2 - gaussian);

  TIME("exact", for (int i = 0; i < 100; i++) blur(gaussian, 22););
  TIME("approx", for (int i = 0; i < 100; i++) fastBlur(gaussian, 22, true););

  }

	*/

	const double dSigma_box = dSigma*2.6086;

	const int nTotalSize = cvRound(dSigma_box);
	const int n1 = nTotalSize / 3;
	const int n2 = (nTotalSize - n1) / 2;
	const int n3 = (nTotalSize - (n1+n2));

	const std::vector<int> aSizes = { n1, n2, n3 };

	for (const int nRad : aSizes)
	{
		const int nDiam = 2 * nRad + 1;
		cv::boxFilter(im, im, im.depth(), cv::Size(nDiam, nDiam));
	}
}
