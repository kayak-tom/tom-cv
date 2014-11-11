/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

#include "util/opencv.h"
#include "util/convert.h"
#include <algorithm>
#include <iostream>

template<int H_BINS, int V_BINS, int MAX_POINTS>
class CPointBin
{
	int aanBinCount[H_BINS][V_BINS];
	struct CPoint{
		int x, y;
		CPoint(int x, int y) : x(x), y(y) {}
		CPoint() {}
	};
	CPoint aaPoints[H_BINS][V_BINS][MAX_POINTS];
	const int nWidth, nHeight;
	bool bInit;

	//Assume images < 1024x1024. Each bin will end up with less than(1+1024/H_BINS)*(1+1024/V_BINS)/(MIN_DIST*MIN_DIST) points
	//Say MAX_POINTS==32 usually plenty
public:
	CPointBin(int nWidth, int nHeight) : nWidth(nWidth), nHeight(nHeight), bInit(false) {}

	void reset()
	{
		DEBUGONLY(int nMaxCount = 0, nFull=0);
		for(int i=0; i<H_BINS; i++)
		{
			int * pBin = aanBinCount[i];
			for(int j=V_BINS; j>0; j--)
			{
				DEBUGONLY(
						if(bInit && *pBin > nMaxCount)
							nMaxCount = *pBin;
						if(bInit && *pBin == MAX_POINTS)
							nFull++;
				);

				*pBin = 0; pBin++;
			}
		}
		DEBUGONLY(
				if(bInit)
				{
					std::cout << "POINT BIN USE: ";
					if(nMaxCount < MAX_POINTS-2)
					{
						std::cout << nMaxCount << " = max count, " << MAX_POINTS << " max" << std::endl;
					}
					else
					{
						double dPropFull = nFull/(double)H_BINS*V_BINS;
						if(dPropFull > 0.05)
							std::cout << dPropFull << " of bins full\n";
						else
							std::cout << "OK\n";
					}
				}
		);
		bInit = true;
	}

	//For measuring spatial distn.:
	int countOccupancy()
	{
		int nCount = 0;
		for(int i=0; i<H_BINS; i++)
		{
			int * pBin = aanBinCount[i];
			for(int j=V_BINS; j>0; j--)
			{
				if(*pBin) nCount++;
				pBin++;
			}
		}
		return nCount;
	}

	bool isTooClose(int x, int y, int rad)
	{
		if(IS_DEBUG) CHECK(!bInit, "Point bin not yet init");

		int radSq = sqr(rad);

		int hBin = (H_BINS * x)/nWidth;
		int vBin = (V_BINS * y)/nHeight;

		int hBinMax = (H_BINS * x+rad)/nWidth;
		int vBinMax = (V_BINS * y+rad)/nHeight;
		int hBinMin = (H_BINS * x-rad)/nWidth;
		int vBinMin = (V_BINS * y-rad)/nHeight;

		const int nHBinStart = std::max<int>(hBinMin, 0); 
		const int nHBinEnd = std::min<int>(hBinMax, H_BINS-1); 
		const int nVBinStart = std::max<int>(vBinMin, 0); 
		const int nVBinEnd = std::min<int>(vBinMax, V_BINS-1); 

		for(int nHSearchBin = nHBinStart; nHSearchBin <= nHBinEnd; nHSearchBin++)
		{
			for(int nVSearchBin = nVBinStart; nVSearchBin <= nVBinEnd; nVSearchBin++)
			{
				int nBinCount = aanBinCount[nHSearchBin][nVSearchBin];
		
				CPoint * p = aaPoints[nHSearchBin][nVSearchBin];
				for(int i=nBinCount; i>0; i--)
				{
					int nDist = sqr(p->x-x)+sqr(p->y-y);
					if(nDist < radSq) return true;
					p++;
				}
			}
		}

		//Add to bin
		int nAlready = aanBinCount[hBin][vBin];
		if(nAlready<MAX_POINTS)
		{
			aaPoints[hBin][vBin][nAlready] = CPoint(x, y);
			aanBinCount[hBin][vBin] = nAlready+1;
		}

		return false;
	}

};

void fastCvGoodFeaturesToTrack( const void* image, void* eigImage, void* tempImage,
                       CvPoint2D64f* corners, int *corner_count,
                       double quality_level, double min_distance,
                       const void* maskImage, int block_size,
                       int use_harris, double harris_k, CPointBin<32, 32, 16> & pointBin );

void fastCvCornerMinEigenVal( const CvMat* srcarr, CvMat* eigenvarr,
                     int block_size, int aperture_size );
void
fastIcvCornerEigenValsVecs( const CvMat* src, CvMat* eigenv, int block_size,
                        int aperture_size, int op_type, double k=0. );

void fastCvFindCornerSubPix( const void* srcarr, CvPoint2D64f* corners,
                    int count, CvSize win, CvSize zeroZone,
                    CvTermCriteria criteria );

