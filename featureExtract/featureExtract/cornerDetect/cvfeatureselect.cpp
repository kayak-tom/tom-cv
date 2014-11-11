/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

/*M///////////////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                        Intel License Agreement
//                For Open Source Computer Vision Library
//
// Copyright (C) 2000, Intel Corporation, all rights reserved.
// Third party copyrights are property of their respective owners.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
//   * Redistribution's of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistribution's in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//   * The name of Intel Corporation may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// This software is provided by the copyright holders and contributors "as is" and
// any express or implied warranties, including, but not limited to, the implied
// warranties of merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the Intel Corporation or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused
// and on any theory of liability, whether in contract, strict liability,
// or tort (including negligence or otherwise) arising in any way out of
// the use of this software, even if advised of the possibility of such damage.
//
//M*/
#ifdef DEPRECATED_WITH_NEW_OCV
#include "../src/_cv.h"
#include "goodFeaturesFast.h"

#define  cmp_features( f1, f2 )  (*(f1) > *(f2))

static CV_IMPLEMENT_QSORT( icvSortFeatures, int *, cmp_features ) //Yuk--this was being used to sort floats

void
fastCvThresholdToZero64( CvMat * pEigIm, double thresh )
{
	double * pdVal = pEigIm->data.db;
	for(int i = pEigIm->height*pEigIm->width; i > 0; i-- )
	{
		if(*pdVal < thresh)
			*pdVal = 0;
		pdVal++;
	}
}

//ALSO casts to INT-- a bit hacky but should improve rest of fn
void
fastCvThresholdToZero32( CvMat * pEigIm, double thresh )
{
	float * pdVal = pEigIm->data.fl;
	for(int i = pEigIm->height*pEigIm->width; i > 0; i-- )
	{
		double dVal = *pdVal;
		if(dVal < thresh)
			*((int *)(void *)pdVal) = 0;
		else
			*((int *)(void *)pdVal) = doubleToInt(dVal * 4096.0); //Scale fairly arbitrary

		pdVal++;
	}
}

//Already converted to int:
void
fastCvThresholdToZero32( CvMat * pEigIm, int thresh )
{
	int * pdVal = (int *)pEigIm->data.ptr;
	for(int i = pEigIm->height*pEigIm->width; i > 0; i-- )
	{
		int dVal = *pdVal;
		if(dVal < thresh)
			*pdVal = 0;

		pdVal++;
	}
}

int fastMax32AndToInt( CvMat * pEigIm )
{
	int dMax = 0;
	volatile float * pdVal = pEigIm->data.fl;
	for(volatile int i = pEigIm->height*pEigIm->width; i > 0; i-- )
	{
		volatile int val = doubleToInt(4096.0 * *pdVal);
		*((int *)(void *)pdVal) = val;
		if(val > dMax)
			dMax = val;
		pdVal++;
	}
	return dMax;
}

double fastMax32( CvMat * pEigIm )
{
	double dMax = 0;
	volatile float * pdVal = pEigIm->data.fl;
	for(volatile int i = pEigIm->height*pEigIm->width; i > 0; i-- )
	{
		volatile double val = (double)*pdVal;
		if(val > dMax)
			dMax = val;
		pdVal++;
	}
	return dMax;
}
double fastMax64( CvMat * pEigIm )
{
	double dMax = 0;
	volatile double * pdVal = pEigIm->data.db;
	for(volatile int i = pEigIm->height*pEigIm->width; i > 0; i-- )
	{
		volatile double val = *pdVal;
		if(val > dMax)
			dMax = val;
		pdVal++;
	}
	return dMax;
}
//CV_IMPL
void
fastCvGoodFeaturesToTrack( const void* image, void* eigImage, void* tempImage,
                       CvPoint2D64f* corners, int *corner_count,
                       double quality_level, double min_distance,
                       const void* maskImage, int block_size,
                       int use_harris, double harris_k, CPointBin<32, 32, 16> & pointBin )
{
    CvMat* _eigImg = 0;
    CvMat* _tmpImg = 0;

    CV_FUNCNAME( "cvGoodFeaturesToTrack" );

    __BEGIN__;

    double max_val = 0;
    int max_count = 0;
    int count = 0;
    int x, y, i, k = 0;
    int min_dist;
    int eig_step;//, tmp_step;

    /* when selecting points, use integer coordinates */
    //CvPoint *ptr = (CvPoint *) corners;

    /* process floating-point images using integer arithmetics */
    int *eig_data = 0;
    //int *tmp_data = 0;
    int **ptr_data = 0;
    uchar *mask_data = 0;
    int  mask_step = 0;
    CvSize size;

    int * aRowBefore_data = 0;
    int * aRowAfter_data = 0;

    int coi1 = 0, coi2 = 0, coi3 = 0;
    CvMat  stub, *img = (CvMat*)image;
    CvMat  eig_stub, *eig = (CvMat*)eigImage;
    CvMat  tmp_stub, *tmp = (CvMat*)tempImage;
    CvMat  mask_stub, *mask = (CvMat*)maskImage;

	pointBin.reset();

    if( corner_count )
    {
        max_count = *corner_count;
        *corner_count = 0;
    }

    CV_CALL( img = cvGetMat( img, &stub, &coi1 ));
    if( eig )
    {
        CV_CALL( eig = cvGetMat( eig, &eig_stub, &coi2 ));
    }
    else
    {
        CV_CALL( _eigImg = cvCreateMat( img->rows, img->cols, CV_32FC1 ));
        eig = _eigImg;
    }

    if( tmp )
    {
        CV_CALL( tmp = cvGetMat( tmp, &tmp_stub, &coi3 ));
    }
    else
    {
        CV_CALL( _tmpImg = cvCreateMat( img->rows, img->cols, CV_32FC1 ));
        tmp = _tmpImg;
    }

    if( mask )
    {
        CV_CALL( mask = cvGetMat( mask, &mask_stub ));
        if( !CV_IS_MASK_ARR( mask ))
        {
            CV_ERROR( CV_StsBadMask, "" );
        }
    }

    if( coi1 != 0 || coi2 != 0 || coi3 != 0 )
        CV_ERROR( CV_BadCOI, "" );

    if( CV_MAT_CN(img->type) != 1 ||
        CV_MAT_CN(eig->type) != 1 ||
        CV_MAT_CN(tmp->type) != 1 )
        CV_ERROR( CV_BadNumChannels, cvUnsupportedFormat );

    if( CV_MAT_DEPTH(tmp->type) != CV_32F ||
        CV_MAT_DEPTH(eig->type) != CV_32F )
        CV_ERROR( CV_BadDepth, cvUnsupportedFormat );

    if( !corners || !corner_count )
        CV_ERROR( CV_StsNullPtr, "" );

    if( max_count <= 0 )
        CV_ERROR( CV_StsBadArg, "maximal corners number is non positive" );

    if( quality_level <= 0 || min_distance < 0 )
        CV_ERROR( CV_StsBadArg, "quality level or min distance are non positive" );

    if( use_harris )
    {
        CV_CALL( cvCornerHarris( img, eig, block_size, 3, harris_k ));
    }
    else
    {
        CV_CALL( fastCvCornerMinEigenVal( img, eig, block_size, 3 ));
		//First slow call (66%) --end up with eig image full of evals and evecs
    }
    //CV_CALL( cvMinMaxLoc( eig, 0, &max_val, 0, 0, mask )); //5% --as floats (?) Speedup by avoiding mask?
	max_val = fastMax32AndToInt(eig); // ALSO ROUNDS TO INT

    CV_CALL( fastCvThresholdToZero32( eig, cvRound(max_val * quality_level))); //inlined? fast
    //CV_CALL( cvDilate( eig, tmp ));//19% -- sets each px of temp to max of 3x3 window

    min_dist = cvRound( min_distance /* * min_distance*/ );

    eig_step = eig->step / sizeof(eig_data[0]);
    //tmp_step = tmp->step / sizeof(tmp_data[0]);

	size = cvGetMatSize( img );
    ptr_data = (int**)(tmp->data.ptr);
    eig_data = (int*)(eig->data.ptr);
    aRowBefore_data = eig_data;
    aRowAfter_data = eig_data + 2*eig_step;
    //tmp_data = (int*)(tmp->data.ptr);
    if( mask )
    {
        mask_data = (uchar*)(mask->data.ptr);
        mask_step = mask->step;
    }

    /* collect list of pointers to features - put them into temporary image
	   Was using ints before--only for faster float== though
	*/
    for( y = 1, k = 0; y < size.height - 1; y++ )
    {
        eig_data += eig_step;
        //tmp_data += tmp_step;
        mask_data += mask_step;

		int * pEigData = eig_data;
		int * pRowBeforeData = aRowBefore_data;
		int * pRowAfterData = aRowAfter_data;
		int nLastVal = *pEigData; pEigData++;
		int nThisVal = *pEigData; pEigData++;
		int nLastValAbove = *pEigData; pRowBeforeData++;
		int nThisValAbove = *pEigData; pRowBeforeData++;
        for( x = 1; x < size.width - 1; x++ )
        {
			pRowAfterData++;

			int nNextVal = *pEigData;
			int nNextValAbove = *pRowBeforeData;
			if(nThisVal > 0 && nThisVal >= nLastVal && nThisVal >= nNextVal)
			{
				if(nThisVal >= nLastValAbove && nThisVal >= nThisValAbove && nThisVal >= nNextValAbove)
					if(nThisVal >= *pRowAfterData && nThisVal >= pRowAfterData[-1] && nThisVal >= pRowAfterData[1])
						ptr_data[k++] = eig_data + x;
			}

			pEigData++;
			nLastVal = nThisVal;
			nThisVal = nNextVal;

			pRowBeforeData++;
			nLastValAbove = nThisValAbove;
			nThisValAbove = nNextValAbove;

			//Check is max of 3x3 neighborhood

            /*if( val != 0 && val == tmp_data[x] && (!mask || mask_data[x]) )
                ptr_data[k++] = eig_data + x;*/
        }
		
		aRowBefore_data += eig_step;
		aRowAfter_data += eig_step;
    }

    icvSortFeatures( ptr_data, k, 0 );

    /* select the strongest features */
    for( i = 0; i < k; i++ )
    {
        int ofs = (int)((uchar*)(ptr_data[i]) - eig->data.ptr);
        y = ofs / eig->step;
        x = (ofs - y * eig->step)/sizeof(float);

		bool bTooClose = pointBin.isTooClose(x, y, min_dist);
        /*{
            for( j = 0; j < count; j++ )
            {
                int dx = x - ptr[j].x;
                int dy = y - ptr[j].y;
                int dist = dx * dx + dy * dy;

                if( dist < min_dist )
                    break;
            }
        }*/

        if( !bTooClose )
        {
        	corners[count].x = x;
        	corners[count].y = y;
            if( ++count >= max_count )
                break;
        }
    }

    /* convert points to floating-point format
    for( i = 0; i < count; i++ )
    {
        assert( (unsigned)ptr[i].x < (unsigned)size.width &&
                (unsigned)ptr[i].y < (unsigned)size.height );

        corners[i].x = ptr[i].x; //could do lookup 0..1024? not much in it
        corners[i].y = ptr[i].y;
    }*/

    *corner_count = count;

    __END__;

    cvReleaseMat( &_eigImg );
    cvReleaseMat( &_tmpImg );
}

/* End of file. */
#endif
