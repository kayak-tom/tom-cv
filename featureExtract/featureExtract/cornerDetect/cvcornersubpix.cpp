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
#include "util/exception.h"
#include "util/convert.h"

//in cvSamplers.cpp --rarely actually called
static const void*
icvAdjustRect( const void* srcptr, int src_step, int pix_size,
               CvSize src_size, CvSize win_size,
               CvPoint ip, CvRect* pRect )
{
    CvRect rect;
    const char* src = (const char*)srcptr;

    if( ip.x >= 0 )
    {
        src += ip.x*pix_size;
        rect.x = 0;
    }
    else
    {
        rect.x = -ip.x;
        if( rect.x > win_size.width )
            rect.x = win_size.width;
    }

    if( ip.x + win_size.width < src_size.width )
        rect.width = win_size.width;
    else
    {
        rect.width = src_size.width - ip.x - 1;
        if( rect.width < 0 )
        {
            src += rect.width*pix_size;
            rect.width = 0;
        }
        assert( rect.width <= win_size.width );
    }

    if( ip.y >= 0 )
    {
        src += ip.y * src_step;
        rect.y = 0;
    }
    else
        rect.y = -ip.y;

    if( ip.y + win_size.height < src_size.height )
        rect.height = win_size.height;
    else
    {
        rect.height = src_size.height - ip.y - 1;
        if( rect.height < 0 )
        {
            src += rect.height*src_step;
            rect.height = 0;
        }
    }

    *pRect = rect;
    return src - rect.x*pix_size;
}


CvStatus fasticvGetRectSubPix_8u32f_C1R
( const uchar* src, int src_step, CvSize src_size,
  double* dst, int dst_step, CvSize win_size, CvPoint2D64f center )
{
    CvPoint ip;
    double  a12, a22, b1, b2;
    double a, b;
    double s = 0;
    int i, j;

    center.x -= (win_size.width-1)*0.5;
    center.y -= (win_size.height-1)*0.5;

    ip.x = cvFloor( center.x );
    ip.y = cvFloor( center.y );

    if( win_size.width <= 0 || win_size.height <= 0 )
        return CV_BADRANGE_ERR;

    a = center.x - ip.x;
    b = center.y - ip.y;
    a = MAX(a,0.0001);
    a12 = a*(1.-b);
    a22 = a*b;
    b1 = 1. - b;
    b2 = b;
    s = (1. - a)/a;

    src_step /= sizeof(src[0]);
    dst_step /= sizeof(dst[0]);

    if( 0 <= ip.x && ip.x + win_size.width < src_size.width &&
        0 <= ip.y && ip.y + win_size.height < src_size.height )
    {
        // extracted rectangle is totally inside the image
        src += ip.y * src_step + ip.x;

        for( ; win_size.height--; src += src_step, dst += dst_step )
        {
            double prev = (1 - a)*(b1*intLookup::uchar2double(src[0]) + b2*intLookup::uchar2double(src[src_step]));
            for( j = 0; j < win_size.width; j++ )
            {
                volatile double t = a12*intLookup::uchar2double(src[j+1]) + a22*intLookup::uchar2double(src[j+1+src_step]);
                dst[j] = (prev + t);
                prev = (t*s);
            }
        }
    }
    else
//    	THROW( "Not implemented--SP window should be well within image")
    {
    	//rarely used--probably indicates an error
        CvRect r;

        src = (const uchar*)icvAdjustRect( src, src_step*sizeof(*src),
                               sizeof(*src), src_size, win_size,ip, &r);

        for( i = 0; i < win_size.height; i++, dst += dst_step )
        {
            const uchar *src2 = src + src_step;

            if( i < r.y || i >= r.height )
                src2 -= src_step;

            for( j = 0; j < r.x; j++ )
            {
                double s0 = CV_8TO32F(src[r.x])*b1 +
                           CV_8TO32F(src2[r.x])*b2;

                dst[j] = s0;
            }

            if( j < r.width )
            {
                double prev = (1 - a)*(b1*intLookup::uchar2double(src[j]) + b2*intLookup::uchar2double(src2[j]));

                for( ; j < r.width; j++ )
                {
                    double t = a12*intLookup::uchar2double(src[j+1]) + a22*intLookup::uchar2double(src2[j+1]);
                    dst[j] = prev + t;
                    prev = (t*s);
                }
            }

            for( ; j < win_size.width; j++ )
            {
                double s0 = intLookup::uchar2double(src[r.width])*b1 +
                		intLookup::uchar2double(src2[r.width])*b2;

                dst[j] = s0;
            }

            if( i < r.height )
                src = src2;
        }
    }

    return CV_OK;
}

/* lightweight convolution with 3x3 kernel */
void fasticvSepConvSmall3_32f( double* src, int src_step, double* dst, int dst_step,
            CvSize src_size, const double* kx, const double* ky, double* buffer )
{
    int  dst_width, buffer_step = 0;
    int  x, y;

    assert( src && dst && src_size.width > 2 && src_size.height > 2 &&
            (src_step & 3) == 0 && (dst_step & 3) == 0 &&
            (kx || ky) && (buffer || !kx || !ky));

    src_step /= sizeof(src[0]);
    dst_step /= sizeof(dst[0]);

    dst_width = src_size.width - 2;

    if(IS_DEBUG) CHECK( !kx, "Need kx cos otherwise buffer is of wrong type" );
    {

        /* set vars, so that vertical convolution
           will write results into destination ROI and
           horizontal convolution won't run */
        /*src_size.width = dst_width;
        buffer_step = dst_step;
        buffer = dst;
        dst_width = 0;*/
    }

    assert( src_step >= src_size.width && dst_step >= dst_width );

    src_size.height -= 3;
    if( !ky )
    {
        /* set vars, so that vertical convolution won't run and
           horizontal convolution will write results into destination ROI */
        src_size.height += 3;
        buffer_step = src_step;
        buffer = src;
        src_size.width = 0;
    }

    for( y = 0; y <= src_size.height; y++, src += src_step,
                                           dst += dst_step,
                                           buffer += buffer_step )
    {
    	double* src2 = src + src_step;
    	double* src3 = src + src_step*2;
        for( x = 0; x < src_size.width; x++ )
        {
            buffer[x] = (ky[0]*src[x] + ky[1]*src2[x] + ky[2]*src3[x]);
        }

        for( x = 0; x < dst_width; x++ )
        {
            dst[x] = (kx[0]*buffer[x] + kx[1]*buffer[x+1] + kx[2]*buffer[x+2]);
        }
    }
}

void
fastCvFindCornerSubPix( const void* srcarr, CvPoint2D64f* corners,
                    int count, CvSize win, CvSize zeroZone,
                    CvTermCriteria criteria )
{
    double* buffer = 0;

    CV_FUNCNAME( "cvFindCornerSubPix" );

    __BEGIN__;

    const int MAX_ITERS = 100;
    const double drv_x[] = { -1.f, 0.f, 1.f };
    const double drv_y[] = { 0.f, 0.5f, 0.f };
    double *maskX;
    double *maskY;
    double *mask;
    double *src_buffer;
    double *gx_buffer;
    double *gy_buffer;
    int win_w = win.width * 2 + 1, win_h = win.height * 2 + 1;
    int win_rect_size = (win_w + 4) * (win_h + 4);
    double coeff;
    CvSize size, src_buf_size;
    int i, j, k, pt_i;
    int max_iters, buffer_size;
    double eps;

    CvMat stub, *src = (CvMat*)srcarr;
    CV_CALL( src = cvGetMat( srcarr, &stub ));

    if( CV_MAT_TYPE( src->type ) != CV_8UC1 )
        CV_ERROR( CV_StsBadMask, "" );

    if( !corners )
        CV_ERROR( CV_StsNullPtr, "" );

    if( count < 0 )
        CV_ERROR( CV_StsBadSize, "" );

    if( count == 0 )
        EXIT;

    if( win.width <= 0 || win.height <= 0 )
        CV_ERROR( CV_StsBadSize, "" );

    size = cvGetMatSize( src );

    if( size.width < win_w + 4 || size.height < win_h + 4 )
        CV_ERROR( CV_StsBadSize, "" );

    /* initialize variables, controlling loop termination */
    switch( criteria.type )
    {
    case CV_TERMCRIT_ITER:
        eps = 0.f;
        max_iters = criteria.max_iter;
        break;
    case CV_TERMCRIT_EPS:
        eps = criteria.epsilon;
        max_iters = MAX_ITERS;
        break;
    case CV_TERMCRIT_ITER | CV_TERMCRIT_EPS:
        eps = criteria.epsilon;
        max_iters = criteria.max_iter;
        break;
    default:
        assert( 0 );
        CV_ERROR( CV_StsBadFlag, "" );
    }

    eps = MAX( eps, 0 );
    eps *= eps;                 /* use square of error in comparsion operations. */

    max_iters = MAX( max_iters, 1 );
    max_iters = MIN( max_iters, MAX_ITERS );

    /* setup buffer */
    buffer_size = (win_rect_size * 5 + win_w + win_h + 32) * sizeof(double);
    buffer = (double*)cvAlloc( buffer_size );

    /* assign pointers */
    maskX = buffer;
    maskY = maskX + win_w + 4;
    mask = maskY + win_h + 4;
    src_buffer = mask + win_w * win_h;
    gx_buffer = src_buffer + win_rect_size;
    gy_buffer = gx_buffer + win_rect_size;

    coeff = 1. / (win.width * win.width);

    /* calculate mask */
    for( i = -win.width, k = 0; i <= win.width; i++, k++ )
    {
        maskX[k] = exp( -i * i * coeff );
    }

    if( win.width == win.height )
    {
        maskY = maskX;
    }
    else
    {
        coeff = 1. / (win.height * win.height);
        for( i = -win.height, k = 0; i <= win.height; i++, k++ )
        {
            maskY[k] = exp( -i * i * coeff );
        }
    }

    for( i = 0; i < win_h; i++ )
    {
        for( j = 0; j < win_w; j++ )
        {
            mask[i * win_w + j] = maskX[j] * maskY[i];
        }
    }


    /* make zero_zone */
    if( zeroZone.width >= 0 && zeroZone.height >= 0 &&
        zeroZone.width * 2 + 1 < win_w && zeroZone.height * 2 + 1 < win_h )
    {
        for( i = win.height - zeroZone.height; i <= win.height + zeroZone.height; i++ )
        {
            for( j = win.width - zeroZone.width; j <= win.width + zeroZone.width; j++ )
            {
                mask[i * win_w + j] = 0;
            }
        }
    }

    /* set sizes of image rectangles, used in convolutions */
    src_buf_size.width = win_w + 2;
    src_buf_size.height = win_h + 2;

    /* do optimization loop for all the points */
    for( pt_i = 0; pt_i < count; pt_i++ )
    {
        CvPoint2D64f cT = corners[pt_i], cI = cT;
        int iter = 0;
        double err;

        do
        {
        	CvPoint2D64f cI2;
            double a, b, c, bb1, bb2;

            IPPI_CALL( fasticvGetRectSubPix_8u32f_C1R( (uchar*)src->data.ptr, src->step, size,
                                        src_buffer, (win_w + 2) * sizeof( src_buffer[0] ),
                                        cvSize( win_w + 2, win_h + 2 ), cI ));

            /* calc derivatives */
            fasticvSepConvSmall3_32f( src_buffer, src_buf_size.width * sizeof(src_buffer[0]),
                                  gx_buffer, win_w * sizeof(gx_buffer[0]),
                                  src_buf_size, drv_x, drv_y, buffer ); //have doubled buffer size earlier

            fasticvSepConvSmall3_32f( src_buffer, src_buf_size.width * sizeof(src_buffer[0]),
                                  gy_buffer, win_w * sizeof(gy_buffer[0]),
                                  src_buf_size, drv_y, drv_x, buffer );

            a = b = c = bb1 = bb2 = 0;

            /* process gradient */
            for( i = 0, k = 0; i < win_h; i++ )
            {
                double py = i - win.height;

                for( j = 0; j < win_w; j++, k++ )
                {
                    double m = mask[k];
                    double tgx = gx_buffer[k];
                    double tgy = gy_buffer[k];
                    double gxx = tgx * tgx * m;
                    double gxy = tgx * tgy * m;
                    double gyy = tgy * tgy * m;
                    double px = j - win.width;

                    a += gxx;
                    b += gxy;
                    c += gyy;

                    bb1 += gxx * px + gxy * py;
                    bb2 += gxy * px + gyy * py;
                }
            }

            {
                //double A[4];
                /*double InvA[4];
                CvMat matA, matInvA;

                A[0] = a;
                A[1] = A[2] = b;
                A[3] = c;

                cvInitMatHeader( &matA, 2, 2, CV_64F, A );
                cvInitMatHeader( &matInvA, 2, 2, CV_64FC1, InvA );

                cvInvert( &matA, &matInvA, CV_SVD );*/
                volatile double det = (a*c-b*b);
                if(fabs(det) < 0.0001)
                	break;
                volatile double det_inv = 1.0/det;
                /*InvA[0] = c * det_inv;
                InvA[1] = InvA[2] = d * det_inv;
                InvA[3] = a * det_inv;*/
                volatile double InvA_diag0 = c * det_inv;
                volatile double InvA_offDiag = -b * det_inv;
                volatile double InvA_diag1 = a * det_inv;
                cI2.x = (cI.x + InvA_diag0*bb1 + InvA_offDiag*bb2);
                cI2.y = (cI.y + InvA_offDiag*bb1 + InvA_diag1*bb2);
            }

            err = (cI2.x - cI.x) * (cI2.x - cI.x) + (cI2.y - cI.y) * (cI2.y - cI.y);
            cI = cI2;
        }
        while( ++iter < max_iters && err > eps );

        /* if new point is too far from initial, it means poor convergence.
           leave initial point as the result */
        if( fabs( cI.x - cT.x ) > win.width || fabs( cI.y - cT.y ) > win.height )
        {
            cI = cT;
        }

        corners[pt_i] = cI;     /* store result */
    }

    __CLEANUP__;
    __END__;

    cvFree( &buffer );
}

/* End of file. */
#endif
