#pragma once
#ifndef _PATCH_H
#define _PATCH_H
//The SIFT Rot. Inv. method is described here: http://en.wikipedia.org/wiki/Scale-invariant_feature_transform
#include "image2.h"
#include "HSI.h"

namespace grc
{

//! Control norm used to compare patches
enum ePatchCompMethod {ePatchEuclid, ePatchEuclidSquared, ePatchL1, ePatchMaxDist, /*ePatchCos, ePatchChiSquared, ePatchJeffereys,*/ ePatchCorrel};

//! Avoid reproducing patch-template parameter lists
#define P_TEMPLATE template<int PATCH_RADIUS, ePatchCompMethod PATCH_COMP_METHOD, ImageDepth ImType, typename imDataType, bool PATCH_ORIENT, bool NORMALISE_PATCH, int PATCH_SCALE>
#define cTPatch cPatch<PATCH_RADIUS, PATCH_COMP_METHOD, ImType,imDataType, PATCH_ORIENT, NORMALISE_PATCH, PATCH_SCALE>

#define SQR(x) ((x)*(x))

//#define DEBUG_PATCH (not threadsafe)
#ifdef DEBUG_PATCH
#include "highgui.h"
#endif

//! Implements a descriptor representing a patch of an image.
/*! Patches may be oriented with the steepest gradient direction to provide rotation invariance,
 *  or may have their intensity normalised for illumination invariance */
P_TEMPLATE
class cPatch
{
	static const int DIAMETER = 2*PATCH_RADIUS+1;
	static const int ROT_BINS = 16;
public:
	static const int SIZE = DIAMETER*DIAMETER*ImType;
private:
    //! Patch data vector:
	unsigned char acPatch[SIZE];

    //! Returns the index of the bin (grouped by angle) that vector (x,y) falls into
    int getRotBin(int x, int y)
	{
		//Todo: lookup table
		int nBin = ROT_BINS/2 + (int)floor(ROT_BINS*atan2((double)x,(double)y)/(2*M_PI));
		if(nBin==ROT_BINS) nBin = 0;
		if(IS_DEBUG) CHECK(nBin<0 || nBin >= ROT_BINS, "getRotBin: Bin OOB")
		return nBin;
	}

    //! Conversions not currently in use in this application:
    inline static double DOUBLE_TO_UCHAR() {return  255.999; }; //MSVC doesn't like constants defined here 
	
    inline unsigned char toUChar(unsigned char v) {return (unsigned char)v;};
	inline unsigned char toUChar(float v) {return (unsigned char)(DOUBLE_TO_UCHAR()*v);};
	inline int toInt(unsigned char v) {return (int)v;};
	inline int toInt(float v) {return (int)(DOUBLE_TO_UCHAR()*v);};
	inline double toDouble(unsigned char v) {return (double)v;};
	inline double toDouble(float v) {return (double)(DOUBLE_TO_UCHAR()*v);};
	/*inline int getI(unsigned char * pcPx)
	{
		int nI = (int)(*pcPx);
		if(ImType==RGB)
	    {
	    	nI += 6 * (int)(*(pcPx+1));//G
	    	nI += 3 * (int)(*(pcPx+2));//R
	    }
		return nI;
	};*/
public:
    //! Computes the patch descriptor around point in pIplIm
	cPatch(const IplImage * pIplIm, CvPoint point) //: A(0), B(0), C(0)
	{
        if(!pIplIm) return;

#ifdef DEBUG_PATCH
		cvNamedWindow("OriPatch");
#endif
        Image<ImType, imDataType> Im(pIplIm, true); //only holds a copy

		int x_min=point.x-PATCH_SCALE*PATCH_RADIUS;
		int y_min=point.y-PATCH_SCALE*PATCH_RADIUS;
		int x_max=point.x+PATCH_SCALE*PATCH_RADIUS;
		int y_max=point.y+PATCH_SCALE*PATCH_RADIUS;
		int step = PATCH_SCALE;

		if(!PATCH_ORIENT)
		{
			int idx=0;

			for(int x=x_min; x <= x_max; x += step)
			{
				for(int y=y_min; y <= y_max; y += step)
				{
					acPatch[idx++] = toUChar(Im.b(x,y));
					if(ImType==RGB)
					{
						acPatch[idx++] = toUChar(Im.g(x,y));
						acPatch[idx++] = toUChar(Im.r(x,y));
					}
				}
            }
        }
		else
		{
			//Oriented patch
			int anSum[ROT_BINS], anCount[ROT_BINS], anResponse[ROT_BINS], anResponse2[ROT_BINS];
			for(int i=0; i<ROT_BINS; i++) { anSum[i] = 0; anCount[i] = 0; }

			for(int x=-PATCH_RADIUS; x <= PATCH_RADIUS; x += 1) //could += PATCH_SCALE??
			{
				for(int y=-PATCH_RADIUS; y <= PATCH_RADIUS; y += 1)
				{
					if(x*x+y*y > SQR(PATCH_RADIUS+1)) continue;
					int nRotBin = getRotBin(x, y);
					if(IS_DEBUG) CHECK(0 > nRotBin || nRotBin >= ROT_BINS, "cPatch: Bin OOB")
					double im_x = x + point.x;
					double im_y = y + point.y;
					int nim_x = (int)im_x;
					int nim_y = (int)im_y;
					if(ImType==RGB)
						anSum[nRotBin] +=  toInt(Im.r(nim_x,nim_y)) + toInt(Im.g(nim_x,nim_y)) ;
					else
						anSum[nRotBin] += toInt(Im.b(nim_x,nim_y));

					anCount[nRotBin] += 1; //to enforce even coverage
				}
			}
			//Hack to compensate for 1st bin missing a few px: or anSum[0] = (anCount[1] * anSum[0]) / anCount[0];
			for(int i=0; i<ROT_BINS; i++)
				anSum[i] = (anCount[1] * anSum[i]) / anCount[i];

			int nMaxBin = 0, nMaxBinVal = -1;
			for(int i=0; i<ROT_BINS; i++)
			{
				int nOppositeBin = (i+ROT_BINS/2) % ROT_BINS;
				if(IS_DEBUG) CHECK(nOppositeBin<0 || nOppositeBin >= ROT_BINS, "cPatch: Idx OOB");
				anResponse[i] = anSum[i] - anSum[nOppositeBin];
				if(nMaxBinVal < anResponse[i])
				{
					nMaxBinVal = anResponse[i];
					nMaxBin = i;
				}
			}

			//////// Optional: look at global resp. rather than just 1/4 of the patch
			nMaxBin = 0; nMaxBinVal = -1;
			for(int i=0;i<ROT_BINS;i++)
			{
				int nPrevBin = (i==0) ? (ROT_BINS-1) : (i-1);
				int nNextBin = (i+1) % ROT_BINS;
				if(IS_DEBUG) CHECK(nNextBin<0 || nNextBin >= ROT_BINS || nPrevBin<0 || nPrevBin >= ROT_BINS, "cPatch: Bin Idx OOB");
				anResponse2[i] = anResponse[nPrevBin] + 2*anResponse[i] + anResponse[nNextBin];
				if(nMaxBinVal < anResponse2[i])
				{
					nMaxBinVal = anResponse2[i];
					nMaxBin = i;
				}
			}
			//////////////

			//Todo: Extrapolate to find peak
			double dAngle = nMaxBin * (-2*M_PI/ROT_BINS); //todo: check sign

			int idx=0;
			for(double x=-PATCH_RADIUS; x <= PATCH_RADIUS; x += 1)
			{
				for(double y=-PATCH_RADIUS; y <= PATCH_RADIUS; y += 1)
				{
					//Rotate backwards to a point
					double s=sin(dAngle), c=cos(dAngle);
					double x_im = c*x - s*y;
					double y_im = s*x + c*y;
					x_im *= PATCH_SCALE;
					y_im *= PATCH_SCALE;

					//Find 4 points about this point and their weight
					double x_lo=floor(x_im);
					double y_lo=floor(y_im);
					double  dB = 0, dG = 0, dR = 0;
					double dx=x_im-x_lo, dy=y_im-y_lo;
					int im_x_lo = (int)x_lo + point.x;
					int im_y_lo = (int)y_lo + point.y;
					int im_x_hi = im_x_lo+1;
					int im_y_hi = im_y_lo+1;
					double coeff_ll = (1-dx)*(1-dy);
					double coeff_lh = (1-dx)*dy;
					double coeff_hl = dx*(1-dy);
					double coeff_hh = dx*dy;
#define BILIN(CHANNEL) Im.CHANNEL(im_x_lo,im_y_lo)*coeff_ll + Im.CHANNEL(im_x_hi,im_y_lo)*coeff_hl + Im.CHANNEL(im_x_lo,im_y_hi)*coeff_lh + Im.CHANNEL(im_x_hi,im_y_hi)*coeff_hh
					dB = BILIN(b);

					acPatch[idx++] = (unsigned char)(dB );
					if(ImType==RGB)
					{
						dG = BILIN(g);
						dR = BILIN(r);
						acPatch[idx++] = (unsigned char)(dG );
						acPatch[idx++] = (unsigned char)(dR );
					}
				}
			}
        }
		
        if(NORMALISE_PATCH)
		{
	        static const int AREA=DIAMETER*DIAMETER;

	        int nTotal = 0;
	        int acPatchHSI[SIZE]; //Copy ISH to here
	        int * pnPatchHSI = acPatchHSI;
	        unsigned char * pcPatch = acPatch;

		    for(int nBin=AREA; nBin>0; nBin--)
	        {
		    	if(ImType==RGB)
		    	{
		    		//Convert to HSI so we can normalise intensity (should we convert back to avoid cyclic stuff?)
		        	int nHue,nSaturation,nIntensity;
		            CHsiInt::HSI(*(pcPatch+2),*(pcPatch+1),*pcPatch,nHue,nSaturation,nIntensity);
		            *pnPatchHSI = nIntensity;
		            *(pnPatchHSI+1) = nSaturation;
		            *(pnPatchHSI+2) = nHue;
		    	}
		    	else
		    	{
		    		*pnPatchHSI = (int)*pcPatch;
		    	}
				nTotal += *pnPatchHSI;

				pcPatch += ImType;
	    		pnPatchHSI += ImType;
	        }

	        static const double PATCH_AREA_INV = 1.0/((double)(AREA));

	        double dAvI = (double)nTotal*PATCH_AREA_INV;

	        double dVarTotal = 0;
	        pnPatchHSI = acPatchHSI;
	        for(int nBin=AREA; nBin>0; nBin--)
	        {
	            double dDif = *pnPatchHSI - dAvI;
			    dVarTotal += dDif*dDif;
			    pnPatchHSI += ImType;
	        }

	        double dSD = sqrt(dVarTotal*PATCH_AREA_INV);

	        const double dTargetMean = 128;
	        const double dTargetSD = 40;
	        const double dScale = dTargetSD/dSD;
	        const double dShift = dTargetMean-dScale*dAvI;

	        //Scale intensities
	        pnPatchHSI = acPatchHSI;
	        pcPatch = acPatch;
		    for(int idx=AREA; idx > 0; idx--)
		    {
	            //Normalise
			    static const int PATCH_MAX_VAL = 256;
	            int nIntensityCorrected = (int)(dShift+dScale*(*pnPatchHSI));
	            if(nIntensityCorrected>=PATCH_MAX_VAL) nIntensityCorrected = PATCH_MAX_VAL-1;
	            if(nIntensityCorrected<0) nIntensityCorrected = 0;

	            // *pnPatchHSI = nIntensityCorrected;
	            //Now copy back, or convert back to RGB

		    	if(ImType==RGB)
		    	{
		    		//convert back to RGB to avoid cyclic stuff?)
		        	int nR,nG,nB;
		            CHsiInt::RGB (*(pnPatchHSI+2),*(pnPatchHSI+1),*pnPatchHSI,nR,nG,nB);
		            *pcPatch = (unsigned char)nB;
		            *(pcPatch+1) = (unsigned char)nG;
		            *(pcPatch+2) = (unsigned char)nR;
		    	}
		    	else
		    	{
		    		*pcPatch = (unsigned char)*pnPatchHSI;
		    	}

				pcPatch += ImType;
	    		pnPatchHSI += ImType;
		    }
		}
			///////////
#ifdef DEBUG_PATCH
			Image<ImType, unsigned char> patch(DIAMETER, DIAMETER);
			idx = 0;
			for(int x=0; x < DIAMETER; x ++)
			{
				for(int y=0; y < DIAMETER; y ++)
				{
					patch.b(x,y) = acPatch[idx++];
					if(ImType==RGB)
					{
						patch.g(x,y) = acPatch[idx++];
						patch.r(x,y) = acPatch[idx++];
					}

				}
			}
			cvShowImage("OriPatch", patch);
			cvWaitKey(0);
#endif
			//////////////
		/*unsigned char * pcPatch = acPatch;
		//for(int i=SIZE; i>0; i--)
		//{
		//	int nBinVal = ((int)*pcPatch);
		//	A += nBinVal;
		//	B += nBinVal * nBinVal;
		//	pcPatch++;
		//}
		//C = intLookup::Sqrt(SIZE*B - A*A)/100; OVERFLOW here
		//double dB = 0;
		for(int i=SIZE; i>0; i--)
		{
			int nBinVal = ((int)*pcPatch);
			A += nBinVal;
			B += nBinVal * nBinVal;
			pcPatch++;
		}
		double dA = (double)A;
		C = sqrt(SIZE*B - dA*dA)/100.0;*/

	};

private:
    //! Calculate distance between 2 patched
    inline static double euclidDistSquared(const double * ad1, const double * ad2, int nLen)
    {
        double d=0;
        for(int i=0;i<nLen;i++)
        {
            double temp = ad1[i] - ad2[i];
            d += temp*temp;
        }

        return d;
    };

    //! Calculate distance between 2 patched
    template<class CHAR>
    inline static int euclidDistSquared(const CHAR * ac1, const CHAR * ac2, int nLen)
    {
        int d=0;
        for(int i=nLen; i>0; i--)
        {
            int temp = (int)(*ac1) - (int)(*ac2);
            d += temp*temp;
            ac1++;
            ac2++;
            //std::cout << temp << ' ';
        }

        return d;
    };

    //! Fast approximation to distance between 2 patches
    template<class CHAR>
    inline static int fastEuclidDistSquared(const CHAR * ac1, const CHAR * ac2, int nLen)
    {
        int d=0;
        //ac1+=3;//Look at the chars that are already aligned
        //ac2+=3;
        for(int i=nLen/4; i>0; i--)
        {
            int temp = (int)(*ac1) - (int)(*ac2);
            d += temp*temp;
            ac1+=4;
            ac2+=4;
            //std::cout << temp << ' ';
        }

        return d;
    };

    //! Fast approximation to distance between 2 patches
    inline static int fastEuclidDistSquared2(const void * ac1, const void * ac2, int nLen)
    {
        int d=0;
        //ac1+=3;//Look at the chars that are already aligned
        //ac2+=3;
        const unsigned int * an1 = (const unsigned int *)ac1;
        const unsigned int * an2 = (const unsigned int *)ac2;
        for(int i=nLen/4; i>0; i--)
        {
            unsigned int temp = ((*an1) >> 24) - ((*an2) >> 24);
            d += temp*temp;
            an1++;
            an2++;
            //std::cout << temp << ' ';
        }

        return d;
    };

public:
    //! Distance between 2 patches
	int distance(const cPatch * pPatch) const
	{
		switch (PATCH_COMP_METHOD)
		{
		case ePatchEuclidSquared:
			return euclidDistSquared(acPatch, pPatch->acPatch, SIZE);
		/*case ePatchMaxDist:
			return MaxDist(acPatch, pPatch->acPatch, SIZE);
		case ePatchL1:
			return L1Dist(acPatch, pPatch->acPatch, SIZE);*/
		/*case ePatchCos:
			return CompareCosine(pH);
	    case ePatchChiSquared:
			return CompareChiSquared(pH);
	    case ePatchJeffereys:
			return CompareJeffreys(pH);*/
	    //case ePatchCorrel:
		//	return correl(pPatch);
		}

	    throw new GRCException("cTPatch::Compare: invalid comparison selected");
	};

    //! Fast approximation to distance between 2 patches
	int fastDistance(const cPatch * pPatch) const
	{
		switch (PATCH_COMP_METHOD)
		{
		case ePatchEuclidSquared:
//			return fastEuclidDistSquared3(acPatch, pPatch->acPatch, SIZE);
			return euclidDistSquared(acPatch, pPatch->acPatch, SIZE/4);
		}

	    throw new GRCException("cTPatch::fastDistance: invalid comparison selected");
	};

	/*int correl(const cPatch * pPatch) const
	{
		double D=0;
		const unsigned char * pcPatch1 = acPatch;
		const unsigned char * pcPatch2 = pPatch->acPatch;
		for(int i=SIZE; i>0; i--)
		{
			D += ((int)*pcPatch1) * ((int)*pcPatch2);
			pcPatch1++;
			pcPatch2++;
		}
		return (SIZE*D - A*pPatch->A) / (C*pPatch->C);
	};*/
};
}

#endif //_PATCH_H
