/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */
#pragma once
//The SIFT Rot. Inv. method is described here: http://en.wikipedia.org/wiki/Scale-invariant_feature_transform
#ifndef PATCHFACTORY_H_
#pragma message("Don't #include this, its slow. You probably want featureExtract/patchFeatureExtractor.h instead")
#endif

#include "image/imageAccess.h"
#include "util/convert.h"
#include "util/fastnorms.h"
#include "image/HSI.h"

#include "patchParams.h"

#define P_TEMPLATE template<typename dataType, int PATCH_RADIUS, CPatchParams::eColourTypeSize ImType, typename imDataType>
#define CTPatch CPatch<dataType, PATCH_RADIUS, ImType, imDataType>

//#define DEBUG_PATCH
#ifdef DEBUG_PATCH
#include "util/opencv_highgui.h"
#endif

P_TEMPLATE
class CPatch
{
	static const int ROT_BINS = 32; //Better than 24, much better than 64
public:
	static const int DIAMETER = 2*PATCH_RADIUS+1;
	static const int SIZE = DIAMETER*DIAMETER*(int)ImType;
    //typedef dataType elType;
protected:
	dataType acPatch[SIZE];

	int getRotBin(const int x, const int y)
	{
		//Only used to compute lookup table
		int nBin = ROT_BINS/2 + intFloor(ROT_BINS*atan2((double)y,(double)x)/(2*M_PI));
		if(nBin==ROT_BINS) nBin = 0;
		if(IS_DEBUG) CHECK(nBin<0 || nBin >= ROT_BINS, "CIplPx<imDataType>::getRotBin: Bin OOB");
		return nBin;
	}

	int * getBinLookup(int nScale)
	{
		const int AREA_RAD = PATCH_RADIUS * nScale;
		const int AREA_WIDTH = (DIAMETER+1) * nScale;
		int * anLookup = new int[sqr(AREA_WIDTH)];
		//int * pnLookup = anLookup;
		for(int y=-AREA_RAD; y <= AREA_RAD; y++)
		{
			int yi = y+AREA_RAD;

			for(int x=-AREA_RAD; x <= AREA_RAD; x++)
			{
				int xi = x+AREA_RAD;
				int nBin = getRotBin(x, y);
				anLookup[xi + AREA_WIDTH*yi] = nBin;
				//std::cout << nBin << '\t';
			}
			//std::cout << std::endl;
		}

		return anLookup;
	}

	/*double * getThreshLookup()
	{
		double * adLookup = new double[ROT_BINS/2];
		adLookup[0] = -1; //Should never be accessed
		for(int i=1; i<ROT_BINS/2; i++)
		{
			if(i==ROT_BINS/4)
				adLookup[i] = 0;
			else
			{
				double dAngle = i*2.0*M_PI/ROT_BINS;
				adLookup[i] = 1.0/tan(dAngle);
			}
		}
		return adLookup;
	}*/
	inline int getRotBinFast(const int x, const int y, const int nScale)
	{
		static int s_nScale = nScale;
		static int * anBins = getBinLookup(nScale);
		if(s_nScale != nScale)
		{
			std::cout << "Patch scale changing...\n";
			s_nScale = nScale;
			delete anBins; anBins = getBinLookup(nScale);
		}
		const int AREA_RAD = PATCH_RADIUS * nScale;
		const int AREA_WIDTH = (DIAMETER+1) * nScale;
		int yi = y+AREA_RAD;
		int xi = x+AREA_RAD;
		int nBin = anBins[xi + AREA_WIDTH*yi];
		if(IS_DEBUG) CHECK(nBin !=  getRotBin(x, y), "getRotBinFast failed");
		return nBin;
	}
/*	int getRotBinFast(const int x, const int y)
	{
		//return getRotBin(x, y);
		static const double * adThresh = getThreshLookup();
		//const double * adThresh2 = adThresh;
		int nBin=-1;
		if(y == 0)
		{
			nBin = (x>0) ? ROT_BINS/2 : 0;
		}
		else
		{
			for(int i=1; i<ROT_BINS/2; i++)
			{
				if(adThresh[i] * y < x)
				{
					std::cout << adThresh[i] - x/(double)y << std::endl;
					nBin = y<0 ? i : i+ROT_BINS/2;
					break;
				}
			}
			if(nBin == -1)
				nBin = (y < 0) ? 0 : ROT_BINS/2;
		}
#ifdef _DEBUG
		int nCorrectBin = getRotBin(x, y);
		if(nBin != nCorrectBin)
		{
			std::cout << x << ',' << y << " = " << nCorrectBin << " (interp=" << nBin << std::endl;
			//CHECK(1, "getRotBinFast failed");
		}
#endif
		return nBin;
	}*/

#define DOUBLE_TO_UCHAR 255.999
	inline dataType toUChar(int v) {return (dataType)v;};
	inline dataType toUChar(unsigned char v) {return (dataType)v;};
	inline dataType toUChar(float v) {return (dataType)(DOUBLE_TO_UCHAR*v);};
	inline dataType toUChar(double v) {return (dataType)(DOUBLE_TO_UCHAR*v);};
	inline int toInt(unsigned char v) {return (int)v;};
	inline int toInt(int v) {return (int)v;};
	inline int toInt(float v) {return (int)(DOUBLE_TO_UCHAR*v);};
	inline double toDouble(unsigned char v) {return (double)v;};
	inline double toDouble(int v) {return (double)v;};
	inline double toDouble(float v) {return (double)(DOUBLE_TO_UCHAR*v);};
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
	//virtual void setupNorm()=0;


	template<bool USE_PATCH_SCALE>
	double getAngle(const IplImage * pIm, CLocation point, const int PATCH_SCALE)
	{
		//Oriented patch
		int anSum[ROT_BINS], anCount[ROT_BINS], anResponse[ROT_BINS];
		for(int i=0; i<ROT_BINS; i++) { anSum[i] = 0; anCount[i] = 0; }


		const int PATCH_STEP = USE_PATCH_SCALE && (PATCH_SCALE>2) ? (PATCH_SCALE/2) : 1;

		const int PATCH_WIDTH = PATCH_RADIUS*PATCH_SCALE - 1;
		for(int x=-PATCH_WIDTH; x <= PATCH_WIDTH; x += PATCH_STEP) //could += PATCH_SCALE??
		{

			for(int y=-PATCH_WIDTH; y <= PATCH_WIDTH; y += PATCH_STEP)
			{
				if(x*x+y*y > sqr(PATCH_WIDTH)) continue;
				int nRotBin = getRotBinFast(x, y, PATCH_SCALE);
				if(IS_DEBUG) CHECK(0 > nRotBin || nRotBin >= ROT_BINS, "CPatch: Bin OOB");
				int im_x = x + point.x();
				int im_y = y + point.y();
				if(ImType!=CPatchParams::GREY)
					anSum[nRotBin] += toInt(CIplPx<imDataType>::getRed(pIm, im_x,im_y)) + toInt(CIplPx<imDataType>::getGreen(pIm, im_x,im_y)) ;
				else
					anSum[nRotBin] += toInt(CIplPx<imDataType>::getGrey(pIm, im_x,im_y));

				anCount[nRotBin]++; //to enforce even coverage
			}
		}
		//Hack to compensate for 1st bin missing a few px: or anSum[0] = (anCount[1] * anSum[0]) / anCount[0];
		for(int i=0; i<ROT_BINS; i++)
			anSum[i] = (anCount[1] * anSum[i]) / std::max<int>(1, anCount[i]); //anCount[i]==0 occasionally if patch scale==1

		int nMaxBin = 0, nMaxBinVal = -1;
		for(int i=0; i<ROT_BINS; i++)
		{
			int nOppositeBin = (i+ROT_BINS/2) % ROT_BINS;
			if(IS_DEBUG) CHECK(nOppositeBin<0 || nOppositeBin >= ROT_BINS, "CPatch: Idx OOB");
			anResponse[i] = anSum[i] - anSum[nOppositeBin];
			if(nMaxBinVal < anResponse[i])
			{
				nMaxBinVal = anResponse[i];
				nMaxBin = i;
			}
		}

		//////// Optional: look at global resp. rather than just 1/4 of the patch
		int anResponse2[ROT_BINS];
		nMaxBin = 0; nMaxBinVal = -1;
		for(int i=0;i<ROT_BINS;i++)
		{
			int nPrevBin = (i==0) ? (ROT_BINS-1) : (i-1);
			int nNextBin = (i+1) % ROT_BINS;
			if(IS_DEBUG) CHECK(nNextBin<0 || nNextBin >= ROT_BINS || nPrevBin<0 || nPrevBin >= ROT_BINS, "CPatch: Bin Idx OOB");
			anResponse2[i] = anResponse[nPrevBin] + (5*anResponse[i])/3 + anResponse[nNextBin];
			if(nMaxBinVal < anResponse2[i])
			{
				nMaxBinVal = anResponse2[i];
				nMaxBin = i;
			}
		}
		////////////// This extra pass is worth it with 32 bins, another pass is not.
		nMaxBin = 0; nMaxBinVal = -1;
		for(int i=0;i<ROT_BINS;i++)
		{
			int nPrevBin = (i==0) ? (ROT_BINS-1) : (i-1);
			int nNextBin = (i+1) % ROT_BINS;
			if(IS_DEBUG) CHECK(nNextBin<0 || nNextBin >= ROT_BINS || nPrevBin<0 || nPrevBin >= ROT_BINS, "CPatch: Bin Idx OOB");
			anResponse[i] = anResponse2[nPrevBin] + 2*anResponse2[i] + anResponse2[nNextBin];
			if(nMaxBinVal < anResponse[i])
			{
				nMaxBinVal = anResponse[i];
				nMaxBin = i;
			}
		}
		//////////////

		//Todo: Extrapolate to find peak : getmin := (y0, y1, y2) -> 1/2*(-y2+y0)/(-2*y1+y2+y0);
		double dAngle = nMaxBin * (2*M_PI/ROT_BINS); //sign is correct, tested 14-12-09
		return dAngle;
	}

public:
	CPatch(const IplImage * pIm, CLocation point, const CPatchParams & PATCH_PARAMS) HOT;

	//For drawing descriptors
	uchar val(int x, int y, int nChannel = 0) const
	{
		if(IS_DEBUG) CHECK(!(x<DIAMETER && y < DIAMETER && x>=0 && y>=0), "val:idx oob");
		const int idx = x*DIAMETER + y;
		return acPatch[ImType * idx + nChannel];
	}
};

//P_TEMPLATE
//CTPatch::getAngle(const IplImage * pIm, CLocation point)

P_TEMPLATE
CTPatch::CPatch(const IplImage * pIm, CLocation point, const CPatchParams & PATCH_PARAMS)
{
#ifdef DEBUG_PATCH
	cvNamedWindow("OriPatch");
#endif
	const int PATCH_SCALE = PATCH_PARAMS.PATCH_SCALE;

	const int x_min=point.x()-PATCH_SCALE*PATCH_RADIUS;
	const int y_min=point.y()-PATCH_SCALE*PATCH_RADIUS;
	const int x_max=point.x()+PATCH_SCALE*PATCH_RADIUS;
	const int y_max=point.y()+PATCH_SCALE*PATCH_RADIUS;
	const int step = PATCH_SCALE;

	if(!PATCH_PARAMS.ORIENT)
	{
		int idx=0;

		for(int x=x_min; x <= x_max; x += step)
		{
			for(int y=y_min; y <= y_max; y += step)
			{
				if(ImType==CPatchParams::GREY)
					acPatch[idx++] = toUChar(CIplPx<imDataType>::getGrey(pIm, x,y));
				else
				{
					if(ImType==CPatchParams::RGB)
						acPatch[idx++] = toUChar(CIplPx<imDataType>::getBlue(pIm, x,y));
					acPatch[idx++] = toUChar(CIplPx<imDataType>::getGreen(pIm, x,y));
					acPatch[idx++] = toUChar(CIplPx<imDataType>::getRed(pIm, x,y));
				}
			}
		}
	}
	else
	{
		//The same about 90% of time. If we blurred more would probably be closer...
		double dAngle = getAngle<true>(pIm, point, PATCH_SCALE);
		//double dAngle = getAngle<false>(pIm, point, PATCH_SCALE);

		//cout << "PATCH_SCALE angle: " << dAngle2 << "Old angle: " << dAngle << endl;

		double s=sin(dAngle), c=cos(dAngle);

		int idx=0;
		for(double x=-PATCH_RADIUS; x <= PATCH_RADIUS; x += 1)
		{
			for(double y=-PATCH_RADIUS; y <= PATCH_RADIUS; y += 1)
			{
				//Rotate backwards to a point
				double x_im = c*x - s*y; //Forward rotation by dAngle
				double y_im = s*x + c*y;
				x_im *= PATCH_SCALE;
				y_im *= PATCH_SCALE;

				//Find 4 points about this point and their weight
				double x_lo=intFloor(x_im);
				double y_lo=intFloor(y_im);
#ifdef _DEBUG
				if(y_im-y_lo<0 || x_im-x_lo < 0)
				{
			    	std::cout << "Margin = " << -1 << " PATCH_SCALE=" << PATCH_SCALE << " PATCH_RADIUS=" << PATCH_RADIUS << " PATCH_ORIENT=" << PATCH_PARAMS.ORIENT <<  std::endl;
					std::cout << "rotated " << x << ',' << y << " by " << dAngle << " to get " << x_im << ',' << y_im << std::endl;
				}
#endif

				if(IS_DEBUG) CHECK(y_im-y_lo<0 || x_im-x_lo < 0, "Margin error");
				double /*dTotalWeight = 0,*/ dB = 0, dG = 0, dR = 0, dI;
/*					for(double im_x0=x_lo; im_x0<=x_lo+1; im_x0++)
					for(double im_y0=y_lo; im_y0<=y_lo+1; im_y0++)
					{
						double dx=x_im-im_x0, dy=y_im-im_y0;
						double weight = 2 - (dx*dx + dy*dy); //todo: this is a dodge weight calc.

						int im_x = (int)im_x0 + point.x;
						int im_y = (int)im_y0 + point.y;
						dB += weight * toDouble(CIplPx<imDataType>::getBlue(pIm, im_x,im_y));
						if(ImType==RGB)
						{
							dG += weight * toDouble(CIplPx<imDataType>::getGreen(pIm, im_x,im_y));
							dR += weight * toDouble(CIplPx<imDataType>::getRed(pIm, im_x,im_y));
						}
						dTotalWeight += weight;
					}
				double dTotalWeight_inv = 1.0/dTotalWeight;*/
				double dx=x_im-x_lo, dy=y_im-y_lo;
				int im_x_lo = (int)x_lo + point.x();
				int im_y_lo = (int)y_lo + point.y();
				int im_x_hi = im_x_lo+1;
				int im_y_hi = im_y_lo+1;
				double coeff_ll = (1-dx)*(1-dy);
				double coeff_lh = (1-dx)*dy;
				double coeff_hl = dx*(1-dy);
				double coeff_hh = dx*dy;
#define BILIN(CHANNEL) CIplPx<imDataType>::get ## CHANNEL(pIm, im_x_lo,im_y_lo)*coeff_ll + CIplPx<imDataType>::get ## CHANNEL(pIm, im_x_hi,im_y_lo)*coeff_hl + CIplPx<imDataType>::get ## CHANNEL(pIm, im_x_lo,im_y_hi)*coeff_lh + CIplPx<imDataType>::get ## CHANNEL(pIm, im_x_hi,im_y_hi)*coeff_hh

				if(ImType==CPatchParams::GREY)
				{
					dI = BILIN(Grey);
					acPatch[idx++] = (dataType)(dI /* dTotalWeight_inv*/);
					if(IS_DEBUG) CHECK((int)dB > 256, "This will never trip but valgrind will detect if dB is uninit");
				} else
				{
					if(ImType==CPatchParams::RGB)
						dB = BILIN(Blue);
					dG = BILIN(Green);
					dR = BILIN(Red);

					if(ImType==CPatchParams::RGB)
						acPatch[idx++] = (dataType)(dB /* dTotalWeight_inv*/);
					acPatch[idx++] = (dataType)(dG /* dTotalWeight_inv*/); //Will break for a float image--not scaling to 0..255
					acPatch[idx++] = (dataType)(dR /* dTotalWeight_inv*/);
					if(IS_DEBUG) CHECK((int)dG < 0, "This will never trip but valgrind will detect if dG is uninit");
				}
			}
			abs((imDataType)1.0); //should break for float images--they need a scaling 3 lines higher
		}

		/*std::cout << "Marking gradient\n";
		CvPoint p1 = cvPoint(point.x(), point.y()), p2=p1; p2.x += (int)(10*c); p2.y += (int)(10*s);
		cvLine(const_cast<IplImage *>(pIm), p1, p2, CV_RGB(255,255,255)); //This is correctly marking a line in direction of darkest change
		cvCircle(const_cast<IplImage *>(pIm), p1, 2, CV_RGB(0,255,0));*/

		if(IS_DEBUG) CHECK(idx != SIZE, "Check all bins are filled");
	}
	if(PATCH_PARAMS.NORMALISE == 1)
	{
		//static const int x = 1/(ImType-2);

		static const int AREA=DIAMETER*DIAMETER;
		static const int PATCH_MAX_VAL = 256;

		int nTotal = 0;
		dataType * pcPatch = acPatch;

		for(int nBin=AREA; nBin>0; nBin--)
		{
			if(ImType==CPatchParams::GREY)
				nTotal += (int)*pcPatch;
			else if(ImType==CPatchParams::RGB)
			{
				nTotal += (int)*(pcPatch+1);
				nTotal += (int)*(pcPatch+2);
			}
			else if(ImType==CPatchParams::RG)
			{
				nTotal += (int)*(pcPatch);
				nTotal += (int)*(pcPatch+1);
			}
			pcPatch += ImType;
		}

		static const double PATCH_AREA_INV = (ImType != CPatchParams::GREY ? 2.0 : 1.0)/((double)(AREA));

		const int nShift = PATCH_MAX_VAL/2 - (int)(PATCH_AREA_INV*nTotal);

		pcPatch = acPatch;

		if(nShift > 0)
			for(int nBin=AREA * (int)ImType; nBin>0; nBin--)
			{
				int newVal = (int)*pcPatch + nShift;
				if(newVal >= PATCH_MAX_VAL) newVal = PATCH_MAX_VAL-1;
				*pcPatch = newVal;
				pcPatch++;
			}
		else if(nShift < 0)
			for(int nBin=AREA * (int)ImType; nBin>0; nBin--)
			{
				int newVal = (int)*pcPatch + nShift;
				if(newVal < 0) newVal = 0;
				*pcPatch = newVal;
				pcPatch++;
			};
	}
	else if(false)//PATCH_PARAMS.NORMALISE == 2)
	{
		//static const int x = 1/(ImType-2);

		static const int AREA=DIAMETER*DIAMETER;

		int nTotal = 0;
		int acPatchHSI[SIZE]; //Copy ISH to here
		int * pnPatchHSI = acPatchHSI;
		dataType * pcPatch = acPatch;

		for(int nBin=AREA; nBin>0; nBin--)
		{
			if(ImType==CPatchParams::RGB)
			{
				//Convert to HSI so we can normalise intensity (should we convert back to avoid cyclic stuff?)
				int nHue,nSaturation,nIntensity;
				CHsiInt::HSI (*(pcPatch+2),*(pcPatch+1),*pcPatch,nHue,nSaturation,nIntensity);
				*pnPatchHSI = nIntensity;
				*(pnPatchHSI+1) = nSaturation;
				*(pnPatchHSI+2) = nHue;
			}
			else
			{
				*pnPatchHSI = (int)*pcPatch;
				if(ImType==CPatchParams::RG) //A bit hacky--essentially normalising R and G
					*(pnPatchHSI + 1) = (int)*(pcPatch + 1);
			}
			nTotal += *pnPatchHSI;
			if(ImType==CPatchParams::RG)
				nTotal += *(pnPatchHSI + 1);

			pcPatch += ImType;
			pnPatchHSI += ImType;
		}

		static const double PATCH_AREA_INV = (ImType == CPatchParams::RG ? 2.0 : 1.0)/((double)(AREA));

		/*for(int nBin=AREA; nBin>0; nBin--)
		{
			nTotal += (int)(pcVal2);
			pcVal += ImType;
		}*/

		double dAvI = (double)nTotal*PATCH_AREA_INV;

		double dVarTotal = 0;
		pnPatchHSI = acPatchHSI;
		for(int nBin=AREA; nBin>0; nBin--)
		{
			double dDif = *pnPatchHSI - dAvI;
			dVarTotal += dDif*dDif;
			if(ImType==CPatchParams::RG)
			{
				dDif = *(pnPatchHSI+1) - dAvI;
				dVarTotal += dDif*dDif;
			}
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
			int nIntensityCorrected = doubleToInt(dShift+dScale*(*pnPatchHSI));
			if(nIntensityCorrected>=PATCH_MAX_VAL) nIntensityCorrected = PATCH_MAX_VAL-1;
			if(nIntensityCorrected<0) nIntensityCorrected = 0;

			if(ImType == CPatchParams::RG)
			{
				int nGreenCorrected = doubleToInt(dShift+dScale*(*(pnPatchHSI+1)));
				if(nGreenCorrected>=PATCH_MAX_VAL) nGreenCorrected = PATCH_MAX_VAL-1;
				if(nGreenCorrected<0) nGreenCorrected = 0;
			}
			//*pnPatchHSI = nIntensityCorrected;
			//Now copy back, or convert back to RGB

			if(ImType==CPatchParams::RGB)
			{
				//convert back to RGB to avoid cyclic stuff?)
				int nR,nG,nB;
				CHsiInt::RGB (*(pnPatchHSI+2),*(pnPatchHSI+1),*pnPatchHSI,nR,nG,nB);
				*pcPatch = toUChar(nB);
				*(pcPatch+1) = toUChar(nG);
				*(pcPatch+2) = toUChar(nR);
			}
			else
			{
				*pcPatch = toUChar(*pnPatchHSI);
				if(ImType == CPatchParams::RG)
					*(pcPatch + 1) = toUChar(*(pnPatchHSI+1));
			}

			pcPatch += ImType;
			pnPatchHSI += ImType;
		}
	}

	///////////
#ifdef DEBUG_PATCH
	Image<ImType==RG ? RGB : ImType, unsigned char> patch(DIAMETER, DIAMETER);
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


#define PN_TEMPLATE template<typename dataType, int PATCH_RADIUS, CPatchDescriptorParams::ePATCH_COMP_METHOD PATCH_COMP_METHOD, CPatchParams::eColourTypeSize ImType, typename imDataType>
#define CTPatchWithNorm CPatchWithNorm<dataType, PATCH_RADIUS, PATCH_COMP_METHOD, ImType, imDataType>

PN_TEMPLATE
class CPatchWithNorm : public CTPatch
{
protected:

public:
	CPatchWithNorm(const IplImage * pIm, CLocation point, const CPatchParams & PATCH_PARAMS) : CTPatch(pIm, point, PATCH_PARAMS) 
	{
		if(PATCH_COMP_METHOD == CPatchDescriptorParams::ePatchL1Parallel)
		{
			//Quantise patch to 0..127 for faster norms
			const int SIZE = CTPatch::SIZE;
			dataType * pcPatch = CTPatch::acPatch;
			for(int idx=SIZE; idx > 0; idx--, pcPatch++)
			{
				*pcPatch /= 2;
			}
		}
	}

	int distance(const CTPatchWithNorm * pPatch) const
	{
		const unsigned char * acPatch1 = CTPatch::acPatch;
		const unsigned char * acPatch2 = pPatch->acPatch;
		const int SIZE = CTPatch::SIZE;

		switch (PATCH_COMP_METHOD)
		{
		case CPatchDescriptorParams::ePatchEuclid:
			return euclidDist(acPatch1, acPatch2, SIZE);
		case CPatchDescriptorParams::ePatchMaxDist:
			return MaxDist(acPatch1, acPatch2, SIZE);
		case CPatchDescriptorParams::ePatchL1:
			return L1Dist(acPatch1, acPatch2, SIZE);
		case CPatchDescriptorParams::ePatchL1Fast:
		{
			const int relSubVecSize = 4; //Look at 1st 1/4 first
			const int subVecLen = SIZE/relSubVecSize;
			int dist = L1Dist(acPatch1, acPatch2, subVecLen);
			if(dist>70*subVecLen) return relSubVecSize*dist;
			return dist + L1Dist(acPatch1+subVecLen, acPatch2+subVecLen, SIZE-subVecLen);
		}
		case CPatchDescriptorParams::ePatchEuclidFast:
		{
			const int relSubVecSize = 4; //Look at 1st 1/4 first
			const int subVecLen = SIZE/relSubVecSize;
			int dist = euclidDist(acPatch1, acPatch2, subVecLen);
			if(dist>sqr(70)*subVecLen) return relSubVecSize*dist;
			return dist + euclidDist(acPatch1+subVecLen, acPatch2+subVecLen, SIZE-subVecLen);
		}
		case CPatchDescriptorParams::ePatchL1Parallel:
		{
			return L1distParallel(acPatch1, acPatch2, SIZE);
		}
		default:
			THROW("CTPatch::Compare: invalid comparison selected");
		/*case CPatchDescriptorParams::ePatchCos:
			return CompareCosine(pH);
	    case CPatchDescriptorParams::ePatchChiSquared:
			return CompareChiSquared(pH);
	    case CPatchDescriptorParams::ePatchJeffereys:
			return CompareJeffreys(pH);*/
	    //case CPatchDescriptorParams::ePatchCorrel:
		//	return correl(pPatch);
		}
	}

	/*int correl(const CPatch * pPatch) const
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
	uchar val(int x, int y, int nChannel = 0) const { return CTPatch::val(x,y,nChannel) * ((PATCH_COMP_METHOD == CPatchDescriptorParams::ePatchL1Parallel) ? 2 : 1); }
};

P_TEMPLATE
class CPatchWithNorm<dataType, PATCH_RADIUS, CPatchDescriptorParams::ePatchEuclidParallel, ImType, imDataType> : public CTPatch
{
	int nTwoTimesSumSquare;
protected:

public:
	CPatchWithNorm(const IplImage * pIm, CLocation point, const CPatchParams & PATCH_PARAMS) : CTPatch(pIm, point, PATCH_PARAMS) 
	{
		//Quantise patch to 0..127 for faster norms
		const int SIZE = CTPatch::SIZE;
		dataType * pcPatch = CTPatch::acPatch;
		for(int idx=SIZE; idx > 0; idx--, pcPatch++)
		{
			*pcPatch /= 2;
		}
		nTwoTimesSumSquare = 2 * sumSquares(CTPatch::acPatch, SIZE);
	}

	int distance(const CPatchWithNorm<dataType, PATCH_RADIUS, CPatchDescriptorParams::ePatchEuclidParallel, ImType, imDataType> * pPatch) const
	{
		const int SIZE = CTPatch::SIZE;
		const unsigned char * acPatch1 = CTPatch::acPatch;
		const unsigned char * acPatch2 = pPatch->acPatch;
		return L2distParallel(acPatch1, acPatch2, SIZE, nTwoTimesSumSquare, pPatch->nTwoTimesSumSquare);
	}

	uchar val(int x, int y, int nChannel = 0) const { return 2*CTPatch::val(x,y,nChannel); }
};
