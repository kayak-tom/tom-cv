/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

#pragma once

#include "image/imageAccess.h"
#include "image/HSI.h"
#include "util/exception.h"
#include "cv.hpp"
//#define TIMEITH
#ifdef TIMEITH
#include "SpeedTest.h"
#endif
#include "util/opencv.h"

extern int nBinBlur; //so can mess with histogram setting from bow_test
extern int HIST_NORMALISE_I;

pragma_warning (push)
pragma_warning (disable: 4127) //const conditional expression

//#define HIST_BINS_H 12 //Weight by saturation??
//#define HIST_BINS_I 4
//#define HIST_RADIUS 4
//#define HIST_DIAMETER (2*HIST_RADIUS-1)
//#define HIST_ISTEP HIST_MAX_VAL/HIST_BINS_I
//#define HIST_HSTEP HIST_MAX_VAL/HIST_BINS_H
//
//#define HIST_BINS (HIST_BINS_H+HIST_BINS_I)
#define HIST_MAX_VAL 256 //Must be 256 cos that's the max val of a channel
#define HIST_SUM 127 //HIST_MAX_VAL// so can convert to un char without breaking (yuk)
//#define HIST_MAX_I 230
//#define HIST_MIN_I 30
//#define HIST_MIN_S 10
extern int HIST_MAX_I ;
extern int HIST_MIN_I ;
extern int HIST_MIN_S ;
extern int S_VS_H; //0..10 10=equal importance

enum eHistCompMethod {eHistEuclid, eHistL1, eHistMaxDist, eHistCos, eHistChiSquared,eHistJeffereys,eHistBhattacharyya};

//Could do template class
#define H_TEMPLATE template<int HIST_BINS_H, int HIST_BINS_S, int HIST_BINS_I, int HIST_RADIUS, eHistCompMethod HIST_COMP_METHOD, bool HIST_CUMULATIVE>
#define CTHist CHist<HIST_BINS_H, HIST_BINS_S, HIST_BINS_I, HIST_RADIUS, HIST_COMP_METHOD, HIST_CUMULATIVE>

H_TEMPLATE
class CHist
{
protected:
	static const int HIST_BINS = HIST_BINS_H+HIST_BINS_I+HIST_BINS_S;
	// CHist is ONLY AN ARRAY--otherwise constructors need changed
	unsigned char Bins[HIST_BINS];
private:
	static const int HIST_DIAMETER = 2*HIST_RADIUS+1;
	static const int HIST_HSTEP = HIST_MAX_VAL/HIST_BINS_H;
	static const int HIST_SSTEP = (int)(HIST_MAX_VAL / (HIST_BINS_S+0.00001));
    static const int HIST_ISTEP = (int)(HIST_MAX_VAL / (HIST_BINS_I+0.00001)) ; // Avoid div by 0

	//Lookup val and return index of bin. Could fuzzify here rather than in comparison fn?
	static void HBinBlur(int nHue, int nSaturation, int nIntensity, double * adBins, double dWeight);
	static void SBinBlur(int nSaturation, double * adBins, double dWeight);
	static void IBinBlur(int nSaturation, double * adBins, double dWeight);
	static void BinBlur(int nVal, double * adBins, double dWeight, int nBinWidth, int nBins); //For S,I
    static void BinBlurWrap(int nVal, double * adBins, double dWeight, int nBinWidth, int nBins); //for H

	static void Normalise(double * adBins, int nBins);
	static void Accumulate(double * adBins, int nBins);

    static int HBin(int nHue, int nSaturation, int nIntensity);
	static int HBinFromIBin(int nHue, int nSaturation, int IBin);

	static int IBin(int nIntensity);
	static int HBin(int nHue);
	static int SBin(int nSaturation);

	static int Bin(int nVal, const int nStep);
    static const int NO_BIN = -1;

	static double * adBinWeights; //[HIST_DIAMETER*HIST_DIAMETER]; //gaussian
	static int bInit; // just so init is called
	static const int NO_BINS = -1;
	static double* Init();

	//Treat array as as one long histogram ATM
	int CompareEuclid(const CHist  * pH) const ;
	int CompareEuclidFast(const CHist  * pH) const ;
	int CompareMaxDist(const CHist  * pH) const ;
	int CompareL1(const CHist  * pH) const ;
	int CompareCosine(const CHist  * pH) const ;
	int CompareChiSquared(const CHist  * pH) const ;
	int CompareJeffreys(const CHist  * pH) const ;
	int CompareJeffreysNew(const CHist  * pH) const ;
	int CompareBhattacharyya(const CHist  * pH) const ;

    template<bool bWholeIm>
	void Construct(const ImageRGB * pIm, CvPoint point);

public:
	inline CHist(const char * adHist) ; //construct from char array (for dup- ing)
	inline CHist(const double * adHist) ; //construct from double array
	inline CHist(const ImageRGB * pIm) { CvPoint temp; Construct<true>(pIm, temp);  }; //whole image hist
	inline CHist(const ImageRGB * pIm, CvPoint point) ;//{ Construct<false>(pIm, point); }; //hist around patch
	CHist(const CHist & H);

	int Compare(const CHist * pH) const;
    inline int Bin(int n) const {
        return (int)Bins[n];
    };

};

H_TEMPLATE
CTHist::CHist(const ImageRGB * pIm, CvPoint point)
{
#ifdef TIMEITH
    CStopWatch s;
    s.startTimer();
#endif

    Construct<false>(pIm, point);  //hist around patch

#ifdef TIMEITH
    s.stopTimer();
    static int nTimeItCount = 10;
    if(nTimeItCount>0)
    {
        nTimeItCount--;
        double dHistTime = s.getElapsedTime();
        cout << dHistTime << "=time per histogram\n";
    }
#endif
}

H_TEMPLATE
double* CTHist::adBinWeights = CTHist::Init();

// constructor for either whole im or a patch
H_TEMPLATE template<bool bWholeIm>
void CTHist::Construct(const ImageRGB * pIm, CvPoint point)
{
	double adBins[HIST_BINS];
	double * adHueBins = adBins;
	double * adSaturationBins = adBins+HIST_BINS_H;
	double * adIntensityBins = adBins+HIST_BINS_H+HIST_BINS_S;
	memset(adBins, 0, HIST_BINS*sizeof(double));

	int x_min, y_min, x_max, y_max, HIST_PIX_WEIGHT = 0, step;

    if(bWholeIm)
	{
		x_min = y_min = 0;
		x_max = pIm->getWidth() - 1;
		y_max = pIm->getHeight() - 1;
		step = HIST_RADIUS;
		HIST_PIX_WEIGHT = step*step*HIST_SUM / (pIm->getWidth() * pIm->getHeight());
	}
	else
	{
		x_min=point.x-HIST_RADIUS;
		y_min=point.y-HIST_RADIUS;
		x_max=point.x+HIST_RADIUS;
		y_max=point.y+HIST_RADIUS;
		step = 1;
	}

	int idx=0;
    int anIntensity[HIST_DIAMETER*HIST_DIAMETER];

	for(int x=x_min; x <= x_max; x += step)
	{
		for(int y=y_min; y <= y_max; y += step)
		{
        	int nHue,nSaturation,nIntensity;
            CHsiInt::HSI (pIm->r(x,y),pIm->g(x,y),pIm->b(x,y),nHue,nSaturation,nIntensity);

            if(HIST_BINS_I)
                anIntensity[idx] = nIntensity;

			double dToAdd = bWholeIm ? HIST_PIX_WEIGHT : adBinWeights[idx];
			//Now bin, weight by gaussian:

            if(nBinBlur)
            {
			    if(HIST_BINS_H)
					HBinBlur(nHue, nSaturation, nIntensity, adHueBins, dToAdd);

			    if(HIST_BINS_S)
					SBinBlur(nSaturation, adSaturationBins, dToAdd);

				if(HIST_BINS_I && (bWholeIm || !HIST_NORMALISE_I))
					IBinBlur(nIntensity, adIntensityBins, dToAdd);

                if(y_max>480)
                    y_max=480;

				if(IS_DEBUG) CHECK(y_max>480, "ERROR: Write OOB");

            } else { //1 bin, simple
                int nHBin;
                //Find bin:
			    if(HIST_BINS_I)
			    {
			        int nIBin = IBin(nIntensity);
                    if(IS_DEBUG) CHECK(nIBin>=HIST_BINS_I, "ERROR: Index OOB");
				    
                    if(bWholeIm || !HIST_NORMALISE_I)
                        adIntensityBins[nIBin] += dToAdd;

                    nHBin = HBinFromIBin(nHue, nSaturation, nIBin);
			    }
                else
			        nHBin = HBin(nHue, nSaturation, nIntensity);

			    //Add to histogram
			    if (nHBin != NO_BIN)
			    {
				    adHueBins[nHBin] += dToAdd;
                    if(IS_DEBUG) CHECK(nHBin>=HIST_BINS_H, "ERROR: Index OOB");
				    //if(HIST_CUMULATIVE && nHBin>0)
				    //	adIntensityBins[nHBin] += adIntensityBins[nHBin-1];
			    }

			    if(HIST_BINS_S)
                {
			        int nSBin = SBin(nSaturation);
                    if(IS_DEBUG) CHECK(nSBin>=HIST_BINS_S, "ERROR: Index OOB");
				    adSaturationBins[nSBin] += dToAdd;
                }
            }
			idx++;
		}
	}
    if(HIST_BINS_I && !bWholeIm && HIST_NORMALISE_I) //Normalise intensity
    { 
        static const int HIST_AREA=HIST_DIAMETER*HIST_DIAMETER;
        //Find mean and var:
        int nTotal = 0;
        int * pnVal = anIntensity;
        static const double HIST_AREA_INV = 1.0/((double)(HIST_AREA));
    	
	    for(int nBin=HIST_AREA; nBin>0; nBin--)
        {
		    nTotal += *pnVal;
            pnVal++;
        }

        double dAvI = (double)nTotal*HIST_AREA_INV;

        double dVarTotal = 0;
        pnVal = anIntensity;
	    for(int nBin=HIST_AREA; nBin>0; nBin--)
        {
            double dDif = *pnVal - dAvI;
		    dVarTotal += dDif*dDif;
            pnVal++;
        }

        double dSD = sqrt(dVarTotal*HIST_AREA_INV);

        const double dTargetMean = 128;
        const double dTargetSD = 40;
        const double dScale = dTargetSD/dSD;
        const double dShift = dTargetMean-dScale*dAvI;

        //bin, correct intensities this time
	    for(int idx=0; idx < HIST_DIAMETER*HIST_DIAMETER; idx ++)
	    {
		    double dToAdd = adBinWeights[idx];
		    
            //Normalise
            int nIntensityCorrected = doubleToInt(dShift+dScale*anIntensity[idx]);
            if(nIntensityCorrected>=HIST_MAX_VAL) nIntensityCorrected = HIST_MAX_VAL-1;
            if(nIntensityCorrected<0) nIntensityCorrected = 0;

            //Now bin, weight by gaussian:
            if(!nBinBlur)
            {
                int nIBin = IBin(nIntensityCorrected);
                adIntensityBins[nIBin] += dToAdd;
            }
            else
		        IBinBlur(nIntensityCorrected, adIntensityBins, dToAdd);
		}
    }

    //Normalise if cos dist
    if(HIST_COMP_METHOD == eHistCos)
    {
		Normalise(adBins, HIST_BINS);
    }

    //Accumulate if cumulative. NB cos dist won't work if cumulative.
    if(HIST_CUMULATIVE)
    {
		Accumulate(adHueBins, HIST_BINS_H);
		Accumulate(adSaturationBins, HIST_BINS_S);
		Accumulate(adIntensityBins, HIST_BINS_I);
    }

	//Now round bins to byte
	for(int nBin=0; nBin<HIST_BINS; nBin++)
		Bins[nBin] = (unsigned char)cvRound(adBins[nBin]);
}

H_TEMPLATE //not for cosine dist
inline void CTHist::Accumulate(double * adBins, int nBins)
{
//	for(int nBin=1; nBin<nBins; nBin++)
//        adBins[nBin] += adBins[nBin-1]; //todo pointer arith
    double * pdBins = adBins;
	for(int nBin=nBins-1; nBin>0; nBin--)
    {
        double dBinTotal = *pdBins;
        pdBins++;
        *pdBins += dBinTotal; //todo check!
    }
}

H_TEMPLATE //for cosine dist only
inline void CTHist::Normalise(double * adBins, int nBins)
{
    double dTotal = 0;
	
	for(int nBin=0; nBin<nBins; nBin++)
    {
        double dtemp = adBins[nBin];
		dTotal += dtemp*dtemp;
    }
    double dNormalisingFac = HIST_SUM/dTotal;
	for(int nBin=0; nBin<HIST_BINS_H; nBin++)
        adBins[nBin] *= dNormalisingFac;
}

/*H_TEMPLATE //for intensity only
inline void CTHist::NormaliseIntensity(double * adBins)
{
    double dTotal = 0;
    double * pdBin = adBins;
    static const double HIST_BINS_I_INV = 1.0/(double)HIST_BINS_I;
	
	for(int nBin=HIST_BINS_I; nBin>0; nBin--)
    {
		dTotal += *pdBin;
        pdBin++;
    }

    double dAv = dTotal*HIST_BINS_I_INV;

    double dVarTotal = 0;
    pdBin = adBins;
	
	for(int nBin=HIST_BINS_I; nBin>0; nBin--)
    {
        double dDif = *pdBin - dAv;
		dVarTotal += dDif*dDif;
        pdBin++;
    }

    double dVar = dVarTotal*HIST_BINS_I_INV;

    double dNormalisingFac = HIST_SUM/dTotal;
	for(int nBin=0; nBin<HIST_BINS_H; nBin++)
        adBins[nBin] *= dNormalisingFac;
}*/


H_TEMPLATE
inline void CTHist::HBinBlur(int nHue, int nSaturation, int nIntensity, double * adBins, double dWeight)
{
	if (nIntensity > HIST_MAX_I
	 || nIntensity < HIST_MIN_I
	 || nSaturation < HIST_MIN_S) return;

    BinBlurWrap(nHue, adBins, dWeight, HIST_HSTEP, HIST_BINS_H);
}

H_TEMPLATE
inline void CTHist::SBinBlur(int nSaturation, double * adBins, double dWeight)
{
    BinBlur(nSaturation, adBins, 0.1*S_VS_H*dWeight, HIST_SSTEP, HIST_BINS_S);
}

H_TEMPLATE
inline void CTHist::IBinBlur(int nIntensity, double * adBins, double dWeight)
{
    BinBlur(nIntensity, adBins, dWeight, HIST_ISTEP, HIST_BINS_I);
}

H_TEMPLATE
void CTHist::BinBlurWrap(int nVal, double * adBins, double dWeight, int nBinWidth, int nBins)
{
    const int bigStep = nBinWidth - 2*nBinBlur;
    const int smallStep = 2*nBinBlur;
    if(IS_DEBUG) CHECK(bigStep <= 0 , "CTHist::BinBlur: nBinBlur too big");
	int nThresh = HIST_MAX_VAL - nBinBlur;

    int nBin=0;

    for(;;)
    {
		if(nThresh < nVal || nBin == nBins) //Rounding errors mean we might not end up here until nBin>nBins otherwise
        {
            //between 2 bins nBin and nBin-1
            double dLastBinWeight = dWeight * ((nVal - nThresh ) / smallStep);
            double dThisBinWeight = dWeight - dLastBinWeight;

            int nThisBin = nBin % nBins;
            adBins[nThisBin] += dThisBinWeight;

            int nLastBin = nThisBin ? (nThisBin-1) : (nBins-1);
            adBins[nLastBin] += dLastBinWeight;

            return;
        }
		nThresh -= bigStep;

        if(nThresh < nVal)
        {
            //in 1 bin
            adBins[nBin] += dWeight;
            return;
        }
		nThresh -= smallStep;
        nBin++;
        if(IS_DEBUG) CHECK(nBin > nBins, "About to write OOB");
    }
}

H_TEMPLATE
void CTHist::BinBlur(int nVal, double * adBins, double dWeight, int nBinWidth, int nBins)
{
    const int bigStep = nBinWidth - 2*nBinBlur;
    const int smallStep = 2*nBinBlur;
    if(IS_DEBUG) CHECK(bigStep <= 0 , "CTHist::BinBlur: nBinBlur too big");
	int nThresh = HIST_MAX_VAL - (nBinBlur + bigStep); //Start off in a big bin, not halfway between endpoints

    int nBin=0;

    for(;;)
    {
        if(nThresh < nVal || nBin == nBins-1) //Rounding errors mean we might not end up here until nBin>nBins otherwise
        {
            //in 1 bin
            adBins[nBin] += dWeight;
            return;
        }
		nThresh -= smallStep;
        nBin++;

		if(nThresh < nVal) 
        {
            //between 2 bins nBin and nBin-1
            double dLastBinWeight = dWeight * ((nVal - nThresh ) / smallStep);
            double dThisBinWeight = dWeight - dLastBinWeight;

            adBins[nBin] += dThisBinWeight;
            adBins[nBin-1] += dLastBinWeight;
            return;
        }
		nThresh -= bigStep;
    }
}

//copy constructor
H_TEMPLATE
CTHist::CHist(const CHist & H)
{
	memcpy(Bins, H, sizeof(CHist));
}

//constructor from double arr
H_TEMPLATE
CTHist::CHist(const double * adHist)
{
	arrayCopy(adHist, Bins, HIST_BINS);
}

//constructor from char arr (for dup- ing)
H_TEMPLATE
CTHist::CHist(const char * acHist)
{
	arrayCopy((char *)acHist, (char *)Bins, HIST_BINS);
}

H_TEMPLATE
inline int CTHist::HBin(int nHue, int nSaturation, int nIntensity)
{
	if (nIntensity > HIST_MAX_I
	 || nIntensity < HIST_MIN_I
	 || nSaturation < HIST_MIN_S) return NO_BIN;

	return HBin(nHue);
}

H_TEMPLATE
inline int CTHist::HBinFromIBin(int nHue, int nSaturation, int nIBin)
{
	if (nIBin == HIST_BINS_I - 1
	 || nIBin == 0
	 || nSaturation < HIST_MIN_S) return NO_BIN;

	return HBin(nHue);
}

H_TEMPLATE
inline int CTHist::IBin(int nIntensity)
{
	return Bin(nIntensity, HIST_ISTEP);
}

H_TEMPLATE
inline int CTHist::HBin(int nHue)
{
	return Bin(nHue, HIST_HSTEP);
}

H_TEMPLATE
inline int CTHist::SBin(int nSaturation)
{
	return Bin(nSaturation, HIST_SSTEP+1); //Todo: this is a hack for when num of bins doens't divide 256
}

H_TEMPLATE
inline int CTHist::Bin(int nVal, const int nStep)
{
	int nThresh = HIST_MAX_VAL - nStep;

	for(int nBin=0; ; nBin++)
	{
		if(nThresh <= nVal) return nBin;
		nThresh -= nStep;
	}
}

H_TEMPLATE
int CTHist::Compare(const CHist * pH) const
{
	switch (HIST_COMP_METHOD)
	{
	case eHistEuclid:
		return CompareEuclidFast(pH);
	case eHistMaxDist:
		return CompareMaxDist(pH);
	case eHistL1:
		return CompareL1(pH);
	case eHistCos:
		return CompareCosine(pH);
    case eHistChiSquared:
		return CompareChiSquared(pH);
    case eHistJeffereys:
		return CompareJeffreys(pH);
    case eHistBhattacharyya:
		return CompareBhattacharyya(pH);
	}
    throw new CException("CTHist::Compare: invalid comparison selected");
}

H_TEMPLATE
inline int CTHist::CompareJeffreysNew(const CHist * pH) const //NB This is positive, but we aren't really normalised
{
	double dDist = 0;
	for(int i=0; i<HIST_BINS; i++)
	{
        double p=Bin(i);
        double q=pH->Bin(i);
        double dLogp = (p*(p+q)>0) ? log(p/((p+q)/2)) : 0;
        double dLogq = (q*(p+q)>0) ? log(q/((p+q)/2)) : 0;

        dDist += p*dLogp+q*dLogq;
	}
    static int k=10;
    if(k>0)
    {
        k--;
        //char pc[50];
        //sprintf_s(pc,50,"%f=Jeffrey's %d=L1\n", dDist, CompareL1(pH));
        //cout<< pc;
    }
	return doubleToInt(dDist);
}

H_TEMPLATE
inline int CTHist::CompareJeffreys(const CHist * pH) const //NB This is positive, but we aren't really normalised
{
	double dDist = 0;
	for(int i=0; i<HIST_BINS; i++)
	{
        double p=Bin(i);
        double q=pH->Bin(i);
        double dLog = (p*q>0) ? log(p/q) : 0;

        dDist += p*dLog-q*dLog;
	}
    static int k=10;
    if(k>0)
    {
        k--;
        //char pc[50];
        //sprintf_s(pc,50,"%f=Jeffrey's\n", dDist);
        //cout<< pc;
    }
	return doubleToInt(0.5*dDist);
}

H_TEMPLATE
inline int CTHist::CompareBhattacharyya(const CHist * pH) const //NB This is positive, but we aren't really normalised
{
	double dDist = 1;
	for(int i=0; i<HIST_BINS; i++)
	{
        double p=Bin(i)/(double)HIST_SUM;
        double q=pH->Bin(i)/(double)HIST_SUM;
        double dSS = (p*q);
        //cout<< dSS <<' ';

        dDist -= sqrt(dSS);
	}
    dDist = (dDist>0) ? sqrt(dDist) : 0;

    static int k=100;
    if(k>0)
    {
        k--;
        //char pc[50];
        //sprintf_s(pc,50,"%f=Bhattacharyya's %d=L1\n", dDist, CompareL1(pH));
        //cout<< pc;
    }
	return doubleToInt(256*dDist);
}

H_TEMPLATE
inline int CTHist::CompareChiSquared(const CHist * pH) const //NB This is positive, but we aren't really normalised
{
	int nDist = 0;
	for(int i=0; i<HIST_BINS; i++)
	{
        int b1=Bin(i), b2=pH->Bin(i);
		int nRelDist =  b1-b2;
        nDist += nRelDist ? (nRelDist*nRelDist)/(b1+b2) : 0;
	}
	return nDist;
}

H_TEMPLATE
inline int CTHist::CompareCosine(const CHist * pH) const //NB This is positive, but we aren't really normalised
{
	int nDist = 0;//HIST_SUM;
	for(int i=0; i<HIST_BINS; i++)
	{
        nDist += (Bin(i)*Bin(i) + pH->Bin(i)*pH->Bin(i))/2;
		nDist -= Bin(i)*pH->Bin(i);
	}
	return nDist/HIST_SUM;
}

H_TEMPLATE
inline int CTHist::CompareEuclid(const CHist * pH) const
{
	int nDist = 0;
	for(int i=0; i<HIST_BINS; i++)
	{
		int nRelDist = Bin(i) - pH->Bin(i);
		nDist += nRelDist*nRelDist;
	}
    return doubleToInt(sqrt((double)nDist)); // or /HIST_SUM
	//return nDist/HIST_SUM;
}

H_TEMPLATE
inline int CTHist::CompareEuclidFast(const CHist * pH) const
{
	int nDist = 0;
    unsigned char * pcBinH1 = (unsigned char *)Bins;
    unsigned char * pcBinH2 = (unsigned char *)pH->Bins;

	for(int i=HIST_BINS; i>0; i--)
	{
        int nBin1 = (int)*pcBinH1;
        int nBin2 = (int)*pcBinH2;
		int nRelDist = nBin1 - nBin2;
		nDist += nRelDist*nRelDist;
        pcBinH1++;
        pcBinH2++;
	}
    return intLookup::Sqrt(nDist); // or /HIST_SUM
	//return nDist/HIST_SUM;
}

H_TEMPLATE
inline int CTHist::CompareL1(const CTHist * pH) const
{
	int nDist = 0;
	for(int i=0; i<HIST_BINS; i++)
	{
		int nRelDist = Bin(i) - pH->Bin(i);
		nDist += abs(nRelDist);
	}
	return nDist;
}

H_TEMPLATE
inline int CTHist::CompareMaxDist(const CTHist * pH) const
{
	int nMaxDist = 0;
	for(int i=0; i<HIST_BINS; i++)
	{
		int nRelDist = Bin(i) - pH->Bin(i);
		if(nMaxDist < nRelDist) nMaxDist = nRelDist;
	}
	return 2*nMaxDist;
}

//Setup the Gaussian
//nHue and nIntensity hists will then each sum to HIST_SUM
H_TEMPLATE
double* CTHist::Init()
{
    static double adBinWeights_temp[HIST_DIAMETER*HIST_DIAMETER];

    CvMat * pMat = cvCreateMat(1, HIST_DIAMETER, CV_64FC1);
    CvSepFilter::init_gaussian_kernel(pMat);

    int i=0;
	for(int r=0; r<HIST_DIAMETER; r++)
    {
	    for(int c=0; c<HIST_DIAMETER; c++)
        {
            //NB the (pMat->data.db)[r] * (pMat->data.db)[c]'s sum to 1 (!)
            adBinWeights_temp[i] = HIST_SUM * (pMat->data.db)[r] * (pMat->data.db)[c]; //HIST_SUM/(double)(HIST_DIAMETER*HIST_DIAMETER); //todo gaussian
            //cout << adBinWeights_temp[i] << ',';
            i++;
        }
    }


	return (double *)adBinWeights_temp;
}

pragma_warning (pop)
