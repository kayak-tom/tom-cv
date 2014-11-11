/*
 * scales.h
 *
 *  Created on: 14 Jan 2010
 *      Author: tom
 */
#pragma once

#ifndef SCALES_H_
#define SCALES_H_

#include "util/exception.h"
#include "util/convert.h"
#include <iostream>

class CLogNormal {
public:

    static void wikipediaParamEstimation(const double dMean, const double dVar, double & d, double & g_sq) {
        d = log(dMean) - 0.5 * log(1 + dVar / sqr(dMean));
        g_sq = log(1 + dVar / sqr(dMean));
    }

    static void medianParamEstimation(const double dMean, const double dVar, double & d, double & g_sq) {
        d = log(dMean);
        g_sq = log(1 + dVar / sqr(dMean)); //why not...
    }

};

class CScale;

class CRelScale {
    friend class CScale;
    double d, g_sq; //Parameters of the lognormal distribution of a *scale change* \delta_{ijk}
    //G>=0 so G<0 <==> uninit

    static inline const double UNINIT() {
        return -1;
    }
public:

    inline double get_d() const {
        return d;
    }

    inline double get_g_sq() const {
        return g_sq;
    }

    enum eSpecialScales {
        eUninformativeScale, eTooBad, eZeroScale /*==uninit*/
    };

    CRelScale(const double d, const double g_sq) : d(d), g_sq(g_sq) {
    }

    CRelScale(const CRelScale & s) : d(s.d), g_sq(s.g_sq) {
    }

    void operator=(const CRelScale & s) {
        d = s.d, g_sq = s.g_sq;
    }

    static const int ZERO_d = -1000;

    CRelScale(const eSpecialScales eScaleType, double dConstantSpeed) {
        switch (eScaleType) {
            case eUninformativeScale: //We have an orientation but no decent scale estimate
                /*//Something like mean 1, s.d. 10
                //Assuming D=0, var=100 gives G about 1.5
                g_sq=1.5;

                //Now want  exp(D+0.5*1.5) = 1 => D = -0.5 G_sq
                d=-0.5*g_sq;*/

                //CLogNormal::medianParamEstimation(dConstantSpeed, 25, d, g_sq);
                d = log(dConstantSpeed);
                g_sq = log(25.0); //Why not... Want all uninformative scales to be equally bad, and worse than scales from pure rotation
                break;

            case eTooBad:
                g_sq = UNINIT();
                break;

            case eZeroScale:
                g_sq = log(5.0); //bit meaningless, todo
                d = ZERO_d;
                break;

            default:
                THROW("CRelScale: Unhandled enum option")
        }
    }

    bool notTooBad() const {
        return g_sq > UNINIT();
    }

    CRelScale inverse() const {
        if (d <= ZERO_d) {
            THROW("Cannot really invert speed of zero");
            return CRelScale(d, g_sq);
        }
        return CRelScale(-d, g_sq); //Easy! All the right properties. EXCEPT WHEN speed is zero
    }
};

class CScale //Represents scale from SLAM or from OR
{
    friend std::ostream& operator<<(std::ostream& s, const CScale & scale);

    double D, G_sq; //Parameters of the lognormal distribution of a scale s_{ij}

    static const double UNINIT() {
        return -1;
    }

    static const double NewComponentVar() {
        return 5;
    }//Start new components with uninformative idea of scale

    //MLE of mean (biased--Shan-1998)

    double scaleMean() const {
        if(IS_DEBUG) CHECK(!hasScale(), "CScale::scaleMean: Don't have scale")
                double dScaleMean = exp(D + 0.5 * G_sq);
        //if(IS_DEBUG) CHECK(dScaleMean<=0 || dScaleMean > 20, "CScale::scaleMean: Bad scales returned");//Todo drop 20 at some stage...
        return dScaleMean;
    }

    //Equals geometric mean

    double scaleMedian() const {
        if(IS_DEBUG) CHECK(!hasScale(), "CScale::scaleMedian: Don't have scale")
                double dScaleMedian = exp(D);
        if(IS_DEBUG) CHECK(dScaleMedian < 0, "CScale::scaleMedian: Scale <= 0");


        const double MAX_SCALE=100;
        if (dScaleMedian > 100) {
            REPEAT(100, std::cout << "Warning: Scale drift overflow: " << dScaleMedian << " resetting to MAX_SCALE\n");
            dScaleMedian = MAX_SCALE;
        }

        return dScaleMedian;
    }
    static double s_MaxBadness;

public:

    CScale() : D(UNINIT()), G_sq(UNINIT()) {
    }

    CScale(double D, double G) : D(D), G_sq(G) {
        if (D > 10) {
            REPEAT(10, std::cout << "Warning: high log scale: " << D << endl);
        }
    }

    void setOrigin(bool bIsFirstComponent)// For new components allow uninformative origins where scale will come from SCORE instead
    {
        D = 0;
        G_sq = bIsFirstComponent ? 0 : NewComponentVar();
    }

    bool isOrigin() const {
        return (D == 0) && (G_sq == 0 || (G_sq == NewComponentVar()));
    }

    bool isNewComponentOrigin() const {
        return (D == 0) && (G_sq == NewComponentVar());
    }

    void setUninit() {
        D = UNINIT();
        G_sq = UNINIT();
    }

    bool hasScale() const {
        return G_sq > UNINIT();
    }

    CScale operator+(const CRelScale & s) const {
        if(IS_DEBUG) CHECK(!hasScale() || !s.notTooBad(), "CScale operator+: Don't have scale");

        const double dNew = D + s.d;

        /*if(dNew > 2) TODO: Parameterise me TODO: Parameterise me
        {
                REPEAT(10, std::cout << "Hack scale for video");
                dNew=2;
        }
        else if (dNew < -1)
        {
                REPEAT(10, std::cout << "Hack scale for video");
                dNew = -1;
        }*/

        const double gNew = G_sq + s.g_sq;
        const bool LIMIT_DRIFT = false;
        if ((fabs(dNew) < 10 || !LIMIT_DRIFT) || s.d == CRelScale::ZERO_d) {
            return CScale(dNew, gNew); //If s.d == CRelScale::ZERO_d this will essentially be zero.
        } else if (dNew < CRelScale::ZERO_d + 10) {
            std::cout << "Returning unit scale as are moving away from stationary\n";
            return CScale(0, gNew);
        } else {
            std::cout << "Returning invalid scale because clearly an outlier (dNew=" << dNew << ")\n";
            return CScale();
        }
    }

    double scaleEstimate() const {
        return scaleMedian();
    }

    //MLE of variance (biased--Shan-1998)

    double scaleVar() const {
        if(IS_DEBUG) CHECK(!hasScale(), "CScale::scaleVar: Don't have scale")
        return (exp(G_sq) - 1) * exp(2 * D + G_sq);
    }

    //Hacked 'variance of median' estimate

    double scaleVarHacked() const {
        if(IS_DEBUG) CHECK(!hasScale(), "CScale::scaleVar: Don't have scale")
        return (exp(G_sq) - 1) * exp(D);
    }

    //Geometric variance of geometric mean

    double scaleVarGeometric() const {
        if(IS_DEBUG) CHECK(!hasScale(), "CScale::scaleVar: Don't have scale")
        return exp(2 * sqrt(G_sq));
    }

    //G high==bad
    //G low==good
    /*double badness() const
    {
            if(IS_DEBUG) CHECK(!hasScale(), "CScale::badness: Too bad, should never be used")
            return G_sq;
    }*/

    //In [0,1]. For converting to a colour for display

    double relBadness() const {
        if (!hasScale())
            return 1;
        if (G_sq > s_MaxBadness) s_MaxBadness = G_sq;
        return G_sq / s_MaxBadness;
    }

    bool isMoreAccurateThan(const CScale & s) const {
        if (!hasScale())
            return false;
        else if (!s.hasScale())
            return true;
        else
            return G_sq < s.getG_sq();
    }

    /*inline bool operator<(const CScale & s) const {
    }

    inline bool operator>(const CScale & s) const {
        return s < *this;
    }*/

    inline bool operator==(const CScale & s) const {
        if (s.hasScale() && hasScale())
            return /*s.getD() == getD() && because direction not necessarily computed correctly */ s.getG_sq() == getG_sq();
        else
            return s.hasScale() == hasScale();
    }

    inline bool operator!=(const CScale & s) const {
        return !(s == *this);
    }

    inline double getG_sq() const {
        if(IS_DEBUG) CHECK(!hasScale(), "CScale::getG_sq: Too bad, should never be used")
        return G_sq;
    }

    inline double getD() const {
        if(IS_DEBUG) CHECK(!hasScale(), "CScale::getD: Too bad, should never be used")
        return D;
    }
};


std::ostream& operator<<(std::ostream& s, const CScale & scale);
std::ostream& operator<<(std::ostream& s, const CRelScale & relscale);

//class CFullRelPose
//C3dNormalisedPose
#endif /* SCALES_H_ */
