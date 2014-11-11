/* Code by Tom Botterill. Documentation and license at http://www.hilandtom.com/tombotterill/code */

//link less

#pragma once

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include<math.h>
#include<stdlib.h>
#include "exception.h"

#define MACHEPS 1.11e-16
class CRandom
{
public:
	static const bool bernoulli(const double d)
	{
		return Uniform(0.0,1.0) < d;
	}
	
    //! Return val ~ Unif[0, d]
    static double Uniform( double dmin, double dmax )
    {
        return dmin+Uniform(dmax-dmin);
    }

    static int Uniform( int nmin, int nmax )
    {
        return nmin+Uniform(nmax-nmin+1);
    }

    static double Uniform( double d )
    {
        static const double dRAND_MAX_INV = 1.0/(double)FASTRAND_MAX;

        return fastrand() * d * dRAND_MAX_INV;
    }
    /*static double Uniform( double d )
    {
        static const double dRAND_MAX_INV = 1.0/(double)RAND_MAX;

        return rand() * d * dRAND_MAX_INV;
    } */    

    //! Return val ~ Unif[0, n)
    static int Uniform( int n )
    {
#ifdef __GNUC__
        int r = rand() % n;
#else
        int r = (rand() * n) / (RAND_MAX+1);
#endif
        if(IS_DEBUG) CHECK(r<0 || r>=n, "Uniform: Random num OOB");
        return r;
    }

    //! Return val ~ Unif[0, n) LIMITED RANGE!
    static int UniformFast( int n )
    {
#ifdef __GNUC__
        return fastrand() % n;
#else
        return Uniform( n );
#endif
    };

    static int perturb( int n ) //returns 1 or -1 1/n of the time
    {
        if(n==0) return 0;
        int r = rand() % (2*n);
        if(r==0) return -1;
        if(r==1) return 1;
        return 0;
    };

    //! Return val ~ Unif[-2sd, 2sd] has standard dev. sd.
    static inline double Unif0( double sd )
    {
        return Uniform(4*sd) - 2*sd;
    };

    static inline void fast_srand( const unsigned int seed )
    {
        g_seed = seed;
    } //fastrand routine returns one integer, similar output value range as C lib.

private:
    static unsigned int g_seed;
public:
    static const unsigned int FASTRAND_MAX = 0x7FFF;
    //Used to seed the generator.

    static inline int fastrand()
    {
        g_seed = (214013*g_seed+2531011);
        return (g_seed>>16)&FASTRAND_MAX;
    }

//public:
    static inline double Uniform()
    {
        static const double dRAND_MAX_INV = 1.0/(double)FASTRAND_MAX;

        return fastrand() * dRAND_MAX_INV;
    };

    static double s_dNextNormal;
    static float s_fNextNormal;
    static const double NEXT_NORMAL_NOT_SET;
    static const float fNEXT_NORMAL_NOT_SET;

    static inline double Normal(const double dMean, const double dSD)
    {
        if(IS_DEBUG) CHECK(dSD < 0, "Normal: Negative SD")
        return Normal() * dSD + dMean;
    }

    static inline double Normal()
    {
        if(s_dNextNormal == NEXT_NORMAL_NOT_SET)
        {
            double U0=Uniform(1.0 - MACHEPS) + MACHEPS;
            double U1=Uniform(1.0 - MACHEPS) + MACHEPS;
            double S = sqrt(-2*log(U0));
            double Z0 = S*cos(2*M_PI*U1);
            double Z1 = S*sin(2*M_PI*U1);

            s_dNextNormal = Z1;
            return Z0;
        }
        else
        {
            double dNN = s_dNextNormal;
            s_dNextNormal = NEXT_NORMAL_NOT_SET;
            return dNN;
        }
    }
    
    //Use sum of a few uniform ints
    //NOT NORMALISED
    static inline double FasterNormalCLT()
    {
        const int NUM=4;
        int sum=0;
        for(int i=0;i<NUM;i++)
            sum += fastrand();
        
        static const int HALF_MAX = NUM*(FASTRAND_MAX/2);
        sum -= HALF_MAX;
        
        //static const double SUM_TO_NORMAL = 1.0/(double)HALF_MAX;
        return (double)sum; 
    }
    static inline double FasterNormal2()
    {
        const int NUM=4;
        int sum=0;
        for(int i=0;i<NUM;i++)
            sum += fastrand();
        
        static const int HALF_MAX = (NUM/2)*FASTRAND_MAX;
        sum -= HALF_MAX;
        
        //Each has var 1/12 FASTRAND_MAX^2, so sum has sd sqrt(NUM/12)*FASTRAND_MAX
        static const double SUM_TO_NORMAL = 1.0/(sqrt((double)NUM/12.0)*FASTRAND_MAX);
        return (double)sum * SUM_TO_NORMAL; 
    }
    
    static inline float FastNormal()
    {
        static const float fM_PI=(float)M_PI;
        if(s_fNextNormal == fNEXT_NORMAL_NOT_SET)
        {
            float U0=(float)(Uniform(1.0 - MACHEPS) + MACHEPS);
            float U1=(float)(Uniform(1.0 - MACHEPS) + MACHEPS);
            float S = sqrtf(-2.0f*logf(U0));
            float Z0 = S*cosf(2.f*fM_PI*U1);
            float Z1 = S*sinf(2.f*fM_PI*U1);

            s_fNextNormal = Z1;
            return Z0;
        }
        else
        {
            float fNN = s_fNextNormal;
            s_fNextNormal = fNEXT_NORMAL_NOT_SET;
            return fNN;
        }
    }
};

