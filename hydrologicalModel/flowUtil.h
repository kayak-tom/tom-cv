/* Drainage functions, etc.
 */
 
#ifndef _FLOWUTIL_H
#define _FLOWUTIL_H

inline double sqr(const double t) { return t*t; }

inline double pseudoHuber(const double x, const double t=0.01)
{
    const double dHuber = t * (sqrt( 1 + sqr(x/t) ) - 1);
    return dHuber;
}

//Better weight fn: linear above 0, logistic below.
// Logistic 4/(1+exp(-x)) has gradient 1, intercept 2 at 0
inline double logisticPlusLinear(const double x, const double t=0.01)
{
    const double L = 4/t;
    if(x>0)
        return x+0.5*t;
        
    return L/(1+exp(-t*x));
}
inline double logisticPlusLinear_inv(const double q, const double t=0.01)
{
    if(q<=0)
        throw "q OOB";
        
    const double L = 4/t;
    if(q>0.5*t)
        return q-0.5*t;
    
    return -log(L/q-1)/t;
}

inline double estimateToState(const double x, const double t=0.1)
{
    return logisticPlusLinear(x, t);
}

inline double stateToEstimate(const double q, const double t=0.1)
{
    return logisticPlusLinear_inv(q, t);
}


#endif
