/* 
 * Library methods manipulating std::vector
 * */

#ifndef VECTORUTIL_H
#define VECTORUTIL_H

#include <util/exception.h>
#include <vector>

/* Removes an element from a vector (throw if not present) */
template<typename T>
void vectorErase(std::vector<T> & vec, const T & t)
{
    for(int i=0; i<(int)vec.size(); i++) {
        if(vec[i]==t) {
            vec.erase(vec.begin() + i);
            return;
        }
    }
    THROW("Element not found in vector");
}

/*template<typename T>
void vectorErase(std::vector<T> & vec, T & t)
{
    vectorErase(vec, t);
}*/

/* Removes an element from a vector (throw if not present) */
template<typename T>
void vectorErase(std::vector<T *> & vec, const T * t)
{
    for(int i=0; i<(int)vec.size(); i++) {
        if(vec[i]==t) {
            vec.erase(vec.begin() + i);
            return;
        }
    }
    THROW("Element not found in vector");
}

double median(std::vector<double> & adNumbers);
void meanSD(std::vector<double> & adNumbers, double & dMean, double & dSD);

//MODIFIES aNumbers
double percentile(const int nPercentile, std::vector<double> & aNumbers);

#endif // VECTORUTIL_H
