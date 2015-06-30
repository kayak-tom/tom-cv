#ifndef _VECTOR_HASH_H
#define _VECTOR_HASH_H

#include <boost/functional/hash.hpp>

//Return double in -1...1
template<typename TVecType>
double vectorHash(const TVecType & x)
{
    const bool bVerbose = false;

    size_t ulHash = 0;
    for(int i=0;i<x.size();i++)
        boost::hash_combine<double>(ulHash, x(i));
        
    const size_t ulQuantisation = (size_t)12345678901UL; //1UL << 42;// std::numeric_limits<std::size_t>::max();// (size_t)-1;
    
    const double dHalfQuantisation = 0.5*ulQuantisation;
        
    const double dHash = (double)(ulHash % ulQuantisation)/dHalfQuantisation - 1.0; //-1 to 1 (yes it is ok that this will lose precision)

    COUT(x);
    COUT(ulHash);
    COUT(dHash);
    COUT(ulQuantisation);
    return dHash;
}

#endif
