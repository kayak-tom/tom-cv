#include "exception.h"
#include <boost/scoped_array.hpp>


/**
 * @class CAlignedVector
 * @brief Substitute for std::vector for vectors of aligned types, i.e. Eigen::Matrix
 */
template<typename T>
class CAlignedVector : boost::noncopyable
{
    int nSize;
    boost::scoped_array<T> aData;
public:
    CAlignedVector(const int nSize) : nSize(nSize), aData(new T[nSize]) {}
    T * begin() { return aData.get(); }
    T * end() { return aData.get() + nSize; }
    T & operator[](int n) { CHECK(n>=nSize || n < 0, "Idx OOB"); return aData[n]; }
    const T * begin() const { return aData.get(); }
    const T * end() const { return aData.get() + nSize; }
    const T & operator[](int n) const { CHECK(n>=nSize || n < 0, "Idx OOB"); return aData[n]; }
    
    const int size() const { return nSize; }
    
    void resize(const int nNewSize)
    {
        if(nNewSize == nSize)
            return;
            
        if(IS_DEBUG) CHECK(nSize>0, "Shouldn't resize vectors which already have a size");
        nSize=nNewSize;
        aData.reset(new T[nSize]);
    }
};