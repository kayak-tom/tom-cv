/*Define a few concepts on top of the basic graph structures
 * */
#ifndef _GRAPHBASE_H
#define _GRAPHBASE_H

#include "g2o/core/base_unary_edge.h"
#include "g2o/core/base_binary_edge.h"
#include "g2o/core/base_multi_edge.h"
#include "g2o/core/base_vertex.h"

using namespace std;
using namespace g2o;

extern g2o::SparseOptimizer optimizer; //Make optimiser global for now (todo: singleton)


const double EXACT_INFO = 1000; // For relations which we assume are known exactly
Eigen::Matrix<double, 1, 1> EXACT_INFO_1D = Eigen::Matrix<double, 1, 1>::Constant(EXACT_INFO);
Eigen::Matrix<double, 1, 1> UNINFORMATIVE_INFO_1D = Eigen::Matrix<double, 1, 1>::Constant(1.0/EXACT_INFO);

#define CAST dynamic_cast

class CVertex1d : public BaseVertex<1, double>
{
protected:
    CVertex1d()
    {
        _estimate = 1;
        
        static int id=0;
        setId(id++);
        optimizer.addVertex(this);        
    }
    
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    virtual void oplus(double * update)
    {
        //cout << *update << endl;
        _estimate += *update;
    }

    virtual bool read(std::istream& is) {
        string label;
        is >> label;
        is >> _estimate;

        return true;
    }
    
    virtual bool write(std::ostream& os) const {
        os << "VERTEX " << _estimate << endl;
        return true;
    }
    
    virtual void setToOrigin() {
        _estimate=1; //Flow should never go to 0 (or below 'delta' really), as can then go negative easily
    }

    virtual bool setEstimateData(const double* est) {
        _estimate = *est;
        return true;
    }

    virtual bool getEstimateData(double* est) const {
        *est = _estimate;
        return true;
    }

    virtual int estimateDimension() const {
        return 1;
    }

    virtual bool setMinimalEstimateData(const double* est) {
        return setEstimateData(est);
    }

    virtual bool getMinimalEstimateData(double* est) const {
        return getEstimateData(est);
    }

    virtual int minimalEstimateDimension() const {
        return estimateDimension();
    }    
};



#endif