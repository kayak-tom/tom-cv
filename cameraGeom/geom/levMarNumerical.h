//Lev-Mar and Conjugate Gradient Descent optimisation, using numerical derivatives
#pragma once
#ifndef LEVMARNUMERICAL_H
#define LEVMARNUMERICAL_H

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <util/exception.h>
#include <boost/noncopyable.hpp>

void spy(const Eigen::SparseMatrix<double> & s);

#define PPMAT(M) pp(#M, M)
void pp(const std::string label, const Eigen::MatrixXd & J, const int nSubmatSize = 5);

class CLMFunction : boost::noncopyable {
    const bool bUseAnalyticDerivatives;
public:

    enum eLMSuccessStatus { eLMSuccess, eLMFail };

	/**
	 * @brief A penalty term for residuals that should prompt LM to step back (and not fail)
	 * @return 
	 */
	inline static const double HUGE_RESIDUAL() { return 1e+12; }
	
    virtual int inputs() const = 0;
    virtual int values() const = 0;

    /**
     * @brief Objective function to optimise
     * @param x Parameter vector (size inputs())
     * @param resids Residual vector to fill (size values())
     * @param bVerbose If CLevMar has a bVerbose flag set, then calls to 'function' are verbose iff not computing numerical derivatives.
     * @param nParamChanged If only one parameter has changed since this residual vector was calculated, nParamChanged is set to that parameter index (the function only needs to update relevent residuals). Otherwise nParamChanges = -1 
     * @return eLMSuccess (unless parameters are invalid, e.g. eLMFail, which will cause the optimisation to either step back or fail without converging. 
     */
    virtual eLMSuccessStatus function(const Eigen::VectorXd &x, Eigen::VectorXd &resids, bool bVerbose = false, const int nParamChanged = -1) = 0;

    double val(const Eigen::VectorXd &x, bool bVerbose = false) //Fast access to functions with 1 val
    {
        if(IS_DEBUG) CHECK(values() != 1, "This function should only be used when the objective function has a 1D output");
        Eigen::VectorXd resids(1);
        function(x, resids, bVerbose, -1);
        return resids(0);
    }
    
    virtual Eigen::VectorXd init() 
    {
        Eigen::VectorXd initParams = Eigen::VectorXd::Zero(inputs());
        return initParams; 
    }

    /*Overwrites 'residuals' vector*/
    double sumSquare(const Eigen::VectorXd &x, Eigen::VectorXd &residuals, bool bVerbose = false) 
    {
        //Eigen::VectorXd resids(values());
        function(x, residuals, bVerbose, -1);
        return residuals.squaredNorm();
    }

    virtual ~CLMFunction() {
    }

    CLMFunction(const bool bUseAnalyticDerivatives = false) : bUseAnalyticDerivatives(bUseAnalyticDerivatives) {
    }

    const bool useAnalyticDerivatives() const {
        return bUseAnalyticDerivatives;
    }

    virtual void analyticDeriv(const Eigen::VectorXd &/*x*/, Eigen::VectorXd &/*residuals*/, const int /*nParam*/) {
        THROW("Should either use numerical derivatives or overload me");
    }

//	static void setLargeResids(Eigen::VectorXd &resids, const int nParamChanged, const bool bVerbose);
};

//template<int DIMS = Eigen::Dynamic>
class CLevMar : boost::noncopyable {
    CLMFunction & function;

    static const int NUM_PARAMS = Eigen::Dynamic;
    static const int NUM_MODEL_VARS = Eigen::Dynamic;
    typedef Eigen::Matrix<double, NUM_PARAMS, 1 > TParamVector;
    typedef Eigen::Matrix<double, NUM_MODEL_VARS, 1 > TResidVector;
    typedef Eigen::Matrix<double, NUM_MODEL_VARS, NUM_PARAMS> TJMatrix;
    typedef Eigen::Matrix<double, NUM_PARAMS, NUM_PARAMS> TJTJMatrix;

    TParamVector paramUpdateVec;
    int nVerbose; //-1 = Nothing, even warnings, 0 = minimal, 1 = summary of each iteration, 2 = full details
    //bool bSuppressAllWarnings;
    double dMinDelta, dMinStepLength /*for both steps and reduction in error*/, pseudoHuber_t /* -1 == off */;

    int nLMIter; //Used for diagnostic output

public:

    enum eMethod {
        eLM, eLMSparse, eLMDiag, eCGD
    };
private:

    const eMethod method;
    bool bMakeSparseJ;

    //Move outside functions to avoid long param lists, and to avoid copying and reallocating
    TJMatrix J;
    TJTJMatrix JTJ;
    Eigen::VectorXd grad;
    Eigen::SparseMatrix<double> spJTJ;
    Eigen::SparseMatrix<double, Eigen::ColMajor> spJ;
    TResidVector residual;
    TParamVector delXn, gradXn_take_1;

    enum eIterState {
        eDescending, eConverged
    };
    struct SVal {
        double alpha, val;

        SVal(double alpha, double val) : alpha(alpha), val(val) {
        }

        SVal() : alpha(HUGE), val(HUGE) {
        }
    };
    
    void computeDerivatives(const Eigen::VectorXd &x);

    CLMFunction::eLMSuccessStatus computeDerivativesNumerically(const Eigen::VectorXd &x, const bool bForwards = true);

    static CLevMar::eIterState dampMarquardt(const double lambda, TParamVector & JTJ_diag);

    static CLevMar::eIterState dampMarquardt(const double lambda, Eigen::MatrixXd & JTJ);

    static CLevMar::eIterState dampMarquardt(const double lambda, Eigen::SparseMatrix<double> & JTJ);

    double optimiseGD(TParamVector & params, const bool bMinimise);
    
    CLMFunction::eLMSuccessStatus robustFunction(const Eigen::VectorXd & params, Eigen::VectorXd &resids, bool bVerbose, const int nParamChanged = -1);
    eIterState levMarIter(Eigen::VectorXd & params, double & lambda, const bool bMarquardt, double & dErr);

    //Also updates params, and evaluated residual at minimum

    eIterState minLineSearch(TParamVector & params, const TParamVector & delXn, double & alpha, double & dInitErr);

    void checkGrad(const TParamVector & grad) const;

    //Notation from http://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method

    eIterState CGDIter(TParamVector & params, const bool bReset, double & alpha, double & dErr);

    //Also updates params, and evaluated residual at minimum
    //See numerical recipes section 5.7
    void testLineSearchDelta(TParamVector & params, const TParamVector & delXn, double & dErr);
        

    eIterState minLineSearchNR(TParamVector & params, const TParamVector & delXn, double & dErr);

//////////// Debugging/analytics
    void testSmoothness(const TParamVector & params, const double dDelta, const int nParam);
    //Compute delta minimising the difference between forward and backward derivatives
    double normalisedDiff(const TParamVector & params);

public:
    CLevMar(CLMFunction & f, const bool bVerbose = false, const double dMinDelta = 2e-9, const eMethod method = eLM);

    double minimise(TParamVector & params, const int MAX_ITERS = 150);

	void setMinStepLength(const double dMinStepLength_in);
    void setPseudoHuberSD(const double t);
    void setVerbose(const int nVerbose_in) { nVerbose = nVerbose_in; }
    void setSuppressWarnings(const bool bSuppressAllWarnings_in) { nVerbose = (bSuppressAllWarnings_in ? -1 : 0); }
	
	

	//For debugging/profiling changes
	int getNumIters() const { return nLMIter; }

    const Eigen::VectorXd & residuals() const;
    
    const double RMSError() const;

    //Simple gradient descent to find a maximum
    double maximiseGD(TParamVector & params);

    //Simple gradient descent to find a minimum
    double minimiseGD(TParamVector & params);

    void testSmoothness(const TParamVector & params) ;
    //Compute delta minimising the difference between forward and backward derivatives
    double computeOptimalDelta(const TParamVector & params);
    
    static Eigen::VectorXd minimiseFunction(CLMFunction & lmFunction, const bool bVerbose = false);
};

class CLMIterLog
{
	std::string name;
	double dCalls, dIters;
public:
	CLMIterLog(const std::string & name) : name(name), dCalls(0), dIters(0) {}
	void log(const CLevMar & LM)
	{
		dCalls++;
		dIters +=	LM.getNumIters();
		cout << name << " LM calls=" << dCalls << " LM iters=" << dIters << " Iters per call=" << (dIters/dCalls) << endl;
	}
};

#endif
