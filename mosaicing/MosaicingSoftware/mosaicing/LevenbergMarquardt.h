#ifndef _GRC_FUNCTION_
#define _GRC_FUNCTION_

#include "TransformErrFunction.h"

namespace grc {

    //! Implements Levenberg-Marquardt iterative nonlinear optimisation.
    /*! Pre-allocate matrices for speed. */
	class LevenbergMarquardt
	{
        CvMat * J_, //!< Jacobian matrix
            * J_vec_, //!< Vector of partial derivatives
            * J_inv_, //!< Inverse of J_+I
            * delta_, //!< Innovation vector
            * gradient_vec_, //!< Direction to move to improve objective function
            * delta_scaled_;//!< Damped descent vector
        TransformErrFunction & errFn_; //!< Objective function to minimise
        int numParams_;//!< Number of parameters to refine
        const int maxLMIterations_;//!< Max number of iterations (needed to maintain real-time performance, and because exact convergence is unneccesary)

        //! Return biggest absolute value in a matrix
        static double maxAbs(const CvMat * M) ;

        //! Return biggest value on the diagonal of a matrix
        static double maxDiagVal(const CvMat * M) ;

        //! Return sum of matrix elements squared
        static double sumSquare(const CvMat * M) ;
        
        //! Return Frobeneus norm of a matrix (like Euclidian norm)
        static double frobeniusNorm(const CvMat * M) ;

        //! Add a damping factor u to the diagonal of a matrix
        static void addLMDampingFactor(CvMat * M, double u) ;

        //! Compute Jacobian matrix from errFn_
		void computeJacobian();
    public:
        //! Performs Levenberg-Marquardt optimnisation on errFn. errFn's parameters are adjusted as appropriate.
        /* \param maxLMIterations Maximum number of iterations (for maintaining real-time performance). */
        LevenbergMarquardt(TransformErrFunction &errFn, const int maxLMIterations);

        ~LevenbergMarquardt();
	};

}

#endif

