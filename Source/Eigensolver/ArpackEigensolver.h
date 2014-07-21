#ifndef ARPACKEIGENSOLVER_H_
#define ARPACKEIGENSOLVER_H_

#include "IterativeEigensolver.h"

#include "ArpackWrapper.h"

#include <Eigen/SparseCore>
#include <Eigen/UmfPackSupport>
#include <Eigen/PardisoSupport>

namespace Math
{

template<typename _MatrixType>
class ArpackSupport : public IterativeEigensolver
{
    typedef _MatrixType MatrixType;
    typedef typename MatrixType::Scalar Scalar;
    typedef typename Eigen::Matrix<Scalar,Eigen::Dynamic,1> VectorType;
    typedef typename MatrixType::RealScalar RealScalar;
    typedef typename std::complex<RealScalar> ComplexScalar;
    typedef typename Eigen::Matrix<ComplexScalar,Eigen::Dynamic,1> ComplexVector;
    typedef typename Eigen::Matrix<ComplexScalar,Eigen::Dynamic,Eigen::Dynamic> ComplexMatrix;

    const MatrixType* A;
    const MatrixType* M;

    Eigen::Matrix<ComplexScalar, Eigen::Dynamic, Eigen::Dynamic> m_eivec;
    ComplexVector m_eigenvalues;

    bool computeEigenvectors;
    bool m_isInitialized;
    bool m_eigenvectorsOk;

    Eigen::ComputationInfo m_info;
    size_t m_nbrConverged;
    size_t m_nbrIterations;

    std::string convertSearchTypeToString(SearchMode searchMode)
    {
        switch( searchMode)
        {
        case SearchMode::LargestMagnitude:
            return "LM";
        case SearchMode::SmallestMagnitude:
            return "SM";
        case SearchMode::LargestRealPart:
            return "LR";
        case SearchMode::SmallestRealPart:
            return "SR";
        case SearchMode::LargestImaginaryPart:
            return "LI";
        case SearchMode::SmallestImaginaryPart:
            return "SI";
        }
        /// We should never reach this!
        return "XX";
    }

    void extractEigenvalues(const Scalar sigma, const RealScalar* eigenvaluesReal, const RealScalar* eigenvaluesImag)
    {
        m_eigenvalues.resize(m_nbrConverged);

        for(unsigned int i=0; i<m_nbrConverged; i++)
        {
            m_eigenvalues[i] = sigma + 1. / ComplexScalar(eigenvaluesReal[i], eigenvaluesImag[i]);
        }
    }

    void extractEigenvalues(const Scalar sigma, const ComplexScalar* eigenvaluesReal, const ComplexScalar* eigenvaluesImag)
    {
        m_eigenvalues.resize(m_nbrConverged);

        for(unsigned int i=0; i<m_nbrConverged; i++)
        {
            m_eigenvalues[i] = sigma + 1. / eigenvaluesReal[i];
        }
    }

    void extractEigenvectors(const Scalar* Z)
    {
        const unsigned int n = A->rows();
        m_eivec.resize(n, m_nbrConverged);

        for (int i=0; i<m_nbrConverged; i++)
          for (int j=0; j<n; j++)
            m_eivec(j,i) = Z[i*n+j];
    }


public:

    /**
     * Constructor for a standard eigenvalue problem
     */
    ArpackSupport(const MatrixType& in_A, bool in_computeEigenvectors = false) : A(&in_A), M(NULL), computeEigenvectors(in_computeEigenvectors)
    {
      eigen_assert(A->rows()==A->cols() && "ArpackSupport: matrix A must be square");
    };

    /**
     * Constructor for a generalized eigenvalue problem
     */
    ArpackSupport(const MatrixType& in_A, const MatrixType& in_M, bool in_computeEigenvectors = false) : A(&in_A), M(&in_M), computeEigenvectors(in_computeEigenvectors)
    {
      eigen_assert(A->rows()==A->cols() && "ArpackSupport: matrix A must be square");
      eigen_assert(M->rows()==M->cols() && "ArpackSupport: matrix M must be square");
      eigen_assert(A->rows()==M->rows() && "ArpackSupport: The matrices A and M must be of the same size");
    };

    /**
     * Find eigenvalues of a generalized eigenvalue problem A*x = lambda*M*x which are closest to a specified target sigma
     */
    void compute(ComplexScalar sigma, int nev, Scalar* initialGuess = NULL)
    {
	RealScalar tol = 0.0;
	unsigned int maxIterations = 100;
        int n = A->rows(); 				///< Size of the eigenproblem.
        char bmat = 'I';				///< So far, we don't support M-orthogonalization!

        /// We use at least 20 Arnoldi vectors, but never more than n.
        int ncv = std::min( std::max(2*nev,20), n);
        maxIterations     = std::max(maxIterations, (unsigned int)2*n/std::max(ncv,1));

        int  lworkl = 3*ncv*ncv+5*ncv+1;    ///< Dimension of array workl.

        Scalar* resid;		// Initial residual vector.
        int     info = 0;      		

	if( initialGuess != NULL )
	{
            resid = initialGuess;
            info = 1;
	}
	else
	{
            resid = new Scalar[n];		/// No initial residual vector provided, so we allocate memory for one
            memset(resid,0,n*sizeof(Scalar));
            info = 0;				///< Should be zero inititally, to start with a random residual vector
	}

        Scalar* V     = new Scalar[n*ncv];	// Arnoldi basis / Schur vectors.
        Scalar* workd = new Scalar[3*n];    	// Original ARPACK internal vector.
        Scalar* workl = new Scalar[lworkl];   // Original ARPACK internal vector.
        Scalar* rwork = new Scalar[ncv];	// Only needed for complex matrices, but for simplicity we allocate it every time!

        int     ipntr[15] = {};  	///< Vector that handles original ARPACK pointers.

        memset(V,0,(n*ncv)*sizeof(Scalar));
        memset(workd,0,3*n*sizeof(Scalar));
        memset(workl,0,lworkl*sizeof(Scalar));
        memset(rwork,0,ncv*sizeof(Scalar));

        /// Setting
        int iparam[11] = {}; 		///< Vector that containse the original ARPACK parameters as specified in the following:
        iparam[0] = 1;			///< Shifting mode for the implicit restart (0:"User provided", 1:"Automatic"). "0" is not supported in this wrapper yet!
        iparam[2] = maxIterations;	///< Maximum number of iterations
        iparam[6] = 1;			///< Indicates the type of the eigenproblem (1: regular, 2: generalized, 3: shift-invert).

        char which[3] = "LM";		///< For the shift-invert mode, only the largest magnitude (i.e. closest to the target) makes sense!
//        which[0] = convertSearchTypeToString(searchMode).c_str()[0];
//        which[1] = convertSearchTypeToString(searchMode).c_str()[1];
//        which[2] = convertSearchTypeToString(searchMode).c_str()[2];

        /// Now, we build and factorize the compound matrix A - sigma*M
        Eigen::SparseMatrix<Scalar> Op = *A - sigma*(*M);

//        std::cout << "Op: " << Op << std::endl;
        Eigen::UmfPackLU<Eigen::SparseMatrix<Scalar>> solver;
//        Eigen::PardisoLU<Eigen::SparseMatrix<Scalar>> solver;

        solver.compute(Op);
        eigen_assert(solver.info()==Eigen::Success && "ArpackSupport: Factorization of combined matrix A-sigma*M failed! ");

        int ido = 0;				///< The communication flag used to indicate the status of the calculation. Must be zero initially to indicate first iteration!
        do
        {
            /// Call the Arpack routine '?naupd' to perform a solving step
            naupd(ido, bmat, n, which, nev, tol, resid, ncv, V, n, iparam, ipntr, workd, workl, lworkl, rwork, info);

            /// Wrap some of the pointers for convenience as Eigen::Vector
            Eigen::Map<const Eigen::Matrix<Scalar,Eigen::Dynamic,1>> b(&workd[ipntr[0]-1],n); // b (rhs)
            Eigen::Map<Eigen::Matrix<Scalar,Eigen::Dynamic,1>> x(&workd[ipntr[1]-1],n); // x (solution target)

            /// Evaluate, which action ARPACK requested
            switch( ido )
            {
                case -1:
                case  1:
                {
                    // Now it depends on which mode we are in:
                    ComplexVector tmp = (*M)*b;
                    x = solver.solve(tmp);
                    eigen_assert(solver.info()==Eigen::Success && "ArpackSupport: Solving of combined matrix A-sigma*M failed! ");
//                    std::cout << "Testing solution: " << (Op*x-tmp).norm() << std::endl;

                    break;
                }

//		case 2:
//                    break;

            }
        }
        while( ido != 99 );

        if (info == 1)
          std::cout << "Problem: No Convergence!" << std::endl;
	    else if (info == 3)
		  std::cout << "Problem: NumericalIssue!" << std::endl;
	    else if (info < 0)
		  std::cout << "Problem: InvalidInput!" << std::endl;
	    else if (info != 0)
		  std::cout << "Problem: Unknown ARPACK return value!" << std::endl;


        m_isInitialized = true;
        m_nbrConverged = iparam[4]; // Extract the number of converged eigenvalues!

        Scalar* evReal = new Scalar[nev+1];
        Scalar* evImag = new Scalar[nev+1];
        char howmny[2] = "A"; 

        /// Combine the converged eigenvectors into a complex matrix
        std::cout << "Computing eigenvectors" << std::endl;
        if( computeEigenvectors )
        {
            m_eivec.resize(n, ncv);
       	     /// Extract the eigenvalues and eigenvalues
            std::cout << "Running neupd: " << computeEigenvectors << " " << howmny << " " << n << " " << tol << " " << resid[0] << " " << ncv << std::endl;
            neupd(computeEigenvectors, howmny, evReal, evImag, m_eivec.data(), n, 0.0, 0.0, bmat, n, which, nev, tol, resid, ncv, V, n, iparam, ipntr, workd, workl, lworkl, rwork, info);
            std::cout << "Finished neupd" << std::endl;
            m_eigenvectorsOk = true;
        }
	else
	{
            /// Extract the eigenvalues alone
        std::cout << "Running neupd: " << computeEigenvectors << " " << howmny << " " << n << " " << tol << " " << resid[0] << " " << ncv << std::endl;
            neupd(computeEigenvectors, howmny, evReal, evImag, NULL, n, 0.0, 0.0, bmat, n, which, nev, tol, resid, ncv, V, n, iparam, ipntr, workd, workl, lworkl, rwork, info);
	}

        /// Combine the converged eigenvalues into a complex vector
        extractEigenvalues(sigma, evReal, evImag);

        /// Free the allocated memory
        delete[] evReal;
        delete[] evImag;
      	if( initialGuess == NULL ) { delete[] resid; }
        delete[] V;
        delete[] workd;
        delete[] workl;
        delete[] rwork;
    }


  /** \brief Returns the converged eigenvalues of given matrix.
   *
   * \returns A const reference to the column vector containing the eigenvalues.
   *
   * \pre The eigenvalues have been computed before.
   *
   * The eigenvalues are repeated according to their algebraic multiplicity,
   * so there are as many eigenvalues as rows in the matrix. The eigenvalues
   * are sorted in increasing order.
   *
   * Example: \include SelfAdjointEigenSolver_eigenvalues.cpp
   * Output: \verbinclude SelfAdjointEigenSolver_eigenvalues.out
   *
   * \sa eigenvectors(), MatrixBase::eigenvalues()
   */
  const ComplexVector& eigenvalues() const
  {
    eigen_assert(m_isInitialized && "ArpackGeneralizedSelfAdjointEigenSolver is not initialized.");
    return m_eigenvalues;
  }

  /** \brief Returns the converged eigenvectors of given matrix.
   *
   * \returns  A const reference to the matrix whose columns are the converged eigenvectors.
   *
   * \pre The eigenvectors have been computed before.
   *
   * Column \f$ k \f$ of the returned matrix is an eigenvector corresponding
   * to eigenvalue number \f$ k \f$ as returned by eigenvalues().  The
   * eigenvectors are normalized to have (Euclidean) norm equal to one. If
   * this object was used to solve the eigenproblem for the selfadjoint
   * matrix \f$ A \f$, then the matrix returned by this function is the
   * matrix \f$ V \f$ in the eigendecomposition \f$ A V = D V \f$.
   * For the generalized eigenproblem, the matrix returned is the solution \f$ A V = D B V \f$
   *
   * Example: \include SelfAdjointEigenSolver_eigenvectors.cpp
   * Output: \verbinclude SelfAdjointEigenSolver_eigenvectors.out
   *
   * \sa eigenvalues()
   */
  const ComplexMatrix& eigenvectors() const
  {
    eigen_assert(m_isInitialized && "ArpackGeneralizedSelfAdjointEigenSolver is not initialized.");
    eigen_assert(m_eigenvectorsOk && "The eigenvectors have not been computed together with the eigenvalues.");
    return m_eivec;
  }

  /** \brief Reports whether previous computation was successful.
   *
   * \returns \c Success if computation was succesful, \c NoConvergence otherwise.
   */
  Eigen::ComputationInfo info() const
  {
    eigen_assert(m_isInitialized && "ArpackGeneralizedSelfAdjointEigenSolver is not initialized.");
    return m_info;
  }

  size_t getNbrConvergedEigenValues() const
  { return m_nbrConverged; }

  size_t getNbrIterations() const
  { return m_nbrIterations; }


};

} // namespace Math
#endif /*ARPACKEIGENSOLVER_H_*/
