#ifndef SPARSEEIGENSOLVER_H_
#define SPARSEEIGENSOLVER_H_

#include <vector>
#include <string>

namespace Math
{

/**
 * Abstract class for a sparse eigensolver
 */
class IterativeEigensolver
{
protected:

public:

    /**
     * Specifies which eigenvalues should be calculated
     */
    enum class SearchMode
    {
        LargestMagnitude,
        SmallestMagnitude,
        LargestRealPart,
        SmallestRealPart,
        LargestImaginaryPart,
        SmallestImaginaryPart
    };

    /**
     * Destructor
     */
    virtual ~SparseEigensolver() {};

    /**
     * Find extremal eigenvalues
     */
    virtual Eigen::VectorXc findEigenvalues(const MatrixType& A, const MatrixType& B, int numEigenvalues, SearchMode searchMode) = 0;

    /**
     * Find inner eigenvalues
     */
    virtual Eigen::VectorXc findEigenvalues(const MatrixType& A, const MatrixType& B, CplxDouble sigma, int numEigenvalues) = 0;
};

} //namespace Math


#endif /*SPARSEEIGENSOLVER_H_*/
