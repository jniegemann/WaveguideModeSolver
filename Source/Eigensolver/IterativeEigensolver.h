#ifndef ITERATIVEEIGENSOLVER_H_
#define ITERATIVEEIGENSOLVER_H_

#include <vector>
#include <string>

namespace Math
{

/**
 * Abstract class for an iterative eigensolver
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
    ~IterativeEigensolver() {};

};

} //namespace Math


#endif /*ITERATIVEEIGENSOLVER_H_*/
