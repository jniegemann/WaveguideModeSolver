#ifndef ARPACKWRAPPER_H_
#define ARPACKWRAPPER_H_

#include "Arpack.h"

/**
  C++ version of ARPACK routine dnaupd that implements a variant of the Arnoldi method. This routine computes approximations to a
  few eigenpairs of a linear operator "OP" with respect to a semi-inner product defined by a symmetric positive semi-definite real
  matrix B. B may be the identity matrix. NOTE: If the linear operator "OP" is real and symmetric with respect to the real positive
  semi-definite symmetric matrix B, i.e. B*OP = (OP')*B, then subroutine saupp should be used instead.

  The computed approximate eigenvalues are called Ritz values and the corresponding approximate eigenvectors are called Ritz vectors.

  naupp is usually called iteratively to solve one of the following problems:

  Mode 1:  A*x = lambda*x.
           ===> OP = A  and  B = I.

  Mode 2:  A*x = lambda*M*x, M symmetric positive definite
           ===> OP = inv[M]*A  and  B = M.
           ===> (If M can be factored see remark 3 below)

  Mode 3:  A*x = lambda*M*x, M symmetric semi-definite
           ===> OP = Real_Part{ inv[A - sigma*M]*M }  and  B = M.
           ===> shift-and-invert mode (in real arithmetic)
           If OP*x = amu*x, then
           amu = 1/2 * [ 1/(lambda-sigma) + 1/(lambda-conjg(sigma)) ].
           Note: If sigma is real, i.e. imaginary part of sigma is zero;
                 Real_Part{ inv[A - sigma*M]*M } == inv[A - sigma*M]*M
                 amu == 1/(lambda-sigma).

  Mode 4:  A*x = lambda*M*x, M symmetric semi-definite
           ===> OP = Imaginary_Part{ inv[A - sigma*M]*M }  and  B = M.
           ===> shift-and-invert mode (in real arithmetic)
           If OP*x = amu*x, then
           amu = 1/2i * [ 1/(lambda-sigma) - 1/(lambda-conjg(sigma)) ].

  Both mode 3 and 4 give the same enhancement to eigenvalues close to the (complex) shift sigma.  However, as lambda goes to infinity,
  the operator OP in mode 4 dampens the eigenvalues more strongly than does OP defined in mode 3.

  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v should
        be accomplished either by a direct method using a sparse matrix
        factorization and solving

                   [A - sigma*M]*w = v  or M*w = v,

        or through an iterative method for solving these systems. If an
        iterative method is used, the convergence test must be more
        stringent than the accuracy requirements for the eigenvalue
        approximations.

  Parameters:

    ido     (Input / Output) Reverse communication flag.  ido must be
            zero on the first call to naupp.  ido will be set
            internally to indicate the type of operation to be
            performed.  Control is then given back to the calling
            routine which has the responsibility to carry out the
            requested operation and call naupp with the result. The
            operand is given in workd[ipntr[1]], the result must be
            put in workd[ipntr[2]].
            ido =  0: first call to the reverse communication interface.
            ido = -1: compute  Y = OP * X  where
                      ipntr[1] is the pointer into workd for X,
                      ipntr[2] is the pointer into workd for Y.
                      This is for the initialization phase to force the
                      starting vector into the range of OP.
            ido =  1: compute  Y = OP * X where
                      ipntr[1] is the pointer into workd for X,
                      ipntr[2] is the pointer into workd for Y.
                      In mode 3 and 4, the vector B * X is already
                      available in workd[ipntr[3]].  It does not
                      need to be recomputed in forming OP * X.
            ido =  2: compute  Y = B * X  where
                      ipntr[1] is the pointer into workd for X,
                      ipntr[2] is the pointer into workd for Y.
            ido =  3: compute the iparam[8] real and imaginary parts
                      of the shifts where inptr[14] is the pointer
                      into workl for placing the shifts. See Remark
                      5 below.
            ido = 99: done.
    bmat    (Input) bmat specifies the type of the matrix B that defines
            the semi-inner product for the operator OP.
            bmat = 'I' -> standard eigenvalue problem A*x = lambda*x;
            bmat = 'G' -> generalized eigenvalue problem A*x = lambda*M*x.
    n       (Input) Dimension of the eigenproblem.
    which   (Input) Specify which of the Ritz values of OP to compute.
            'LM' - compute the NEV eigenvalues of largest magnitude.
            'SM' - compute the NEV eigenvalues of smallest magnitude.
            'LR' - compute the NEV eigenvalues of largest real part.
            'SR' - compute the NEV eigenvalues of smallest real part.
            'LI' - compute the NEV eigenvalues of largest imaginary part.
            'SI' - compute the NEV eigenvalues of smallest imaginary part.
    nev     (Input) Number of eigenvalues to be computed. 0 < nev < n-1.
    tol     (Input) Stopping criterion: the relative accuracy of the
            Ritz value is considered acceptable if BOUNDS[i] <=
            tol*abs(RITZ[i]),where ABS(RITZ[i]) is the magnitude when
            RITZ[i] is complex. If tol<=0.0 is passed, the machine
            precision as computed by the LAPACK auxiliary subroutine
            _LAMCH is used.
    resid   (Input / Output) Array of length n.
            On input:
            If info==0, a random initial residual vector is used.
            If info!=0, resid contains the initial residual vector,
                        possibly from a previous run.
            On output:
            resid contains the final residual vector.
    ncv     (Input) Number of Arnoldi vectors that are generated at each
            iteration. After the startup phase in which nev Arnoldi
            vectors are generated, the algorithm generates ncv-nev
            Arnoldi vectors at each subsequent update iteration. Most of
            the cost in generating each Arnoldi vector is in the
            matrix-vector product OP*x.
            NOTE: 2 <= NCV-NEV in order that complex conjugate pairs of
            Ritz values are kept together (see remark 4 below).
    V       (Output) Double precision array of length ncv*n+1. V contains
            the ncv Arnoldi basis vectors. The first element V[0] is never
            referenced.
    ldv     (Input) Dimension of the basis vectors contained in V. This
            parameter MUST be set to n.
    iparam  (Input / Output) Array of length 12.
            iparam[1]  = ISHIFT: method for selecting the implicit shifts.
            The shifts selected at each iteration are used to restart
            the Arnoldi iteration in an implicit fashion.
            -------------------------------------------------------------
            ISHIFT = 0: the shifts are provided by the user via
                        reverse communication. The real and imaginary
                        parts of the NCV eigenvalues of the Hessenberg
                        matrix H are returned in the part of the WORKL
                        array corresponding to RITZR and RITZI. See remark
                        5 below.
            ISHIFT = 1: exact shifts with respect to the current
                        Hessenberg matrix H.  This is equivalent to
                        restarting the iteration with a starting vector
                        that is a linear combination of approximate Schur
                        vectors associated with the "wanted" Ritz values.
            -------------------------------------------------------------
            iparam[2] is no longer referenced.
            iparam[3]  = MXITER
            On INPUT:  maximum number of Arnoldi update iterations allowed.
            On OUTPUT: actual number of Arnoldi update iterations taken.
            iparam[4]  = NB: blocksize to be used in the recurrence.
            The code currently works only for NB = 1.
            iparam[5]  = NCONV: number of "converged" Ritz values.
            This represents the number of Ritz values that satisfy
            the convergence criterion.
            iparam[6] is no longer referenced.
            iparam[7]  = MODE. On INPUT determines what type of
            eigenproblem is being solved. Must be 1,2,3,4.
            iparam[8]  = NP. When ido = 3 and the user provides shifts
            through reverse communication (iparam[1]=0), naupp returns
            NP, the number of shifts the user is to provide.
            0 < NP <=ncv-nev. See Remark 5 below.
            iparam[9]  =  total number of OP*x operations.
            iparam[10] = total number of B*x operations if bmat='G'.
            iparam[11] = total number of steps of re-orthogonalization.
    ipntr   (Output) Array of length 14. Pointer to mark the starting
            locations in the workd and workl arrays for matrices/vectors
            used by the Arnoldi iteration.
            ipntr[1] : pointer to the current operand vector X in workd.
            ipntr[2] : pointer to the current result vector Y in workd.
            ipntr[3] : pointer to the vector B * X in workd when used in
                       the shift-and-invert mode.
            ipntr[4] : pointer to the next available location in workl
                       that is untouched by the program.
            ipntr[5] : pointer to the ncv by ncv upper Hessenberg matrix
                       H in workl.
            ipntr[6] : pointer to the real part of the ritz value array
                       RITZR in workl.
            ipntr[7] : pointer to the imaginary part of the ritz value
                       array RITZI in workl.
            ipntr[8] : pointer to the Ritz estimates in array workl
                       associated with the Ritz values located in RITZR
                       and RITZI in workl.
            ipntr[14]: pointer to the np shifts in workl. See Remark 6.
            Note: ipntr[9:13] is only referenced by neupp. See Remark 2.
            ipntr[9] : pointer to the real part of the ncv RITZ values of
                       the original system.
            ipntr[10]: pointer to the imaginary part of the ncv RITZ values
                       of the original system.
            ipntr[11]: pointer to the ncv corresponding error bounds.
            ipntr[12]: pointer to the ncv by ncv upper quasi-triangular
                       Schur matrix for H.
            ipntr[13]: pointer to the ncv by ncv matrix of eigenvectors
                       of the upper Hessenberg matrix H. Only referenced by
                       neupp if rvec == TRUE. See Remark 2 below.
    workd   (Input / Output) Array of length 3*N+1.
            Distributed array to be used in the basic Arnoldi iteration
            for reverse communication.  The user should not use workd as
            temporary workspace during the iteration. Upon termination
            workd[1:n] contains B*resid[1:n]. If the Ritz vectors are
            desired subroutine neupp uses this output.
    workl   (Output) Array of length lworkl+1. Private (replicated) array
            on each PE or array allocated on the front end.
    lworkl  (Input) lworkl must be at least 3*ncv*(ncv+2).
    info    (Input / Output) On input, if info = 0, a randomly initial
            residual vector is used, otherwise resid contains the initial
            residual vector, possibly from a previous run.
            On output, info works as a error flag:
            =  0   : Normal exit.
            =  1   : Maximum number of iterations taken. All possible
                     eigenvalues of OP has been found. iparam[5]
                     returns the number of wanted converged Ritz values.
            =  3   : No shifts could be applied during a cycle of the
                     Implicitly restarted Arnoldi iteration. One
                     possibility is to increase the size of NCV relative
                     to nev. See remark 4 below.
            = -1   : n must be positive.
            = -2   : nev must be positive.
            = -3   : ncv must satisfy nev+2 <= ncv <= n.
            = -4   : The maximum number of Arnoldi update iterations
                     allowed must be greater than zero.
            = -5   : which must be one of 'LM','SM','LR','SR','LI','SI'.
            = -6   : bmat must be one of 'I' or 'G'.
            = -7   : Length of private work array workl is not sufficient.
            = -8   : Error return from LAPACK eigenvalue calculation.
            = -9   : Starting vector is zero.
            = -10  : iparam[7] must be 1,2,3,4.
            = -11  : iparam[7] = 1 and bmat = 'G' are incompatible.
            = -12  : iparam[1] must be equal to 0 or 1.
            = -13  : nev and which = 'BE' are incompatible.
            = -9999: Could not build an Arnoldi factorization. iparam[5]
                     returns the size of the current Arnoldi factorization.
                     The user is advised to check that enough workspace
                     and array storage has been allocated.

  Remarks:
   1. The computed Ritz values are approximate eigenvalues of OP. The
      selection of "which" should be made with this in mind when
      Mode = 3 and 4.  After convergence, approximate eigenvalues of the
      original problem may be obtained with the ARPACK subroutine neupp.
   2. If a basis for the invariant subspace corresponding to the converged
      Ritz values is needed, the user must call neupp immediately following
      completion of naupp. This is new starting with release 2 of ARPACK.
   3. If M can be factored into a Cholesky factorization M = LL'
      then Mode = 2 should not be selected.  Instead one should use
      Mode = 1 with  OP = inv(L)*A*inv(L').  Appropriate triangular
      linear systems should be solved with L and L' rather
      than computing inverses.  After convergence, an approximate
      eigenvector z of the original problem is recovered by solving
      L'z = x  where x is a Ritz vector of OP.
   4. At present there is no a-priori analysis to guide the selection
      of ncv relative to nev.  The only formal requrement is that ncv
      >= nev+2. However, it is recommended that ncv >= 2*nev+1. If many
      problems of the same type are to be solved, one should experiment
      with increasing ncv while keeping ncv fixed for a given test
      problem. This will usually decrease the required number of OP*x
      operations but it also increases the work and storage required to
      maintain the orthogonal basis vectors.   The optimal "cross-over"
      with respect to CPU time is problem dependent and must be
      determined empirically.
   5. When iparam[1] = 0, and ido = 3, the user needs to provide the
      NP = iparam[8] real and imaginary parts of the shifts in locations
          real part                  imaginary part
          -----------------------    --------------
      1   workl[ipntr[14]]           workl[ipntr[14]+NP]
      2   workl[ipntr[14]+1]         workl[ipntr[14]+NP+1]
                         .                          .
                         .                          .
                         .                          .
      NP  workl[ipntr[14]+NP-1]      workl[ipntr[14]+2*NP-1].

      Only complex conjugate pairs of shifts may be applied and the pairs
      must be placed in consecutive locations. The real part of the
      eigenvalues of the current upper Hessenberg matrix are located in
      workl[ipntr[6]] through workl[ipntr[6]+ncv-1] and the imaginary part
      in workl[ipntr[7]] through workl[ipntr[7]+ncv-1]. They are ordered
      according to the order defined by which. The complex conjugate pairs
      are kept together and the associated Ritz estimates are located in
      workl[ipntr[8]], workl[ipntr[8]+1], ... , workl[ipntr[8]+ncv-1].
*/

namespace Math
{

inline void naupd(int& ido, char bmat, int n, char* which, int nev, float tol, float* resid, int ncv, float* V, int ldv,
                  int* iparam, int* ipntr, float* workd, float* workl, int lworkl, float* rwork, int& info)
{
    snaupd_(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv, V, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
}

inline void naupd(int& ido, char bmat, int n, char* which, int nev, double tol, double* resid, int ncv, double* V, int ldv,
                  int* iparam, int* ipntr, double* workd, double* workl, int lworkl, double* rwork, int& info)
{
    dnaupd_(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv, V, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);
}

/// Actually, rwork should only be of type float*, but making it of type std::complex<float>* does not do any harm and simplifies the implementation
inline void naupd(int& ido, char bmat, int n, char* which, int nev, float tol, std::complex<float>* resid, int ncv, std::complex<float>* V, int ldv,
                  int* iparam, int* ipntr, std::complex<float>* workd, std::complex<float>* workl, int lworkl, std::complex<float>* rwork, int& info)
{
    cnaupd_(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv, V, &ldv, iparam, ipntr, workd, workl, &lworkl, (float*)rwork, &info);
}

/// Actually, rwork should only be of type double*, but making it of type std::complex<double>* does not do any harm and simplifies the implementation
inline void naupd(int& ido, char bmat, int n, char* which, int nev, double& tol, std::complex<double>* resid, int ncv, std::complex<double>* V, int ldv,
                  int* iparam, int* ipntr, std::complex<double>* workd, std::complex<double>* workl, int lworkl, std::complex<double>* rwork, int& info)
{
    znaupd_(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv, V, &ldv, iparam, ipntr, workd, workl, &lworkl, (double*)rwork, &info);
}


/*
  c++ version of ARPACK routine dneupd.
  This subroutine returns the converged approximations to eigenvalues
  of A*z = lambda*B*z and (optionally):

  (1) the corresponding approximate eigenvectors,
  (2) an orthonormal basis for the associated approximate
      invariant subspace,

  There is negligible additional cost to obtain eigenvectors. An
  orthonormal basis is always computed.  There is an additional storage cost
  of n*nev if both are requested (in this case a separate array Z must be
  supplied).
  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
  are derived from approximate eigenvalues and eigenvectors of
  of the linear operator OP prescribed by the MODE selection in the
  call to naupp. naupp must be called before this routine is called.
  These approximate eigenvalues and vectors are commonly called Ritz
  values and Ritz vectors respectively.  They are referred to as such
  in the comments that follow.  The computed orthonormal basis for the
  invariant subspace corresponding to these Ritz values is referred to
  as a Schur basis.
  See documentation in the header of the subroutine naupp for
  definition of OP as well as other terms and the relation of computed
  Ritz values and Ritz vectors of OP with respect to the given problem
  A*z = lambda*B*z. For a brief description, see definitions of
  iparam[7], MODE and which in the documentation of naupp.

  Parameters:

    rvec    (Input) Specifies whether Ritz vectors corresponding to the
            Ritz value approximations to the eigenproblem A*z = lambda*B*z
            are computed.
            rvec = false: Compute Ritz values only.
            rvec = true : Compute the Ritz vectors or Schur vectors.
                          See Remarks below.
    HowMny  (Input) Specifies the form of the basis for the invariant
            subspace corresponding to the converged Ritz values that
            is to be computed.
            = 'A': Compute nev Ritz vectors;
            = 'P': Compute nev Schur vectors;
    dr      (Output) Array of dimension nev+1.
            If iparam[7] = 1,2 or 3 and sigmai=0.0  then on exit: dr
            contains the real part of the Ritz  approximations to the
            eigenvalues of A*z = lambda*B*z.
            If iparam[7] = 3, 4 and sigmai is not equal to zero, then on
            exit: dr contains the real part of the Ritz values of OP
            computed by naupp. A further computation must be performed by
            the user to transform the Ritz values computed for OP by naupp
            to those of the original system A*z = lambda*B*z. See remark 3.
    di      (Output) Array of dimension nev+1.
            On exit, di contains the imaginary part of the Ritz value
            approximations to the eigenvalues of A*z = lambda*B*z
            associated with dr.
            NOTE: When Ritz values are complex, they will come in complex
                  conjugate pairs.  If eigenvectors are requested, the
                  corresponding Ritz vectors will also come in conjugate
                  pairs and the real and imaginary parts of these are
                  represented in two consecutive columns of the array Z
                  (see below).
    Z       (Output) Array of dimension nev*n if rvec = TRUE and HowMny =
            'A'.  if rvec = TRUE and HowMny = 'A', then Z contains the
            approximate eigenvectors (Ritz vectors) corresponding to the
            NCONV=iparam[5] Ritz values for eigensystem A*z = lambda*B*z.
            The complex Ritz vector associated with the Ritz value
            with positive imaginary part is stored in two consecutive
            columns.  The first column holds the real part of the Ritz
            vector and the second column holds the imaginary part.  The
            Ritz vector associated with the Ritz value with negative
            imaginary part is simply the complex conjugate of the Ritz
            vector associated with the positive imaginary part.
            If rvec = .FALSE. or HowMny = 'P', then Z is not referenced.
            NOTE: If rvec = TRUE and a Schur basis is not required,
                  the array Z may be set equal to first nev+1 columns of
                  the Arnoldi basis array V computed by naupp.  In this
                  case the Arnoldi basis will be destroyed and overwritten
                  with the eigenvector basis.
    ldz     (Input) Dimension of the vectors contained in Z. This
            parameter MUST be set to n.
    sigmar  (Input) If iparam[7] = 3 or 4, represents the real part of
            the shift. Not referenced if iparam[7] = 1 or 2.
    sigmai  (Input) If iparam[7] = 3 or 4, represents the imaginary part
            of the shift. Not referenced if iparam[7] = 1 or 2. See
            remark 3 below.
    workv   (Workspace) Array of dimension 3*ncv.
    V       (Input/Output) Array of dimension n*ncv+1.
            Upon Input: V contains the ncv vectors of the Arnoldi basis
                        for OP as constructed by naupp.
            Upon Output: If rvec = TRUE the first NCONV=iparam[5] columns
                        contain approximate Schur vectors that span the
                        desired invariant subspace.  See Remark 2 below.
            NOTE: If the array Z has been set equal to first nev+1 columns
                  of the array V and rvec = TRUE and HowMny = 'A', then
                  the Arnoldi basis held by V has been overwritten by the
                  desired Ritz vectors.  If a separate array Z has been
                  passed then the first NCONV=iparam[5] columns of V will
                  contain approximate Schur vectors that span the desired
                  invariant subspace.
    workl   (Input / Output) Array of length lworkl+1.
            workl[1:ncv*ncv+3*ncv] contains information obtained in
            naupp. They are not changed by neupp.
            workl[ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv] holds the real and
            imaginary part of the untransformed Ritz values, the upper
            quasi-triangular matrix for H, and the associated matrix
            representation of the invariant subspace for H.
    ipntr   (Input / Output) Array of length 14. Pointer to mark the
            starting locations in the workl array for matrices/vectors
            used by naupp and neupp.
            ipntr[9]:  pointer to the real part of the ncv RITZ values
                       of the original system.
            ipntr[10]: pointer to the imaginary part of the ncv RITZ
                       values of the original system.
            ipntr[11]: pointer to the ncv corresponding error bounds.
            ipntr[12]: pointer to the ncv by ncv upper quasi-triangular
                       Schur matrix for H.
            ipntr[13]: pointer to the ncv by ncv matrix of eigenvectors
                       of the upper Hessenberg matrix H. Only referenced
                       by neupp if rvec = TRUE. See Remark 2 below.
    info    (Output) Error flag.
            =  0 : Normal exit.
            =  1 : The Schur form computed by LAPACK routine dlahqr
                   could not be reordered by LAPACK routine dtrsen.
                   Re-enter subroutine neupp with iparam[5] = ncv and
                   increase the size of the arrays DR and DI to have
                   dimension at least dimension ncv and allocate at least
                   ncv columns for Z. NOTE: Not necessary if Z and V share
                   the same space. Please notify the authors if this error
                   occurs.
            = -1 : n must be positive.
            = -2 : nev must be positive.
            = -3 : ncv must satisfy nev+2 <= ncv <= n.
            = -5 : which must be one of 'LM','SM','LR','SR','LI','SI'.
            = -6 : bmat must be one of 'I' or 'G'.
            = -7 : Length of private work workl array is not sufficient.
            = -8 : Error return from calculation of a real Schur form.
                   Informational error from LAPACK routine dlahqr.
            = -9 : Error return from calculation of eigenvectors.
                   Informational error from LAPACK routine dtrevc.
            = -10: iparam[7] must be 1,2,3,4.
            = -11: iparam[7] = 1 and bmat = 'G' are incompatible.
            = -12: HowMny = 'S' not yet implemented
            = -13: HowMny must be one of 'A' or 'P' if rvec = TRUE.
            = -14: naupp did not find any eigenvalues to sufficient
                   accuracy.

  NOTE:     The following arguments

            bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam,
            ipntr, workd, workl, lworkl, info

            must be passed directly to neupp following the last call
            to naupp.  These arguments MUST NOT BE MODIFIED between
            the the last call to naupp and the call to neupp.

  Remarks
    1. Currently only HowMny = 'A' and 'P' are implemented.
    2. Schur vectors are an orthogonal representation for the basis of
       Ritz vectors. Thus, their numerical properties are often superior.
       Let X' denote the transpose of X. If rvec = .TRUE. then the
       relationship A * V[:,1:iparam[5]] = V[:,1:iparam[5]] * T, and
       V[:,1:iparam[5]]' * V[:,1:iparam[5]] = I are approximately satisfied.
       Here T is the leading submatrix of order iparam[5] of the real
       upper quasi-triangular matrix stored workl[ipntr[12]]. That is,
       T is block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
       each 2-by-2 diagonal block has its diagonal elements equal and its
       off-diagonal elements of opposite sign.  Corresponding to each
       2-by-2 diagonal block is a complex conjugate pair of Ritz values.
       The real Ritz values are stored on the diagonal of T.
    3. If iparam[7] = 3 or 4 and sigmai is not equal zero, then the user
       must form the iparam[5] Rayleigh quotients in order to transform the
       Ritz values computed by naupp for OP to those of A*z = lambda*B*z.
       Set rvec = TRUE. and HowMny = 'A', and compute
       Z[:,I]' * A * Z[:,I] if di[I] = 0.
       If di[I] is not equal to zero and di[I+1] = - D[I],
       then the desired real and imaginary parts of the Ritz value are
       Z[:,I]' * A * Z[:,I] +  Z[:,I+1]' * A * Z[:,I+1],
       Z[:,I]' * A * Z[:,I+1] -  Z[:,I+1]' * A * Z[:,I], respectively.
       Another possibility is to set rvec = .true. and HowMny = 'P' and
       compute V[:,1:iparam[5]]' * A * V[:,1:iparam[5]] and then an upper
       quasi-triangular matrix of order iparam[5] is computed. See remark
       2 above.
*/
inline void neupd(bool rvec, char* howMny, float* dr, float* di, float* Z, int ldz, float sigmar, float sigmai,
                  char bmat, int n, char* which, int nev, float tol, float* resid, int ncv, float* V, int ldv,
                  int* iparam, int* ipntr, float* workd, float* workl, int lworkl, float* dummy, int& info)
{
    int    irvec   = rvec ? 1 : 0;
    int*   iselect = new int[ncv];
    float* workv   = new float[3*ncv];
    float* iZ      = (Z == NULL) ? V : Z;

    std::cout << "calling sneupd_" << std::endl;
    sneupd_(&irvec, howMny, iselect, dr, di, iZ, &ldz, &sigmar, &sigmai, workv,
            &bmat, &n, which, &nev, &tol, resid, &ncv, V, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);

    delete[] iselect;
    delete[] workv;
}

inline void neupd(bool rvec, char* howMny, double* dr, double* di, double* Z, int ldz, double sigmar, double sigmai,
                  char bmat, int n, char* which, int nev, double tol, double* resid, int ncv, double* V, int ldv,
                  int* iparam, int* ipntr, double* workd, double* workl, int lworkl, double* dummy, int& info)
{
    int     irvec   = rvec ? 1 : 0;
    int*    iselect = new int[ncv];
    double* workv   = new double[3*ncv];
    double* iZ      = (Z == NULL) ? V : Z;

    std::cout << "calling dneupd_" << std::endl;
    dneupd_(&irvec, howMny, iselect, dr, di, iZ, &ldz, &sigmar, &sigmai, workv,
            &bmat, &n, which, &nev, &tol, resid, &ncv, V, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);

    delete[] iselect;
    delete[] workv;
}

inline void neupd(bool rvec, char* howMny, std::complex<float>* d, std::complex<float>* dummy, std::complex<float>* Z, int ldz, std::complex<float> sigma, std::complex<float> dummy2,
                  char bmat, int n, char* which, int nev, float tol, std::complex<float>* resid, int ncv, std::complex<float>* V, int ldv,
                  int* iparam, int* ipntr, std::complex<float>* workd, std::complex<float>* workl, int lworkl, std::complex<float>* rwork, int& info)
{
    int        irvec   = rvec ? 1 : 0;
    int*       iselect = new int[ncv];
    std::complex<float>* workev = new std::complex<float>[2*ncv];
    std::complex<float>* iZ = (Z == NULL) ? V : Z;
    std::cout << "calling cneupd_" << std::endl;

    cneupd_(&irvec, howMny, iselect, d, iZ, &ldz, &sigma, workev, &bmat, &n, which, &nev, &tol, resid,
            &ncv, V, &ldv, iparam, ipntr, workd, workl, &lworkl, (float*)rwork, &info);

    delete[] iselect;
    delete[] workev;
}

inline void neupd(bool rvec, char* howMny, std::complex<double>* d, std::complex<double>* dummy, std::complex<double>* Z, int ldz, std::complex<double> sigma, std::complex<double> dummy2, 
char bmat, int n, char* which, int nev, double tol, std::complex<double>* resid, int ncv, std::complex<double>* V, int ldv,
                  int* iparam, int* ipntr, std::complex<double>* workd, std::complex<double>* workl, int lworkl, std::complex<double>* rwork, int& info)
{
    int        irvec   = rvec ? 1 : 0;
    int*       iselect = new int[ncv];
    std::complex<double>* workev = new std::complex<double>[2*ncv];
    std::complex<double>* iZ = (Z == NULL) ? V : Z;
    std::cout << "calling zneupd_" << std::endl;
    zneupd_(&irvec, howMny, iselect, d, iZ, &ldz, &sigma, workev, &bmat, &n, which, &nev, &tol, resid, &ncv, V, &ldv, iparam, ipntr, workd, workl, &lworkl, (double*)rwork, &info);

    delete[] iselect;
    delete[] workev;
}

} // namespace Math

#endif /*ARPACKWRAPPER_H_*/
