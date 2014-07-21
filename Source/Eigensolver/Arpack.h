#ifndef ARPACK_H_
#define ARPACK_H_

#include <complex> 

#ifdef __cplusplus
extern "C" {
#endif

    /**
     * Single precision symmetric routines
     */

	// Performs the actual iterative solving procedure
    void ssaupd_(int *ido, char *bmat, int *n, char *which, int *nev, float *tol, float *resid, int *ncv, float *V, int *ldv,
                 int *iparam, int *ipntr, float *workd, float *workl, int *lworkl, int *info);

    // Performs the postprocessing, i.e. the extraction of the eigenvalues
    void sseupd_(int *rvec, char *HowMny, int *select, float *d, float *Z, int *ldz, float *sigma, char *bmat, int *n, char *which,
                 int *nev, float *tol, float *resid, int *ncv, float *V, int *ldv, int *iparam, int *ipntr, float *workd, float *workl,
                 int *lworkl, int *info);

    /**
     * Single precision nonsymmetric routines
     */
    void snaupd_(int *ido, char *bmat, int *n, char *which, int *nev, float *tol, float *resid, int *ncv, float *V, int *ldv,
                 int *iparam, int *ipntr, float *workd, float *workl, int *lworkl, int *info);

    void sneupd_(int *rvec, char *HowMny, int *select, float *dr, float *di, float *Z, int *ldz, float *sigmar, float *sigmai,
                 float *workev, char *bmat, int *n, char *which, int *nev, float *tol, float *resid, int *ncv, float *V, int *ldv,
                 int *iparam, int *ipntr, float *workd, float *workl, int *lworkl, int *info);

    /**
     * Double precision symmetric routines
     */
    void dsaupd_(int *ido, char *bmat, int *n, char *which, int *nev, double *tol, double *resid, int *ncv, double *V, int *ldv,
                 int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *info);

    void dseupd_(int *rvec, char *HowMny, int *select, double *d, double *Z, int *ldz, double *sigma, char *bmat, int *n, char *which,
                 int *nev, double *tol, double *resid, int *ncv, double *V, int *ldv, int *iparam, int *ipntr, double *workd,
                 double *workl, int *lworkl, int *info);

    /**
     * Double precision nonsymmetric routines.
     */
    void dnaupd_(int *ido, char *bmat, int *n, char *which, int *nev, double *tol, double *resid, int *ncv, double *V, int *ldv,
                 int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *info);

    void dneupd_(int *rvec, char *HowMny, int *select, double *dr, double *di, double *Z, int *ldz, double *sigmar, double *sigmai,
                 double *workev, char *bmat, int *n, char *which, int *nev, double *tol, double *resid, int *ncv, double *V, int *ldv,
                 int *iparam, int *ipntr, double *workd, double *workl, int *lworkl, int *info);


    /**
     * Single precision complex routines.
     */
    void cnaupd_(int *ido, char *bmat, int *n, char *which, int *nev, float *tol, std::complex<float> *resid, int *ncv, std::complex<float> *V, int *ldv,
                 int *iparam, int *ipntr, std::complex<float> *workd, std::complex<float> *workl, int *lworkl, float *rwork, int *info);

    void cneupd_(int *rvec, char *HowMny, int *select, std::complex<float> *d, std::complex<float> *Z, int *ldz, std::complex<float> *sigma, std::complex<float> *workev,
                 char *bmat, int *n, char *which, int *nev, float *tol, std::complex<float> *resid, int *ncv, std::complex<float> *V, int *ldv, int *iparam,
                 int *ipntr, std::complex<float> *workd, std::complex<float> *workl, int *lworkl, float *rwork, int *info);

    /**
     * Double precision complex routines
     */
    void znaupd_(int *ido, char *bmat, int *n, char *which, int *nev, double *tol, std::complex<double> *resid, int *ncv, std::complex<double> *V, int *ldv,
                 int *iparam, int *ipntr, std::complex<double> *workd, std::complex<double> *workl, int *lworkl, double *rwork, int *info);

    void zneupd_(int *rvec, char *HowMny, int *select, std::complex<double> *d, std::complex<double> *Z, int *ldz, std::complex<double> *sigma, std::complex<double> *workev, 
                 char *bmat, int *n, char *which, int *nev, double *tol, std::complex<double> *resid, int *ncv, std::complex<double> *V, int *ldv,
                 int *iparam, int *ipntr, std::complex<double> *workd, std::complex<double> *workl, int *lworkl, double *rwork, int *info);

#ifdef __cplusplus
}
#endif

#endif /*ARPACK_H_*/
