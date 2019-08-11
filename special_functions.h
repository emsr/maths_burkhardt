#ifndef BURKHARDT_SPECIAL_FUNCTIONS_H
#define BURKHARDT_SPECIAL_FUNCTIONS_H 1

#include <stddef.h>
#ifdef __cplusplus
#include <complex>
#define __GFORTRAN_FLOAT_COMPLEX std::complex<float>
#define __GFORTRAN_DOUBLE_COMPLEX std::complex<double>
#define __GFORTRAN_LONG_DOUBLE_COMPLEX std::complex<long double>
extern "C" {
#else
#define __GFORTRAN_FLOAT_COMPLEX float _Complex
#define __GFORTRAN_DOUBLE_COMPLEX double _Complex
#define __GFORTRAN_LONG_DOUBLE_COMPLEX long double _Complex
#endif

void airyab (double x, double *ai, double *bi, double *ad, double *bd);
void confhyp (double a, double b, double x, double *chg);
void hyperg2f1 (double a, double b, double c, double x, double *hf);
void struveh (double nu, double x, double *sh);
void struvel (double nu, double x, double *sl);
void tricomiu (double a, double b, double x, double *hu, int *md);

#ifdef __cplusplus
}
#endif

#endif // BURKHARDT_SPECIAL_FUNCTIONS_H
