#ifndef MATH_FUNCTIONS_H
#define MATH_FUNCTIONS_H

/* Uniformly distrubuted random numbers */
static inline Real RandomReal();

/* Normally distributed random numbers */
static Real RandomNormal(Real, Real);

/* Vector analysis functions */
static inline Real Norm (Real vec[3]);
static inline Real Dot (Real v1[3], Real v2[3]);

/* Linear interpolation of array values */
static int locate(Real *array, Real target, const int n);
static Real interpolate(Real *xvals, Real *yvals, Real x_target, const int n);

#include "math_functions.c"

#endif /* MATH_FUNCTIONS_H */
