#include <stdlib.h>
#include "math_functions.h"

/*==============================================================================
 * RANDOM NUMBERS
 *
 * Real RandomNormal()       - Normally distributed random numbers
 * Real RandomReal()         - Uniformly distributed random numbers
 *
 *----------------------------------------------------------------------------*/

static inline Real RandomReal()
/* Returns a uniformly distributed random number in the interval [0,1].       */
{
  return (Real) rand()/RAND_MAX;
}

static Real RandomNormal(Real mu, Real sigma)
/* Implements the box-muller routine.  Gives a mean of mu, and a
   standard deviation sigma.                                                  */
{
  Real x1, x2, w, y1;
  static Real y2;
  static int use_last = 0;

  if (use_last){ /* use value from previous call */
    y1 = y2;
    use_last = 0;
  }
  else {
    do {
      x1 = 2.0 * RandomReal() - 1.0;
      x2 = 2.0 * RandomReal() - 1.0;
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = sqrt((-2.0 * log(w)) / w);
    y1 = x1 * w;
    y2 = x2 * w;
    use_last = 1;
  }

  return (mu + y1 * sigma);
}



/*==============================================================================
 * VECTOR ANALYSIS
 *
 * Real Dot()                - dot product and norm of vectors
 * Real Norm()
 *
 *----------------------------------------------------------------------------*/

static inline Real Norm (Real vec[3])
{
  return sqrt(Dot(vec,vec));
}

static inline Real Dot (Real vec1[3], Real vec2[3])
{
  Real sum = 0.0;
  int i;
  for (i=0; i<3; i++){
    sum += vec1[i] * vec2[i];
  }
  return sum;
}


/* Linear interpolation of array values */
static int locate (Real *array, Real x, const int n)
/* Given an array array[0..n-1], and given a value x, returns the
   index j, such that x is between array[j] and array[j+1].  array
   must be monotonic.  Returns -1 if x is out of range.  */
{
  /* Coordinates */
  int ju = n, jl = 0;
  int jm = (ju + jl)/2;

  int ascending = (array[n-1] >= array[0]);

  /* Check to see whether x is in the array */
  if (ascending){
    if ((x > array[n-1]) || (x < array[0]))
      return -1;
  }
  else{
    if ((x < array[n-1]) || (x > array[0]))
      return -1;
  }


  /* Squeeze jl and ju around the point x */
  while (ju-jl > 1){
    jm = (ju+jl)/2;

    if ((ascending && x >= array[jm])
        || (!ascending && x <= array[jm])) jl = jm;
    else ju = jm;
  }

  return jl;
}

static Real interpolate(Real *xvals, Real *yvals, Real x_target, const int n)
{
  Real m;
  int j = locate(xvals, x_target, n);
  int ascending = (xvals[n-1] >= xvals[0]);

  if (j == -1 || j == 0 || j == n-1){
    if ((ascending && x_target <= xvals[0])
        || (!ascending && x_target >= xvals[0]))
      return yvals[0];
    else
      return yvals[n-1];
  }

  m = (yvals[j+1]-yvals[j])/(xvals[j+1]-xvals[j]);
  return yvals[j] + m * (x_target-xvals[j]);
}
