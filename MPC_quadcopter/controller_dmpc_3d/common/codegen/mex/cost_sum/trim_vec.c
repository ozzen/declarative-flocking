/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * trim_vec.c
 *
 * Code generation for function 'trim_vec'
 *
 */

/* Include files */
#include "trim_vec.h"
#include "cost_sum.h"
#include "mwmathutil.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void trim_vec(real_T v[3], real_T mag)
{
  real_T scale;
  real_T absxk;
  real_T t;
  real_T d;

  /*  Usama Mehmood - Oct 2019 */
  scale = 3.3121686421112381E-170;
  absxk = muDoubleScalarAbs(v[0]);
  if (absxk > 3.3121686421112381E-170) {
    d = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    d = t * t;
  }

  absxk = muDoubleScalarAbs(v[1]);
  if (absxk > scale) {
    t = scale / absxk;
    d = d * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    d += t * t;
  }

  absxk = muDoubleScalarAbs(v[2]);
  if (absxk > scale) {
    t = scale / absxk;
    d = d * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    d += t * t;
  }

  d = scale * muDoubleScalarSqrt(d);
  if (d > mag) {
    scale = mag / d;
    v[0] *= scale;
    v[1] *= scale;
    v[2] *= scale;
  }
}

/* End of code generation (trim_vec.c) */
