/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * constraints_types.h
 *
 * Code generation for function 'constraints_types'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef typedef_struct1_T
#define typedef_struct1_T

typedef struct {
  real_T cx;
  real_T cy;
  real_T w;
  real_T l;
  real_T theta;
} struct1_T;

#endif                                 /*typedef_struct1_T*/

#ifndef typedef_struct0_T
#define typedef_struct0_T

typedef struct {
  real_T predator;
  real_T pFactor;
  real_T pred_radius;
  real_T wp;
  real_T steps;
  real_T n;
  real_T h;
  real_T dt;
  real_T ct;
  real_T amax;
  real_T vmax;
  real_T minP;
  real_T maxP;
  real_T minV;
  real_T maxV;
  real_T Dmax;
  real_T dmin;
  real_T Ds;
  real_T Dc;
  real_T gamma;
  real_T rs;
  real_T ws;
  real_T wc;
  real_T wo;
  real_T wt;
  struct1_T rects[5];
} struct0_T;

#endif                                 /*typedef_struct0_T*/

/* End of code generation (constraints_types.h) */
