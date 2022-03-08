/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * dynamics_quadcopter.h
 *
 * Code generation for function 'dynamics_quadcopter'
 *
 */

#pragma once

/* Include files */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include "emlrt.h"
#include "rtwtypes.h"
#include "cost_sum_types.h"

/* Function Declarations */
void dynamics_quadcopter(const emlrtStack *sp, real_T s[12], const real_T u[4],
  const struct0_T *params, real_T omegas[4]);

/* End of code generation (dynamics_quadcopter.h) */
