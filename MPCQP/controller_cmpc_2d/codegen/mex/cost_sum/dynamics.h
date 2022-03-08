/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * dynamics.h
 *
 * Code generation for function 'dynamics'
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
void dynamics(const emlrtStack *sp, real_T pos[30], real_T vel[30], const real_T
              acc_data[], const int32_T acc_size[2], real_T params_n, real_T
              params_dt, real_T params_vmax);

/* End of code generation (dynamics.h) */
