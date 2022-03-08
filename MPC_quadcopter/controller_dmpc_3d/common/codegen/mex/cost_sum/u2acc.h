/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * u2acc.h
 *
 * Code generation for function 'u2acc'
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
void u2acc(const emlrtStack *sp, const real_T u[18], real_T params_h, real_T
           acc_data[], int32_T acc_size[2]);

/* End of code generation (u2acc.h) */
