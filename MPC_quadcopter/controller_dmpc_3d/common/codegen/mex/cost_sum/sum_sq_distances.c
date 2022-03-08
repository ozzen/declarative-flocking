/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sum_sq_distances.c
 *
 * Code generation for function 'sum_sq_distances'
 *
 */

/* Include files */
#include "sum_sq_distances.h"
#include "cost_sum.h"
#include "cost_sum_data.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtBCInfo d_emlrtBCI = { 1,   /* iFirst */
  8,                                   /* iLast */
  10,                                  /* lineNo */
  43,                                  /* colNo */
  "pos",                               /* aName */
  "sum_sq_distances",                  /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/cost_functions/sum_sq_distances.m",/* pName */
  0                                    /* checkKind */
};

/* Function Definitions */
real_T sum_sq_distances(const emlrtStack *sp, const real_T s[96])
{
  real_T res;
  int32_T i;
  int32_T b_i;
  int32_T j;
  int32_T c_i;
  real_T d;
  real_T z1_idx_0;
  real_T z1_idx_1;

  /*  Usama Mehmood - Oct 2019 */
  res = 0.0;
  for (i = 0; i < 8; i++) {
    b_i = 6 - i;
    for (j = 0; j <= b_i; j++) {
      c_i = (i + j) + 2;
      if (c_i > 8) {
        emlrtDynamicBoundsCheckR2012b(c_i, 1, 8, &d_emlrtBCI, sp);
      }

      c_i = 12 * (c_i - 1);
      d = s[12 * i] - s[c_i];
      z1_idx_0 = d * d;
      d = s[12 * i + 1] - s[c_i + 1];
      z1_idx_1 = d * d;
      d = s[12 * i + 2] - s[c_i + 2];
      res += (z1_idx_0 + z1_idx_1) + d * d;
      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b(sp);
      }
    }

    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  return res;
}

/* End of code generation (sum_sq_distances.c) */
