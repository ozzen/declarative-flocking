/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * separation.c
 *
 * Code generation for function 'separation'
 *
 */

/* Include files */
#include "separation.h"
#include "cost_sum.h"
#include "cost_sum_data.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtBCInfo e_emlrtBCI = { 1,   /* iFirst */
  8,                                   /* iLast */
  10,                                  /* lineNo */
  45,                                  /* colNo */
  "pos",                               /* aName */
  "separation",                        /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/cost_functions/separation.m",/* pName */
  0                                    /* checkKind */
};

/* Function Definitions */
real_T separation(const emlrtStack *sp, const real_T s[96])
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
        emlrtDynamicBoundsCheckR2012b(c_i, 1, 8, &e_emlrtBCI, sp);
      }

      c_i = 12 * (c_i - 1);
      d = s[12 * i] - s[c_i];
      z1_idx_0 = d * d;
      d = s[12 * i + 1] - s[c_i + 1];
      z1_idx_1 = d * d;
      d = s[12 * i + 2] - s[c_i + 2];
      res += 1.0 / ((z1_idx_0 + z1_idx_1) + d * d);
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

/* End of code generation (separation.c) */
