/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * constraints.c
 *
 * Code generation for function 'constraints'
 *
 */

/* Include files */
#include "constraints.h"
#include "mwmathutil.h"
#include "power.h"
#include "rt_nonfinite.h"
#include "sum.h"
#include "u2acc.h"
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 4,     /* lineNo */
  "constraints",                       /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_cmpc_2d/common/constraints.m"/* pathName */
};

static emlrtRSInfo b_emlrtRSI = { 5,   /* lineNo */
  "constraints",                       /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_cmpc_2d/common/constraints.m"/* pathName */
};

static emlrtRTEInfo emlrtRTEI = { 13,  /* lineNo */
  9,                                   /* colNo */
  "sqrt",                              /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/elfun/sqrt.m"/* pName */
};

/* Function Definitions */
void constraints(const emlrtStack *sp, const real_T u[24], const struct0_T
                 *params, real_T c_data[], int32_T c_size[1])
{
  real_T acc_data[1152];
  int32_T acc_size[3];
  real_T tmp_data[1152];
  int32_T tmp_size[3];
  real_T res_data[576];
  boolean_T p;
  int32_T loop_ub;
  int32_T k;
  int32_T i;
  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;

  /*  constraints - non-linear constraints to enforce magnitude of acceleration */
  /*  */
  st.site = &emlrtRSI;
  u2acc(&st, u, params->n, params->h, acc_data, acc_size);
  st.site = &b_emlrtRSI;
  power(acc_data, acc_size, tmp_data, tmp_size);
  b_st.site = &b_emlrtRSI;
  sum(&b_st, tmp_data, tmp_size, res_data, acc_size);
  p = false;
  loop_ub = acc_size[0] * acc_size[2];
  for (k = 0; k < loop_ub; k++) {
    if (p || (res_data[k] < 0.0)) {
      p = true;
    }
  }

  if (p) {
    emlrtErrorWithMessageIdR2018a(&st, &emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  for (k = 0; k < loop_ub; k++) {
    res_data[k] = muDoubleScalarSqrt(res_data[k]);
  }

  k = acc_size[0] * acc_size[1] * acc_size[2];
  for (i = 0; i < k; i++) {
    res_data[i] -= params->amax;
  }

  c_size[0] = loop_ub;
  if (0 <= loop_ub - 1) {
    memcpy(&c_data[0], &res_data[0], loop_ub * sizeof(real_T));
  }
}

/* End of code generation (constraints.c) */
