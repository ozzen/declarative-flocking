/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * cost_sum.c
 *
 * Code generation for function 'cost_sum'
 *
 */

/* Include files */
#include "cost_sum.h"
#include "cost_sum_data.h"
#include "cost_sum_emxutil.h"
#include "dynamics.h"
#include "fitness.h"
#include "rt_nonfinite.h"
#include "u2acc.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 3,     /* lineNo */
  "cost_sum",                          /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_cmpc_2d/common/cost_sum.m"/* pathName */
};

static emlrtRSInfo b_emlrtRSI = { 9,   /* lineNo */
  "cost_sum",                          /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_cmpc_2d/common/cost_sum.m"/* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 10,  /* lineNo */
  "cost_sum",                          /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_cmpc_2d/common/cost_sum.m"/* pathName */
};

static emlrtRTEInfo emlrtRTEI = { 6,   /* lineNo */
  13,                                  /* colNo */
  "cost_sum",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_cmpc_2d/common/cost_sum.m"/* pName */
};

static emlrtBCInfo emlrtBCI = { -1,    /* iFirst */
  -1,                                  /* iLast */
  7,                                   /* lineNo */
  21,                                  /* colNo */
  "acc",                               /* aName */
  "cost_sum",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_cmpc_2d/common/cost_sum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtRTEInfo b_emlrtRTEI = { 8, /* lineNo */
  17,                                  /* colNo */
  "cost_sum",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_cmpc_2d/common/cost_sum.m"/* pName */
};

static emlrtRTEInfo k_emlrtRTEI = { 3, /* lineNo */
  5,                                   /* colNo */
  "cost_sum",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_cmpc_2d/common/cost_sum.m"/* pName */
};

/* Function Definitions */
real_T cost_sum(const emlrtStack *sp, const real_T u[180], real_T pos[30],
                real_T vel[30], const struct0_T *params)
{
  real_T ret;
  emxArray_real_T *acc;
  real_T control_steps;
  real_T res;
  int32_T i;
  int32_T h;
  int32_T loop_ub;
  int32_T i1;
  real_T x[180];
  int32_T a_size[2];
  real_T a_data[360];
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInit_real_T(sp, &acc, 3, &k_emlrtRTEI, true);
  st.site = &emlrtRSI;
  u2acc(&st, u, params->n, params->h, acc);
  control_steps = params->ct / params->dt;
  res = 0.0;
  i = (int32_T)params->h;
  emlrtForLoopVectorCheckR2012b(1.0, 1.0, params->h, mxDOUBLE_CLASS, (int32_T)
    params->h, &emlrtRTEI, sp);
  for (h = 0; h < i; h++) {
    loop_ub = acc->size[0];
    i1 = (int32_T)(h + 1U);
    if ((i1 < 1) || (i1 > acc->size[2])) {
      emlrtDynamicBoundsCheckR2012b(i1, 1, acc->size[2], &emlrtBCI, sp);
    }

    a_size[0] = 2;
    a_size[1] = acc->size[0];
    for (i1 = 0; i1 < loop_ub; i1++) {
      a_data[2 * i1] = acc->data[i1 + acc->size[0] * 2 * h];
      a_data[2 * i1 + 1] = acc->data[(i1 + acc->size[0]) + acc->size[0] * 2 * h];
    }

    i1 = (int32_T)control_steps;
    emlrtForLoopVectorCheckR2012b(1.0, 1.0, control_steps, mxDOUBLE_CLASS,
      (int32_T)control_steps, &b_emlrtRTEI, sp);
    for (loop_ub = 0; loop_ub < i1; loop_ub++) {
      st.site = &b_emlrtRSI;
      dynamics(&st, pos, vel, a_data, a_size, params->n, params->dt,
               params->vmax);
      st.site = &c_emlrtRSI;
      res += fitness(&st, pos, vel, params->n);
      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b(sp);
      }
    }

    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  emxFree_real_T(&acc);
  for (loop_ub = 0; loop_ub < 180; loop_ub++) {
    x[loop_ub] = u[loop_ub] * u[loop_ub];
  }

  control_steps = x[0];
  for (loop_ub = 0; loop_ub < 179; loop_ub++) {
    control_steps += x[loop_ub + 1];
  }

  ret = res / params->h + control_steps / params->h;
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
  return ret;
}

/* End of code generation (cost_sum.c) */
