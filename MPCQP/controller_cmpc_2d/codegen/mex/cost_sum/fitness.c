/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fitness.c
 *
 * Code generation for function 'fitness'
 *
 */

/* Include files */
#include "fitness.h"
#include "cost_sum.h"
#include "cost_sum_data.h"
#include "mwmathutil.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo f_emlrtRSI = { 2,   /* lineNo */
  "fitness",                           /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_cmpc_2d/fitness.m"/* pathName */
};

static emlrtRSInfo g_emlrtRSI = { 4,   /* lineNo */
  "fitness",                           /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_cmpc_2d/fitness.m"/* pathName */
};

static emlrtRSInfo h_emlrtRSI = { 5,   /* lineNo */
  "fitness",                           /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_cmpc_2d/fitness.m"/* pathName */
};

static emlrtBCInfo b_emlrtBCI = { 1,   /* iFirst */
  15,                                  /* iLast */
  8,                                   /* lineNo */
  49,                                  /* colNo */
  "pos",                               /* aName */
  "separation_polynomial",             /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_cmpc_2d/cost_functions/separation_polynomial.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = { 1,   /* iFirst */
  15,                                  /* iLast */
  8,                                   /* lineNo */
  47,                                  /* colNo */
  "pos",                               /* aName */
  "sum_sq_distances",                  /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_cmpc_2d/cost_functions/sum_sq_distances.m",/* pName */
  0                                    /* checkKind */
};

static emlrtRTEInfo h_emlrtRTEI = { 4, /* lineNo */
  9,                                   /* colNo */
  "velocity_matching",                 /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_cmpc_2d/cost_functions/velocity_matching.m"/* pName */
};

static emlrtRTEInfo i_emlrtRTEI = { 5, /* lineNo */
  13,                                  /* colNo */
  "velocity_matching",                 /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_cmpc_2d/cost_functions/velocity_matching.m"/* pName */
};

static emlrtBCInfo d_emlrtBCI = { 1,   /* iFirst */
  15,                                  /* iLast */
  6,                                   /* lineNo */
  32,                                  /* colNo */
  "vel",                               /* aName */
  "velocity_matching",                 /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_cmpc_2d/cost_functions/velocity_matching.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo e_emlrtBCI = { 1,   /* iFirst */
  15,                                  /* iLast */
  6,                                   /* lineNo */
  43,                                  /* colNo */
  "vel",                               /* aName */
  "velocity_matching",                 /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_cmpc_2d/cost_functions/velocity_matching.m",/* pName */
  0                                    /* checkKind */
};

/* Function Definitions */
real_T fitness(const emlrtStack *sp, const real_T pos[30], const real_T vel[30],
               real_T params_n)
{
  real_T b_sp;
  int32_T i;
  int32_T b_i;
  int32_T j;
  int32_T c_i;
  real_T asd;
  int32_T i1;
  real_T vm;
  real_T scale;
  int32_T absxk_tmp;
  real_T absxk;
  real_T t;
  real_T y;
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &f_emlrtRSI;

  /*  Given an pos -> 2 x n matrix, res is sum of squared pairwise distances.  */
  /*  */
  b_sp = 0.0;
  for (i = 0; i < 14; i++) {
    b_i = 13 - i;
    for (j = 0; j <= b_i; j++) {
      c_i = (i + j) + 2;
      if (c_i > 15) {
        emlrtDynamicBoundsCheckR2012b(c_i, 1, 15, &b_emlrtBCI, &st);
      }

      i1 = i << 1;
      c_i = (c_i - 1) << 1;
      vm = pos[i1] - pos[c_i];
      scale = vm * vm;
      vm = pos[i1 + 1] - pos[c_i + 1];
      b_sp += 1.0 / (scale + vm * vm);
      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b(&st);
      }
    }

    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(&st);
    }
  }

  b_sp /= 105.0;

  /*  spC = separation_constrained(pos, params); */
  st.site = &g_emlrtRSI;

  /*  Given an pos -> 2 x n matrix, res is sum of squared pairwise distances.  */
  /*  */
  asd = 0.0;
  for (i = 0; i < 14; i++) {
    b_i = 13 - i;
    for (j = 0; j <= b_i; j++) {
      c_i = (i + j) + 2;
      if (c_i > 15) {
        emlrtDynamicBoundsCheckR2012b(c_i, 1, 15, &c_emlrtBCI, &st);
      }

      i1 = i << 1;
      c_i = (c_i - 1) << 1;
      vm = pos[i1] - pos[c_i];
      scale = vm * vm;
      vm = pos[i1 + 1] - pos[c_i + 1];
      asd += scale + vm * vm;
      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b(&st);
      }
    }

    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(&st);
    }
  }

  asd /= 105.0;
  st.site = &h_emlrtRSI;

  /*  Given an vel -> 2 x n matrix. */
  vm = 0.0;
  b_i = (int32_T)(params_n - 1.0);
  emlrtForLoopVectorCheckR2012b(1.0, 1.0, params_n - 1.0, mxDOUBLE_CLASS,
    (int32_T)(params_n - 1.0), &h_emlrtRTEI, &st);
  for (i = 0; i < b_i; i++) {
    i1 = (int32_T)(params_n + (1.0 - (((real_T)i + 1.0) + 1.0)));
    emlrtForLoopVectorCheckR2012b(((real_T)i + 1.0) + 1.0, 1.0, params_n,
      mxDOUBLE_CLASS, i1, &i_emlrtRTEI, &st);
    for (j = 0; j < i1; j++) {
      c_i = (int32_T)(i + 1U);
      if ((c_i < 1) || (c_i > 15)) {
        emlrtDynamicBoundsCheckR2012b(c_i, 1, 15, &d_emlrtBCI, &st);
      }

      c_i = (int32_T)((((real_T)i + 1.0) + 1.0) + (real_T)j);
      if ((c_i < 1) || (c_i > 15)) {
        emlrtDynamicBoundsCheckR2012b(c_i, 1, 15, &e_emlrtBCI, &st);
      }

      scale = 3.3121686421112381E-170;
      absxk_tmp = i << 1;
      c_i = (c_i - 1) << 1;
      absxk = muDoubleScalarAbs(vel[absxk_tmp] - vel[c_i]);
      if (absxk > 3.3121686421112381E-170) {
        y = 1.0;
        scale = absxk;
      } else {
        t = absxk / 3.3121686421112381E-170;
        y = t * t;
      }

      absxk = muDoubleScalarAbs(vel[absxk_tmp + 1] - vel[c_i + 1]);
      if (absxk > scale) {
        t = scale / absxk;
        y = y * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }

      y = scale * muDoubleScalarSqrt(y);
      vm += y;
      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b(&st);
      }
    }

    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(&st);
    }
  }

  vm /= params_n * (params_n - 1.0) / 2.0;

  /*  res = 10 * asd + 1000000 * spC; */
  return (2.0 * asd + 5000.0 * b_sp) + 1.0E+7 * vm;
}

/* End of code generation (fitness.c) */
