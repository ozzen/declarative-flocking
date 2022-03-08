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
#include "eml_int_forloop_overflow_check.h"
#include "mwmathutil.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 17,    /* lineNo */
  "cost_sum",                          /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/cost_sum.m"/* pathName */
};

static emlrtRSInfo b_emlrtRSI = { 22,  /* lineNo */
  "cost_sum",                          /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/cost_sum.m"/* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 25,  /* lineNo */
  "cost_sum",                          /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/cost_sum.m"/* pathName */
};

static emlrtRSInfo d_emlrtRSI = { 27,  /* lineNo */
  "cost_sum",                          /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/cost_sum.m"/* pathName */
};

static emlrtRSInfo e_emlrtRSI = { 6,   /* lineNo */
  "u2acc",                             /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/u2acc.m"/* pathName */
};

static emlrtRSInfo f_emlrtRSI = { 18,  /* lineNo */
  "reshapeSizeChecks",                 /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/reshapeSizeChecks.m"/* pathName */
};

static emlrtRSInfo g_emlrtRSI = { 21,  /* lineNo */
  "eml_int_forloop_overflow_check",    /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"/* pathName */
};

static emlrtRSInfo h_emlrtRSI = { 5,   /* lineNo */
  "dynamics",                          /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/dynamics.m"/* pathName */
};

static emlrtRSInfo i_emlrtRSI = { 24,  /* lineNo */
  "fitness",                           /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/fitness.m"/* pathName */
};

static emlrtDCInfo emlrtDCI = { 14,    /* lineNo */
  19,                                  /* colNo */
  "cost_sum",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/cost_sum.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo emlrtBCI = { 1,     /* iFirst */
  15,                                  /* iLast */
  14,                                  /* lineNo */
  19,                                  /* colNo */
  "pos",                               /* aName */
  "cost_sum",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/cost_sum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo b_emlrtDCI = { 15,  /* lineNo */
  19,                                  /* colNo */
  "cost_sum",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/cost_sum.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = { 1,   /* iFirst */
  15,                                  /* iLast */
  15,                                  /* lineNo */
  19,                                  /* colNo */
  "vel",                               /* aName */
  "cost_sum",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/cost_sum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtRTEInfo emlrtRTEI = { 20,  /* lineNo */
  13,                                  /* colNo */
  "cost_sum",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/cost_sum.m"/* pName */
};

static emlrtBCInfo c_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  21,                                  /* lineNo */
  19,                                  /* colNo */
  "acc",                               /* aName */
  "cost_sum",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/cost_sum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  24,                                  /* lineNo */
  28,                                  /* colNo */
  "acceleration",                      /* aName */
  "cost_sum",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/cost_sum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtECInfo emlrtECI = { -1,    /* nDims */
  24,                                  /* lineNo */
  13,                                  /* colNo */
  "cost_sum",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/cost_sum.m"/* pName */
};

static emlrtRTEInfo b_emlrtRTEI = { 59,/* lineNo */
  23,                                  /* colNo */
  "reshapeSizeChecks",                 /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/reshapeSizeChecks.m"/* pName */
};

static emlrtRTEInfo c_emlrtRTEI = { 57,/* lineNo */
  23,                                  /* colNo */
  "reshapeSizeChecks",                 /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/reshapeSizeChecks.m"/* pName */
};

static emlrtRTEInfo d_emlrtRTEI = { 52,/* lineNo */
  13,                                  /* colNo */
  "reshapeSizeChecks",                 /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/reshapeSizeChecks.m"/* pName */
};

static emlrtECInfo b_emlrtECI = { 2,   /* nDims */
  3,                                   /* lineNo */
  7,                                   /* colNo */
  "dynamics",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/dynamics.m"/* pName */
};

static emlrtBCInfo e_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  6,                                   /* lineNo */
  22,                                  /* colNo */
  "vel",                               /* aName */
  "dynamics",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/dynamics.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo f_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  7,                                   /* lineNo */
  68,                                  /* colNo */
  "vel",                               /* aName */
  "dynamics",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/dynamics.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo g_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  7,                                   /* lineNo */
  82,                                  /* colNo */
  "vel",                               /* aName */
  "dynamics",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/dynamics.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo h_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  11,                                  /* lineNo */
  19,                                  /* colNo */
  "vel",                               /* aName */
  "dynamics",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/dynamics.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo i_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  12,                                  /* lineNo */
  45,                                  /* colNo */
  "vel",                               /* aName */
  "dynamics",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/dynamics.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo j_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  12,                                  /* lineNo */
  59,                                  /* colNo */
  "vel",                               /* aName */
  "dynamics",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/dynamics.m",/* pName */
  0                                    /* checkKind */
};

static emlrtECInfo c_emlrtECI = { 2,   /* nDims */
  15,                                  /* lineNo */
  7,                                   /* colNo */
  "dynamics",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/dynamics.m"/* pName */
};

static emlrtBCInfo k_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  12,                                  /* lineNo */
  15,                                  /* colNo */
  "vel",                               /* aName */
  "dynamics",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/dynamics.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo l_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  7,                                   /* lineNo */
  19,                                  /* colNo */
  "vel",                               /* aName */
  "dynamics",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/dynamics.m",/* pName */
  0                                    /* checkKind */
};

static emlrtRTEInfo e_emlrtRTEI = { 49,/* lineNo */
  19,                                  /* colNo */
  "assertValidSizeArg",                /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/assertValidSizeArg.m"/* pName */
};

static emlrtRTEInfo f_emlrtRTEI = { 64,/* lineNo */
  15,                                  /* colNo */
  "assertValidSizeArg",                /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/assertValidSizeArg.m"/* pName */
};

static emlrtBCInfo m_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  6,                                   /* lineNo */
  22,                                  /* colNo */
  "pos",                               /* aName */
  "murmuration_cost",                  /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/cost_functions/murmuration_cost.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo n_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  6,                                   /* lineNo */
  31,                                  /* colNo */
  "pos",                               /* aName */
  "murmuration_cost",                  /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/cost_functions/murmuration_cost.m",/* pName */
  0                                    /* checkKind */
};

static emlrtRTEInfo g_emlrtRTEI = { 13,/* lineNo */
  13,                                  /* colNo */
  "toLogicalCheck",                    /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/toLogicalCheck.m"/* pName */
};

/* Function Definitions */
real_T cost_sum(const emlrtStack *sp, const real_T u[6], const real_T pos[30],
                const real_T vel[30], const struct0_T *params)
{
  real_T ret;
  boolean_T b;
  int32_T loop_ub;
  int32_T pos_size[2];
  int32_T i;
  int32_T pos_tmp;
  int32_T b_loop_ub;
  real_T pos_data[30];
  int32_T vel_size[2];
  real_T vel_data[30];
  real_T varargin_1[2];
  int32_T k;
  boolean_T guard1 = false;
  int32_T exitg2;
  boolean_T exitg1;
  real_T b_params;
  int32_T num_idx_0_tmp;
  real_T acc_data[12];
  uint8_T b_u;
  int32_T h;
  int32_T acceleration_size[2];
  real_T acceleration_data[30];
  int32_T iv[1];
  int32_T iv1[1];
  int32_T i1;
  real_T t;
  real_T scale;
  real_T d;
  real_T absxk_tmp;
  real_T y;
  real_T b_absxk_tmp;
  int32_T i2;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  b = false;

  /*  cost_sum() - dmpc */
  /*  Input: */
  /*    - u            % 2.h x 1 */
  /*    - pos          % 2 x n - state of the flock  */
  /*    - vel          % 2 x n - state of the flock  */
  /*    - params       % parameters */
  /*  Output: */
  /*    - a            % 2 x n - The sequence of control action over the horizon  */
  /*    - fit_val      % cost for the optimal sequence of accelerations. */
  /*  Usama Mehmood - Feb 2020 */
  /*  */
  if (1.0 > params->nn) {
    loop_ub = 0;
  } else {
    if (params->nn != (int32_T)muDoubleScalarFloor(params->nn)) {
      emlrtIntegerCheckR2012b(params->nn, &emlrtDCI, sp);
    }

    loop_ub = (int32_T)params->nn;
    if ((loop_ub < 1) || (loop_ub > 15)) {
      emlrtDynamicBoundsCheckR2012b(loop_ub, 1, 15, &emlrtBCI, sp);
    }
  }

  pos_size[0] = 2;
  pos_size[1] = loop_ub;
  for (i = 0; i < loop_ub; i++) {
    pos_tmp = i << 1;
    pos_data[2 * i] = pos[pos_tmp];
    pos_data[2 * i + 1] = pos[pos_tmp + 1];
  }

  if (1.0 > params->nn) {
    b_loop_ub = 0;
  } else {
    if (params->nn != (int32_T)muDoubleScalarFloor(params->nn)) {
      emlrtIntegerCheckR2012b(params->nn, &b_emlrtDCI, sp);
    }

    b_loop_ub = (int32_T)params->nn;
    if ((b_loop_ub < 1) || (b_loop_ub > 15)) {
      emlrtDynamicBoundsCheckR2012b(b_loop_ub, 1, 15, &b_emlrtBCI, sp);
    }
  }

  vel_size[0] = 2;
  vel_size[1] = b_loop_ub;
  for (i = 0; i < b_loop_ub; i++) {
    pos_tmp = i << 1;
    vel_data[2 * i] = vel[pos_tmp];
    vel_data[2 * i + 1] = vel[pos_tmp + 1];
  }

  st.site = &emlrtRSI;

  /*  u2acc - dmpc - convert column vector u to acc matric */
  /*  Output: */
  /*    - acc            % 2 x h - The sequence of control action over the horizon */
  /*  */
  b_st.site = &e_emlrtRSI;
  varargin_1[0] = params->h;
  varargin_1[1] = 2.0;
  c_st.site = &f_emlrtRSI;
  k = 0;
  guard1 = false;
  do {
    exitg2 = 0;
    if (k < 2) {
      if ((varargin_1[k] != muDoubleScalarFloor(varargin_1[k])) ||
          muDoubleScalarIsInf(varargin_1[k])) {
        guard1 = true;
        exitg2 = 1;
      } else {
        k++;
        guard1 = false;
      }
    } else {
      k = 0;
      exitg2 = 2;
    }
  } while (exitg2 == 0);

  if (exitg2 != 1) {
    exitg1 = false;
    while ((!exitg1) && (k < 2)) {
      if ((varargin_1[k] < -2.147483648E+9) || (varargin_1[k] > 2.147483647E+9))
      {
        guard1 = true;
        exitg1 = true;
      } else {
        k++;
      }
    }
  }

  if (guard1) {
    emlrtErrorWithMessageIdR2018a(&c_st, &e_emlrtRTEI,
      "Coder:toolbox:eml_assert_valid_size_arg_invalidSizeVector",
      "Coder:toolbox:eml_assert_valid_size_arg_invalidSizeVector", 4, 12,
      MIN_int32_T, 12, MAX_int32_T);
  }

  if (params->h <= 0.0) {
    b_params = 0.0;
  } else {
    b_params = params->h;
  }

  if (!(b_params * 2.0 <= 2.147483647E+9)) {
    emlrtErrorWithMessageIdR2018a(&c_st, &f_emlrtRTEI, "Coder:MATLAB:pmaxsize",
      "Coder:MATLAB:pmaxsize", 0);
  }

  num_idx_0_tmp = (int32_T)params->h;
  if (num_idx_0_tmp > 6) {
    emlrtErrorWithMessageIdR2018a(&b_st, &d_emlrtRTEI,
      "Coder:toolbox:reshape_emptyReshapeLimit",
      "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }

  if (num_idx_0_tmp < 0) {
    emlrtErrorWithMessageIdR2018a(&b_st, &c_emlrtRTEI,
      "MATLAB:checkDimCommon:nonnegativeSize",
      "MATLAB:checkDimCommon:nonnegativeSize", 0);
  }

  if (num_idx_0_tmp << 1 != 6) {
    emlrtErrorWithMessageIdR2018a(&b_st, &b_emlrtRTEI,
      "Coder:MATLAB:getReshapeDims_notSameNumel",
      "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }

  for (i = 0; i < num_idx_0_tmp; i++) {
    acc_data[2 * i] = u[i];
    acc_data[2 * i + 1] = u[i + num_idx_0_tmp];
  }

  /* TODO */
  b_params = muDoubleScalarRound(params->ct / params->dt);
  if (b_params < 256.0) {
    if (b_params >= 0.0) {
      b_u = (uint8_T)b_params;
    } else {
      b_u = 0U;
    }
  } else if (b_params >= 256.0) {
    b_u = MAX_uint8_T;
  } else {
    b_u = 0U;
  }

  ret = 0.0;
  emlrtForLoopVectorCheckR2012b(1.0, 1.0, params->h, mxDOUBLE_CLASS, (int32_T)
    params->h, &emlrtRTEI, sp);
  for (h = 0; h < num_idx_0_tmp; h++) {
    i = h + 1;
    if (i > num_idx_0_tmp) {
      emlrtDynamicBoundsCheckR2012b(i, 1, num_idx_0_tmp, &c_emlrtBCI, sp);
    }

    st.site = &b_emlrtRSI;
    if ((1 <= b_u) && (b_u > 254)) {
      b_st.site = &g_emlrtRSI;
      check_forloop_overflow_error(&b_st);
    }

    i = b_u;
    for (k = 0; k < i; k++) {
      acceleration_size[1] = loop_ub;
      pos_tmp = loop_ub << 1;
      if (0 <= pos_tmp - 1) {
        memset(&acceleration_data[0], 0, pos_tmp * sizeof(real_T));
      }

      if (1 > loop_ub) {
        emlrtDynamicBoundsCheckR2012b(1, 1, 0, &d_emlrtBCI, sp);
      }

      if (!b) {
        iv[0] = 2;
        iv1[0] = 2;
        b = true;
      }

      emlrtSubAssignSizeCheckR2012b(&iv[0], 1, &iv1[0], 1, &emlrtECI, sp);
      acceleration_data[0] = acc_data[2 * h];
      acceleration_data[1] = acc_data[2 * h + 1];
      st.site = &c_emlrtRSI;

      /*  dynamics() - dmpc */
      pos_tmp = 2 * loop_ub;
      acceleration_size[0] = 2;
      for (i1 = 0; i1 < pos_tmp; i1++) {
        acceleration_data[i1] *= params->dt;
      }

      emlrtSizeEqCheckNDR2012b(vel_size, acceleration_size, &b_emlrtECI, &st);
      pos_tmp = 2 * b_loop_ub;
      vel_size[0] = 2;
      for (i1 = 0; i1 < pos_tmp; i1++) {
        vel_data[i1] += acceleration_data[i1];
      }

      for (pos_tmp = 0; pos_tmp < loop_ub; pos_tmp++) {
        b_st.site = &h_emlrtRSI;
        if (muDoubleScalarIsNaN(params->predator)) {
          emlrtErrorWithMessageIdR2018a(&b_st, &g_emlrtRTEI,
            "MATLAB:nologicalnan", "MATLAB:nologicalnan", 0);
        }

        if ((params->predator != 0.0) && ((real_T)pos_tmp + 1.0 == params->n)) {
          i1 = pos_tmp + 1;
          if (i1 > b_loop_ub) {
            emlrtDynamicBoundsCheckR2012b(i1, 1, b_loop_ub, &e_emlrtBCI, &st);
          }

          scale = 3.3121686421112381E-170;
          b_params = vel_data[2 * pos_tmp];
          absxk_tmp = muDoubleScalarAbs(b_params);
          if (absxk_tmp > 3.3121686421112381E-170) {
            y = 1.0;
            scale = absxk_tmp;
          } else {
            t = absxk_tmp / 3.3121686421112381E-170;
            y = t * t;
          }

          i1 = 2 * pos_tmp + 1;
          b_absxk_tmp = muDoubleScalarAbs(vel_data[i1]);
          if (b_absxk_tmp > scale) {
            t = scale / b_absxk_tmp;
            y = y * t * t + 1.0;
            scale = b_absxk_tmp;
          } else {
            t = b_absxk_tmp / scale;
            y += t * t;
          }

          y = scale * muDoubleScalarSqrt(y);
          d = params->pFactor * params->vmax;
          if (y > d) {
            i2 = pos_tmp + 1;
            if (i2 > b_loop_ub) {
              emlrtDynamicBoundsCheckR2012b(i2, 1, b_loop_ub, &f_emlrtBCI, &st);
            }

            scale = 3.3121686421112381E-170;
            if (absxk_tmp > 3.3121686421112381E-170) {
              y = 1.0;
              scale = absxk_tmp;
            } else {
              t = absxk_tmp / 3.3121686421112381E-170;
              y = t * t;
            }

            if (b_absxk_tmp > scale) {
              t = scale / b_absxk_tmp;
              y = y * t * t + 1.0;
              scale = b_absxk_tmp;
            } else {
              t = b_absxk_tmp / scale;
              y += t * t;
            }

            y = scale * muDoubleScalarSqrt(y);
            scale = d / y;
            i2 = pos_tmp + 1;
            if (i2 > b_loop_ub) {
              emlrtDynamicBoundsCheckR2012b(i2, 1, b_loop_ub, &g_emlrtBCI, &st);
            }

            absxk_tmp = scale * vel_data[i1];
            i2 = pos_tmp + 1;
            if (i2 > b_loop_ub) {
              emlrtDynamicBoundsCheckR2012b(i2, 1, b_loop_ub, &l_emlrtBCI, &st);
            }

            b_params *= scale;
            vel_data[2 * pos_tmp] = b_params;
            vel_data[i1] = absxk_tmp;
          }
        } else {
          i1 = pos_tmp + 1;
          if (i1 > b_loop_ub) {
            emlrtDynamicBoundsCheckR2012b(i1, 1, b_loop_ub, &h_emlrtBCI, &st);
          }

          scale = 3.3121686421112381E-170;
          b_params = vel_data[2 * pos_tmp];
          absxk_tmp = muDoubleScalarAbs(b_params);
          if (absxk_tmp > 3.3121686421112381E-170) {
            y = 1.0;
            scale = absxk_tmp;
          } else {
            t = absxk_tmp / 3.3121686421112381E-170;
            y = t * t;
          }

          i1 = 2 * pos_tmp + 1;
          b_absxk_tmp = muDoubleScalarAbs(vel_data[i1]);
          if (b_absxk_tmp > scale) {
            t = scale / b_absxk_tmp;
            y = y * t * t + 1.0;
            scale = b_absxk_tmp;
          } else {
            t = b_absxk_tmp / scale;
            y += t * t;
          }

          y = scale * muDoubleScalarSqrt(y);
          if (y > params->vmax) {
            i2 = pos_tmp + 1;
            if (i2 > b_loop_ub) {
              emlrtDynamicBoundsCheckR2012b(i2, 1, b_loop_ub, &i_emlrtBCI, &st);
            }

            scale = 3.3121686421112381E-170;
            if (absxk_tmp > 3.3121686421112381E-170) {
              y = 1.0;
              scale = absxk_tmp;
            } else {
              t = absxk_tmp / 3.3121686421112381E-170;
              y = t * t;
            }

            if (b_absxk_tmp > scale) {
              t = scale / b_absxk_tmp;
              y = y * t * t + 1.0;
              scale = b_absxk_tmp;
            } else {
              t = b_absxk_tmp / scale;
              y += t * t;
            }

            y = scale * muDoubleScalarSqrt(y);
            scale = params->vmax / y;
            i2 = pos_tmp + 1;
            if (i2 > b_loop_ub) {
              emlrtDynamicBoundsCheckR2012b(i2, 1, b_loop_ub, &j_emlrtBCI, &st);
            }

            absxk_tmp = scale * vel_data[i1];
            i2 = pos_tmp + 1;
            if (i2 > b_loop_ub) {
              emlrtDynamicBoundsCheckR2012b(i2, 1, b_loop_ub, &k_emlrtBCI, &st);
            }

            b_params *= scale;
            vel_data[2 * pos_tmp] = b_params;
            vel_data[i1] = absxk_tmp;
          }
        }

        if (*emlrtBreakCheckR2012bFlagVar != 0) {
          emlrtBreakCheckR2012b(&st);
        }
      }

      acceleration_size[0] = 2;
      acceleration_size[1] = b_loop_ub;
      pos_tmp = 2 * b_loop_ub;
      for (i1 = 0; i1 < pos_tmp; i1++) {
        acceleration_data[i1] = params->dt * vel_data[i1];
      }

      emlrtSizeEqCheckNDR2012b(pos_size, acceleration_size, &c_emlrtECI, &st);
      pos_tmp = 2 * loop_ub;
      pos_size[0] = 2;
      for (i1 = 0; i1 < pos_tmp; i1++) {
        pos_data[i1] += acceleration_data[i1];
      }

      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b(sp);
      }
    }

    st.site = &d_emlrtRSI;

    /*  dmpc - fitness */
    /*  asd = average_squared_distance(pos); */
    /*  sp = separation(pos); */
    /*  Obstacle avoidance */
    /*  rects = params.rects; */
    /*  obst = obstacle_cost(pos, rects, params); */
    /*  Target seeking */
    /*  tar = [60; 50]; */
    /*  tgt = target(pos, tar); */
    /*  Predator Avoidance */
    /*  res = predator_avoidance(pos, params); */
    /*  weighted sum */
    /*  res = params.wc * asd ... */
    /*      + params.ws * sp; */
    /*      + params.wo * obst ... */
    /*      + params.wt * tgt; */
    /*  murmuration */
    b_st.site = &i_emlrtRSI;
    t = 0.0;
    for (pos_tmp = 0; pos_tmp <= loop_ub - 2; pos_tmp++) {
      if (1 > loop_ub) {
        emlrtDynamicBoundsCheckR2012b(1, 1, 0, &m_emlrtBCI, &b_st);
      }

      i = pos_tmp + 2;
      if (i > loop_ub) {
        emlrtDynamicBoundsCheckR2012b(i, 1, loop_ub, &n_emlrtBCI, &b_st);
      }

      /*  codegen */
      /*  - average_squared_distance() - dmpc */
      /*  codegen */
      /*  - separation() - dmpc */
      i = 2 * (pos_tmp + 1);
      b_params = pos_data[0] - pos_data[i];
      d = pos_data[1] - pos_data[i + 1];
      scale = b_params * b_params + d * d;
      if (scale != 0.0) {
        absxk_tmp = 1.0 / scale;
      } else {
        absxk_tmp = 1.0E+7;
      }

      t = (t + params->w_m[pos_tmp] * (params->wc * scale)) + params->ws *
        absxk_tmp;
      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b(&b_st);
      }
    }

    ret += t;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  /*      ret = res/params.h + sum(u.^2)/params.h; */
  return ret;
}

/* End of code generation (cost_sum.c) */
