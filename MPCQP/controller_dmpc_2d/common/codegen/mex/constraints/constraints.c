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
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 4,     /* lineNo */
  "constraints",                       /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/constraints.m"/* pathName */
};

static emlrtRSInfo b_emlrtRSI = { 5,   /* lineNo */
  "constraints",                       /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/constraints.m"/* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 6,   /* lineNo */
  "u2acc",                             /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPCQP/controller_dmpc_2d/common/u2acc.m"/* pathName */
};

static emlrtRSInfo d_emlrtRSI = { 18,  /* lineNo */
  "reshapeSizeChecks",                 /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/reshapeSizeChecks.m"/* pathName */
};

static emlrtRTEInfo emlrtRTEI = { 59,  /* lineNo */
  23,                                  /* colNo */
  "reshapeSizeChecks",                 /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/reshapeSizeChecks.m"/* pName */
};

static emlrtRTEInfo b_emlrtRTEI = { 57,/* lineNo */
  23,                                  /* colNo */
  "reshapeSizeChecks",                 /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/reshapeSizeChecks.m"/* pName */
};

static emlrtRTEInfo c_emlrtRTEI = { 52,/* lineNo */
  13,                                  /* colNo */
  "reshapeSizeChecks",                 /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/reshapeSizeChecks.m"/* pName */
};

static emlrtRTEInfo d_emlrtRTEI = { 13,/* lineNo */
  9,                                   /* colNo */
  "sqrt",                              /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/elfun/sqrt.m"/* pName */
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

/* Function Definitions */
void constraints(const emlrtStack *sp, const real_T u[6], const struct0_T
                 *params, real_T c_data[], int32_T c_size[1])
{
  real_T varargin_1[2];
  int32_T k;
  boolean_T guard1 = false;
  int32_T exitg2;
  boolean_T exitg1;
  real_T b_params;
  int32_T xoffset;
  int32_T nx;
  real_T acc_data[12];
  int32_T y_size_idx_1;
  real_T y_data[12];
  boolean_T p;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;

  /*  constraints - non-linear constraints to enforce magnitude of acceleration */
  /*  */
  st.site = &emlrtRSI;

  /*  u2acc - dmpc - convert column vector u to acc matric */
  /*  Output: */
  /*    - acc            % 2 x h - The sequence of control action over the horizon */
  /*  */
  b_st.site = &c_emlrtRSI;
  varargin_1[0] = params->h;
  varargin_1[1] = 2.0;
  c_st.site = &d_emlrtRSI;
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

  xoffset = (int32_T)params->h;
  if (xoffset > 6) {
    emlrtErrorWithMessageIdR2018a(&b_st, &c_emlrtRTEI,
      "Coder:toolbox:reshape_emptyReshapeLimit",
      "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }

  if (xoffset < 0) {
    emlrtErrorWithMessageIdR2018a(&b_st, &b_emlrtRTEI,
      "MATLAB:checkDimCommon:nonnegativeSize",
      "MATLAB:checkDimCommon:nonnegativeSize", 0);
  }

  if (xoffset << 1 != 6) {
    emlrtErrorWithMessageIdR2018a(&b_st, &emlrtRTEI,
      "Coder:MATLAB:getReshapeDims_notSameNumel",
      "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }

  for (nx = 0; nx < xoffset; nx++) {
    acc_data[2 * nx] = u[nx];
    acc_data[2 * nx + 1] = u[nx + xoffset];
  }

  y_size_idx_1 = (int8_T)xoffset;
  nx = (int8_T)xoffset << 1;
  for (k = 0; k < nx; k++) {
    y_data[k] = acc_data[k] * acc_data[k];
  }

  if (y_size_idx_1 == 0) {
    c_data[0] = 0.0;
    c_data[1] = 0.0;
  } else {
    c_data[0] = y_data[0];
    c_data[1] = y_data[1];
    for (k = 2; k <= y_size_idx_1; k++) {
      xoffset = (k - 1) << 1;
      c_data[0] += y_data[xoffset];
      c_data[1] += y_data[xoffset + 1];
    }
  }

  st.site = &b_emlrtRSI;
  p = false;
  if ((c_data[0] < 0.0) || (c_data[1] < 0.0)) {
    p = true;
  }

  if (p) {
    emlrtErrorWithMessageIdR2018a(&st, &d_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  c_data[0] = muDoubleScalarSqrt(c_data[0]);
  c_data[1] = muDoubleScalarSqrt(c_data[1]);
  c_size[0] = 2;
  c_data[0] -= params->amax;
  c_data[1] -= params->amax;
}

/* End of code generation (constraints.c) */
