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
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 12,    /* lineNo */
  "constraints",                       /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/common/constraints.m"/* pathName */
};

static emlrtRSInfo b_emlrtRSI = { 13,  /* lineNo */
  "constraints",                       /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/common/constraints.m"/* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 9,   /* lineNo */
  "u2acc",                             /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/common/u2acc.m"/* pathName */
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
void constraints(const emlrtStack *sp, const real_T u[18], const struct0_T
                 *params, real_T c_data[], int32_T c_size[1])
{
  real_T sz[2];
  int32_T k;
  boolean_T guard1 = false;
  int32_T exitg2;
  boolean_T exitg1;
  real_T b_params;
  int32_T y_size_idx_1;
  int32_T nx;
  real_T y_data[54];
  real_T res_data[18];
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

  /*  contraints - for mpc  */
  /*  1. bound on acceleration magnitude. */
  /*  2. bound on maximum deviation in direction of acceleration.  */
  /*  Input: */
  /*     u    - 3.n.h x 1 vector */
  /*     a    - 3 x n vector */
  /*  Output */
  /*     c    - ~x1 vector. */
  /*  Usama Mehmood - Oct 2019 */
  /*  bound on acceleration magnitude. */
  st.site = &emlrtRSI;

  /*  u2acc - 3D DMPC - convert col vector representation of acc to matrix form.  */
  /*  Input: */
  /*     u    - 3.reh x 1 vector */
  /*  Output */
  /*     acc    - 3xnxh Matric. */
  /*  Usama Mehmood - Oct 2019 */
  b_st.site = &c_emlrtRSI;
  sz[0] = 3.0;
  sz[1] = params->h;
  c_st.site = &d_emlrtRSI;
  k = 0;
  guard1 = false;
  do {
    exitg2 = 0;
    if (k < 2) {
      if ((sz[k] != muDoubleScalarFloor(sz[k])) || muDoubleScalarIsInf(sz[k])) {
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
      if ((sz[k] < -2.147483648E+9) || (sz[k] > 2.147483647E+9)) {
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
    b_params = 3.0 * params->h;
  }

  if (!(b_params <= 2.147483647E+9)) {
    emlrtErrorWithMessageIdR2018a(&c_st, &f_emlrtRTEI, "Coder:MATLAB:pmaxsize",
      "Coder:MATLAB:pmaxsize", 0);
  }

  k = (int32_T)params->h;
  if (k > 18) {
    emlrtErrorWithMessageIdR2018a(&b_st, &c_emlrtRTEI,
      "Coder:toolbox:reshape_emptyReshapeLimit",
      "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }

  if (k < 0) {
    emlrtErrorWithMessageIdR2018a(&b_st, &b_emlrtRTEI,
      "MATLAB:checkDimCommon:nonnegativeSize",
      "MATLAB:checkDimCommon:nonnegativeSize", 0);
  }

  if (3 * k != 18) {
    emlrtErrorWithMessageIdR2018a(&b_st, &emlrtRTEI,
      "Coder:MATLAB:getReshapeDims_notSameNumel",
      "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }

  y_size_idx_1 = (int8_T)k;
  nx = 3 * y_size_idx_1;
  for (k = 0; k < nx; k++) {
    y_data[k] = u[k] * u[k];
  }

  if (y_size_idx_1 == 0) {
    y_size_idx_1 = 0;
  } else {
    for (k = 0; k < y_size_idx_1; k++) {
      nx = k * 3;
      res_data[k] = (y_data[nx] + y_data[nx + 1]) + y_data[nx + 2];
    }
  }

  st.site = &b_emlrtRSI;
  p = false;
  for (k = 0; k < y_size_idx_1; k++) {
    if (p || (res_data[k] < 0.0)) {
      p = true;
    }
  }

  if (p) {
    emlrtErrorWithMessageIdR2018a(&st, &d_emlrtRTEI,
      "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
      4, "sqrt");
  }

  for (k = 0; k < y_size_idx_1; k++) {
    res_data[k] = muDoubleScalarSqrt(res_data[k]);
  }

  k = y_size_idx_1 - 1;
  for (nx = 0; nx <= k; nx++) {
    res_data[nx] -= params->amax;
  }

  c_size[0] = y_size_idx_1;
  if (0 <= y_size_idx_1 - 1) {
    memcpy(&c_data[0], &res_data[0], y_size_idx_1 * sizeof(real_T));
  }

  /*  bound on maximum deviation in direction of acceleration - prediction horizon */
  /*      c_2 = zeros(params.h-1, params.n); */
  /*      for i = 1:params.n */
  /*          for h = 1:params.h-1 */
  /*              u = acc(:,i,h); */
  /*              v = acc(:,i,h+1); */
  /*              c_2(h, i) = angle_vectors(u,v) - params.delta_angle; */
  /*          end */
  /*      end */
  /*  %% bound on maximum deviation in direction of acceleration - first step */
  /*      c_3 = zeros(params.n, 1); */
  /*      for i = 1:params.n */
  /*          u = a(:,i); */
  /*          v = acc(:,i,1); */
  /*          c_3(i) = angle_vectors(u,v) - params.delta_angle; */
  /*      end   */
  /*      c = cat(1, c_1, c_2(:), c_3); */
}

/* End of code generation (constraints.c) */
