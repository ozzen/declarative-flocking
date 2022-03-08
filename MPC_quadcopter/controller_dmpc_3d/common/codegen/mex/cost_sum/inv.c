/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * inv.c
 *
 * Code generation for function 'inv'
 *
 */

/* Include files */
#include "inv.h"
#include "cost_sum.h"
#include "cost_sum_data.h"
#include "eml_int_forloop_overflow_check.h"
#include "mwmathutil.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo x_emlrtRSI = { 173, /* lineNo */
  "invNxN",                            /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/inv.m"/* pathName */
};

static emlrtRSInfo y_emlrtRSI = { 190, /* lineNo */
  "invNxN",                            /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/inv.m"/* pathName */
};

static emlrtRSInfo ab_emlrtRSI = { 30, /* lineNo */
  "xgetrf",                            /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+lapack/xgetrf.m"/* pathName */
};

static emlrtRSInfo bb_emlrtRSI = { 50, /* lineNo */
  "xzgetrf",                           /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+reflapack/xzgetrf.m"/* pathName */
};

static emlrtRSInfo cb_emlrtRSI = { 58, /* lineNo */
  "xzgetrf",                           /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+reflapack/xzgetrf.m"/* pathName */
};

static emlrtRSInfo db_emlrtRSI = { 45, /* lineNo */
  "xgeru",                             /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+blas/xgeru.m"/* pathName */
};

static emlrtRSInfo eb_emlrtRSI = { 45, /* lineNo */
  "xger",                              /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+blas/xger.m"/* pathName */
};

static emlrtRSInfo fb_emlrtRSI = { 15, /* lineNo */
  "xger",                              /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+refblas/xger.m"/* pathName */
};

static emlrtRSInfo gb_emlrtRSI = { 54, /* lineNo */
  "xgerx",                             /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+refblas/xgerx.m"/* pathName */
};

static emlrtRSInfo hb_emlrtRSI = { 41, /* lineNo */
  "xgerx",                             /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+refblas/xgerx.m"/* pathName */
};

static emlrtRSInfo ib_emlrtRSI = { 59, /* lineNo */
  "xtrsm",                             /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+blas/xtrsm.m"/* pathName */
};

static emlrtRSInfo jb_emlrtRSI = { 51, /* lineNo */
  "xtrsm",                             /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/+refblas/xtrsm.m"/* pathName */
};

/* Function Definitions */
void invNxN(const emlrtStack *sp, const real_T x[16], real_T y[16])
{
  real_T b_x[16];
  int8_T ipiv[4];
  int32_T j;
  int32_T b;
  int32_T jj;
  int8_T p[4];
  int32_T jp1j;
  int32_T n;
  int32_T jy;
  int32_T ix;
  real_T smax;
  int32_T iy;
  real_T s;
  int8_T i;
  int32_T b_i;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack h_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  f_st.prev = &e_st;
  f_st.tls = e_st.tls;
  g_st.prev = &f_st;
  g_st.tls = f_st.tls;
  h_st.prev = &g_st;
  h_st.tls = g_st.tls;
  memset(&y[0], 0, 16U * sizeof(real_T));
  st.site = &x_emlrtRSI;
  memcpy(&b_x[0], &x[0], 16U * sizeof(real_T));
  b_st.site = &ab_emlrtRSI;
  ipiv[0] = 1;
  ipiv[1] = 2;
  ipiv[2] = 3;
  for (j = 0; j < 3; j++) {
    b = j * 5;
    jj = j * 5;
    jp1j = b + 2;
    n = 4 - j;
    jy = 0;
    ix = b;
    smax = muDoubleScalarAbs(b_x[jj]);
    for (iy = 2; iy <= n; iy++) {
      ix++;
      s = muDoubleScalarAbs(b_x[ix]);
      if (s > smax) {
        jy = iy - 1;
        smax = s;
      }
    }

    if (b_x[jj + jy] != 0.0) {
      if (jy != 0) {
        iy = j + jy;
        ipiv[j] = (int8_T)(iy + 1);
        smax = b_x[j];
        b_x[j] = b_x[iy];
        b_x[iy] = smax;
        ix = j + 4;
        iy += 4;
        smax = b_x[ix];
        b_x[ix] = b_x[iy];
        b_x[iy] = smax;
        ix += 4;
        iy += 4;
        smax = b_x[ix];
        b_x[ix] = b_x[iy];
        b_x[iy] = smax;
        ix += 4;
        iy += 4;
        smax = b_x[ix];
        b_x[ix] = b_x[iy];
        b_x[iy] = smax;
      }

      jy = (jj - j) + 4;
      c_st.site = &bb_emlrtRSI;
      for (b_i = jp1j; b_i <= jy; b_i++) {
        b_x[b_i - 1] /= b_x[jj];
      }
    }

    n = 2 - j;
    jy = b + 4;
    c_st.site = &cb_emlrtRSI;
    d_st.site = &db_emlrtRSI;
    e_st.site = &eb_emlrtRSI;
    f_st.site = &fb_emlrtRSI;
    iy = jj + 6;
    g_st.site = &hb_emlrtRSI;
    for (jp1j = 0; jp1j <= n; jp1j++) {
      smax = b_x[jy];
      if (b_x[jy] != 0.0) {
        ix = jj + 1;
        b = (iy - j) + 2;
        g_st.site = &gb_emlrtRSI;
        if ((iy <= b) && (b > 2147483646)) {
          h_st.site = &p_emlrtRSI;
          check_forloop_overflow_error(&h_st);
        }

        for (b_i = iy; b_i <= b; b_i++) {
          b_x[b_i - 1] += b_x[ix] * -smax;
          ix++;
        }
      }

      jy += 4;
      iy += 4;
    }
  }

  p[0] = 1;
  p[1] = 2;
  p[2] = 3;
  p[3] = 4;
  if (ipiv[0] > 1) {
    jy = ipiv[0] - 1;
    iy = p[jy];
    p[jy] = 1;
    p[0] = (int8_T)iy;
  }

  if (ipiv[1] > 2) {
    jy = ipiv[1] - 1;
    iy = p[jy];
    p[jy] = p[1];
    p[1] = (int8_T)iy;
  }

  if (ipiv[2] > 3) {
    jy = ipiv[2] - 1;
    iy = p[jy];
    p[jy] = p[2];
    p[2] = (int8_T)iy;
  }

  i = p[0];
  y[(p[0] - 1) << 2] = 1.0;
  for (j = 1; j < 5; j++) {
    if (y[(j + ((i - 1) << 2)) - 1] != 0.0) {
      iy = j + 1;
      for (b_i = iy; b_i < 5; b_i++) {
        jp1j = (b_i + ((i - 1) << 2)) - 1;
        y[jp1j] -= y[(j + ((i - 1) << 2)) - 1] * b_x[(b_i + ((j - 1) << 2)) - 1];
      }
    }
  }

  i = p[1];
  y[((p[1] - 1) << 2) + 1] = 1.0;
  for (j = 2; j < 5; j++) {
    if (y[(j + ((i - 1) << 2)) - 1] != 0.0) {
      iy = j + 1;
      for (b_i = iy; b_i < 5; b_i++) {
        jp1j = (b_i + ((i - 1) << 2)) - 1;
        y[jp1j] -= y[(j + ((i - 1) << 2)) - 1] * b_x[(b_i + ((j - 1) << 2)) - 1];
      }
    }
  }

  i = p[2];
  y[((p[2] - 1) << 2) + 2] = 1.0;
  for (j = 3; j < 5; j++) {
    jy = (j + ((i - 1) << 2)) - 1;
    if (y[jy] != 0.0) {
      iy = j + 1;
      for (b_i = iy; b_i < 5; b_i++) {
        jp1j = ((i - 1) << 2) + 3;
        y[jp1j] -= y[jy] * b_x[((j - 1) << 2) + 3];
      }
    }
  }

  y[((p[3] - 1) << 2) + 3] = 1.0;
  st.site = &y_emlrtRSI;
  b_st.site = &ib_emlrtRSI;
  for (j = 0; j < 4; j++) {
    jy = j << 2;
    smax = y[jy + 3];
    if (smax != 0.0) {
      y[jy + 3] = smax / b_x[15];
      c_st.site = &jb_emlrtRSI;
      for (b_i = 0; b_i < 3; b_i++) {
        jp1j = b_i + jy;
        y[jp1j] -= y[jy + 3] * b_x[b_i + 12];
      }
    }

    smax = y[jy + 2];
    if (smax != 0.0) {
      y[jy + 2] = smax / b_x[10];
      c_st.site = &jb_emlrtRSI;
      for (b_i = 0; b_i < 2; b_i++) {
        jp1j = b_i + jy;
        y[jp1j] -= y[jy + 2] * b_x[b_i + 8];
      }
    }

    smax = y[jy + 1];
    if (smax != 0.0) {
      y[jy + 1] = smax / b_x[5];
      c_st.site = &jb_emlrtRSI;
      for (b_i = 0; b_i < 1; b_i++) {
        y[jy] -= y[jy + 1] * b_x[4];
      }
    }

    if (y[jy] != 0.0) {
      y[jy] /= b_x[0];
      c_st.site = &jb_emlrtRSI;
    }
  }
}

/* End of code generation (inv.c) */
