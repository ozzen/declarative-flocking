/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * dynamics_quadcopter.c
 *
 * Code generation for function 'dynamics_quadcopter'
 *
 */

/* Include files */
#include "dynamics_quadcopter.h"
#include "cost_sum.h"
#include "cost_sum_data.h"
#include "cost_sum_mexutil.h"
#include "inv.h"
#include "mwmathutil.h"
#include "rt_nonfinite.h"
#include "trim_vec.h"
#include "warning.h"

/* Variable Definitions */
static emlrtRSInfo r_emlrtRSI = { 45,  /* lineNo */
  "mpower",                            /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/ops/mpower.m"/* pathName */
};

static emlrtRSInfo t_emlrtRSI = { 9,   /* lineNo */
  "dynamics_quadcopter",               /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/Quadcopter_model_new/dynamics_quadcopter.m"/* pathName */
};

static emlrtRSInfo u_emlrtRSI = { 12,  /* lineNo */
  "dynamics_quadcopter",               /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/Quadcopter_model_new/dynamics_quadcopter.m"/* pathName */
};

static emlrtRSInfo v_emlrtRSI = { 21,  /* lineNo */
  "inv",                               /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/inv.m"/* pathName */
};

static emlrtRSInfo w_emlrtRSI = { 22,  /* lineNo */
  "inv",                               /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/inv.m"/* pathName */
};

static emlrtRSInfo kb_emlrtRSI = { 42, /* lineNo */
  "checkcond",                         /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/inv.m"/* pathName */
};

static emlrtRSInfo lb_emlrtRSI = { 46, /* lineNo */
  "checkcond",                         /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/matfun/inv.m"/* pathName */
};

static emlrtRSInfo mb_emlrtRSI = { 6,  /* lineNo */
  "rotor_satuation",                   /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/Quadcopter_model_new/rotor_satuation.m"/* pathName */
};

static emlrtRSInfo nb_emlrtRSI = { 7,  /* lineNo */
  "rotor_satuation",                   /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/Quadcopter_model_new/rotor_satuation.m"/* pathName */
};

static emlrtMCInfo c_emlrtMCI = { 53,  /* lineNo */
  19,                                  /* colNo */
  "flt2str",                           /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/flt2str.m"/* pName */
};

static emlrtRSInfo pb_emlrtRSI = { 53, /* lineNo */
  "flt2str",                           /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/flt2str.m"/* pathName */
};

/* Function Declarations */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[14]);
static const mxArray *b_sprintf(const emlrtStack *sp, const mxArray *b, const
  mxArray *c, emlrtMCInfo *location);
static void emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *a__output_of_sprintf_, const char_T *identifier, char_T y[14]);
static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[14]);

/* Function Definitions */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[14])
{
  l_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static const mxArray *b_sprintf(const emlrtStack *sp, const mxArray *b, const
  mxArray *c, emlrtMCInfo *location)
{
  const mxArray *pArrays[2];
  const mxArray *m;
  pArrays[0] = b;
  pArrays[1] = c;
  return emlrtCallMATLABR2012b(sp, 1, &m, 2, pArrays, "sprintf", true, location);
}

static void emlrt_marshallIn(const emlrtStack *sp, const mxArray
  *a__output_of_sprintf_, const char_T *identifier, char_T y[14])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  b_emlrt_marshallIn(sp, emlrtAlias(a__output_of_sprintf_), &thisId, y);
  emlrtDestroyArray(&a__output_of_sprintf_);
}

static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[14])
{
  static const int32_T dims[2] = { 1, 14 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "char", false, 2U, dims);
  emlrtImportCharArrayR2015b(sp, src, &ret[0], 14);
  emlrtDestroyArray(&src);
}

void dynamics_quadcopter(const emlrtStack *sp, real_T s[12], const real_T u[4],
  const struct0_T *params, real_T omegas[4])
{
  real_T M[16];
  real_T b_s;
  real_T rc;
  real_T a[16];
  real_T n1x;
  int32_T j;
  boolean_T exitg1;
  real_T n1xinv;
  int32_T s_tmp;
  const mxArray *y;
  const mxArray *m;
  static const int32_T iv[2] = { 1, 6 };

  static const char_T rfmt[6] = { '%', '1', '4', '.', '6', 'e' };

  char_T str[14];
  real_T omegas_sq[4];
  real_T s_new_tmp;
  real_T s_new_idx_3;
  real_T b_s_new_tmp;
  real_T c_s_new_tmp;
  real_T dds[6];
  real_T ds[6];
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;

  /*  Usama Mehmood - Oct 2019 */
  M[0] = params->k_f;
  M[4] = params->k_f;
  M[8] = params->k_f;
  M[12] = params->k_f;
  M[1] = 0.0;
  b_s = params->k_f * params->L;
  M[5] = b_s;
  M[9] = 0.0;
  rc = -params->k_f * params->L;
  M[13] = rc;
  M[2] = rc;
  M[6] = 0.0;
  M[10] = b_s;
  M[14] = 0.0;
  M[3] = params->k_m;
  M[7] = -params->k_m;
  M[11] = params->k_m;
  M[15] = -params->k_m;
  st.site = &t_emlrtRSI;
  b_st.site = &v_emlrtRSI;
  invNxN(&b_st, M, a);
  b_st.site = &w_emlrtRSI;
  n1x = 0.0;
  j = 0;
  exitg1 = false;
  while ((!exitg1) && (j < 4)) {
    s_tmp = j << 2;
    b_s = ((muDoubleScalarAbs(M[s_tmp]) + muDoubleScalarAbs(M[s_tmp + 1])) +
           muDoubleScalarAbs(M[s_tmp + 2])) + muDoubleScalarAbs(M[s_tmp + 3]);
    if (muDoubleScalarIsNaN(b_s)) {
      n1x = rtNaN;
      exitg1 = true;
    } else {
      if (b_s > n1x) {
        n1x = b_s;
      }

      j++;
    }
  }

  n1xinv = 0.0;
  j = 0;
  exitg1 = false;
  while ((!exitg1) && (j < 4)) {
    s_tmp = j << 2;
    b_s = ((muDoubleScalarAbs(a[s_tmp]) + muDoubleScalarAbs(a[s_tmp + 1])) +
           muDoubleScalarAbs(a[s_tmp + 2])) + muDoubleScalarAbs(a[s_tmp + 3]);
    if (muDoubleScalarIsNaN(b_s)) {
      n1xinv = rtNaN;
      exitg1 = true;
    } else {
      if (b_s > n1xinv) {
        n1xinv = b_s;
      }

      j++;
    }
  }

  rc = 1.0 / (n1x * n1xinv);
  if ((n1x == 0.0) || (n1xinv == 0.0) || (rc == 0.0)) {
    c_st.site = &kb_emlrtRSI;
    warning(&c_st);
  } else {
    if (muDoubleScalarIsNaN(rc) || (rc < 2.2204460492503131E-16)) {
      c_st.site = &lb_emlrtRSI;
      y = NULL;
      m = emlrtCreateCharArray(2, iv);
      emlrtInitCharArrayR2013a(&c_st, 6, m, &rfmt[0]);
      emlrtAssign(&y, m);
      d_st.site = &pb_emlrtRSI;
      emlrt_marshallIn(&d_st, b_sprintf(&d_st, y, emlrt_marshallOut(rc),
        &c_emlrtMCI), "<output of sprintf>", str);
      c_st.site = &lb_emlrtRSI;
      b_warning(&c_st, str);
    }
  }

  for (j = 0; j < 4; j++) {
    omegas_sq[j] = ((a[j] * u[0] + a[j + 4] * u[1]) + a[j + 8] * u[2]) + a[j +
      12] * u[3];
  }

  /*      %% Actuator saturation. */
  st.site = &u_emlrtRSI;
  b_st.site = &mb_emlrtRSI;
  omegas[0] = muDoubleScalarSign(omegas_sq[0]) * muDoubleScalarMin
    (muDoubleScalarSqrt(muDoubleScalarAbs(omegas_sq[0])),
     params->max_rotor_speed);
  b_st.site = &nb_emlrtRSI;
  c_st.site = &r_emlrtRSI;
  if (*emlrtBreakCheckR2012bFlagVar != 0) {
    emlrtBreakCheckR2012b(&st);
  }

  b_st.site = &mb_emlrtRSI;
  omegas[1] = muDoubleScalarSign(omegas_sq[1]) * muDoubleScalarMin
    (muDoubleScalarSqrt(muDoubleScalarAbs(omegas_sq[1])),
     params->max_rotor_speed);
  b_st.site = &nb_emlrtRSI;
  c_st.site = &r_emlrtRSI;
  if (*emlrtBreakCheckR2012bFlagVar != 0) {
    emlrtBreakCheckR2012b(&st);
  }

  b_st.site = &mb_emlrtRSI;
  omegas[2] = muDoubleScalarSign(omegas_sq[2]) * muDoubleScalarMin
    (muDoubleScalarSqrt(muDoubleScalarAbs(omegas_sq[2])),
     params->max_rotor_speed);
  b_st.site = &nb_emlrtRSI;
  c_st.site = &r_emlrtRSI;
  if (*emlrtBreakCheckR2012bFlagVar != 0) {
    emlrtBreakCheckR2012b(&st);
  }

  b_st.site = &mb_emlrtRSI;
  omegas[3] = muDoubleScalarSign(omegas_sq[3]) * muDoubleScalarMin
    (muDoubleScalarSqrt(muDoubleScalarAbs(omegas_sq[3])),
     params->max_rotor_speed);
  b_st.site = &nb_emlrtRSI;
  c_st.site = &r_emlrtRSI;
  if (*emlrtBreakCheckR2012bFlagVar != 0) {
    emlrtBreakCheckR2012b(&st);
  }

  /*       u = M * omegas_sq; */
  /*       disp(omegas) */
  /*      %% */
  /*  acceleration - angles are in radians, omega in rad/s, alpha in rad/s^2 */
  /*  Usama Mehmood - Oct 2019 */
  /* roll */
  /* pitch */
  /* yaw */
  /*  */
  b_s = muDoubleScalarSin(s[8]);
  rc = muDoubleScalarSin(s[6]);
  n1x = muDoubleScalarCos(s[8]);
  n1xinv = muDoubleScalarCos(s[6]);
  s_new_tmp = u[0] / params->m;
  s_new_idx_3 = (s[10] * s[11] * ((params->Iyy - params->Izz) / params->Ixx) +
                 params->j_r * s[10] * (((omegas[0] + -omegas[1]) + omegas[2]) +
    -omegas[3]) / params->Ixx) + u[1] / params->Ixx;
  b_s_new_tmp = (s[11] * s[9] * ((params->Izz - params->Ixx) / params->Iyy) +
                 params->j_r * s[9] * (((-omegas[0] + omegas[1]) + -omegas[2]) +
    omegas[3]) / params->Iyy) + u[2] / params->Iyy;
  c_s_new_tmp = s[9] * s[10] * ((params->Ixx - params->Iyy) / params->Izz) + u[3]
    / params->Izz;

  /*  */
  /*      dds = [(cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi)) * (u(1)/params.m);... */
  /*             (cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi)) * (u(1)/params.m);... */
  /*             -params.g + cos(phi)*cos(theta)*(u(1)/params.m);... */
  /*             u(2)/params.Ixx;... */
  /*             u(3)/params.Iyy;... */
  /*             u(4)/params.Izz]; */
  dds[0] = (n1xinv * muDoubleScalarSin(s[7]) * n1x + rc * b_s) * s_new_tmp;
  dds[1] = (muDoubleScalarCos(s[6]) * muDoubleScalarSin(s[7]) * b_s - rc * n1x) *
    s_new_tmp;
  dds[2] = -params->g + n1xinv * muDoubleScalarCos(s[7]) * s_new_tmp;
  dds[3] = muDoubleScalarSign(s_new_idx_3) * muDoubleScalarMin(muDoubleScalarAbs
    (s_new_idx_3), params->max_angular_acc);
  dds[4] = muDoubleScalarSign(b_s_new_tmp) * muDoubleScalarMin(muDoubleScalarAbs
    (b_s_new_tmp), params->max_angular_acc);
  dds[5] = muDoubleScalarSign(c_s_new_tmp) * muDoubleScalarMin(muDoubleScalarAbs
    (c_s_new_tmp), params->max_angular_acc);

  /* update first-derivative [vx, vy, ...] */
  for (j = 0; j < 6; j++) {
    dds[j] *= params->dt;
  }

  b_s = s[3] + dds[0];
  ds[0] = b_s;
  ds[3] = s[9] + dds[3];
  s[3] = b_s;
  b_s = s[4] + dds[1];
  ds[1] = b_s;
  ds[4] = s[10] + dds[4];
  s[4] = b_s;
  b_s = s[5] + dds[2];
  ds[2] = b_s;
  ds[5] = s[11] + dds[5];
  s[5] = b_s;

  /*       s(4:6) = sign(s(4:6)).*min(abs(s(4:6)),params.vmax); */
  trim_vec(*(real_T (*)[3])&s[3], params->vmax);
  s[9] = ds[3];
  s[9] = muDoubleScalarSign(s[9]) * muDoubleScalarMin(muDoubleScalarAbs(s[9]),
    params->max_angular_speed);
  s[10] = ds[4];
  s[10] = muDoubleScalarSign(s[10]) * muDoubleScalarMin(muDoubleScalarAbs(s[10]),
    params->max_angular_speed);
  s[11] = ds[5];
  s[11] = muDoubleScalarSign(s[11]) * muDoubleScalarMin(muDoubleScalarAbs(s[11]),
    params->max_angular_speed);

  /* update zeroth-derivative [x, y, ...] */
  for (j = 0; j < 6; j++) {
    ds[j] *= params->dt;
  }

  b_s = s[6] + ds[3];
  s[0] += ds[0];
  s[6] = b_s;
  b_s = s[7] + ds[4];
  s[1] += ds[1];
  s[7] = b_s;
  b_s = s[8] + ds[5];
  s[2] += ds[2];
  s[8] = b_s;

  /* reset angles */
  /*       s(7:9) = mod(s(7:9), 2*pi); */
}

/* End of code generation (dynamics_quadcopter.c) */
