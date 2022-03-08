/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_cost_sum_api.c
 *
 * Code generation for function '_coder_cost_sum_api'
 *
 */

/* Include files */
#include "_coder_cost_sum_api.h"
#include "cost_sum.h"
#include "cost_sum_data.h"
#include "cost_sum_mexutil.h"
#include "dynamics_quadcopter.h"
#include "rt_nonfinite.h"
#include <string.h>

/* Function Declarations */
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  char_T *identifier))[18];
static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[18];
static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *s, const
  char_T *identifier))[96];
static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[96];
static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *params,
  const char_T *identifier, struct0_T *y);
static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, struct0_T *y);
static real_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static void j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[2]);
static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[20]);
static real_T (*m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[18];
static real_T (*n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[96];
static real_T o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);
static void p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[2]);
static void q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[20]);

/* Function Definitions */
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  char_T *identifier))[18]
{
  real_T (*y)[18];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(u), &thisId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[18]
{
  real_T (*y)[18];
  y = m_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *s, const
  char_T *identifier))[96]
{
  real_T (*y)[96];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(s), &thisId);
  emlrtDestroyArray(&s);
  return y;
}
  static real_T (*f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId))[96]
{
  real_T (*y)[96];
  y = n_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *params,
  const char_T *identifier, struct0_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  h_emlrt_marshallIn(sp, emlrtAlias(params), &thisId, y);
  emlrtDestroyArray(&params);
}

static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, struct0_T *y)
{
  emlrtMsgIdentifier thisId;
  static const char * fieldNames[37] = { "n", "h", "dt", "ct", "vmax", "amax",
    "ipos", "ivel", "delta_angle", "t_end", "steps", "m", "L", "dmin", "g",
    "j_r", "Ixx", "Iyy", "Izz", "max_rotor_speed", "max_angular_speed",
    "max_angular_acc", "k_f", "k_m", "wc", "ws", "quad", "cmpc_prediction_model",
    "knn", "turn", "eta", "w_m", "start_turn", "num_leaders", "t_fix",
    "turn_angle", "turn_acc" };

  static const int32_T dims = 0;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b(sp, parentId, u, 37, fieldNames, 0U, &dims);
  thisId.fIdentifier = "n";
  y->n = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 0, "n")),
    &thisId);
  thisId.fIdentifier = "h";
  y->h = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 1, "h")),
    &thisId);
  thisId.fIdentifier = "dt";
  y->dt = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 2,
    "dt")), &thisId);
  thisId.fIdentifier = "ct";
  y->ct = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 3,
    "ct")), &thisId);
  thisId.fIdentifier = "vmax";
  y->vmax = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 4,
    "vmax")), &thisId);
  thisId.fIdentifier = "amax";
  y->amax = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 5,
    "amax")), &thisId);
  thisId.fIdentifier = "ipos";
  j_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 6, "ipos")),
                     &thisId, y->ipos);
  thisId.fIdentifier = "ivel";
  j_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 7, "ivel")),
                     &thisId, y->ivel);
  thisId.fIdentifier = "delta_angle";
  y->delta_angle = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 8, "delta_angle")), &thisId);
  thisId.fIdentifier = "t_end";
  y->t_end = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 9,
    "t_end")), &thisId);
  thisId.fIdentifier = "steps";
  y->steps = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 10,
    "steps")), &thisId);
  thisId.fIdentifier = "m";
  y->m = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 11, "m")),
    &thisId);
  thisId.fIdentifier = "L";
  y->L = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 12, "L")),
    &thisId);
  thisId.fIdentifier = "dmin";
  y->dmin = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 13,
    "dmin")), &thisId);
  thisId.fIdentifier = "g";
  y->g = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 14, "g")),
    &thisId);
  thisId.fIdentifier = "j_r";
  y->j_r = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 15,
    "j_r")), &thisId);
  thisId.fIdentifier = "Ixx";
  y->Ixx = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 16,
    "Ixx")), &thisId);
  thisId.fIdentifier = "Iyy";
  y->Iyy = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 17,
    "Iyy")), &thisId);
  thisId.fIdentifier = "Izz";
  y->Izz = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 18,
    "Izz")), &thisId);
  thisId.fIdentifier = "max_rotor_speed";
  y->max_rotor_speed = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp,
    u, 0, 19, "max_rotor_speed")), &thisId);
  thisId.fIdentifier = "max_angular_speed";
  y->max_angular_speed = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b
    (sp, u, 0, 20, "max_angular_speed")), &thisId);
  thisId.fIdentifier = "max_angular_acc";
  y->max_angular_acc = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp,
    u, 0, 21, "max_angular_acc")), &thisId);
  thisId.fIdentifier = "k_f";
  y->k_f = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 22,
    "k_f")), &thisId);
  thisId.fIdentifier = "k_m";
  y->k_m = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 23,
    "k_m")), &thisId);
  thisId.fIdentifier = "wc";
  y->wc = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 24,
    "wc")), &thisId);
  thisId.fIdentifier = "ws";
  y->ws = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 25,
    "ws")), &thisId);
  thisId.fIdentifier = "quad";
  y->quad = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 26,
    "quad")), &thisId);
  thisId.fIdentifier = "cmpc_prediction_model";
  y->cmpc_prediction_model = i_emlrt_marshallIn(sp, emlrtAlias
    (emlrtGetFieldR2017b(sp, u, 0, 27, "cmpc_prediction_model")), &thisId);
  thisId.fIdentifier = "knn";
  y->knn = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 28,
    "knn")), &thisId);
  thisId.fIdentifier = "turn";
  y->turn = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 29,
    "turn")), &thisId);
  thisId.fIdentifier = "eta";
  y->eta = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 30,
    "eta")), &thisId);
  thisId.fIdentifier = "w_m";
  k_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 31, "w_m")),
                     &thisId, y->w_m);
  thisId.fIdentifier = "start_turn";
  y->start_turn = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    32, "start_turn")), &thisId);
  thisId.fIdentifier = "num_leaders";
  y->num_leaders = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 33, "num_leaders")), &thisId);
  thisId.fIdentifier = "t_fix";
  y->t_fix = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 34,
    "t_fix")), &thisId);
  thisId.fIdentifier = "turn_angle";
  y->turn_angle = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    35, "turn_angle")), &thisId);
  thisId.fIdentifier = "turn_acc";
  y->turn_acc = i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    36, "turn_acc")), &thisId);
  emlrtDestroyArray(&u);
}

static real_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = o_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[2])
{
  p_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[20])
{
  q_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static real_T (*m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[18]
{
  real_T (*ret)[18];
  static const int32_T dims[1] = { 18 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 1U, dims);
  ret = (real_T (*)[18])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static real_T (*n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[96]
{
  real_T (*ret)[96];
  static const int32_T dims[3] = { 1, 12, 8 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 3U, dims);
  ret = (real_T (*)[96])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static void p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[2])
{
  static const int32_T dims[2] = { 1, 2 };

  real_T (*r)[2];
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  r = (real_T (*)[2])emlrtMxGetData(src);
  ret[0] = (*r)[0];
  ret[1] = (*r)[1];
  emlrtDestroyArray(&src);
}

static void q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[20])
{
  static const int32_T dims[2] = { 1, 20 };

  real_T (*r)[20];
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  r = (real_T (*)[20])emlrtMxGetData(src);
  memcpy(&ret[0], &(*r)[0], 20U * sizeof(real_T));
  emlrtDestroyArray(&src);
}

void cost_sum_api(const mxArray * const prhs[3], int32_T nlhs, const mxArray
                  *plhs[1])
{
  const mxArray *prhs_copy_idx_1;
  real_T (*u)[18];
  real_T (*s)[96];
  struct0_T params;
  real_T cost;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  (void)nlhs;
  st.tls = emlrtRootTLSGlobal;
  prhs_copy_idx_1 = emlrtProtectR2012b(prhs[1], 1, false, -1);

  /* Marshall function inputs */
  u = c_emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "u");
  s = e_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_1), "s");
  g_emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "params", &params);

  /* Invoke the target function */
  cost = cost_sum(&st, *u, *s, &params);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(cost);
}

/* End of code generation (_coder_cost_sum_api.c) */
