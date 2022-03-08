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
#include "rt_nonfinite.h"
#include <string.h>

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[6];
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *pos,
  const char_T *identifier))[30];
static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[30];
static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *params,
  const char_T *identifier, struct0_T *y);
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  char_T *identifier))[6];
static const mxArray *emlrt_marshallOut(const real_T u);
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, struct0_T *y);
static real_T g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[15]);
static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, struct1_T y[5]);
static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[6];
static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[30];
static real_T l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);
static void m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[15]);

/* Function Definitions */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[6]
{
  real_T (*y)[6];
  y = j_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *pos,
  const char_T *identifier))[30]
{
  real_T (*y)[30];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(pos), &thisId);
  emlrtDestroyArray(&pos);
  return y;
}

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[30]
{
  real_T (*y)[30];
  y = k_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *params,
  const char_T *identifier, struct0_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  f_emlrt_marshallIn(sp, emlrtAlias(params), &thisId, y);
  emlrtDestroyArray(&params);
}

static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  char_T *identifier))[6]
{
  real_T (*y)[6];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(u), &thisId);
  emlrtDestroyArray(&u);
  return y;
}
  static const mxArray *emlrt_marshallOut(const real_T u)
{
  const mxArray *y;
  const mxArray *m;
  y = NULL;
  m = emlrtCreateDoubleScalar(u);
  emlrtAssign(&y, m);
  return y;
}

static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, struct0_T *y)
{
  emlrtMsgIdentifier thisId;
  static const char * fieldNames[29] = { "predator", "pFactor", "pred_radius",
    "wp", "steps", "n", "h", "dt", "ct", "amax", "vmax", "minP", "maxP", "minV",
    "maxV", "Dmax", "dmin", "Ds", "Dc", "gamma", "rs", "ws", "wc", "wo", "wt",
    "nn", "knn", "w_m", "rects" };

  static const int32_T dims = 0;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b(sp, parentId, u, 29, fieldNames, 0U, &dims);
  thisId.fIdentifier = "predator";
  y->predator = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    0, "predator")), &thisId);
  thisId.fIdentifier = "pFactor";
  y->pFactor = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 1,
    "pFactor")), &thisId);
  thisId.fIdentifier = "pred_radius";
  y->pred_radius = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 2, "pred_radius")), &thisId);
  thisId.fIdentifier = "wp";
  y->wp = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 3,
    "wp")), &thisId);
  thisId.fIdentifier = "steps";
  y->steps = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 4,
    "steps")), &thisId);
  thisId.fIdentifier = "n";
  y->n = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 5, "n")),
    &thisId);
  thisId.fIdentifier = "h";
  y->h = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 6, "h")),
    &thisId);
  thisId.fIdentifier = "dt";
  y->dt = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 7,
    "dt")), &thisId);
  thisId.fIdentifier = "ct";
  y->ct = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 8,
    "ct")), &thisId);
  thisId.fIdentifier = "amax";
  y->amax = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 9,
    "amax")), &thisId);
  thisId.fIdentifier = "vmax";
  y->vmax = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 10,
    "vmax")), &thisId);
  thisId.fIdentifier = "minP";
  y->minP = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 11,
    "minP")), &thisId);
  thisId.fIdentifier = "maxP";
  y->maxP = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 12,
    "maxP")), &thisId);
  thisId.fIdentifier = "minV";
  y->minV = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 13,
    "minV")), &thisId);
  thisId.fIdentifier = "maxV";
  y->maxV = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 14,
    "maxV")), &thisId);
  thisId.fIdentifier = "Dmax";
  y->Dmax = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 15,
    "Dmax")), &thisId);
  thisId.fIdentifier = "dmin";
  y->dmin = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 16,
    "dmin")), &thisId);
  thisId.fIdentifier = "Ds";
  y->Ds = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 17,
    "Ds")), &thisId);
  thisId.fIdentifier = "Dc";
  y->Dc = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 18,
    "Dc")), &thisId);
  thisId.fIdentifier = "gamma";
  y->gamma = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 19,
    "gamma")), &thisId);
  thisId.fIdentifier = "rs";
  y->rs = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 20,
    "rs")), &thisId);
  thisId.fIdentifier = "ws";
  y->ws = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 21,
    "ws")), &thisId);
  thisId.fIdentifier = "wc";
  y->wc = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 22,
    "wc")), &thisId);
  thisId.fIdentifier = "wo";
  y->wo = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 23,
    "wo")), &thisId);
  thisId.fIdentifier = "wt";
  y->wt = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 24,
    "wt")), &thisId);
  thisId.fIdentifier = "nn";
  y->nn = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 25,
    "nn")), &thisId);
  thisId.fIdentifier = "knn";
  y->knn = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 26,
    "knn")), &thisId);
  thisId.fIdentifier = "w_m";
  h_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 27, "w_m")),
                     &thisId, y->w_m);
  thisId.fIdentifier = "rects";
  i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 28, "rects")),
                     &thisId, y->rects);
  emlrtDestroyArray(&u);
}

static real_T g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = l_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[15])
{
  m_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, struct1_T y[5])
{
  emlrtMsgIdentifier thisId;
  static const char * fieldNames[5] = { "cx", "cy", "w", "l", "theta" };

  static const int32_T dims[2] = { 1, 5 };

  int32_T i;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b(sp, parentId, u, 5, fieldNames, 2U, dims);
  for (i = 0; i < 5; i++) {
    thisId.fIdentifier = "cx";
    y[i].cx = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, i, 0,
      "cx")), &thisId);
    thisId.fIdentifier = "cy";
    y[i].cy = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, i, 1,
      "cy")), &thisId);
    thisId.fIdentifier = "w";
    y[i].w = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, i, 2,
      "w")), &thisId);
    thisId.fIdentifier = "l";
    y[i].l = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, i, 3,
      "l")), &thisId);
    thisId.fIdentifier = "theta";
    y[i].theta = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, i,
      4, "theta")), &thisId);
  }

  emlrtDestroyArray(&u);
}

static real_T (*j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[6]
{
  real_T (*ret)[6];
  static const int32_T dims[1] = { 6 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 1U, dims);
  ret = (real_T (*)[6])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static real_T (*k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[30]
{
  real_T (*ret)[30];
  static const int32_T dims[2] = { 2, 15 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[30])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static void m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[15])
{
  static const int32_T dims[2] = { 1, 15 };

  real_T (*r)[15];
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  r = (real_T (*)[15])emlrtMxGetData(src);
  memcpy(&ret[0], &(*r)[0], 15U * sizeof(real_T));
  emlrtDestroyArray(&src);
}

void cost_sum_api(const mxArray * const prhs[4], int32_T nlhs, const mxArray
                  *plhs[1])
{
  real_T (*u)[6];
  real_T (*pos)[30];
  real_T (*vel)[30];
  struct0_T params;
  real_T ret;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  (void)nlhs;
  st.tls = emlrtRootTLSGlobal;

  /* Marshall function inputs */
  u = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "u");
  pos = c_emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "pos");
  vel = c_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "vel");
  e_emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "params", &params);

  /* Invoke the target function */
  ret = cost_sum(&st, *u, *pos, *vel, &params);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(ret);
}

/* End of code generation (_coder_cost_sum_api.c) */
