/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_constraints_api.c
 *
 * Code generation for function '_coder_constraints_api'
 *
 */

/* Include files */
#include "_coder_constraints_api.h"
#include "constraints.h"
#include "constraints_data.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[24];
static const mxArray *b_emlrt_marshallOut(void);
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *params,
  const char_T *identifier, struct0_T *y);
static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, struct0_T *y);
static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  char_T *identifier))[24];
static const mxArray *emlrt_marshallOut(const real_T u_data[], const int32_T
  u_size[1]);
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, struct1_T y[5]);
static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[24];
static real_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);

/* Function Definitions */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[24]
{
  real_T (*y)[24];
  y = g_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static const mxArray *b_emlrt_marshallOut(void)
{
  const mxArray *y;
  const mxArray *m;
  static const int32_T iv[2] = { 0, 0 };

  static const int32_T iv1[2] = { 0, 0 };

  y = NULL;
  m = emlrtCreateNumericArray(2, iv, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, NULL);
  emlrtSetDimensions((mxArray *)m, iv1, 2);
  emlrtAssign(&y, m);
  return y;
}

static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *params,
  const char_T *identifier, struct0_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  d_emlrt_marshallIn(sp, emlrtAlias(params), &thisId, y);
  emlrtDestroyArray(&params);
}

static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, struct0_T *y)
{
  emlrtMsgIdentifier thisId;
  static const char * fieldNames[26] = { "predator", "pFactor", "pred_radius",
    "wp", "steps", "n", "h", "dt", "ct", "amax", "vmax", "minP", "maxP", "minV",
    "maxV", "Dmax", "dmin", "Ds", "Dc", "gamma", "rs", "ws", "wc", "wo", "wt",
    "rects" };

  static const int32_T dims = 0;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b(sp, parentId, u, 26, fieldNames, 0U, &dims);
  thisId.fIdentifier = "predator";
  y->predator = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0,
    0, "predator")), &thisId);
  thisId.fIdentifier = "pFactor";
  y->pFactor = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 1,
    "pFactor")), &thisId);
  thisId.fIdentifier = "pred_radius";
  y->pred_radius = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u,
    0, 2, "pred_radius")), &thisId);
  thisId.fIdentifier = "wp";
  y->wp = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 3,
    "wp")), &thisId);
  thisId.fIdentifier = "steps";
  y->steps = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 4,
    "steps")), &thisId);
  thisId.fIdentifier = "n";
  y->n = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 5, "n")),
    &thisId);
  thisId.fIdentifier = "h";
  y->h = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 6, "h")),
    &thisId);
  thisId.fIdentifier = "dt";
  y->dt = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 7,
    "dt")), &thisId);
  thisId.fIdentifier = "ct";
  y->ct = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 8,
    "ct")), &thisId);
  thisId.fIdentifier = "amax";
  y->amax = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 9,
    "amax")), &thisId);
  thisId.fIdentifier = "vmax";
  y->vmax = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 10,
    "vmax")), &thisId);
  thisId.fIdentifier = "minP";
  y->minP = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 11,
    "minP")), &thisId);
  thisId.fIdentifier = "maxP";
  y->maxP = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 12,
    "maxP")), &thisId);
  thisId.fIdentifier = "minV";
  y->minV = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 13,
    "minV")), &thisId);
  thisId.fIdentifier = "maxV";
  y->maxV = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 14,
    "maxV")), &thisId);
  thisId.fIdentifier = "Dmax";
  y->Dmax = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 15,
    "Dmax")), &thisId);
  thisId.fIdentifier = "dmin";
  y->dmin = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 16,
    "dmin")), &thisId);
  thisId.fIdentifier = "Ds";
  y->Ds = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 17,
    "Ds")), &thisId);
  thisId.fIdentifier = "Dc";
  y->Dc = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 18,
    "Dc")), &thisId);
  thisId.fIdentifier = "gamma";
  y->gamma = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 19,
    "gamma")), &thisId);
  thisId.fIdentifier = "rs";
  y->rs = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 20,
    "rs")), &thisId);
  thisId.fIdentifier = "ws";
  y->ws = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 21,
    "ws")), &thisId);
  thisId.fIdentifier = "wc";
  y->wc = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 22,
    "wc")), &thisId);
  thisId.fIdentifier = "wo";
  y->wo = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 23,
    "wo")), &thisId);
  thisId.fIdentifier = "wt";
  y->wt = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 24,
    "wt")), &thisId);
  thisId.fIdentifier = "rects";
  f_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 25, "rects")),
                     &thisId, y->rects);
  emlrtDestroyArray(&u);
}

static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = h_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  char_T *identifier))[24]
{
  real_T (*y)[24];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(u), &thisId);
  emlrtDestroyArray(&u);
  return y;
}
  static const mxArray *emlrt_marshallOut(const real_T u_data[], const int32_T
  u_size[1])
{
  const mxArray *y;
  const mxArray *m;
  static const int32_T iv[1] = { 0 };

  y = NULL;
  m = emlrtCreateNumericArray(1, iv, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, (void *)&u_data[0]);
  emlrtSetDimensions((mxArray *)m, u_size, 1);
  emlrtAssign(&y, m);
  return y;
}

static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
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
    y[i].cx = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, i, 0,
      "cx")), &thisId);
    thisId.fIdentifier = "cy";
    y[i].cy = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, i, 1,
      "cy")), &thisId);
    thisId.fIdentifier = "w";
    y[i].w = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, i, 2,
      "w")), &thisId);
    thisId.fIdentifier = "l";
    y[i].l = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, i, 3,
      "l")), &thisId);
    thisId.fIdentifier = "theta";
    y[i].theta = e_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, i,
      4, "theta")), &thisId);
  }

  emlrtDestroyArray(&u);
}

static real_T (*g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[24]
{
  real_T (*ret)[24];
  static const int32_T dims[1] = { 24 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 1U, dims);
  ret = (real_T (*)[24])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static real_T h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

void constraints_api(const mxArray * const prhs[2], int32_T nlhs, const mxArray *
                     plhs[2])
{
  real_T (*c_data)[576];
  real_T (*u)[24];
  struct0_T params;
  int32_T c_size[1];
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  c_data = (real_T (*)[576])mxMalloc(sizeof(real_T [576]));
  mxMalloc(0U);

  /* Marshall function inputs */
  u = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "u");
  c_emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "params", &params);

  /* Invoke the target function */
  constraints(&st, *u, &params, *c_data, c_size);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(*c_data, c_size);
  if (nlhs > 1) {
    plhs[1] = b_emlrt_marshallOut();
  }
}

/* End of code generation (_coder_constraints_api.c) */
