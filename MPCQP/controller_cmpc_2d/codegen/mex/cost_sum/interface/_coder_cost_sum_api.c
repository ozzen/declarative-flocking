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

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[180];
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *pos,
  const char_T *identifier))[30];
static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[30];
static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *params,
  const char_T *identifier, struct0_T *y);
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  char_T *identifier))[180];
static const mxArray *emlrt_marshallOut(const real_T u);
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, struct0_T *y);
static real_T g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[180];
static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[30];
static real_T j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);

/* Function Definitions */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[180]
{
  real_T (*y)[180];
  y = h_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
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
  y = i_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
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
  char_T *identifier))[180]
{
  real_T (*y)[180];
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
  static const char * fieldNames[17] = { "steps", "n", "h", "dt", "ct", "amax",
    "vmax", "minP", "maxP", "minV", "maxV", "Dmax", "dmin", "Ds", "Dc", "gamma",
    "rs" };

  static const int32_T dims = 0;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b(sp, parentId, u, 17, fieldNames, 0U, &dims);
  thisId.fIdentifier = "steps";
  y->steps = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 0,
    "steps")), &thisId);
  thisId.fIdentifier = "n";
  y->n = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 1, "n")),
    &thisId);
  thisId.fIdentifier = "h";
  y->h = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 2, "h")),
    &thisId);
  thisId.fIdentifier = "dt";
  y->dt = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 3,
    "dt")), &thisId);
  thisId.fIdentifier = "ct";
  y->ct = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 4,
    "ct")), &thisId);
  thisId.fIdentifier = "amax";
  y->amax = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 5,
    "amax")), &thisId);
  thisId.fIdentifier = "vmax";
  y->vmax = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 6,
    "vmax")), &thisId);
  thisId.fIdentifier = "minP";
  y->minP = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 7,
    "minP")), &thisId);
  thisId.fIdentifier = "maxP";
  y->maxP = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 8,
    "maxP")), &thisId);
  thisId.fIdentifier = "minV";
  y->minV = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 9,
    "minV")), &thisId);
  thisId.fIdentifier = "maxV";
  y->maxV = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 10,
    "maxV")), &thisId);
  thisId.fIdentifier = "Dmax";
  y->Dmax = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 11,
    "Dmax")), &thisId);
  thisId.fIdentifier = "dmin";
  y->dmin = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 12,
    "dmin")), &thisId);
  thisId.fIdentifier = "Ds";
  y->Ds = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 13,
    "Ds")), &thisId);
  thisId.fIdentifier = "Dc";
  y->Dc = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 14,
    "Dc")), &thisId);
  thisId.fIdentifier = "gamma";
  y->gamma = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 15,
    "gamma")), &thisId);
  thisId.fIdentifier = "rs";
  y->rs = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b(sp, u, 0, 16,
    "rs")), &thisId);
  emlrtDestroyArray(&u);
}

static real_T g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = j_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T (*h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[180]
{
  real_T (*ret)[180];
  static const int32_T dims[1] = { 180 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 1U, dims);
  ret = (real_T (*)[180])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static real_T (*i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[30]
{
  real_T (*ret)[30];
  static const int32_T dims[2] = { 2, 15 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[30])emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

void cost_sum_api(const mxArray * const prhs[4], int32_T nlhs, const mxArray
                  *plhs[1])
{
  const mxArray *prhs_copy_idx_1;
  const mxArray *prhs_copy_idx_2;
  real_T (*u)[180];
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
  prhs_copy_idx_1 = emlrtProtectR2012b(prhs[1], 1, false, -1);
  prhs_copy_idx_2 = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  u = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "u");
  pos = c_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_1), "pos");
  vel = c_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_2), "vel");
  e_emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "params", &params);

  /* Invoke the target function */
  ret = cost_sum(&st, *u, *pos, *vel, &params);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(ret);
}

/* End of code generation (_coder_cost_sum_api.c) */
