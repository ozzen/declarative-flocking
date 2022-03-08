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
#include "cost_sum_emxutil.h"
#include "eml_int_forloop_overflow_check.h"
#include "mwmathutil.h"
#include "rt_nonfinite.h"
#include "separation.h"
#include "sum_sq_distances.h"

/* Variable Definitions */
static emlrtRSInfo h_emlrtRSI = { 3,   /* lineNo */
  "fitness",                           /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/fitness.m"/* pathName */
};

static emlrtRSInfo i_emlrtRSI = { 4,   /* lineNo */
  "fitness",                           /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/fitness.m"/* pathName */
};

static emlrtRSInfo j_emlrtRSI = { 5,   /* lineNo */
  "fitness",                           /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/fitness.m"/* pathName */
};

static emlrtRSInfo k_emlrtRSI = { 6,   /* lineNo */
  "fitness",                           /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/fitness.m"/* pathName */
};

static emlrtRSInfo l_emlrtRSI = { 8,   /* lineNo */
  "fitness",                           /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/fitness.m"/* pathName */
};

static emlrtRSInfo m_emlrtRSI = { 3,   /* lineNo */
  "vm",                                /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/cost_functions/vm.m"/* pathName */
};

static emlrtRSInfo n_emlrtRSI = { 28,  /* lineNo */
  "repmat",                            /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/elmat/repmat.m"/* pathName */
};

static emlrtRSInfo o_emlrtRSI = { 64,  /* lineNo */
  "repmat",                            /* fcnName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/elmat/repmat.m"/* pathName */
};

static emlrtRTEInfo k_emlrtRTEI = { 58,/* lineNo */
  23,                                  /* colNo */
  "assertValidSizeArg",                /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/assertValidSizeArg.m"/* pName */
};

static emlrtDCInfo d_emlrtDCI = { 31,  /* lineNo */
  14,                                  /* colNo */
  "repmat",                            /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/lib/matlab/elmat/repmat.m",/* pName */
  4                                    /* checkKind */
};

static emlrtECInfo b_emlrtECI = { 2,   /* nDims */
  6,                                   /* lineNo */
  11,                                  /* colNo */
  "vm",                                /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/cost_functions/vm.m"/* pName */
};

static emlrtBCInfo c_emlrtBCI = { 1,   /* iFirst */
  20,                                  /* iLast */
  6,                                   /* lineNo */
  24,                                  /* colNo */
  "params.w_m",                        /* aName */
  "vm",                                /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/cost_functions/vm.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo e_emlrtDCI = { 6,   /* lineNo */
  24,                                  /* colNo */
  "vm",                                /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/cost_functions/vm.m",/* pName */
  1                                    /* checkKind */
};

static emlrtECInfo c_emlrtECI = { 2,   /* nDims */
  3,                                   /* lineNo */
  23,                                  /* colNo */
  "vm",                                /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/cost_functions/vm.m"/* pName */
};

static emlrtRTEInfo l_emlrtRTEI = { 13,/* lineNo */
  13,                                  /* colNo */
  "toLogicalCheck",                    /* fName */
  "/Applications/MATLAB_R2019b.app/toolbox/eml/eml/+coder/+internal/toLogicalCheck.m"/* pName */
};

static emlrtRTEInfo o_emlrtRTEI = { 3, /* lineNo */
  23,                                  /* colNo */
  "vm",                                /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/cost_functions/vm.m"/* pName */
};

static emlrtRTEInfo p_emlrtRTEI = { 1, /* lineNo */
  18,                                  /* colNo */
  "fitness",                           /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/fitness.m"/* pName */
};

/* Function Definitions */
real_T fitness(const emlrtStack *sp, const real_T s[96], real_T params_wc,
               real_T params_ws, real_T params_knn, real_T params_turn, const
               real_T params_w_m[20])
{
  real_T res;
  real_T d;
  real_T b_params_knn;
  emxArray_real_T *r;
  int32_T jtilecol;
  int32_T loop_ub;
  int32_T ibtile;
  int32_T iv[2];
  real_T a[21];
  int32_T a_tmp;
  real_T z1[21];
  boolean_T p;
  real_T rel_speed[7];
  int32_T iv1[2];
  real_T params_w_m_data[20];
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = sp;
  e_st.tls = sp->tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);

  /*  Usama Mehmood - Oct 2019 */
  /*  fitness 3D DMPC Usama Mehmood - Oct 2019 */
  st.site = &h_emlrtRSI;
  if (muDoubleScalarIsNaN(params_turn)) {
    emlrtErrorWithMessageIdR2018a(&st, &l_emlrtRTEI, "MATLAB:nologicalnan",
      "MATLAB:nologicalnan", 0);
  }

  if (params_turn != 0.0) {
    st.site = &k_emlrtRSI;
    b_st.site = &m_emlrtRSI;
    c_st.site = &n_emlrtRSI;
    d = muDoubleScalarFloor(params_knn);
    if ((params_knn != d) || muDoubleScalarIsInf(params_knn) || (params_knn <
         -2.147483648E+9) || (params_knn > 2.147483647E+9)) {
      emlrtErrorWithMessageIdR2018a(&c_st, &k_emlrtRTEI,
        "Coder:MATLAB:NonIntegerInput", "Coder:MATLAB:NonIntegerInput", 4, 12,
        MIN_int32_T, 12, MAX_int32_T);
    }

    if (params_knn <= 0.0) {
      b_params_knn = 0.0;
    } else {
      b_params_knn = params_knn;
    }

    if (!(b_params_knn <= 2.147483647E+9)) {
      emlrtErrorWithMessageIdR2018a(&c_st, &j_emlrtRTEI, "Coder:MATLAB:pmaxsize",
        "Coder:MATLAB:pmaxsize", 0);
    }

    if (!(params_knn >= 0.0)) {
      emlrtNonNegativeCheckR2012b(params_knn, &d_emlrtDCI, &b_st);
    }

    emxInit_real_T(&b_st, &r, 2, &p_emlrtRTEI, true);
    jtilecol = r->size[0] * r->size[1];
    r->size[0] = 3;
    loop_ub = (int32_T)params_knn;
    r->size[1] = loop_ub;
    emxEnsureCapacity_real_T(&b_st, r, jtilecol, &o_emlrtRTEI);
    c_st.site = &o_emlrtRSI;
    if ((1 <= loop_ub) && (loop_ub > 2147483646)) {
      d_st.site = &p_emlrtRSI;
      check_forloop_overflow_error(&d_st);
    }

    for (jtilecol = 0; jtilecol < loop_ub; jtilecol++) {
      ibtile = jtilecol * 3;
      r->data[ibtile] = s[3];
      r->data[ibtile + 1] = s[4];
      r->data[ibtile + 2] = s[5];
    }

    iv[0] = 3;
    iv[1] = 7;
    emlrtSizeEqCheckNDR2012b(*(int32_T (*)[2])r->size, iv, &c_emlrtECI, &st);
    for (jtilecol = 0; jtilecol < 7; jtilecol++) {
      ibtile = 12 * (jtilecol + 1);
      a[3 * jtilecol] = r->data[3 * jtilecol] - s[ibtile + 3];
      a_tmp = 3 * jtilecol + 1;
      a[a_tmp] = r->data[a_tmp] - s[ibtile + 4];
      a_tmp = 3 * jtilecol + 2;
      a[a_tmp] = r->data[a_tmp] - s[ibtile + 5];
    }

    emxFree_real_T(&r);
    for (jtilecol = 0; jtilecol < 21; jtilecol++) {
      z1[jtilecol] = a[jtilecol] * a[jtilecol];
    }

    for (jtilecol = 0; jtilecol < 7; jtilecol++) {
      ibtile = jtilecol * 3;
      rel_speed[jtilecol] = (z1[ibtile] + z1[ibtile + 1]) + z1[ibtile + 2];
    }

    b_st.site = &m_emlrtRSI;
    p = false;
    for (jtilecol = 0; jtilecol < 7; jtilecol++) {
      if (p || (rel_speed[jtilecol] < 0.0)) {
        p = true;
      }
    }

    if (p) {
      emlrtErrorWithMessageIdR2018a(&b_st, &e_emlrtRTEI,
        "Coder:toolbox:ElFunDomainError", "Coder:toolbox:ElFunDomainError", 3, 4,
        4, "sqrt");
    }

    for (jtilecol = 0; jtilecol < 7; jtilecol++) {
      rel_speed[jtilecol] = muDoubleScalarSqrt(rel_speed[jtilecol]);
    }

    /*  AWNing */
    if (1.0 > params_knn) {
      loop_ub = 0;
    } else {
      if (params_knn != (int32_T)d) {
        emlrtIntegerCheckR2012b(params_knn, &e_emlrtDCI, &st);
      }

      if ((loop_ub < 1) || (loop_ub > 20)) {
        emlrtDynamicBoundsCheckR2012b(loop_ub, 1, 20, &c_emlrtBCI, &st);
      }
    }

    iv[0] = 1;
    iv[1] = loop_ub;
    iv1[0] = 1;
    iv1[1] = 7;
    if (loop_ub != 7) {
      emlrtSizeEqCheckNDR2012b(&iv[0], &iv1[0], &b_emlrtECI, &st);
    }

    for (jtilecol = 0; jtilecol < 7; jtilecol++) {
      params_w_m_data[jtilecol] = params_w_m[jtilecol];
    }

    for (jtilecol = 0; jtilecol < 7; jtilecol++) {
      rel_speed[jtilecol] *= params_w_m_data[jtilecol];
    }

    res = rel_speed[0];
    for (jtilecol = 0; jtilecol < 6; jtilecol++) {
      res += rel_speed[jtilecol + 1];
    }

    st.site = &i_emlrtRSI;
    e_st.site = &j_emlrtRSI;
    res += params_wc * sum_sq_distances(&st, s) + params_ws * separation(&e_st,
      s);
  } else {
    st.site = &l_emlrtRSI;
    res = params_wc * sum_sq_distances(&st, s) + params_ws * separation(&st, s);
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
  return res;
}

/* End of code generation (fitness.c) */
