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
#include "cost_sum_emxutil.h"
#include "dynamics_quadcopter.h"
#include "fitness.h"
#include "mwmathutil.h"
#include "rt_nonfinite.h"
#include "trim_vec.h"
#include "u2acc.h"
#include <string.h>

/* Type Definitions */
#ifndef typedef_struct_T
#define typedef_struct_T

typedef struct {
  real_T prev_e_R[3];
  real_T sum_e_R[3];
  real_T prev_e_w[3];
  real_T sum_e_w[3];
} struct_T;

#endif                                 /*typedef_struct_T*/

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 16,    /* lineNo */
  "cost_sum",                          /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/common/cost_sum.m"/* pathName */
};

static emlrtRSInfo b_emlrtRSI = { 34,  /* lineNo */
  "cost_sum",                          /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/common/cost_sum.m"/* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 41,  /* lineNo */
  "cost_sum",                          /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/common/cost_sum.m"/* pathName */
};

static emlrtRSInfo d_emlrtRSI = { 44,  /* lineNo */
  "cost_sum",                          /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/common/cost_sum.m"/* pathName */
};

static emlrtRSInfo e_emlrtRSI = { 47,  /* lineNo */
  "cost_sum",                          /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/common/cost_sum.m"/* pathName */
};

static emlrtRSInfo q_emlrtRSI = { 23,  /* lineNo */
  "acc2thrust",                        /* fcnName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/Quadcopter_model_new/acc2thrust.m"/* pathName */
};

static emlrtRTEInfo emlrtRTEI = { 30,  /* lineNo */
  13,                                  /* colNo */
  "cost_sum",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/common/cost_sum.m"/* pName */
};

static emlrtRTEInfo b_emlrtRTEI = { 31,/* lineNo */
  17,                                  /* colNo */
  "cost_sum",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/common/cost_sum.m"/* pName */
};

static emlrtECInfo emlrtECI = { -1,    /* nDims */
  18,                                  /* lineNo */
  1,                                   /* colNo */
  "cost_sum",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/common/cost_sum.m"/* pName */
};

static emlrtBCInfo emlrtBCI = { -1,    /* iFirst */
  -1,                                  /* iLast */
  32,                                  /* lineNo */
  41,                                  /* colNo */
  "a",                                 /* aName */
  "cost_sum",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/common/cost_sum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtRTEInfo c_emlrtRTEI = { 38,/* lineNo */
  13,                                  /* colNo */
  "cost_sum",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/common/cost_sum.m"/* pName */
};

static emlrtRTEInfo d_emlrtRTEI = { 40,/* lineNo */
  21,                                  /* colNo */
  "cost_sum",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/common/cost_sum.m"/* pName */
};

static emlrtBCInfo b_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  41,                                  /* lineNo */
  43,                                  /* colNo */
  "a",                                 /* aName */
  "cost_sum",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/common/cost_sum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo emlrtDCI = { 17,    /* lineNo */
  26,                                  /* colNo */
  "cost_sum",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/common/cost_sum.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo b_emlrtDCI = { 17,  /* lineNo */
  26,                                  /* colNo */
  "cost_sum",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/common/cost_sum.m",/* pName */
  4                                    /* checkKind */
};

static emlrtDCInfo c_emlrtDCI = { 17,  /* lineNo */
  1,                                   /* colNo */
  "cost_sum",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/common/cost_sum.m",/* pName */
  1                                    /* checkKind */
};

static emlrtRTEInfo n_emlrtRTEI = { 17,/* lineNo */
  1,                                   /* colNo */
  "cost_sum",                          /* fName */
  "/Users/Usama/UmehmoodGoogle/Work/Code/Swarms New/MPC_quadcopter/controller_dmpc_3d/common/cost_sum.m"/* pName */
};

/* Function Definitions */
real_T cost_sum(const emlrtStack *sp, const real_T u[18], real_T s[96], const
                struct0_T *params)
{
  real_T cost;
  emxArray_real_T *a;
  real_T a_i_data[54];
  int32_T a_i_size[2];
  int32_T i;
  real_T d;
  int32_T loop_ub;
  int32_T b_loop_ub;
  int32_T iv[3];
  real_T control_step;
  struct_T e[8];
  real_T u_dash_prev[32];
  int32_T h;
  int32_T k;
  int32_T i1;
  real_T acc_total;
  real_T z1[18];
  int32_T b_i;
  int32_T absxk_tmp;
  int32_T b_absxk_tmp;
  int32_T s_tmp;
  real_T scale;
  real_T absxk;
  real_T t;
  real_T thrust;
  real_T d1;
  real_T c;
  real_T b_s;
  real_T R_des[9];
  real_T b_c[9];
  real_T dv[9];
  static const int8_T iv1[3] = { 1, 0, 0 };

  real_T R[9];
  real_T u_quad[4];
  int32_T i2;
  real_T b_R_des[9];
  static const int8_T R_z_des[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

  int32_T i3;
  real_T unusedExpr[4];
  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInit_real_T(sp, &a, 3, &n_emlrtRTEI, true);

  /*  cost_sum - 3D DMPC */
  /*  Non-liear constraints passed to fmincon as function handle.  */
  /*  like this:  */
  /*      x = fmincon(@myfun,x0,A,b,Aeq,beq,lb,ub,@mycon) */
  /*  Input: */
  /*    - s      % 1x12xn - state of the flock  */
  /*    - u      % 3.h x 1 - The sequence of control action over the horizon  */
  /*    - params % parameters */
  /*  Output: */
  /*    - cost      % combined cost */
  /*  Usama Mehmood - Oct 2019 */
  /*  Add zero acc of neighbors */
  st.site = &emlrtRSI;
  u2acc(&st, u, params->h, a_i_data, a_i_size);
  i = a->size[0] * a->size[1] * a->size[2];
  a->size[0] = 3;
  a->size[1] = 8;
  emxEnsureCapacity_real_T(sp, a, i, &n_emlrtRTEI);
  if (!(params->h >= 0.0)) {
    emlrtNonNegativeCheckR2012b(params->h, &b_emlrtDCI, sp);
  }

  d = (int32_T)muDoubleScalarFloor(params->h);
  if (params->h != d) {
    emlrtIntegerCheckR2012b(params->h, &emlrtDCI, sp);
  }

  i = a->size[0] * a->size[1] * a->size[2];
  loop_ub = (int32_T)params->h;
  a->size[2] = loop_ub;
  emxEnsureCapacity_real_T(sp, a, i, &n_emlrtRTEI);
  if (params->h != d) {
    emlrtIntegerCheckR2012b(params->h, &c_emlrtDCI, sp);
  }

  b_loop_ub = 24 * loop_ub;
  for (i = 0; i < b_loop_ub; i++) {
    a->data[i] = 0.0;
  }

  iv[0] = 3;
  iv[1] = 1;
  iv[2] = loop_ub;
  emlrtSubAssignSizeCheckR2012b(&iv[0], 3, &a_i_size[0], 2, &emlrtECI, sp);
  for (i = 0; i < loop_ub; i++) {
    a->data[24 * i] = a_i_data[3 * i];
    a->data[24 * i + 1] = a_i_data[3 * i + 1];
    a->data[24 * i + 2] = a_i_data[3 * i + 2];
  }

  /*  */
  control_step = params->ct / params->dt;
  cost = 0.0;
  for (b_loop_ub = 0; b_loop_ub < 8; b_loop_ub++) {
    e[b_loop_ub].prev_e_R[0] = 0.0;
    e[b_loop_ub].prev_e_R[1] = 0.0;
    e[b_loop_ub].prev_e_R[2] = 0.0;
    e[b_loop_ub].sum_e_R[0] = 0.0;
    e[b_loop_ub].sum_e_R[1] = 0.0;
    e[b_loop_ub].sum_e_R[2] = 0.0;
    e[b_loop_ub].prev_e_w[0] = 0.0;
    e[b_loop_ub].prev_e_w[1] = 0.0;
    e[b_loop_ub].prev_e_w[2] = 0.0;
    e[b_loop_ub].sum_e_w[0] = 0.0;
    e[b_loop_ub].sum_e_w[1] = 0.0;
    e[b_loop_ub].sum_e_w[2] = 0.0;
  }

  if (params->cmpc_prediction_model == 1.0) {
    /* point-model */
    emlrtForLoopVectorCheckR2012b(1.0, 1.0, params->h, mxDOUBLE_CLASS, (int32_T)
      params->h, &emlrtRTEI, sp);
    for (h = 0; h < loop_ub; h++) {
      i = (int32_T)control_step;
      emlrtForLoopVectorCheckR2012b(1.0, 1.0, control_step, mxDOUBLE_CLASS,
        (int32_T)control_step, &b_emlrtRTEI, sp);
      for (k = 0; k < i; k++) {
        i1 = h + 1;
        if ((i1 < 1) || (i1 > a->size[2])) {
          emlrtDynamicBoundsCheckR2012b(i1, 1, a->size[2], &emlrtBCI, sp);
        }

        /*  dynamics_point  */
        /*  Input: */
        /*    - s      % 1x12xn - state of the flock  */
        /*    - a      % 3xn - The sequence of control action over the horizon  */
        /*    - params % parameters */
        /*  Output: */
        /*    - s_new  % next state */
        /*  Usama Mehmood - Oct 2019 */
        d = params->dt;
        for (b_i = 0; b_i < 8; b_i++) {
          /* update vel */
          b_loop_ub = 12 * b_i + 3;
          absxk_tmp = 3 * b_i + 24 * h;
          s[b_loop_ub] += d * a->data[absxk_tmp];
          b_absxk_tmp = 12 * b_i + 4;
          s[b_absxk_tmp] += d * a->data[absxk_tmp + 1];
          s_tmp = 12 * b_i + 5;
          s[s_tmp] += d * a->data[absxk_tmp + 2];
          trim_vec(*(real_T (*)[3])&s[12 * b_i + 3], params->vmax);

          /* update pos */
          s[12 * b_i] += d * s[b_loop_ub];
          b_loop_ub = 12 * b_i + 1;
          s[b_loop_ub] += d * s[b_absxk_tmp];
          b_loop_ub = 12 * b_i + 2;
          s[b_loop_ub] += d * s[s_tmp];
          if (*emlrtBreakCheckR2012bFlagVar != 0) {
            emlrtBreakCheckR2012b(sp);
          }
        }

        if (*emlrtBreakCheckR2012bFlagVar != 0) {
          emlrtBreakCheckR2012b(sp);
        }
      }

      st.site = &b_emlrtRSI;
      cost += fitness(&st, s, params->wc, params->ws, params->knn, params->turn,
                      params->w_m);
      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b(sp);
      }
    }
  } else {
    if (params->cmpc_prediction_model == 2.0) {
      /* quadcopter model */
      memset(&u_dash_prev[0], 0, 32U * sizeof(real_T));
      emlrtForLoopVectorCheckR2012b(1.0, 1.0, params->h, mxDOUBLE_CLASS,
        (int32_T)params->h, &c_emlrtRTEI, sp);
      if (0 <= loop_ub - 1) {
        i1 = (int32_T)control_step;
      }

      for (h = 0; h < loop_ub; h++) {
        for (b_i = 0; b_i < 8; b_i++) {
          emlrtForLoopVectorCheckR2012b(1.0, 1.0, control_step, mxDOUBLE_CLASS,
            (int32_T)control_step, &d_emlrtRTEI, sp);
          for (k = 0; k < i1; k++) {
            st.site = &c_emlrtRSI;
            i = h + 1;
            if ((i < 1) || (i > a->size[2])) {
              emlrtDynamicBoundsCheckR2012b(i, 1, a->size[2], &b_emlrtBCI, &st);
            }

            /*  Usama Mehmood - Oct 2019 */
            /*  roll = 0; */
            /*  if all( a(1:2) == 0 ) */
            /*      pitch = 0; */
            /*  else */
            /*      pitch = atand(a(2)/ sqrt((a(1)^2) + (a(2)^2))); */
            /*  end */
            /*  yaw = 0; */
            /*  yaw = atan360(a(2), a(1)); */
            /*  if norm(a) ~= 0  */
            /*      thrust = norm(a) * params.m;  */
            /*      R = vec2rot(a); */
            /*      angles = rotm2eul(R); %Order is ZYX */
            /*  end */
            /*  Small angle approximation */
            scale = 3.3121686421112381E-170;
            b_loop_ub = 3 * b_i + 24 * h;
            absxk = muDoubleScalarAbs(a->data[b_loop_ub]);
            if (absxk > 3.3121686421112381E-170) {
              acc_total = 1.0;
              scale = absxk;
            } else {
              t = absxk / 3.3121686421112381E-170;
              acc_total = t * t;
            }

            absxk_tmp = b_loop_ub + 1;
            absxk = muDoubleScalarAbs(a->data[absxk_tmp]);
            if (absxk > scale) {
              t = scale / absxk;
              acc_total = acc_total * t * t + 1.0;
              scale = absxk;
            } else {
              t = absxk / scale;
              acc_total += t * t;
            }

            b_absxk_tmp = b_loop_ub + 2;
            absxk = muDoubleScalarAbs(a->data[b_absxk_tmp]);
            if (absxk > scale) {
              t = scale / absxk;
              acc_total = acc_total * t * t + 1.0;
              scale = absxk;
            } else {
              t = absxk / scale;
              acc_total += t * t;
            }

            acc_total = scale * muDoubleScalarSqrt(acc_total);
            if (acc_total != 0.0) {
              /*  net acceleration along body-frame z-axis. */
              scale = a->data[b_absxk_tmp] + params->g;
              b_st.site = &q_emlrtRSI;
              acc_total = (a->data[b_loop_ub] * a->data[b_loop_ub] + a->
                           data[absxk_tmp] * a->data[absxk_tmp]) + scale * scale;
              if (acc_total < 0.0) {
                emlrtErrorWithMessageIdR2018a(&b_st, &e_emlrtRTEI,
                  "Coder:toolbox:ElFunDomainError",
                  "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
              }

              acc_total = muDoubleScalarSqrt(acc_total);
              scale = 1.0 / acc_total;
              absxk = a->data[b_loop_ub];
              t = a->data[absxk_tmp];
              d1 = scale * 0.0 * absxk + -scale * t;
              d = scale * absxk + scale * 0.0 * t;
              thrust = acc_total * params->m;
            } else {
              thrust = 0.0;
              i = b_i << 2;
              d = u_dash_prev[i + 2];
              d1 = u_dash_prev[i + 1];
            }

            /*  exact calculation */
            /*  if norm(a) ~= 0  */
            /*  psi = 0; */
            /*      % net acceleration along body-frame z-axis. */
            /*      acc_total = sqrt(a(1)^2 + a(2)^2 + (a(3)+params.g)^2); */
            /*      M = [sin(psi), -cos(psi);... */
            /*          cos(psi), sin(psi)]; */
            /*      trigs = 1/acc_total * M * [a(1); a(2)]; */
            /*      phi = asin(trigs(1)); */
            /*      theta = asin(trigs(2)/cos(phi)); */
            /*      angles(3) = phi; */
            /*      angles(2) = theta; */
            /*      angles(1) = psi; */
            /*      thrust = acc_total * params.m; */
            /*  end */
            /*  */
            b_loop_ub = b_i << 2;
            u_dash_prev[b_loop_ub] = thrust;
            u_dash_prev[(b_i << 2) + 1] = d1;
            u_dash_prev[(b_i << 2) + 2] = d;
            u_dash_prev[b_loop_ub + 3] = 0.0;

            /*  Usama Mehmood - Oct 2019 */
            /* 25]; */
            /* 250]; */
            /* 25]; */
            /* 250]; */
            /* roll */
            /* pitch */
            /* yaw */
            /*  Computing the error in orientation, e_R        */
            /*  R_dim - angles are in radians */
            /*   */
            b_loop_ub = 12 * b_i + 6;
            scale = muDoubleScalarCos(s[b_loop_ub]);
            absxk = muDoubleScalarSin(s[b_loop_ub]);

            /*  R_dim - angles are in radians */
            /*   */
            b_loop_ub = 12 * b_i + 7;
            t = muDoubleScalarCos(s[b_loop_ub]);
            acc_total = muDoubleScalarSin(s[b_loop_ub]);

            /*  R_dim - angles are in radians */
            /*   */
            b_loop_ub = 12 * b_i + 8;
            c = muDoubleScalarCos(s[b_loop_ub]);
            b_s = muDoubleScalarSin(s[b_loop_ub]);
            R_des[1] = 0.0;
            R_des[4] = scale;
            R_des[7] = absxk;
            R_des[2] = 0.0;
            R_des[5] = -absxk;
            R_des[8] = scale;
            b_c[0] = t;
            b_c[3] = 0.0;
            b_c[6] = -acc_total;
            R_des[0] = 1.0;
            b_c[1] = 0.0;
            R_des[3] = 0.0;
            b_c[4] = 1.0;
            R_des[6] = 0.0;
            b_c[7] = 0.0;
            b_c[2] = acc_total;
            b_c[5] = 0.0;
            b_c[8] = t;
            for (i = 0; i < 3; i++) {
              scale = R_des[i + 3];
              absxk = R_des[i + 6];
              for (absxk_tmp = 0; absxk_tmp < 3; absxk_tmp++) {
                dv[i + 3 * absxk_tmp] = ((real_T)(int32_T)R_des[i] * b_c[3 *
                  absxk_tmp] + scale * b_c[3 * absxk_tmp + 1]) + absxk * b_c[3 *
                  absxk_tmp + 2];
              }
            }

            b_c[0] = c;
            b_c[3] = b_s;
            b_c[6] = 0.0;
            b_c[1] = -b_s;
            b_c[4] = c;
            b_c[7] = 0.0;
            b_c[2] = 0.0;
            b_c[5] = 0.0;
            b_c[8] = 1.0;

            /*  R_dim - angles are in radians */
            /*   */
            scale = muDoubleScalarCos(d1);
            absxk = muDoubleScalarSin(d1);

            /*  R_dim - angles are in radians */
            /*   */
            t = muDoubleScalarCos(d);
            acc_total = muDoubleScalarSin(d);

            /*  R_dim - angles are in radians */
            /*   */
            for (i = 0; i < 3; i++) {
              d = dv[i + 3];
              d1 = dv[i + 6];
              for (absxk_tmp = 0; absxk_tmp < 3; absxk_tmp++) {
                R[i + 3 * absxk_tmp] = (dv[i] * b_c[3 * absxk_tmp] + d * b_c[3 *
                  absxk_tmp + 1]) + d1 * b_c[3 * absxk_tmp + 2];
              }

              R_des[3 * i] = iv1[i];
            }

            R_des[1] = 0.0;
            R_des[4] = scale;
            R_des[7] = absxk;
            R_des[2] = 0.0;
            R_des[5] = -absxk;
            R_des[8] = scale;
            b_c[0] = t;
            b_c[3] = 0.0;
            b_c[6] = -acc_total;
            b_c[1] = 0.0;
            b_c[4] = 1.0;
            b_c[7] = 0.0;
            b_c[2] = acc_total;
            b_c[5] = 0.0;
            b_c[8] = t;
            for (i = 0; i < 3; i++) {
              d = R_des[i + 3];
              d1 = R_des[i + 6];
              for (absxk_tmp = 0; absxk_tmp < 3; absxk_tmp++) {
                dv[i + 3 * absxk_tmp] = ((real_T)(int32_T)R_des[i] * b_c[3 *
                  absxk_tmp] + d * b_c[3 * absxk_tmp + 1]) + d1 * b_c[3 *
                  absxk_tmp + 2];
              }

              d = dv[i + 3];
              d1 = dv[i + 6];
              for (absxk_tmp = 0; absxk_tmp < 3; absxk_tmp++) {
                b_R_des[i + 3 * absxk_tmp] = (dv[i] * (real_T)R_z_des[3 *
                  absxk_tmp] + d * (real_T)R_z_des[3 * absxk_tmp + 1]) + d1 *
                  (real_T)R_z_des[3 * absxk_tmp + 2];
              }
            }

            for (i = 0; i < 3; i++) {
              absxk_tmp = 3 * i + 1;
              b_absxk_tmp = 3 * i + 2;
              for (s_tmp = 0; s_tmp < 3; s_tmp++) {
                b_loop_ub = i + 3 * s_tmp;
                b_c[b_loop_ub] = 0.0;
                i2 = 3 * s_tmp + 1;
                i3 = 3 * s_tmp + 2;
                b_c[b_loop_ub] = (R[3 * i] * b_R_des[3 * s_tmp] + R[absxk_tmp] *
                                  b_R_des[i2]) + R[b_absxk_tmp] * b_R_des[i3];
                R_des[b_loop_ub] = (b_R_des[3 * i] * R[3 * s_tmp] +
                                    b_R_des[absxk_tmp] * R[i2]) +
                  b_R_des[b_absxk_tmp] * R[i3];
              }
            }

            for (i = 0; i < 9; i++) {
              R_des[i] = 0.5 * (R_des[i] - b_c[i]);
            }

            /*  e_R = error_angle_difference(phi, phi_des, theta, theta_des, psi, psi_des); */
            /*  Computing error in the angular velocity e_w */
            /*  PID */
            u_quad[0] = thrust;
            d = 0.0 - s[12 * b_i + 9];
            d1 = e[b_i].sum_e_R[0] + R_des[5];
            u_quad[1] = ((((0.24434609527920614 * R_des[5] + 1.7453292519943295 *
                            (R_des[5] - e[b_i].prev_e_R[0])) + 0.0 * d1) + -0.0 *
                          d) + -0.0 * (d - e[b_i].prev_e_w[0])) + 0.0 * (e[b_i].
              sum_e_w[0] + d);
            e[b_i].sum_e_R[0] = d1;
            e[b_i].prev_e_R[0] = R_des[5];
            e[b_i].sum_e_w[0] += R_des[5];
            e[b_i].prev_e_w[0] = R_des[5];
            d = 0.0 - s[12 * b_i + 10];
            d1 = e[b_i].sum_e_R[1] + R_des[6];
            u_quad[2] = ((((0.24434609527920614 * R_des[6] + 1.7453292519943295 *
                            (R_des[6] - e[b_i].prev_e_R[1])) + 0.0 * d1) + -0.0 *
                          d) + -0.0 * (d - e[b_i].prev_e_w[1])) + 0.0 * (e[b_i].
              sum_e_w[1] + d);
            e[b_i].sum_e_R[1] = d1;
            e[b_i].prev_e_R[1] = R_des[6];
            e[b_i].sum_e_w[1] += R_des[6];
            e[b_i].prev_e_w[1] = R_des[6];
            d = 0.0 - s[12 * b_i + 11];
            d1 = e[b_i].sum_e_R[2] + R_des[1];
            u_quad[3] = ((((0.0024434609527920616 * R_des[1] +
                            0.017453292519943295 * (R_des[1] - e[b_i].prev_e_R[2]))
                           + 0.0 * d1) + -0.0 * d) + -0.0 * (d - e[b_i]
              .prev_e_w[2])) + 0.0 * (e[b_i].sum_e_w[2] + d);
            e[b_i].sum_e_R[2] = d1;
            e[b_i].prev_e_R[2] = R_des[1];
            e[b_i].sum_e_w[2] += R_des[1];
            e[b_i].prev_e_w[2] = R_des[1];
            st.site = &d_emlrtRSI;
            dynamics_quadcopter(&st, *(real_T (*)[12])&s[12 * b_i], u_quad,
                                params, unusedExpr);
            if (*emlrtBreakCheckR2012bFlagVar != 0) {
              emlrtBreakCheckR2012b(sp);
            }
          }

          if (*emlrtBreakCheckR2012bFlagVar != 0) {
            emlrtBreakCheckR2012b(sp);
          }
        }

        st.site = &e_emlrtRSI;
        cost += fitness(&st, s, params->wc, params->ws, params->knn,
                        params->turn, params->w_m);
        if (*emlrtBreakCheckR2012bFlagVar != 0) {
          emlrtBreakCheckR2012b(sp);
        }
      }
    }
  }

  emxFree_real_T(&a);
  for (k = 0; k < 18; k++) {
    z1[k] = u[k] * u[k];
  }

  acc_total = z1[0];
  for (k = 0; k < 17; k++) {
    acc_total += z1[k + 1];
  }

  cost = cost / params->h + 0.5 * (acc_total / params->h);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
  return cost;
}

/* End of code generation (cost_sum.c) */
