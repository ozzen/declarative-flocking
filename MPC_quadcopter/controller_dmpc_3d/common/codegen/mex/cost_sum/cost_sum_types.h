/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * cost_sum_types.h
 *
 * Code generation for function 'cost_sum_types'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  real_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_real_T*/

#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T

typedef struct emxArray_real_T emxArray_real_T;

#endif                                 /*typedef_emxArray_real_T*/

#ifndef typedef_struct0_T
#define typedef_struct0_T

typedef struct {
  real_T n;
  real_T h;
  real_T dt;
  real_T ct;
  real_T vmax;
  real_T amax;
  real_T ipos[2];
  real_T ivel[2];
  real_T delta_angle;
  real_T t_end;
  real_T steps;
  real_T m;
  real_T L;
  real_T dmin;
  real_T g;
  real_T j_r;
  real_T Ixx;
  real_T Iyy;
  real_T Izz;
  real_T max_rotor_speed;
  real_T max_angular_speed;
  real_T max_angular_acc;
  real_T k_f;
  real_T k_m;
  real_T wc;
  real_T ws;
  real_T quad;
  real_T cmpc_prediction_model;
  real_T knn;
  real_T turn;
  real_T eta;
  real_T w_m[20];
  real_T start_turn;
  real_T num_leaders;
  real_T t_fix;
  real_T turn_angle;
  real_T turn_acc;
} struct0_T;

#endif                                 /*typedef_struct0_T*/

/* End of code generation (cost_sum_types.h) */
