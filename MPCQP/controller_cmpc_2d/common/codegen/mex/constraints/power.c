/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * power.c
 *
 * Code generation for function 'power'
 *
 */

/* Include files */
#include "power.h"
#include "constraints.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void power(const real_T a_data[], const int32_T a_size[3], real_T y_data[],
           int32_T y_size[3])
{
  int8_T unnamed_idx_0;
  int8_T unnamed_idx_2;
  int32_T nx;
  int32_T k;
  unnamed_idx_0 = (int8_T)a_size[0];
  unnamed_idx_2 = (int8_T)a_size[2];
  y_size[0] = unnamed_idx_0;
  y_size[1] = 2;
  y_size[2] = unnamed_idx_2;
  nx = (unnamed_idx_0 << 1) * unnamed_idx_2;
  for (k = 0; k < nx; k++) {
    y_data[k] = a_data[k] * a_data[k];
  }
}

/* End of code generation (power.c) */
