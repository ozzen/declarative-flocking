/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * fitness.h
 *
 * Code generation for function 'fitness'
 *
 */

#pragma once

/* Include files */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mex.h"
#include "emlrt.h"
#include "rtwtypes.h"
#include "cost_sum_types.h"

/* Function Declarations */
real_T fitness(const emlrtStack *sp, const real_T pos[30], const real_T vel[30],
               real_T params_n);

/* End of code generation (fitness.h) */
