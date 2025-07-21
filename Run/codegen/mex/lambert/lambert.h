/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * lambert.h
 *
 * Code generation for function 'lambert'
 *
 */

#pragma once

/* Include files */
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void lambert(const emlrtStack *sp, real_T r1vec[3], real_T r2vec[3], real_T tf,
             real_T m, real_T muC, real_T V1[3], real_T V2[3],
             real_T extremal_distances[2], real_T *exitflag);

void sigmax_init(void);

/* End of code generation (lambert.h) */
