/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * propagateFromKeplerians.h
 *
 * Code generation for function 'propagateFromKeplerians'
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
void propagateFromKeplerians(const emlrtStack *sp, const real_T Ki[6],
                             real_T mu, const real_T etR_data[],
                             const int32_T etR_size[2], real_T Xs_data[],
                             int32_T Xs_size[2]);

/* End of code generation (propagateFromKeplerians.h) */
