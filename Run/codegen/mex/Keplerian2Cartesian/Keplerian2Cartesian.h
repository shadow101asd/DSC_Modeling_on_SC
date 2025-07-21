/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Keplerian2Cartesian.h
 *
 * Code generation for function 'Keplerian2Cartesian'
 *
 */

#pragma once

/* Include files */
#include "Keplerian2Cartesian_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void Keplerian2Cartesian(const emlrtStack *sp, const emxArray_real_T *a,
                         const emxArray_real_T *e, const emxArray_real_T *i,
                         emxArray_real_T *Om, emxArray_real_T *w,
                         const emxArray_real_T *f0, real_T mu,
                         emxArray_real_T *y);

/* End of code generation (Keplerian2Cartesian.h) */
