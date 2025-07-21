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
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void Keplerian2Cartesian(const emlrtStack *sp, real_T a, real_T e, real_T i,
                         real_T Om, real_T w, real_T f0, real_T mu,
                         real_T y[6]);

/* End of code generation (Keplerian2Cartesian.h) */
