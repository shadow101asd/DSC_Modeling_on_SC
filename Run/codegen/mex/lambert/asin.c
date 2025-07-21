/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * asin.c
 *
 * Code generation for function 'asin'
 *
 */

/* Include files */
#include "asin.h"
#include "lambert_data.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Function Definitions */
void b_asin(const emlrtStack *sp, real_T x[2])
{
  int32_T k;
  boolean_T p;
  p = false;
  for (k = 0; k < 2; k++) {
    if (p) {
      p = true;
    } else {
      real_T d;
      d = x[k];
      if ((d < -1.0) || (d > 1.0)) {
        p = true;
      }
    }
  }
  if (p) {
    emlrtErrorWithMessageIdR2018a(
        sp, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "asin");
  }
  x[0] = muDoubleScalarAsin(x[0]);
  x[1] = muDoubleScalarAsin(x[1]);
}

/* End of code generation (asin.c) */
