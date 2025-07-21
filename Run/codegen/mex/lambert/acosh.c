/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * acosh.c
 *
 * Code generation for function 'acosh'
 *
 */

/* Include files */
#include "acosh.h"
#include "lambert_data.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo nb_emlrtRSI = {
    16,      /* lineNo */
    "acosh", /* fcnName */
    "/Applications/MATLAB_R2024a.app/toolbox/eml/lib/matlab/elfun/acosh.m" /* pathName
                                                                            */
};

static emlrtRSInfo ob_emlrtRSI = {
    17,                           /* lineNo */
    "applyScalarFunctionInPlace", /* fcnName */
    "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/"
    "applyScalarFunctionInPlace.m" /* pathName */
};

static emlrtRSInfo pb_emlrtRSI = {
    11,      /* lineNo */
    "acosh", /* fcnName */
    "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/+scalar/"
    "acosh.m" /* pathName */
};

static emlrtRSInfo qb_emlrtRSI = {
    39,          /* lineNo */
    "eml_acosh", /* fcnName */
    "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/+scalar/"
    "acosh.m" /* pathName */
};

/* Function Definitions */
void b_acosh(const emlrtStack *sp, real_T *x)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &nb_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  b_st.site = &ob_emlrtRSI;
  c_st.site = &pb_emlrtRSI;
  if (*x >= 2.68435456E+8) {
    *x = muDoubleScalarLog(*x) + 0.69314718055994529;
  } else if (*x > 2.0) {
    *x = muDoubleScalarLog(*x + muDoubleScalarSqrt(*x * *x - 1.0));
  } else {
    real_T absz;
    (*x)--;
    d_st.site = &qb_emlrtRSI;
    absz = 2.0 * *x + *x * *x;
    if (absz < 0.0) {
      emlrtErrorWithMessageIdR2018a(
          &d_st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
          "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
    }
    absz = muDoubleScalarSqrt(absz);
    *x += absz;
    absz = muDoubleScalarAbs(*x);
    if ((absz > 4.503599627370496E+15) ||
        (muDoubleScalarIsInf(*x) || muDoubleScalarIsNaN(*x))) {
      (*x)++;
      *x = muDoubleScalarLog(*x);
    } else if (!(absz < 2.2204460492503131E-16)) {
      *x = muDoubleScalarLog(*x + 1.0) * (*x / ((*x + 1.0) - 1.0));
    }
  }
}

/* End of code generation (acosh.c) */
