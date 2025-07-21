/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * lambert_initialize.c
 *
 * Code generation for function 'lambert_initialize'
 *
 */

/* Include files */
#include "lambert_initialize.h"
#include "_coder_lambert_mex.h"
#include "lambert.h"
#include "lambert_data.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void lambert_once(void);

/* Function Definitions */
static void lambert_once(void)
{
  mex_InitInfAndNan();
  sigmax_init();
}

void lambert_initialize(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2022b(&st);
  emlrtClearAllocCountR2012b(&st, false, 0U, NULL);
  emlrtEnterRtStackR2012b(&st);
  if (emlrtFirstTimeR2012b(emlrtRootTLSGlobal)) {
    lambert_once();
  }
}

/* End of code generation (lambert_initialize.c) */
