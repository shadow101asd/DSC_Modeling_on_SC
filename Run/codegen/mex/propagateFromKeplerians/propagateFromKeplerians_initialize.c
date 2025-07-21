/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * propagateFromKeplerians_initialize.c
 *
 * Code generation for function 'propagateFromKeplerians_initialize'
 *
 */

/* Include files */
#include "propagateFromKeplerians_initialize.h"
#include "_coder_propagateFromKeplerians_mex.h"
#include "propagateFromKeplerians_data.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void propagateFromKeplerians_once(void);

/* Function Definitions */
static void propagateFromKeplerians_once(void)
{
  mex_InitInfAndNan();
}

void propagateFromKeplerians_initialize(void)
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
    propagateFromKeplerians_once();
  }
}

/* End of code generation (propagateFromKeplerians_initialize.c) */
