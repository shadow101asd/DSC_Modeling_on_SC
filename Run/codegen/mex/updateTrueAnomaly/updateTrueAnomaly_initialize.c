/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * updateTrueAnomaly_initialize.c
 *
 * Code generation for function 'updateTrueAnomaly_initialize'
 *
 */

/* Include files */
#include "updateTrueAnomaly_initialize.h"
#include "_coder_updateTrueAnomaly_mex.h"
#include "rt_nonfinite.h"
#include "updateTrueAnomaly_data.h"

/* Function Declarations */
static void updateTrueAnomaly_once(void);

/* Function Definitions */
static void updateTrueAnomaly_once(void)
{
  mex_InitInfAndNan();
}

void updateTrueAnomaly_initialize(void)
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
    updateTrueAnomaly_once();
  }
}

/* End of code generation (updateTrueAnomaly_initialize.c) */
