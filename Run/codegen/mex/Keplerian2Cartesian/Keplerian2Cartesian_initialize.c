/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Keplerian2Cartesian_initialize.c
 *
 * Code generation for function 'Keplerian2Cartesian_initialize'
 *
 */

/* Include files */
#include "Keplerian2Cartesian_initialize.h"
#include "Keplerian2Cartesian_data.h"
#include "_coder_Keplerian2Cartesian_mex.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void Keplerian2Cartesian_once(void);

/* Function Definitions */
static void Keplerian2Cartesian_once(void)
{
  mex_InitInfAndNan();
}

void Keplerian2Cartesian_initialize(void)
{
  static const volatile char_T *emlrtBreakCheckR2012bFlagVar = NULL;
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
    Keplerian2Cartesian_once();
  }
}

/* End of code generation (Keplerian2Cartesian_initialize.c) */
