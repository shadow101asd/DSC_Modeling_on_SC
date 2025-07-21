/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_lambert_mex.c
 *
 * Code generation for function '_coder_lambert_mex'
 *
 */

/* Include files */
#include "_coder_lambert_mex.h"
#include "_coder_lambert_api.h"
#include "lambert_data.h"
#include "lambert_initialize.h"
#include "lambert_terminate.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void lambert_mexFunction(int32_T nlhs, mxArray *plhs[4], int32_T nrhs,
                         const mxArray *prhs[5])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *outputs[4];
  int32_T i;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 5) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 5, 4,
                        7, "lambert");
  }
  if (nlhs > 4) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 7,
                        "lambert");
  }
  /* Call the function. */
  lambert_api(prhs, nlhs, outputs);
  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    i = 1;
  } else {
    i = nlhs;
  }
  emlrtReturnArrays(i, &plhs[0], &outputs[0]);
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&lambert_atexit);
  /* Module initialization. */
  lambert_initialize();
  /* Dispatch the entry-point. */
  lambert_mexFunction(nlhs, plhs, nrhs, prhs);
  /* Module termination. */
  lambert_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL, "UTF-8", true);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_lambert_mex.c) */
