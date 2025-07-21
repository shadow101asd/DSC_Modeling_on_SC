/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_propagateFromKeplerians_mex.c
 *
 * Code generation for function '_coder_propagateFromKeplerians_mex'
 *
 */

/* Include files */
#include "_coder_propagateFromKeplerians_mex.h"
#include "_coder_propagateFromKeplerians_api.h"
#include "propagateFromKeplerians_data.h"
#include "propagateFromKeplerians_initialize.h"
#include "propagateFromKeplerians_terminate.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&propagateFromKeplerians_atexit);
  /* Module initialization. */
  propagateFromKeplerians_initialize();
  /* Dispatch the entry-point. */
  propagateFromKeplerians_mexFunction(nlhs, plhs, nrhs, prhs);
  /* Module termination. */
  propagateFromKeplerians_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL, "UTF-8", true);
  return emlrtRootTLSGlobal;
}

void propagateFromKeplerians_mexFunction(int32_T nlhs, mxArray *plhs[1],
                                         int32_T nrhs, const mxArray *prhs[3])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *outputs;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 3, 4,
                        23, "propagateFromKeplerians");
  }
  if (nlhs > 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 23,
                        "propagateFromKeplerians");
  }
  /* Call the function. */
  propagateFromKeplerians_api(prhs, &outputs);
  /* Copy over outputs to the caller. */
  emlrtReturnArrays(1, &plhs[0], &outputs);
}

/* End of code generation (_coder_propagateFromKeplerians_mex.c) */
