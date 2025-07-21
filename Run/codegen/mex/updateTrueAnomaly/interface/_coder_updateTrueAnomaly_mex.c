/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_updateTrueAnomaly_mex.c
 *
 * Code generation for function '_coder_updateTrueAnomaly_mex'
 *
 */

/* Include files */
#include "_coder_updateTrueAnomaly_mex.h"
#include "_coder_updateTrueAnomaly_api.h"
#include "rt_nonfinite.h"
#include "updateTrueAnomaly_data.h"
#include "updateTrueAnomaly_initialize.h"
#include "updateTrueAnomaly_terminate.h"

/* Function Definitions */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&updateTrueAnomaly_atexit);
  /* Module initialization. */
  updateTrueAnomaly_initialize();
  /* Dispatch the entry-point. */
  updateTrueAnomaly_mexFunction(nlhs, plhs, nrhs, prhs);
  /* Module termination. */
  updateTrueAnomaly_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL, "UTF-8", true);
  return emlrtRootTLSGlobal;
}

void updateTrueAnomaly_mexFunction(int32_T nlhs, mxArray *plhs[1], int32_T nrhs,
                                   const mxArray *prhs[8])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *outputs;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 8) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 8, 4,
                        17, "updateTrueAnomaly");
  }
  if (nlhs > 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 17,
                        "updateTrueAnomaly");
  }
  /* Call the function. */
  updateTrueAnomaly_api(prhs, &outputs);
  /* Copy over outputs to the caller. */
  emlrtReturnArrays(1, &plhs[0], &outputs);
}

/* End of code generation (_coder_updateTrueAnomaly_mex.c) */
