/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_Keplerian2Cartesian_mex.c
 *
 * Code generation for function '_coder_Keplerian2Cartesian_mex'
 *
 */

/* Include files */
#include "_coder_Keplerian2Cartesian_mex.h"
#include "Keplerian2Cartesian_data.h"
#include "Keplerian2Cartesian_initialize.h"
#include "Keplerian2Cartesian_terminate.h"
#include "_coder_Keplerian2Cartesian_api.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void Keplerian2Cartesian_mexFunction(int32_T nlhs, mxArray *plhs[1],
                                     int32_T nrhs, const mxArray *prhs[7])
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  const mxArray *outputs;
  st.tls = emlrtRootTLSGlobal;
  /* Check for proper number of arguments. */
  if (nrhs != 7) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 7, 4,
                        19, "Keplerian2Cartesian");
  }
  if (nlhs > 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 19,
                        "Keplerian2Cartesian");
  }
  /* Call the function. */
  Keplerian2Cartesian_api(prhs, &outputs);
  /* Copy over outputs to the caller. */
  emlrtReturnArrays(1, &plhs[0], &outputs);
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&Keplerian2Cartesian_atexit);
  /* Module initialization. */
  Keplerian2Cartesian_initialize();
  /* Dispatch the entry-point. */
  Keplerian2Cartesian_mexFunction(nlhs, plhs, nrhs, prhs);
  /* Module termination. */
  Keplerian2Cartesian_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL, "UTF-8", true);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_Keplerian2Cartesian_mex.c) */
