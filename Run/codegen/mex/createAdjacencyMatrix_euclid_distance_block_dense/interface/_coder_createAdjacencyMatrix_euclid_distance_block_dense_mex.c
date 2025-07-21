/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_createAdjacencyMatrix_euclid_distance_block_dense_mex.c
 *
 * Code generation for function
 * '_coder_createAdjacencyMatrix_euclid_distance_block_dense_mex'
 *
 */

/* Include files */
#include "_coder_createAdjacencyMatrix_euclid_distance_block_dense_mex.h"
#include "_coder_createAdjacencyMatrix_euclid_distance_block_dense_api.h"
#include "createAdjacencyMatrix_euclid_distance_block_dense_data.h"
#include "createAdjacencyMatrix_euclid_distance_block_dense_initialize.h"
#include "createAdjacencyMatrix_euclid_distance_block_dense_terminate.h"
#include "rt_nonfinite.h"

/* Function Definitions */
void createAdjacencyMatrix_euclid_distance_block_dense_mexFunction(
    int32_T nlhs, mxArray *plhs[1], int32_T nrhs, const mxArray *prhs[3])
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
    emlrtErrMsgIdAndTxt(
        &st, "EMLRT:runTime:WrongNumberOfInputsFAVDefaultValues", 5, 12, 3, 4,
        49, "createAdjacencyMatrix_euclid_distance_block_dense");
  }
  if (nlhs > 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 49,
                        "createAdjacencyMatrix_euclid_distance_block_dense");
  }
  /* Call the function. */
  c_createAdjacencyMatrix_euclid_(prhs, &outputs);
  /* Copy over outputs to the caller. */
  emlrtReturnArrays(1, &plhs[0], &outputs);
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs,
                 const mxArray *prhs[])
{
  mexAtExit(&createAdjacencyMatrix_euclid_distance_block_dense_atexit);
  /* Module initialization. */
  createAdjacencyMatrix_euclid_distance_block_dense_initialize();
  /* Dispatch the entry-point. */
  createAdjacencyMatrix_euclid_distance_block_dense_mexFunction(nlhs, plhs,
                                                                nrhs, prhs);
  /* Module termination. */
  createAdjacencyMatrix_euclid_distance_block_dense_terminate();
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1,
                           NULL, "UTF-8", true);
  return emlrtRootTLSGlobal;
}

/* End of code generation
 * (_coder_createAdjacencyMatrix_euclid_distance_block_dense_mex.c) */
