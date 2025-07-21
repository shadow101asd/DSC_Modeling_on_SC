/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * createAdjacencyMatrix_euclid_distance_block_dense_initialize.c
 *
 * Code generation for function
 * 'createAdjacencyMatrix_euclid_distance_block_dense_initialize'
 *
 */

/* Include files */
#include "createAdjacencyMatrix_euclid_distance_block_dense_initialize.h"
#include "_coder_createAdjacencyMatrix_euclid_distance_block_dense_mex.h"
#include "createAdjacencyMatrix_euclid_distance_block_dense_data.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void createAdjacencyMatrix_euclid_distance_block_dense_once(void);

/* Function Definitions */
static void createAdjacencyMatrix_euclid_distance_block_dense_once(void)
{
  mex_InitInfAndNan();
}

void createAdjacencyMatrix_euclid_distance_block_dense_initialize(void)
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
    createAdjacencyMatrix_euclid_distance_block_dense_once();
  }
}

/* End of code generation
 * (createAdjacencyMatrix_euclid_distance_block_dense_initialize.c) */
