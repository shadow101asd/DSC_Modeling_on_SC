/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * createAdjacencyMatrix_euclid_distance_block_dense_terminate.c
 *
 * Code generation for function
 * 'createAdjacencyMatrix_euclid_distance_block_dense_terminate'
 *
 */

/* Include files */
#include "createAdjacencyMatrix_euclid_distance_block_dense_terminate.h"
#include "_coder_createAdjacencyMatrix_euclid_distance_block_dense_mex.h"
#include "createAdjacencyMatrix_euclid_distance_block_dense_data.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void emlrtExitTimeCleanupDtorFcn(const void *r);

/* Function Definitions */
static void emlrtExitTimeCleanupDtorFcn(const void *r)
{
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void createAdjacencyMatrix_euclid_distance_block_dense_atexit(void)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtPushHeapReferenceStackR2021a(
      &st, false, NULL, (void *)&emlrtExitTimeCleanupDtorFcn, NULL, NULL, NULL);
  emlrtEnterRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
  emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void createAdjacencyMatrix_euclid_distance_block_dense_terminate(void)
{
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation
 * (createAdjacencyMatrix_euclid_distance_block_dense_terminate.c) */
