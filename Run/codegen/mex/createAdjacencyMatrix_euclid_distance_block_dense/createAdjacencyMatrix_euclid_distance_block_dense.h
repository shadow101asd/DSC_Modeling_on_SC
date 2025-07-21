/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * createAdjacencyMatrix_euclid_distance_block_dense.h
 *
 * Code generation for function
 * 'createAdjacencyMatrix_euclid_distance_block_dense'
 *
 */

#pragma once

/* Include files */
#include "createAdjacencyMatrix_euclid_distance_block_dense_types.h"
#include "rtwtypes.h"
#include "emlrt.h"
#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Function Declarations */
void createAdjacencyMatrix_euclid_distance_block_dense(
    const emlrtStack *sp, const emxArray_real_T *Xs, const real_T cutoff_data[],
    const int32_T cutoff_size[2], real_T blockSize, emxArray_real_T *A);

/* End of code generation (createAdjacencyMatrix_euclid_distance_block_dense.h)
 */
