/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * createAdjacencyMatrix_euclid_distance_block_dense_emxutil.h
 *
 * Code generation for function
 * 'createAdjacencyMatrix_euclid_distance_block_dense_emxutil'
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
void emxEnsureCapacity_real_T(const emlrtStack *sp, emxArray_real_T *emxArray,
                              int32_T oldNumel,
                              const emlrtRTEInfo *srcLocation);

void emxFree_real_T(const emlrtStack *sp, emxArray_real_T **pEmxArray);

void emxInit_real_T(const emlrtStack *sp, emxArray_real_T **pEmxArray,
                    int32_T numDimensions, const emlrtRTEInfo *srcLocation);

/* End of code generation
 * (createAdjacencyMatrix_euclid_distance_block_dense_emxutil.h) */
