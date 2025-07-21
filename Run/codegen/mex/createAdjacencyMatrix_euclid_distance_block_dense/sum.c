/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * sum.c
 *
 * Code generation for function 'sum'
 *
 */

/* Include files */
#include "sum.h"
#include "createAdjacencyMatrix_euclid_distance_block_dense_emxutil.h"
#include "createAdjacencyMatrix_euclid_distance_block_dense_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRTEInfo n_emlrtRTEI = {
    20,    /* lineNo */
    1,     /* colNo */
    "sum", /* fName */
    "/Applications/MATLAB_R2024a.app/toolbox/eml/lib/matlab/datafun/sum.m" /* pName
                                                                            */
};

static emlrtRTEInfo o_emlrtRTEI = {
    146,                /* lineNo */
    24,                 /* colNo */
    "blockedSummation", /* fName */
    "/Applications/MATLAB_R2024a.app/toolbox/eml/lib/matlab/datafun/private/"
    "blockedSummation.m" /* pName */
};

/* Function Definitions */
void sum(const emlrtStack *sp, const emxArray_real_T *x, emxArray_real_T *y)
{
  const real_T *x_data;
  real_T *y_data;
  int32_T k;
  int32_T xj;
  x_data = x->data;
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    int32_T xoffset;
    xj = y->size[0] * y->size[1];
    y->size[0] = (int16_T)x->size[0];
    y->size[1] = (int16_T)x->size[1];
    emxEnsureCapacity_real_T(sp, y, xj, &n_emlrtRTEI);
    y_data = y->data;
    xoffset = (int16_T)x->size[0] * (int16_T)x->size[1];
    for (xj = 0; xj < xoffset; xj++) {
      y_data[xj] = 0.0;
    }
  } else {
    int32_T vstride_tmp;
    vstride_tmp = x->size[0] * x->size[1];
    xj = y->size[0] * y->size[1];
    y->size[0] = (int16_T)x->size[0];
    y->size[1] = (int16_T)x->size[1];
    emxEnsureCapacity_real_T(sp, y, xj, &o_emlrtRTEI);
    y_data = y->data;
    for (xj = 0; xj < vstride_tmp; xj++) {
      y_data[xj] = x_data[xj];
    }
    for (k = 0; k < 2; k++) {
      int32_T xoffset;
      xoffset = (k + 1) * vstride_tmp;
      for (xj = 0; xj < vstride_tmp; xj++) {
        y_data[xj] += x_data[xoffset + xj];
      }
    }
  }
}

/* End of code generation (sum.c) */
