/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * squeeze.c
 *
 * Code generation for function 'squeeze'
 *
 */

/* Include files */
#include "squeeze.h"
#include "createAdjacencyMatrix_euclid_distance_block_dense_data.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo
    e_emlrtRSI =
        {
            38,        /* lineNo */
            "squeeze", /* fcnName */
            "/Applications/MATLAB_R2024a.app/toolbox/eml/lib/matlab/elmat/"
            "squeeze.m" /* pathName */
};

/* Function Definitions */
void squeeze(const emlrtStack *sp, const real_T a_data[],
             const int32_T a_size[3], real_T b_data[], int32_T b_size[2])
{
  emlrtStack st;
  int32_T n;
  int32_T nx;
  int16_T szb_idx_1;
  st.prev = sp;
  st.tls = sp->tls;
  szb_idx_1 = 1;
  if (a_size[2] != 1) {
    szb_idx_1 = (int16_T)a_size[2];
  }
  st.site = &e_emlrtRSI;
  nx = 3 * a_size[2];
  n = 3;
  if (a_size[2] > 3) {
    n = a_size[2];
  }
  if (szb_idx_1 > muIntScalarMax_sint32(nx, n)) {
    emlrtErrorWithMessageIdR2018a(&st, &e_emlrtRTEI,
                                  "Coder:toolbox:reshape_emptyReshapeLimit",
                                  "Coder:toolbox:reshape_emptyReshapeLimit", 0);
  }
  n = 3 * szb_idx_1;
  if (n != nx) {
    emlrtErrorWithMessageIdR2018a(
        &st, &c_emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
        "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }
  b_size[0] = 3;
  b_size[1] = szb_idx_1;
  if (n - 1 >= 0) {
    memcpy(&b_data[0], &a_data[0], (uint32_T)n * sizeof(real_T));
  }
}

/* End of code generation (squeeze.c) */
