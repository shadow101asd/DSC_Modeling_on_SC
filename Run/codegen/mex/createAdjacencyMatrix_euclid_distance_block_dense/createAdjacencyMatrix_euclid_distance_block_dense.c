/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * createAdjacencyMatrix_euclid_distance_block_dense.c
 *
 * Code generation for function
 * 'createAdjacencyMatrix_euclid_distance_block_dense'
 *
 */

/* Include files */
#include "createAdjacencyMatrix_euclid_distance_block_dense.h"
#include "assertValidSizeArg.h"
#include "createAdjacencyMatrix_euclid_distance_block_dense_data.h"
#include "createAdjacencyMatrix_euclid_distance_block_dense_emxutil.h"
#include "createAdjacencyMatrix_euclid_distance_block_dense_types.h"
#include "rt_nonfinite.h"
#include "squeeze.h"
#include "sum.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = {
    27,                                                  /* lineNo */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m" /* pathName */
};

static emlrtRSInfo b_emlrtRSI = {
    44,                                                  /* lineNo */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m" /* pathName */
};

static emlrtRSInfo c_emlrtRSI = {
    45,                                                  /* lineNo */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m" /* pathName */
};

static emlrtRSInfo d_emlrtRSI = {
    48,                                                  /* lineNo */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m" /* pathName */
};

static emlrtRSInfo h_emlrtRSI = {
    40,                  /* lineNo */
    "reshapeSizeChecks", /* fcnName */
    "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/"
    "reshapeSizeChecks.m" /* pathName */
};

static emlrtBCInfo emlrtBCI = {
    -1,                                                  /* iFirst */
    -1,                                                  /* iLast */
    27,                                                  /* lineNo */
    28,                                                  /* colNo */
    "Xs",                                                /* aName */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    0                                                      /* checkKind */
};

static emlrtRTEInfo emlrtRTEI = {
    33,                                                  /* lineNo */
    18,                                                  /* colNo */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m" /* pName */
};

static emlrtDCInfo emlrtDCI = {
    35,                                                  /* lineNo */
    21,                                                  /* colNo */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    1                                                      /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = {
    -1,                                                  /* iFirst */
    -1,                                                  /* iLast */
    35,                                                  /* lineNo */
    21,                                                  /* colNo */
    "Xt",                                                /* aName */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    0                                                      /* checkKind */
};

static emlrtDCInfo b_emlrtDCI = {
    35,                                                  /* lineNo */
    24,                                                  /* colNo */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    1                                                      /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = {
    -1,                                                  /* iFirst */
    -1,                                                  /* iLast */
    35,                                                  /* lineNo */
    24,                                                  /* colNo */
    "Xt",                                                /* aName */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    0                                                      /* checkKind */
};

static emlrtRTEInfo b_emlrtRTEI = {
    37,                                                  /* lineNo */
    22,                                                  /* colNo */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m" /* pName */
};

static emlrtDCInfo c_emlrtDCI = {
    39,                                                  /* lineNo */
    25,                                                  /* colNo */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    1                                                      /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = {
    -1,                                                  /* iFirst */
    -1,                                                  /* iLast */
    39,                                                  /* lineNo */
    25,                                                  /* colNo */
    "Xt",                                                /* aName */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    0                                                      /* checkKind */
};

static emlrtDCInfo d_emlrtDCI = {
    39,                                                  /* lineNo */
    28,                                                  /* colNo */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    1                                                      /* checkKind */
};

static emlrtBCInfo e_emlrtBCI = {
    -1,                                                  /* iFirst */
    -1,                                                  /* iLast */
    39,                                                  /* lineNo */
    28,                                                  /* colNo */
    "Xt",                                                /* aName */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    0                                                      /* checkKind */
};

static emlrtDCInfo e_emlrtDCI = {
    56,                                                  /* lineNo */
    19,                                                  /* colNo */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    1                                                      /* checkKind */
};

static emlrtBCInfo f_emlrtBCI = {
    -1,                                                  /* iFirst */
    -1,                                                  /* iLast */
    56,                                                  /* lineNo */
    19,                                                  /* colNo */
    "A",                                                 /* aName */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    0                                                      /* checkKind */
};

static emlrtDCInfo f_emlrtDCI = {
    56,                                                  /* lineNo */
    22,                                                  /* colNo */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    1                                                      /* checkKind */
};

static emlrtBCInfo g_emlrtBCI = {
    -1,                                                  /* iFirst */
    -1,                                                  /* iLast */
    56,                                                  /* lineNo */
    22,                                                  /* colNo */
    "A",                                                 /* aName */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    0                                                      /* checkKind */
};

static emlrtDCInfo g_emlrtDCI = {
    56,                                                  /* lineNo */
    26,                                                  /* colNo */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    1                                                      /* checkKind */
};

static emlrtBCInfo h_emlrtBCI = {
    -1,                                                  /* iFirst */
    -1,                                                  /* iLast */
    56,                                                  /* lineNo */
    26,                                                  /* colNo */
    "A",                                                 /* aName */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    0                                                      /* checkKind */
};

static emlrtDCInfo h_emlrtDCI = {
    56,                                                  /* lineNo */
    29,                                                  /* colNo */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    1                                                      /* checkKind */
};

static emlrtBCInfo i_emlrtBCI = {
    -1,                                                  /* iFirst */
    -1,                                                  /* iLast */
    56,                                                  /* lineNo */
    29,                                                  /* colNo */
    "A",                                                 /* aName */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    0                                                      /* checkKind */
};

static emlrtBCInfo j_emlrtBCI = {
    -1,                                                  /* iFirst */
    -1,                                                  /* iLast */
    56,                                                  /* lineNo */
    33,                                                  /* colNo */
    "A",                                                 /* aName */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    0                                                      /* checkKind */
};

static emlrtECInfo emlrtECI = {
    -1,                                                  /* nDims */
    56,                                                  /* lineNo */
    17,                                                  /* colNo */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m" /* pName */
};

static emlrtDCInfo i_emlrtDCI = {
    58,                                                  /* lineNo */
    23,                                                  /* colNo */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    1                                                      /* checkKind */
};

static emlrtBCInfo k_emlrtBCI = {
    -1,                                                  /* iFirst */
    -1,                                                  /* iLast */
    58,                                                  /* lineNo */
    23,                                                  /* colNo */
    "A",                                                 /* aName */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    0                                                      /* checkKind */
};

static emlrtDCInfo j_emlrtDCI = {
    58,                                                  /* lineNo */
    26,                                                  /* colNo */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    1                                                      /* checkKind */
};

static emlrtBCInfo l_emlrtBCI = {
    -1,                                                  /* iFirst */
    -1,                                                  /* iLast */
    58,                                                  /* lineNo */
    26,                                                  /* colNo */
    "A",                                                 /* aName */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    0                                                      /* checkKind */
};

static emlrtDCInfo k_emlrtDCI = {
    58,                                                  /* lineNo */
    30,                                                  /* colNo */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    1                                                      /* checkKind */
};

static emlrtBCInfo m_emlrtBCI = {
    -1,                                                  /* iFirst */
    -1,                                                  /* iLast */
    58,                                                  /* lineNo */
    30,                                                  /* colNo */
    "A",                                                 /* aName */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    0                                                      /* checkKind */
};

static emlrtDCInfo l_emlrtDCI = {
    58,                                                  /* lineNo */
    33,                                                  /* colNo */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    1                                                      /* checkKind */
};

static emlrtBCInfo n_emlrtBCI = {
    -1,                                                  /* iFirst */
    -1,                                                  /* iLast */
    58,                                                  /* lineNo */
    33,                                                  /* colNo */
    "A",                                                 /* aName */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    0                                                      /* checkKind */
};

static emlrtBCInfo o_emlrtBCI = {
    -1,                                                  /* iFirst */
    -1,                                                  /* iLast */
    58,                                                  /* lineNo */
    37,                                                  /* colNo */
    "A",                                                 /* aName */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    0                                                      /* checkKind */
};

static emlrtECInfo b_emlrtECI = {
    -1,                                                  /* nDims */
    58,                                                  /* lineNo */
    21,                                                  /* colNo */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m" /* pName */
};

static emlrtRTEInfo d_emlrtRTEI = {
    79,                  /* lineNo */
    23,                  /* colNo */
    "reshapeSizeChecks", /* fName */
    "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/"
    "reshapeSizeChecks.m" /* pName */
};

static emlrtRTEInfo f_emlrtRTEI = {
    13,     /* lineNo */
    9,      /* colNo */
    "sqrt", /* fName */
    "/Applications/MATLAB_R2024a.app/toolbox/eml/lib/matlab/elfun/sqrt.m" /* pName
                                                                           */
};

static emlrtBCInfo p_emlrtBCI = {
    -1,                                                  /* iFirst */
    -1,                                                  /* iLast */
    30,                                                  /* lineNo */
    31,                                                  /* colNo */
    "cutoff",                                            /* aName */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    0                                                      /* checkKind */
};

static emlrtBCInfo q_emlrtBCI = {
    -1,                                                  /* iFirst */
    -1,                                                  /* iLast */
    52,                                                  /* lineNo */
    27,                                                  /* colNo */
    "dists",                                             /* aName */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m", /* pName */
    0                                                      /* checkKind */
};

static emlrtRTEInfo i_emlrtRTEI = {
    24,                                                  /* lineNo */
    5,                                                   /* colNo */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m" /* pName */
};

static emlrtRTEInfo j_emlrtRTEI = {
    47,                                                  /* lineNo */
    25,                                                  /* colNo */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m" /* pName */
};

static emlrtRTEInfo k_emlrtRTEI = {
    31,            /* lineNo */
    30,            /* colNo */
    "unsafeSxfun", /* fName */
    "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/"
    "unsafeSxfun.m" /* pName */
};

static emlrtRTEInfo l_emlrtRTEI = {
    58,                                                  /* lineNo */
    42,                                                  /* colNo */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m" /* pName */
};

static emlrtRTEInfo m_emlrtRTEI = {
    48,                                                  /* lineNo */
    17,                                                  /* colNo */
    "createAdjacencyMatrix_euclid_distance_block_dense", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "createAdjacencyMatrix_euclid_distance_block_dense.m" /* pName */
};

/* Function Definitions */
void createAdjacencyMatrix_euclid_distance_block_dense(
    const emlrtStack *sp, const emxArray_real_T *Xs, const real_T cutoff_data[],
    const int32_T cutoff_size[2], real_T blockSize, emxArray_real_T *A)
{
  emlrtStack b_st;
  emlrtStack st;
  emxArray_real_T *Xt;
  emxArray_real_T *b_dists;
  emxArray_real_T *dists;
  emxArray_real_T *r;
  real_T tmp_data[4506];
  real_T dv[3];
  real_T dv1[3];
  const real_T *Xs_data;
  real_T *A_data;
  real_T *b_dists_data;
  real_T *dists_data;
  int32_T Xs_size[3];
  int32_T iv[2];
  int32_T tmp_size[2];
  int32_T N_tmp;
  int32_T b_i1;
  int32_T b_j1;
  int32_T b_nx;
  int32_T i;
  int32_T i1;
  int32_T i2;
  int32_T i3;
  int32_T i5;
  int32_T i6;
  int32_T i7;
  int32_T loop_ub;
  int32_T n;
  int32_T n_tmp;
  int32_T nx;
  int32_T t;
  int32_T unnamed_idx_0;
  boolean_T out;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  Xs_data = Xs->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  /* CREATEADJACENCYMATRIX_EUCLID_DISTANCE_BLOCK_DENSE */
  /*    Computes symmetric, time-varying adjacency matrix of Euclidean distances
   */
  /*    in blocks, using dense output matrix. Should have an identical output */
  /*    to createAdjacencyMatrix_euclid_distance.  */
  /*  */
  /*    Inputs: */
  /*        Xs        - [3, nT, N] positions of N nodes over nT timesteps */
  /*        cutoff    - [1 x nT] vector; distances above cutoff set to 0 */
  /*        blockSize - scalar; number of nodes to process per block */
  /*  */
  /*    Output: */
  /*        A - [N, N, nT] dense symmetric adjacency matrix */
  /*  default value */
  N_tmp = Xs->size[2];
  i = A->size[0] * A->size[1] * A->size[2];
  A->size[0] = Xs->size[2];
  A->size[1] = Xs->size[2];
  i1 = Xs->size[1];
  A->size[2] = Xs->size[1];
  emxEnsureCapacity_real_T(sp, A, i, &i_emlrtRTEI);
  A_data = A->data;
  loop_ub = Xs->size[2] * Xs->size[2] * Xs->size[1];
  for (i = 0; i < loop_ub; i++) {
    A_data[i] = 0.0;
  }
  /*  Preallocate full result */
  emxInit_real_T(sp, &dists, 2, &m_emlrtRTEI);
  emxInit_real_T(sp, &Xt, 3, &j_emlrtRTEI);
  emxInit_real_T(sp, &r, 3, &k_emlrtRTEI);
  emxInit_real_T(sp, &b_dists, 2, &l_emlrtRTEI);
  if (i1 - 1 >= 0) {
    Xs_size[0] = 3;
    Xs_size[1] = 1;
    Xs_size[2] = N_tmp;
    i2 = (int32_T)(((real_T)N_tmp + (blockSize - 1.0)) / blockSize);
  }
  for (t = 0; t < i1; t++) {
    real_T Xt_data[4506];
    real_T cutoff_t;
    if (t + 1 > i1) {
      emlrtDynamicBoundsCheckR2012b(t + 1, 1, i1, &emlrtBCI, (emlrtConstCTX)sp);
    }
    for (i = 0; i < N_tmp; i++) {
      Xt_data[3 * i] = Xs_data[3 * t + 3 * Xs->size[1] * i];
      Xt_data[3 * i + 1] = Xs_data[(3 * t + 3 * Xs->size[1] * i) + 1];
      Xt_data[3 * i + 2] = Xs_data[(3 * t + 3 * Xs->size[1] * i) + 2];
    }
    st.site = &emlrtRSI;
    squeeze(&st, Xt_data, Xs_size, tmp_data, tmp_size);
    loop_ub = tmp_size[1];
    for (i = 0; i < 3; i++) {
      for (i3 = 0; i3 < loop_ub; i3++) {
        Xt_data[i3 + loop_ub * i] = tmp_data[i + 3 * i3];
      }
    }
    /*  [N x 3] */
    cutoff_t = rtInf;
    if (cutoff_size[1] != 0) {
      if (t + 1 > cutoff_size[1]) {
        emlrtDynamicBoundsCheckR2012b(t + 1, 1, cutoff_size[1], &p_emlrtBCI,
                                      (emlrtConstCTX)sp);
      }
      cutoff_t = cutoff_data[t];
    }
    emlrtForLoopVectorCheckR2021a(1.0, blockSize, N_tmp, mxDOUBLE_CLASS, i2,
                                  &emlrtRTEI, (emlrtConstCTX)sp);
    for (b_i1 = 0; b_i1 < i2; b_i1++) {
      real_T b_i2;
      real_T c_i1;
      int32_T i4;
      c_i1 = (real_T)b_i1 * blockSize + 1.0;
      b_i2 = muDoubleScalarMin((c_i1 + blockSize) - 1.0, N_tmp);
      if (c_i1 > b_i2) {
        i = 0;
        i3 = 0;
      } else {
        if (c_i1 != (int32_T)muDoubleScalarFloor(c_i1)) {
          emlrtIntegerCheckR2012b(c_i1, &emlrtDCI, (emlrtConstCTX)sp);
        }
        if (((int32_T)c_i1 < 1) || ((int32_T)c_i1 > loop_ub)) {
          emlrtDynamicBoundsCheckR2012b((int32_T)c_i1, 1, loop_ub, &b_emlrtBCI,
                                        (emlrtConstCTX)sp);
        }
        i = (int32_T)c_i1 - 1;
        if (b_i2 != (int32_T)muDoubleScalarFloor(b_i2)) {
          emlrtIntegerCheckR2012b(b_i2, &b_emlrtDCI, (emlrtConstCTX)sp);
        }
        if (((int32_T)b_i2 < 1) || ((int32_T)b_i2 > loop_ub)) {
          emlrtDynamicBoundsCheckR2012b((int32_T)b_i2, 1, loop_ub, &c_emlrtBCI,
                                        (emlrtConstCTX)sp);
        }
        i3 = (int32_T)b_i2;
      }
      /*  [Bi x 3] */
      i4 = (int32_T)(((real_T)N_tmp + (blockSize - c_i1)) / blockSize);
      emlrtForLoopVectorCheckR2021a(c_i1, blockSize, N_tmp, mxDOUBLE_CLASS, i4,
                                    &b_emlrtRTEI, (emlrtConstCTX)sp);
      if (i4 - 1 >= 0) {
        dv[0] = (b_i2 - c_i1) + 1.0;
        dv[1] = 1.0;
        dv[2] = 3.0;
        nx = (i3 - i) * 3;
        out = ((int32_T)dv[0] >= 0);
        dv1[0] = 1.0;
        dv1[2] = 3.0;
        unnamed_idx_0 = (int32_T)dv[0];
        n_tmp = i3 - i;
      }
      for (b_j1 = 0; b_j1 < i4; b_j1++) {
        real_T c_j1;
        real_T j2;
        int32_T b_loop_ub;
        boolean_T p;
        c_j1 = c_i1 + (real_T)b_j1 * blockSize;
        j2 = muDoubleScalarMin((c_j1 + blockSize) - 1.0, N_tmp);
        if (c_j1 > j2) {
          i3 = 0;
          i5 = 0;
        } else {
          if (c_j1 != (int32_T)muDoubleScalarFloor(c_j1)) {
            emlrtIntegerCheckR2012b(c_j1, &c_emlrtDCI, (emlrtConstCTX)sp);
          }
          if (((int32_T)c_j1 < 1) || ((int32_T)c_j1 > loop_ub)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)c_j1, 1, loop_ub,
                                          &d_emlrtBCI, (emlrtConstCTX)sp);
          }
          i3 = (int32_T)c_j1 - 1;
          if (j2 != (int32_T)muDoubleScalarFloor(j2)) {
            emlrtIntegerCheckR2012b(j2, &d_emlrtDCI, (emlrtConstCTX)sp);
          }
          if (((int32_T)j2 < 1) || ((int32_T)j2 > loop_ub)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)j2, 1, loop_ub, &e_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          i5 = (int32_T)j2;
        }
        /*  [Bj x 3] */
        st.site = &b_emlrtRSI;
        b_st.site = &h_emlrtRSI;
        assertValidSizeArg(&b_st, dv);
        n = n_tmp;
        if (n_tmp < 3) {
          n = 3;
        }
        if ((int32_T)dv[0] > muIntScalarMax_sint32(nx, n)) {
          emlrtErrorWithMessageIdR2018a(
              &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
              "Coder:toolbox:reshape_emptyReshapeLimit", 0);
        }
        if (!out) {
          emlrtErrorWithMessageIdR2018a(
              &st, &d_emlrtRTEI, "MATLAB:checkDimCommon:nonnegativeSize",
              "MATLAB:checkDimCommon:nonnegativeSize", 0);
        }
        if ((int32_T)dv[0] * 3 != nx) {
          emlrtErrorWithMessageIdR2018a(
              &st, &c_emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
              "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
        }
        dv1[1] = (j2 - c_j1) + 1.0;
        st.site = &c_emlrtRSI;
        n = i5 - i3;
        b_nx = n * 3;
        b_st.site = &h_emlrtRSI;
        assertValidSizeArg(&b_st, dv1);
        if (n < 3) {
          n = 3;
        }
        if ((int32_T)dv1[1] > muIntScalarMax_sint32(b_nx, n)) {
          emlrtErrorWithMessageIdR2018a(
              &st, &e_emlrtRTEI, "Coder:toolbox:reshape_emptyReshapeLimit",
              "Coder:toolbox:reshape_emptyReshapeLimit", 0);
        }
        if ((int32_T)dv1[1] < 0) {
          emlrtErrorWithMessageIdR2018a(
              &st, &d_emlrtRTEI, "MATLAB:checkDimCommon:nonnegativeSize",
              "MATLAB:checkDimCommon:nonnegativeSize", 0);
        }
        if ((int32_T)dv1[1] * 3 != b_nx) {
          emlrtErrorWithMessageIdR2018a(
              &st, &c_emlrtRTEI, "Coder:MATLAB:getReshapeDims_notSameNumel",
              "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
        }
        /*  [Bi x Bj x 3] */
        st.site = &d_emlrtRSI;
        i5 = Xt->size[0] * Xt->size[1] * Xt->size[2];
        Xt->size[0] = unnamed_idx_0;
        b_loop_ub = (int32_T)dv1[1];
        Xt->size[1] = (int32_T)dv1[1];
        Xt->size[2] = 3;
        emxEnsureCapacity_real_T(&st, Xt, i5, &j_emlrtRTEI);
        dists_data = Xt->data;
        for (i5 = 0; i5 < 3; i5++) {
          for (i6 = 0; i6 < b_loop_ub; i6++) {
            for (i7 = 0; i7 < unnamed_idx_0; i7++) {
              n = loop_ub * i5;
              dists_data[(i7 + Xt->size[0] * i6) +
                         Xt->size[0] * Xt->size[1] * i5] =
                  Xt_data[(i + i7) + n] - Xt_data[(i3 + i6) + n];
            }
          }
        }
        i3 = r->size[0] * r->size[1] * r->size[2];
        r->size[0] = Xt->size[0];
        r->size[1] = (int32_T)dv1[1];
        r->size[2] = 3;
        emxEnsureCapacity_real_T(&st, r, i3, &k_emlrtRTEI);
        b_dists_data = r->data;
        b_loop_ub = Xt->size[0] * Xt->size[1] * 3;
        for (i3 = 0; i3 < b_loop_ub; i3++) {
          real_T varargin_1;
          varargin_1 = dists_data[i3];
          b_dists_data[i3] = varargin_1 * varargin_1;
        }
        b_st.site = &d_emlrtRSI;
        sum(&b_st, r, dists);
        dists_data = dists->data;
        p = false;
        i3 = dists->size[0] * dists->size[1];
        for (n = 0; n < i3; n++) {
          if (p || (dists_data[n] < 0.0)) {
            p = true;
          }
        }
        if (p) {
          emlrtErrorWithMessageIdR2018a(
              &st, &f_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
              "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
        }
        for (n = 0; n < i3; n++) {
          dists_data[n] = muDoubleScalarSqrt(dists_data[n]);
        }
        /*  [Bi x Bj] */
        /*  Apply cutoff if needed */
        if ((!muDoubleScalarIsInf(cutoff_t)) &&
            (!muDoubleScalarIsNaN(cutoff_t))) {
          n = i3 - 1;
          for (b_nx = 0; b_nx <= n; b_nx++) {
            if (dists_data[b_nx] > cutoff_t) {
              i3 = dists->size[0] * dists->size[1] - 1;
              if (b_nx > i3) {
                emlrtDynamicBoundsCheckR2012b(b_nx, 0, i3, &q_emlrtBCI,
                                              (emlrtConstCTX)sp);
              }
              dists_data[b_nx] = 0.0;
            }
          }
        }
        /*  Fill full matrix (both halves for symmetry) */
        if (c_i1 > b_i2) {
          i3 = 0;
          i5 = 0;
        } else {
          if (c_i1 != (int32_T)muDoubleScalarFloor(c_i1)) {
            emlrtIntegerCheckR2012b(c_i1, &e_emlrtDCI, (emlrtConstCTX)sp);
          }
          if (((int32_T)c_i1 < 1) || ((int32_T)c_i1 > N_tmp)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)c_i1, 1, N_tmp, &f_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          i3 = (int32_T)c_i1 - 1;
          if (b_i2 != (int32_T)muDoubleScalarFloor(b_i2)) {
            emlrtIntegerCheckR2012b(b_i2, &f_emlrtDCI, (emlrtConstCTX)sp);
          }
          if (((int32_T)b_i2 < 1) || ((int32_T)b_i2 > N_tmp)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)b_i2, 1, N_tmp, &g_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          i5 = (int32_T)b_i2;
        }
        if (c_j1 > j2) {
          i6 = 0;
          i7 = 0;
        } else {
          if (c_j1 != (int32_T)muDoubleScalarFloor(c_j1)) {
            emlrtIntegerCheckR2012b(c_j1, &g_emlrtDCI, (emlrtConstCTX)sp);
          }
          if (((int32_T)c_j1 < 1) || ((int32_T)c_j1 > N_tmp)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)c_j1, 1, N_tmp, &h_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          i6 = (int32_T)c_j1 - 1;
          if (j2 != (int32_T)muDoubleScalarFloor(j2)) {
            emlrtIntegerCheckR2012b(j2, &h_emlrtDCI, (emlrtConstCTX)sp);
          }
          if (((int32_T)j2 < 1) || ((int32_T)j2 > N_tmp)) {
            emlrtDynamicBoundsCheckR2012b((int32_T)j2, 1, N_tmp, &i_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          i7 = (int32_T)j2;
        }
        if (t + 1 > i1) {
          emlrtDynamicBoundsCheckR2012b(t + 1, 1, i1, &j_emlrtBCI,
                                        (emlrtConstCTX)sp);
        }
        b_loop_ub = i5 - i3;
        iv[0] = b_loop_ub;
        n = i7 - i6;
        iv[1] = n;
        emlrtSubAssignSizeCheckR2012b(&iv[0], 2, &dists->size[0], 2, &emlrtECI,
                                      (emlrtCTX)sp);
        for (i5 = 0; i5 < n; i5++) {
          for (i7 = 0; i7 < b_loop_ub; i7++) {
            A_data[((i3 + i7) + A->size[0] * (i6 + i5)) +
                   A->size[0] * A->size[1] * t] =
                dists_data[i7 + b_loop_ub * i5];
          }
        }
        if (c_i1 != c_j1) {
          int32_T dists_tmp;
          if (c_j1 > j2) {
            i3 = 0;
            i5 = 0;
          } else {
            if (c_j1 != (int32_T)muDoubleScalarFloor(c_j1)) {
              emlrtIntegerCheckR2012b(c_j1, &i_emlrtDCI, (emlrtConstCTX)sp);
            }
            if (((int32_T)c_j1 < 1) || ((int32_T)c_j1 > N_tmp)) {
              emlrtDynamicBoundsCheckR2012b((int32_T)c_j1, 1, N_tmp,
                                            &k_emlrtBCI, (emlrtConstCTX)sp);
            }
            i3 = (int32_T)c_j1 - 1;
            if (j2 != (int32_T)muDoubleScalarFloor(j2)) {
              emlrtIntegerCheckR2012b(j2, &j_emlrtDCI, (emlrtConstCTX)sp);
            }
            if (((int32_T)j2 < 1) || ((int32_T)j2 > N_tmp)) {
              emlrtDynamicBoundsCheckR2012b((int32_T)j2, 1, N_tmp, &l_emlrtBCI,
                                            (emlrtConstCTX)sp);
            }
            i5 = (int32_T)j2;
          }
          if (c_i1 > b_i2) {
            i6 = 0;
            i7 = 0;
          } else {
            if (c_i1 != (int32_T)muDoubleScalarFloor(c_i1)) {
              emlrtIntegerCheckR2012b(c_i1, &k_emlrtDCI, (emlrtConstCTX)sp);
            }
            if (((int32_T)c_i1 < 1) || ((int32_T)c_i1 > N_tmp)) {
              emlrtDynamicBoundsCheckR2012b((int32_T)c_i1, 1, N_tmp,
                                            &m_emlrtBCI, (emlrtConstCTX)sp);
            }
            i6 = (int32_T)c_i1 - 1;
            if (b_i2 != (int32_T)muDoubleScalarFloor(b_i2)) {
              emlrtIntegerCheckR2012b(b_i2, &l_emlrtDCI, (emlrtConstCTX)sp);
            }
            if (((int32_T)b_i2 < 1) || ((int32_T)b_i2 > N_tmp)) {
              emlrtDynamicBoundsCheckR2012b((int32_T)b_i2, 1, N_tmp,
                                            &n_emlrtBCI, (emlrtConstCTX)sp);
            }
            i7 = (int32_T)b_i2;
          }
          if (t + 1 > i1) {
            emlrtDynamicBoundsCheckR2012b(t + 1, 1, i1, &o_emlrtBCI,
                                          (emlrtConstCTX)sp);
          }
          b_loop_ub = i5 - i3;
          iv[0] = b_loop_ub;
          n = i7 - i6;
          iv[1] = n;
          b_nx = dists->size[1];
          tmp_size[0] = dists->size[1];
          dists_tmp = dists->size[0];
          tmp_size[1] = dists->size[0];
          emlrtSubAssignSizeCheckR2012b(&iv[0], 2, &tmp_size[0], 2, &b_emlrtECI,
                                        (emlrtCTX)sp);
          i5 = b_dists->size[0] * b_dists->size[1];
          b_dists->size[0] = dists->size[1];
          b_dists->size[1] = dists->size[0];
          emxEnsureCapacity_real_T(sp, b_dists, i5, &l_emlrtRTEI);
          b_dists_data = b_dists->data;
          for (i5 = 0; i5 < dists_tmp; i5++) {
            for (i7 = 0; i7 < b_nx; i7++) {
              b_dists_data[i7 + b_dists->size[0] * i5] =
                  dists_data[i5 + dists->size[0] * i7];
            }
          }
          for (i5 = 0; i5 < n; i5++) {
            for (i7 = 0; i7 < b_loop_ub; i7++) {
              A_data[((i3 + i7) + A->size[0] * (i6 + i5)) +
                     A->size[0] * A->size[1] * t] =
                  b_dists_data[i7 + b_loop_ub * i5];
            }
          }
          /*  symmetry */
        }
        if (*emlrtBreakCheckR2012bFlagVar != 0) {
          emlrtBreakCheckR2012b((emlrtConstCTX)sp);
        }
      }
      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b((emlrtConstCTX)sp);
      }
    }
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b((emlrtConstCTX)sp);
    }
  }
  emxFree_real_T(sp, &b_dists);
  emxFree_real_T(sp, &r);
  emxFree_real_T(sp, &Xt);
  emxFree_real_T(sp, &dists);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

/* End of code generation (createAdjacencyMatrix_euclid_distance_block_dense.c)
 */
