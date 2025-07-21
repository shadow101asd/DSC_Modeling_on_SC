/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_Keplerian2Cartesian_api.c
 *
 * Code generation for function '_coder_Keplerian2Cartesian_api'
 *
 */

/* Include files */
#include "_coder_Keplerian2Cartesian_api.h"
#include "Keplerian2Cartesian.h"
#include "Keplerian2Cartesian_data.h"
#include "Keplerian2Cartesian_emxutil.h"
#include "Keplerian2Cartesian_types.h"
#include "rt_nonfinite.h"

/* Variable Definitions */
static emlrtRTEInfo hd_emlrtRTEI = {
    1,                                /* lineNo */
    1,                                /* colNo */
    "_coder_Keplerian2Cartesian_api", /* fName */
    ""                                /* pName */
};

/* Function Declarations */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y);

static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                 const char_T *identifier);

static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId);

static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret);

static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                             const char_T *identifier, emxArray_real_T *y);

static const mxArray *emlrt_marshallOut(const emxArray_real_T *u);

static real_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId);

/* Function Definitions */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                               const emlrtMsgIdentifier *parentId,
                               emxArray_real_T *y)
{
  e_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                                 const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  real_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static real_T d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = f_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                               const emlrtMsgIdentifier *msgId,
                               emxArray_real_T *ret)
{
  static const int32_T dims[2] = {1, -1};
  int32_T iv[2];
  int32_T i;
  boolean_T bv[2] = {false, true};
  emlrtCheckVsBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 2U,
                            (const void *)&dims[0], &bv[0], &iv[0]);
  ret->allocatedSize = iv[0] * iv[1];
  i = ret->size[0] * ret->size[1];
  ret->size[0] = iv[0];
  ret->size[1] = iv[1];
  emxEnsureCapacity_real_T(sp, ret, i, (emlrtRTEInfo *)NULL);
  ret->data = (real_T *)emlrtMxGetData(src);
  ret->canFreeData = false;
  emlrtDestroyArray(&src);
}

static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                             const char_T *identifier, emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  b_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId, y);
  emlrtDestroyArray(&nullptr);
}

static const mxArray *emlrt_marshallOut(const emxArray_real_T *u)
{
  static const int32_T iv[2] = {0, 0};
  const mxArray *m;
  const mxArray *y;
  const real_T *u_data;
  u_data = u->data;
  y = NULL;
  m = emlrtCreateNumericArray(2, (const void *)&iv[0], mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m, (void *)&u_data[0]);
  emlrtSetDimensions((mxArray *)m, &u->size[0], 2);
  emlrtAssign(&y, m);
  return y;
}

static real_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId)
{
  static const int32_T dims = 0;
  real_T ret;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 0U,
                          (const void *)&dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

void Keplerian2Cartesian_api(const mxArray *const prhs[7], const mxArray **plhs)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  emxArray_real_T *Om;
  emxArray_real_T *a;
  emxArray_real_T *e;
  emxArray_real_T *f0;
  emxArray_real_T *i;
  emxArray_real_T *w;
  emxArray_real_T *y;
  const mxArray *prhs_copy_idx_3;
  const mxArray *prhs_copy_idx_4;
  real_T mu;
  st.tls = emlrtRootTLSGlobal;
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  prhs_copy_idx_3 = emlrtProtectR2012b(prhs[3], 3, false, -1);
  prhs_copy_idx_4 = emlrtProtectR2012b(prhs[4], 4, false, -1);
  /* Marshall function inputs */
  emxInit_real_T(&st, &a, &hd_emlrtRTEI);
  a->canFreeData = false;
  emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "a", a);
  emxInit_real_T(&st, &e, &hd_emlrtRTEI);
  e->canFreeData = false;
  emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "e", e);
  emxInit_real_T(&st, &i, &hd_emlrtRTEI);
  i->canFreeData = false;
  emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "i", i);
  emxInit_real_T(&st, &Om, &hd_emlrtRTEI);
  Om->canFreeData = false;
  emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_3), "Om", Om);
  emxInit_real_T(&st, &w, &hd_emlrtRTEI);
  w->canFreeData = false;
  emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_4), "w", w);
  emxInit_real_T(&st, &f0, &hd_emlrtRTEI);
  f0->canFreeData = false;
  emlrt_marshallIn(&st, emlrtAlias(prhs[5]), "f0", f0);
  mu = c_emlrt_marshallIn(&st, emlrtAliasP(prhs[6]), "mu");
  /* Invoke the target function */
  emxInit_real_T(&st, &y, &hd_emlrtRTEI);
  Keplerian2Cartesian(&st, a, e, i, Om, w, f0, mu, y);
  emxFree_real_T(&st, &f0);
  emxFree_real_T(&st, &w);
  emxFree_real_T(&st, &Om);
  emxFree_real_T(&st, &i);
  emxFree_real_T(&st, &e);
  emxFree_real_T(&st, &a);
  /* Marshall function outputs */
  y->canFreeData = false;
  *plhs = emlrt_marshallOut(y);
  emxFree_real_T(&st, &y);
  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

/* End of code generation (_coder_Keplerian2Cartesian_api.c) */
