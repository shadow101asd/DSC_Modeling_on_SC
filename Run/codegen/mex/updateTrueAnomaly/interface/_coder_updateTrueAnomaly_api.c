/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_updateTrueAnomaly_api.c
 *
 * Code generation for function '_coder_updateTrueAnomaly_api'
 *
 */

/* Include files */
#include "_coder_updateTrueAnomaly_api.h"
#include "rt_nonfinite.h"
#include "updateTrueAnomaly.h"
#include "updateTrueAnomaly_data.h"
#include "updateTrueAnomaly_mexutil.h"

/* Function Declarations */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId);

static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
                                 const emlrtMsgIdentifier *msgId);

static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                               const char_T *identifier);

/* Function Definitions */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
                                 const emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = c_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static real_T c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
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

static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
                               const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  real_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

void updateTrueAnomaly_api(const mxArray *const prhs[8], const mxArray **plhs)
{
  emlrtStack st = {
      NULL, /* site */
      NULL, /* tls */
      NULL  /* prev */
  };
  real_T Om;
  real_T a;
  real_T dt;
  real_T e;
  real_T f0;
  real_T i;
  real_T mu;
  real_T w;
  st.tls = emlrtRootTLSGlobal;
  /* Marshall function inputs */
  a = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "a");
  e = emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "e");
  i = emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "i");
  Om = emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "Om");
  w = emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "w");
  f0 = emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "f0");
  mu = emlrt_marshallIn(&st, emlrtAliasP(prhs[6]), "mu");
  dt = emlrt_marshallIn(&st, emlrtAliasP(prhs[7]), "dt");
  /* Invoke the target function */
  a = updateTrueAnomaly(&st, a, e, i, Om, w, f0, mu, dt);
  /* Marshall function outputs */
  *plhs = emlrt_marshallOut(a);
}

/* End of code generation (_coder_updateTrueAnomaly_api.c) */
