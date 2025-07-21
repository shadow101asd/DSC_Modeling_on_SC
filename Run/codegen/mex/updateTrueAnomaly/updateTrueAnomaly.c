/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * updateTrueAnomaly.c
 *
 * Code generation for function 'updateTrueAnomaly'
 *
 */

/* Include files */
#include "updateTrueAnomaly.h"
#include "rt_nonfinite.h"
#include "updateTrueAnomaly_data.h"
#include "updateTrueAnomaly_mexutil.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = {
    9,                   /* lineNo */
    "updateTrueAnomaly", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "updateTrueAnomaly.m" /* pathName */
};

static emlrtRSInfo b_emlrtRSI = {
    17,                  /* lineNo */
    "updateTrueAnomaly", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "updateTrueAnomaly.m" /* pathName */
};

static emlrtRSInfo c_emlrtRSI = {
    44,                  /* lineNo */
    "updateTrueAnomaly", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "updateTrueAnomaly.m" /* pathName */
};

static emlrtRSInfo d_emlrtRSI = {
    50,                  /* lineNo */
    "updateTrueAnomaly", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "updateTrueAnomaly.m" /* pathName */
};

static emlrtRSInfo e_emlrtRSI = {
    51,                  /* lineNo */
    "updateTrueAnomaly", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "updateTrueAnomaly.m" /* pathName */
};

static emlrtRSInfo
    f_emlrtRSI =
        {
            44,       /* lineNo */
            "mpower", /* fcnName */
            "/Applications/MATLAB_R2024a.app/toolbox/eml/lib/matlab/matfun/"
            "mpower.m" /* pathName */
};

static emlrtRSInfo
    h_emlrtRSI =
        {
            53,        /* lineNo */
            "warning", /* fcnName */
            "/Applications/MATLAB_R2024a.app/toolbox/eml/lib/matlab/lang/"
            "warning.m" /* pathName */
};

static emlrtMCInfo emlrtMCI = {
    84,                         /* lineNo */
    21,                         /* colNo */
    "WarningState/callWarning", /* fName */
    "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/"
    "WarningState.m" /* pName */
};

static emlrtRTEInfo emlrtRTEI = {
    13,     /* lineNo */
    9,      /* colNo */
    "sqrt", /* fName */
    "/Applications/MATLAB_R2024a.app/toolbox/eml/lib/matlab/elfun/sqrt.m" /* pName
                                                                           */
};

static emlrtRSInfo i_emlrtRSI = {
    84,                         /* lineNo */
    "WarningState/callWarning", /* fcnName */
    "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/"
    "WarningState.m" /* pathName */
};

/* Function Declarations */
static void feval(const emlrtStack *sp, const mxArray *m, const mxArray *m1,
                  const mxArray *m2, emlrtMCInfo *location);

/* Function Definitions */
static void feval(const emlrtStack *sp, const mxArray *m, const mxArray *m1,
                  const mxArray *m2, emlrtMCInfo *location)
{
  const mxArray *pArrays[3];
  pArrays[0] = m;
  pArrays[1] = m1;
  pArrays[2] = m2;
  emlrtCallMATLABR2012b((emlrtConstCTX)sp, 0, NULL, 3, &pArrays[0], "feval",
                        true, location);
}

real_T updateTrueAnomaly(const emlrtStack *sp, real_T a, real_T e, real_T i,
                         real_T Om, real_T w, real_T f0, real_T mu, real_T dt)
{
  static const int32_T iv[2] = {1, 7};
  static const int32_T iv1[2] = {1, 44};
  static const char_T varargin_1[44] = {
      'N', 'e', 'w', 't', 'o', 'n', ' ', 'm', 'e', 't', 'h', 'o', 'd', ' ', 'd',
      'i', 'd', ' ', 'n', 'o', 't', ' ', 'c', 'o', 'n', 'v', 'e', 'r', 'g', 'e',
      '.', ' ', 'R', 'e', 's', 'i', 'd', 'u', 'a', 'l', ':', ' ', '%', 'g'};
  static const char_T u[7] = {'w', 'a', 'r', 'n', 'i', 'n', 'g'};
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  const mxArray *b_y;
  const mxArray *m;
  const mxArray *y;
  real_T Ecur;
  real_T Eupdated;
  real_T f;
  real_T f0new;
  real_T n;
  int32_T k;
  (void)i;
  (void)Om;
  (void)w;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  /* UPDATETRUEANOMALY Summary of this function goes here */
  /*    Detailed explanation goes here */
  /*  Compute current eccentric anomaly */
  st.site = &emlrtRSI;
  Ecur = (1.0 - e) / (e + 1.0);
  if (Ecur < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  Ecur = muDoubleScalarSqrt(Ecur);
  Ecur = 2.0 * muDoubleScalarAtan(Ecur * muDoubleScalarTan(f0 / 2.0));
  /*  Compute Current Mean Anomaly from Kepler's Formula */
  /*  Update Mean Anomaly */
  st.site = &b_emlrtRSI;
  b_st.site = &f_emlrtRSI;
  st.site = &b_emlrtRSI;
  n = mu / muDoubleScalarPower(a, 3.0);
  if (n < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  n = muDoubleScalarSqrt(n);
  /*  Mean motion */
  Ecur = (Ecur - e * muDoubleScalarSin(Ecur)) + n * dt;
  /*  Compute Updated Eccentric Anomaly from Meupdated */
  /*  Numerical method: */
  /*  fun = @(E) E - e*sin(E) - Meupdated; */
  /*  Eupdated = fzero(fun, Ecur); % Ecur used as initial guess */
  /*  Use Newton-Raphson: Runs faster than fzero afaict */
  /*  Initial guess (Mupdated is good when e is small) */
  Eupdated = Ecur;
  k = 0;
  int32_T exitg1;
  do {
    exitg1 = 0;
    if (k < 15) {
      f = (Eupdated - e * muDoubleScalarSin(Eupdated)) - Ecur;
      n = -f / (1.0 - e * muDoubleScalarCos(Eupdated));
      Eupdated += n;
      /*  Convergence criteria */
      if (muDoubleScalarAbs(n) < 1.0E-5) {
        /*  Compute Updated True Anomaly from Eupdated */
        st.site = &c_emlrtRSI;
        b_st.site = &f_emlrtRSI;
        Ecur = e * e;
        st.site = &c_emlrtRSI;
        if (1.0 - Ecur < 0.0) {
          emlrtErrorWithMessageIdR2018a(
              &st, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
              "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
        }
        f0new = muDoubleScalarAtan2(muDoubleScalarSin(Eupdated) *
                                        muDoubleScalarSqrt(1.0 - Ecur),
                                    muDoubleScalarCos(Eupdated) - e);
        exitg1 = 1;
      } else {
        k++;
        if (*emlrtBreakCheckR2012bFlagVar != 0) {
          emlrtBreakCheckR2012b((emlrtConstCTX)sp);
        }
      }
    } else {
      /*  Warn if not converged, but still report result */
      st.site = &d_emlrtRSI;
      b_st.site = &f_emlrtRSI;
      Ecur = e * e;
      st.site = &d_emlrtRSI;
      if (1.0 - Ecur < 0.0) {
        emlrtErrorWithMessageIdR2018a(
            &st, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
            "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
      }
      f0new = muDoubleScalarAtan2(muDoubleScalarSin(Eupdated) *
                                      muDoubleScalarSqrt(1.0 - Ecur),
                                  muDoubleScalarCos(Eupdated) - e);
      st.site = &e_emlrtRSI;
      b_st.site = &h_emlrtRSI;
      y = NULL;
      m = emlrtCreateCharArray(2, &iv[0]);
      emlrtInitCharArrayR2013a(&b_st, 7, m, &u[0]);
      emlrtAssign(&y, m);
      b_y = NULL;
      m = emlrtCreateCharArray(2, &iv1[0]);
      emlrtInitCharArrayR2013a(&b_st, 44, m, &varargin_1[0]);
      emlrtAssign(&b_y, m);
      c_st.site = &i_emlrtRSI;
      feval(&c_st, y, b_y, emlrt_marshallOut(muDoubleScalarAbs(f)), &emlrtMCI);
      exitg1 = 1;
    }
  } while (exitg1 == 0);
  return f0new;
}

/* End of code generation (updateTrueAnomaly.c) */
