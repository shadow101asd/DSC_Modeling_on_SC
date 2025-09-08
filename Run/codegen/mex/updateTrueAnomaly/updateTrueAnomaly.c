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
    41,                  /* lineNo */
    "updateTrueAnomaly", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "updateTrueAnomaly.m" /* pathName */
};

static emlrtRSInfo d_emlrtRSI = {
    47,                  /* lineNo */
    "updateTrueAnomaly", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "updateTrueAnomaly.m" /* pathName */
};

static emlrtRSInfo e_emlrtRSI = {
    48,                  /* lineNo */
    "updateTrueAnomaly", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "updateTrueAnomaly.m" /* pathName */
};

static emlrtRSInfo f_emlrtRSI = {
    52,                  /* lineNo */
    "updateTrueAnomaly", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "updateTrueAnomaly.m" /* pathName */
};

static emlrtRSInfo g_emlrtRSI = {
    53,                  /* lineNo */
    "updateTrueAnomaly", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "updateTrueAnomaly.m" /* pathName */
};

static emlrtRSInfo
    h_emlrtRSI =
        {
            44,       /* lineNo */
            "mpower", /* fcnName */
            "/Applications/MATLAB_R2024a.app/toolbox/eml/lib/matlab/matfun/"
            "mpower.m" /* pathName */
};

static emlrtRSInfo
    j_emlrtRSI =
        {
            53,        /* lineNo */
            "warning", /* fcnName */
            "/Applications/MATLAB_R2024a.app/toolbox/eml/lib/matlab/lang/"
            "warning.m" /* pathName */
};

static emlrtRSInfo k_emlrtRSI =
    {
        91,      /* lineNo */
        "fzero", /* fcnName */
        "/Applications/MATLAB_R2024a.app/toolbox/eml/lib/matlab/optimfun/"
        "fzero.m" /* pathName */
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

static emlrtRTEInfo b_emlrtRTEI =
    {
        72,      /* lineNo */
        19,      /* colNo */
        "fzero", /* fName */
        "/Applications/MATLAB_R2024a.app/toolbox/eml/lib/matlab/optimfun/"
        "fzero.m" /* pName */
};

static emlrtRTEInfo c_emlrtRTEI =
    {
        83,      /* lineNo */
        9,       /* colNo */
        "fzero", /* fName */
        "/Applications/MATLAB_R2024a.app/toolbox/eml/lib/matlab/optimfun/"
        "fzero.m" /* pName */
};

static emlrtRSInfo l_emlrtRSI = {
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
  real_T Meupdated;
  real_T f0new;
  real_T fx;
  real_T n;
  real_T r;
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
  n = (1.0 - e) / (e + 1.0);
  if (n < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  n = muDoubleScalarSqrt(n);
  Ecur = 2.0 * muDoubleScalarAtan(n * muDoubleScalarTan(f0 / 2.0));
  /*  Compute Current Mean Anomaly from Kepler's Formula */
  /*  Update Mean Anomaly */
  st.site = &b_emlrtRSI;
  b_st.site = &h_emlrtRSI;
  st.site = &b_emlrtRSI;
  n = mu / muDoubleScalarPower(a, 3.0);
  if (n < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  n = muDoubleScalarSqrt(n);
  /*  Mean motion */
  fx = Ecur - e * muDoubleScalarSin(Ecur);
  Meupdated = fx + n * dt;
  /*  Compute Updated Eccentric Anomaly from Meupdated */
  /*  Numerical method: */
  if (e < 0.9) {
    real_T Eupdated;
    int32_T k;
    /*  Use Newton-Raphson: Runs faster than fzero afaict */
    /*  Initial guess (Mupdated is good when e is small) */
    Eupdated = Meupdated;
    k = 0;
    int32_T exitg1;
    do {
      exitg1 = 0;
      if (k < 20) {
        r = (Eupdated - e * muDoubleScalarSin(Eupdated)) - Meupdated;
        n = -r / (1.0 - e * muDoubleScalarCos(Eupdated));
        Eupdated += n;
        /*  Convergence criteria */
        if (muDoubleScalarAbs(n) < 1.0E-5) {
          real_T c;
          /*  Compute Updated True Anomaly from Eupdated */
          st.site = &c_emlrtRSI;
          b_st.site = &h_emlrtRSI;
          c = e * e;
          st.site = &c_emlrtRSI;
          if (1.0 - c < 0.0) {
            emlrtErrorWithMessageIdR2018a(
                &st, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
                "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
          }
          f0new = muDoubleScalarAtan2(muDoubleScalarSin(Eupdated) *
                                          muDoubleScalarSqrt(1.0 - c),
                                      muDoubleScalarCos(Eupdated) - e);
          exitg1 = 1;
        } else {
          k++;
          if (*emlrtBreakCheckR2012bFlagVar != 0) {
            emlrtBreakCheckR2012b((emlrtConstCTX)sp);
          }
        }
      } else {
        real_T c;
        /*  Warn if not converged, but still report result */
        st.site = &d_emlrtRSI;
        b_st.site = &h_emlrtRSI;
        c = e * e;
        st.site = &d_emlrtRSI;
        if (1.0 - c < 0.0) {
          emlrtErrorWithMessageIdR2018a(
              &st, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
              "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
        }
        f0new = muDoubleScalarAtan2(muDoubleScalarSin(Eupdated) *
                                        muDoubleScalarSqrt(1.0 - c),
                                    muDoubleScalarCos(Eupdated) - e);
        st.site = &e_emlrtRSI;
        b_st.site = &j_emlrtRSI;
        y = NULL;
        m = emlrtCreateCharArray(2, &iv[0]);
        emlrtInitCharArrayR2013a(&b_st, 7, m, &u[0]);
        emlrtAssign(&y, m);
        b_y = NULL;
        m = emlrtCreateCharArray(2, &iv1[0]);
        emlrtInitCharArrayR2013a(&b_st, 44, m, &varargin_1[0]);
        emlrtAssign(&b_y, m);
        c_st.site = &l_emlrtRSI;
        feval(&c_st, y, b_y, emlrt_marshallOut(muDoubleScalarAbs(r)),
              &emlrtMCI);
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  } else {
    real_T Eupdated;
    real_T c;
    /*  fzero is more robust for large eccentricities */
    st.site = &f_emlrtRSI;
    if (muDoubleScalarIsInf(Ecur) || muDoubleScalarIsNaN(Ecur)) {
      emlrtErrorWithMessageIdR2018a(&st, &b_emlrtRTEI,
                                    "MATLAB:optimfun:fzero:Arg2NotFinite",
                                    "MATLAB:optimfun:fzero:Arg2NotFinite", 0);
    }
    fx -= Meupdated;
    if (fx == 0.0) {
      Eupdated = Ecur;
    } else {
      real_T b_a;
      real_T fb;
      int32_T exitg1;
      if (muDoubleScalarIsInf(fx) || muDoubleScalarIsNaN(fx)) {
        emlrtErrorWithMessageIdR2018a(
            &st, &c_emlrtRTEI,
            "MATLAB:optimfun:fzero:ValueAtInitGuessComplexOrNotFinite",
            "MATLAB:optimfun:fzero:ValueAtInitGuessComplexOrNotFinite", 0);
      }
      if (Ecur != 0.0) {
        n = Ecur / 50.0;
      } else {
        n = 0.02;
      }
      b_st.site = &k_emlrtRSI;
      b_a = Ecur;
      Eupdated = Ecur;
      fb = fx;
      int32_T k;
      do {
        exitg1 = 0;
        k = (fb > 0.0);
        if ((fx > 0.0) == k) {
          n *= 1.4142135623730951;
          b_a = Ecur - n;
          fx = (b_a - e * muDoubleScalarSin(b_a)) - Meupdated;
          if (muDoubleScalarIsInf(fx) || muDoubleScalarIsNaN(fx)) {
            Eupdated = rtNaN;
            exitg1 = 1;
          } else if (muDoubleScalarIsInf(b_a) || muDoubleScalarIsNaN(b_a)) {
            Eupdated = rtNaN;
            exitg1 = 1;
          } else if ((fx > 0.0) != k) {
            exitg1 = 2;
          } else {
            Eupdated = Ecur + n;
            fb = (Eupdated - e * muDoubleScalarSin(Eupdated)) - Meupdated;
            if (muDoubleScalarIsInf(fb) || muDoubleScalarIsNaN(fb)) {
              Eupdated = rtNaN;
              exitg1 = 1;
            } else if (muDoubleScalarIsInf(Eupdated) ||
                       muDoubleScalarIsNaN(Eupdated)) {
              Eupdated = rtNaN;
              exitg1 = 1;
            }
          }
        } else {
          exitg1 = 2;
        }
      } while (exitg1 == 0);
      if (exitg1 != 1) {
        real_T b_e;
        real_T d;
        real_T fc;
        boolean_T exitg2;
        fc = fb;
        c = Eupdated;
        b_e = 0.0;
        d = 0.0;
        exitg2 = false;
        while ((!exitg2) && ((fb != 0.0) && (b_a != Eupdated))) {
          real_T b_m;
          real_T toler;
          if ((fb > 0.0) == (fc > 0.0)) {
            c = b_a;
            fc = fx;
            d = Eupdated - b_a;
            b_e = d;
          }
          if (muDoubleScalarAbs(fc) < muDoubleScalarAbs(fb)) {
            b_a = Eupdated;
            Eupdated = c;
            c = b_a;
            fx = fb;
            fb = fc;
            fc = fx;
          }
          b_m = 0.5 * (c - Eupdated);
          toler = 4.4408920985006262E-16 *
                  muDoubleScalarMax(muDoubleScalarAbs(Eupdated), 1.0);
          if ((muDoubleScalarAbs(b_m) <= toler) || (fb == 0.0)) {
            exitg2 = true;
          } else {
            if ((muDoubleScalarAbs(b_e) < toler) ||
                (muDoubleScalarAbs(fx) <= muDoubleScalarAbs(fb))) {
              d = b_m;
              b_e = b_m;
            } else {
              real_T s;
              s = fb / fx;
              if (b_a == c) {
                n = 2.0 * b_m * s;
                Ecur = 1.0 - s;
              } else {
                Ecur = fx / fc;
                r = fb / fc;
                n = s * (2.0 * b_m * Ecur * (Ecur - r) -
                         (Eupdated - b_a) * (r - 1.0));
                Ecur = (Ecur - 1.0) * (r - 1.0) * (s - 1.0);
              }
              if (n > 0.0) {
                Ecur = -Ecur;
              } else {
                n = -n;
              }
              if ((2.0 * n <
                   3.0 * b_m * Ecur - muDoubleScalarAbs(toler * Ecur)) &&
                  (n < muDoubleScalarAbs(0.5 * b_e * Ecur))) {
                b_e = d;
                d = n / Ecur;
              } else {
                d = b_m;
                b_e = b_m;
              }
            }
            b_a = Eupdated;
            fx = fb;
            if (muDoubleScalarAbs(d) > toler) {
              Eupdated += d;
            } else if (Eupdated > c) {
              Eupdated -= toler;
            } else {
              Eupdated += toler;
            }
            fb = (Eupdated - e * muDoubleScalarSin(Eupdated)) - Meupdated;
          }
        }
      }
    }
    /*  Ecur used as initial guess */
    st.site = &g_emlrtRSI;
    b_st.site = &h_emlrtRSI;
    c = e * e;
    st.site = &g_emlrtRSI;
    if (1.0 - c < 0.0) {
      emlrtErrorWithMessageIdR2018a(
          &st, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
          "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
    }
    f0new = muDoubleScalarAtan2(muDoubleScalarSin(Eupdated) *
                                    muDoubleScalarSqrt(1.0 - c),
                                muDoubleScalarCos(Eupdated) - e);
  }
  return f0new;
}

/* End of code generation (updateTrueAnomaly.c) */
