/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * propagateFromKeplerians.c
 *
 * Code generation for function 'propagateFromKeplerians'
 *
 */

/* Include files */
#include "propagateFromKeplerians.h"
#include "Keplerian2Cartesian.h"
#include "propagateFromKeplerians_data.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"
#include <string.h>

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = {
    13,                        /* lineNo */
    "propagateFromKeplerians", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "propagateFromKeplerians.m" /* pathName */
};

static emlrtRSInfo b_emlrtRSI = {
    17,                        /* lineNo */
    "propagateFromKeplerians", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "propagateFromKeplerians.m" /* pathName */
};

static emlrtRSInfo c_emlrtRSI = {
    18,                        /* lineNo */
    "propagateFromKeplerians", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "propagateFromKeplerians.m" /* pathName */
};

static emlrtRSInfo j_emlrtRSI = {
    9,                   /* lineNo */
    "updateTrueAnomaly", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "updateTrueAnomaly.m" /* pathName */
};

static emlrtRSInfo k_emlrtRSI = {
    17,                  /* lineNo */
    "updateTrueAnomaly", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "updateTrueAnomaly.m" /* pathName */
};

static emlrtRSInfo l_emlrtRSI = {
    41,                  /* lineNo */
    "updateTrueAnomaly", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "updateTrueAnomaly.m" /* pathName */
};

static emlrtRSInfo m_emlrtRSI = {
    47,                  /* lineNo */
    "updateTrueAnomaly", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "updateTrueAnomaly.m" /* pathName */
};

static emlrtRSInfo n_emlrtRSI = {
    48,                  /* lineNo */
    "updateTrueAnomaly", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "updateTrueAnomaly.m" /* pathName */
};

static emlrtRSInfo o_emlrtRSI = {
    52,                  /* lineNo */
    "updateTrueAnomaly", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "updateTrueAnomaly.m" /* pathName */
};

static emlrtRSInfo p_emlrtRSI = {
    53,                  /* lineNo */
    "updateTrueAnomaly", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "updateTrueAnomaly.m" /* pathName */
};

static emlrtRSInfo
    q_emlrtRSI =
        {
            44,       /* lineNo */
            "mpower", /* fcnName */
            "/Applications/MATLAB_R2024a.app/toolbox/eml/lib/matlab/matfun/"
            "mpower.m" /* pathName */
};

static emlrtRSInfo
    r_emlrtRSI =
        {
            53,        /* lineNo */
            "warning", /* fcnName */
            "/Applications/MATLAB_R2024a.app/toolbox/eml/lib/matlab/lang/"
            "warning.m" /* pathName */
};

static emlrtRSInfo s_emlrtRSI =
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

static emlrtBCInfo emlrtBCI = {
    -1,                        /* iFirst */
    -1,                        /* iLast */
    9,                         /* lineNo */
    23,                        /* colNo */
    "etR",                     /* aName */
    "propagateFromKeplerians", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "propagateFromKeplerians.m", /* pName */
    0                            /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = {
    -1,                        /* iFirst */
    -1,                        /* iLast */
    18,                        /* lineNo */
    14,                        /* colNo */
    "Xs",                      /* aName */
    "propagateFromKeplerians", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "propagateFromKeplerians.m", /* pName */
    0                            /* checkKind */
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

static emlrtRSInfo t_emlrtRSI = {
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

void propagateFromKeplerians(const emlrtStack *sp, const real_T Ki[6],
                             real_T mu, const real_T etR_data[],
                             const int32_T etR_size[2], real_T Xs_data[],
                             int32_T Xs_size[2])
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
  emlrtStack d_st;
  emlrtStack st;
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *m;
  const mxArray *y;
  real_T a;
  real_T dt;
  real_T dx;
  real_T e;
  real_T f;
  int32_T i;
  int32_T loop_ub;
  int32_T n;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  /* PROPAGATEFROMKEPLERIANS Returns Cartesian Trajectories from initial */
  /* orbital elements. All orbital elements are assumed constant except for the
   */
  /* true anomaly. */
  i = etR_size[1];
  if (etR_size[1] < 1) {
    emlrtDynamicBoundsCheckR2012b(1, 1, etR_size[1], &emlrtBCI,
                                  (emlrtConstCTX)sp);
  }
  dt = (etR_data[etR_size[1] - 1] - etR_data[0]) / ((real_T)etR_size[1] - 1.0);
  /* UNPACKKEPLERIAN Unpacks the 6-vector when called */
  a = Ki[0];
  e = Ki[1];
  Xs_size[0] = 6;
  Xs_size[1] = etR_size[1];
  loop_ub = 6 * etR_size[1];
  if (loop_ub - 1 >= 0) {
    memset(&Xs_data[0], 0, (uint32_T)loop_ub * sizeof(real_T));
  }
  /*  preallocating */
  st.site = &emlrtRSI;
  Keplerian2Cartesian(&st, Ki[0], Ki[1], Ki[2], Ki[3], Ki[4], Ki[5], mu,
                      &Xs_data[0]);
  dx = Ki[5];
  for (n = 0; n <= i - 2; n++) {
    real_T Ecur;
    real_T Meupdated;
    real_T b_n;
    st.site = &b_emlrtRSI;
    /* UPDATETRUEANOMALY Summary of this function goes here */
    /*    Detailed explanation goes here */
    /*  Compute current eccentric anomaly */
    b_st.site = &j_emlrtRSI;
    b_n = (1.0 - e) / (e + 1.0);
    if (b_n < 0.0) {
      emlrtErrorWithMessageIdR2018a(
          &b_st, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
          "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
    }
    b_n = muDoubleScalarSqrt(b_n);
    Ecur = 2.0 * muDoubleScalarAtan(b_n * muDoubleScalarTan(dx / 2.0));
    /*  Compute Current Mean Anomaly from Kepler's Formula */
    /*  Update Mean Anomaly */
    b_st.site = &k_emlrtRSI;
    c_st.site = &q_emlrtRSI;
    b_st.site = &k_emlrtRSI;
    b_n = mu / muDoubleScalarPower(a, 3.0);
    if (b_n < 0.0) {
      emlrtErrorWithMessageIdR2018a(
          &b_st, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
          "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
    }
    b_n = muDoubleScalarSqrt(b_n);
    /*  Mean motion */
    dx = Ecur - e * muDoubleScalarSin(Ecur);
    Meupdated = dx + b_n * dt;
    /*  Compute Updated Eccentric Anomaly from Meupdated */
    /*  Numerical method: */
    if (e < 0.9) {
      real_T Eupdated;
      /*  Use Newton-Raphson: Runs faster than fzero afaict */
      /*  Initial guess (Mupdated is good when e is small) */
      Eupdated = Meupdated;
      loop_ub = 0;
      int32_T exitg1;
      do {
        exitg1 = 0;
        if (loop_ub < 20) {
          f = (Eupdated - e * muDoubleScalarSin(Eupdated)) - Meupdated;
          dx = -f / (1.0 - e * muDoubleScalarCos(Eupdated));
          Eupdated += dx;
          /*  Convergence criteria */
          if (muDoubleScalarAbs(dx) < 1.0E-5) {
            real_T c;
            /*  Compute Updated True Anomaly from Eupdated */
            b_st.site = &l_emlrtRSI;
            c_st.site = &q_emlrtRSI;
            c = e * e;
            b_st.site = &l_emlrtRSI;
            if (1.0 - c < 0.0) {
              emlrtErrorWithMessageIdR2018a(
                  &b_st, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
                  "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
            }
            dx = muDoubleScalarAtan2(muDoubleScalarSin(Eupdated) *
                                         muDoubleScalarSqrt(1.0 - c),
                                     muDoubleScalarCos(Eupdated) - e);
            exitg1 = 1;
          } else {
            loop_ub++;
            if (*emlrtBreakCheckR2012bFlagVar != 0) {
              emlrtBreakCheckR2012b(&st);
            }
          }
        } else {
          real_T c;
          /*  Warn if not converged, but still report result */
          b_st.site = &m_emlrtRSI;
          c_st.site = &q_emlrtRSI;
          c = e * e;
          b_st.site = &m_emlrtRSI;
          if (1.0 - c < 0.0) {
            emlrtErrorWithMessageIdR2018a(
                &b_st, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
                "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
          }
          dx = muDoubleScalarAtan2(muDoubleScalarSin(Eupdated) *
                                       muDoubleScalarSqrt(1.0 - c),
                                   muDoubleScalarCos(Eupdated) - e);
          b_st.site = &n_emlrtRSI;
          c_st.site = &r_emlrtRSI;
          y = NULL;
          m = emlrtCreateCharArray(2, &iv[0]);
          emlrtInitCharArrayR2013a(&c_st, 7, m, &u[0]);
          emlrtAssign(&y, m);
          b_y = NULL;
          m = emlrtCreateCharArray(2, &iv1[0]);
          emlrtInitCharArrayR2013a(&c_st, 44, m, &varargin_1[0]);
          emlrtAssign(&b_y, m);
          c_y = NULL;
          m = emlrtCreateDoubleScalar(muDoubleScalarAbs(f));
          emlrtAssign(&c_y, m);
          d_st.site = &t_emlrtRSI;
          feval(&d_st, y, b_y, c_y, &emlrtMCI);
          exitg1 = 1;
        }
      } while (exitg1 == 0);
    } else {
      real_T Eupdated;
      real_T c;
      /*  fzero is more robust for large eccentricities */
      b_st.site = &o_emlrtRSI;
      if (muDoubleScalarIsInf(Ecur) || muDoubleScalarIsNaN(Ecur)) {
        emlrtErrorWithMessageIdR2018a(&b_st, &b_emlrtRTEI,
                                      "MATLAB:optimfun:fzero:Arg2NotFinite",
                                      "MATLAB:optimfun:fzero:Arg2NotFinite", 0);
      }
      b_n = dx - Meupdated;
      if (b_n == 0.0) {
        Eupdated = Ecur;
      } else {
        real_T b_a;
        real_T fb;
        int32_T exitg1;
        if (muDoubleScalarIsInf(b_n) || muDoubleScalarIsNaN(b_n)) {
          emlrtErrorWithMessageIdR2018a(
              &b_st, &c_emlrtRTEI,
              "MATLAB:optimfun:fzero:ValueAtInitGuessComplexOrNotFinite",
              "MATLAB:optimfun:fzero:ValueAtInitGuessComplexOrNotFinite", 0);
        }
        if (Ecur != 0.0) {
          dx = Ecur / 50.0;
        } else {
          dx = 0.02;
        }
        c_st.site = &s_emlrtRSI;
        b_a = Ecur;
        Eupdated = Ecur;
        fb = b_n;
        do {
          exitg1 = 0;
          loop_ub = (fb > 0.0);
          if ((b_n > 0.0) == loop_ub) {
            dx *= 1.4142135623730951;
            b_a = Ecur - dx;
            b_n = (b_a - e * muDoubleScalarSin(b_a)) - Meupdated;
            if (muDoubleScalarIsInf(b_n) || muDoubleScalarIsNaN(b_n)) {
              Eupdated = rtNaN;
              exitg1 = 1;
            } else if (muDoubleScalarIsInf(b_a) || muDoubleScalarIsNaN(b_a)) {
              Eupdated = rtNaN;
              exitg1 = 1;
            } else if ((b_n > 0.0) != loop_ub) {
              exitg1 = 2;
            } else {
              Eupdated = Ecur + dx;
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
              fc = b_n;
              d = Eupdated - b_a;
              b_e = d;
            }
            if (muDoubleScalarAbs(fc) < muDoubleScalarAbs(fb)) {
              b_a = Eupdated;
              Eupdated = c;
              c = b_a;
              b_n = fb;
              fb = fc;
              fc = b_n;
            }
            b_m = 0.5 * (c - Eupdated);
            toler = 4.4408920985006262E-16 *
                    muDoubleScalarMax(muDoubleScalarAbs(Eupdated), 1.0);
            if ((muDoubleScalarAbs(b_m) <= toler) || (fb == 0.0)) {
              exitg2 = true;
            } else {
              if ((muDoubleScalarAbs(b_e) < toler) ||
                  (muDoubleScalarAbs(b_n) <= muDoubleScalarAbs(fb))) {
                d = b_m;
                b_e = b_m;
              } else {
                real_T s;
                s = fb / b_n;
                if (b_a == c) {
                  b_n = 2.0 * b_m * s;
                  dx = 1.0 - s;
                } else {
                  dx = b_n / fc;
                  Ecur = fb / fc;
                  b_n = s * (2.0 * b_m * dx * (dx - Ecur) -
                             (Eupdated - b_a) * (Ecur - 1.0));
                  dx = (dx - 1.0) * (Ecur - 1.0) * (s - 1.0);
                }
                if (b_n > 0.0) {
                  dx = -dx;
                } else {
                  b_n = -b_n;
                }
                if ((2.0 * b_n <
                     3.0 * b_m * dx - muDoubleScalarAbs(toler * dx)) &&
                    (b_n < muDoubleScalarAbs(0.5 * b_e * dx))) {
                  b_e = d;
                  d = b_n / dx;
                } else {
                  d = b_m;
                  b_e = b_m;
                }
              }
              b_a = Eupdated;
              b_n = fb;
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
      b_st.site = &p_emlrtRSI;
      c_st.site = &q_emlrtRSI;
      c = e * e;
      b_st.site = &p_emlrtRSI;
      if (1.0 - c < 0.0) {
        emlrtErrorWithMessageIdR2018a(
            &b_st, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
            "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
      }
      dx = muDoubleScalarAtan2(muDoubleScalarSin(Eupdated) *
                                   muDoubleScalarSqrt(1.0 - c),
                               muDoubleScalarCos(Eupdated) - e);
    }
    if (n + 2 > Xs_size[1]) {
      emlrtDynamicBoundsCheckR2012b(n + 2, 1, Xs_size[1], &b_emlrtBCI,
                                    (emlrtConstCTX)sp);
    }
    st.site = &c_emlrtRSI;
    Keplerian2Cartesian(&st, a, e, Ki[2], Ki[3], Ki[4], dx, mu,
                        &Xs_data[6 * (n + 1)]);
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b((emlrtConstCTX)sp);
    }
  }
}

/* End of code generation (propagateFromKeplerians.c) */
