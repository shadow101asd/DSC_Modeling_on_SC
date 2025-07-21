/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Keplerian2Cartesian.c
 *
 * Code generation for function 'Keplerian2Cartesian'
 *
 */

/* Include files */
#include "Keplerian2Cartesian.h"
#include "propagateFromKeplerians_data.h"
#include "rt_nonfinite.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo h_emlrtRSI = {
    45,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

/* Function Definitions */
void Keplerian2Cartesian(const emlrtStack *sp, real_T a, real_T e, real_T i,
                         real_T Om, real_T w, real_T f0, real_T mu, real_T y[6])
{
  emlrtStack st;
  real_T IJKOverPQW1;
  real_T IJKOverPQW1_tmp;
  real_T IJKOverPQW1_tmp_tmp;
  real_T IJKOverPQW2;
  real_T IJKOverPQW3;
  real_T IJKOverPQW3_tmp;
  real_T IJKOverPQW4;
  real_T IJKOverPQW5;
  real_T IJKOverPQW6;
  real_T IJKOverPQW7;
  real_T cosnu;
  real_T p;
  real_T rPQW_idx_0;
  real_T rPQW_idx_1;
  real_T rPQW_temp_idx_0;
  real_T rootmup;
  real_T sinnu;
  real_T vPQW_idx_0;
  real_T vPQW_idx_1;
  st.prev = sp;
  st.tls = sp->tls;
  /*  Keplerian2Cartesian - Convert Keplerian elements to Cartesian state
   * vectors */
  /*  */
  /*  Inputs (all 1×n): */
  /*    a   - semi-major axis */
  /*    e   - eccentricity or [3×n] eccentricity vectors */
  /*    i   - inclination [rad] */
  /*    Om  - RAAN [rad] */
  /*    w   - argument of periapsis [rad] */
  /*    f0  - true anomaly [rad] */
  /*    mu  - gravitational parameter (scalar or 1×n) */
  /*  */
  /*  Output: */
  /*    y   - 6×n matrix [r; v] in inertial frame */
  /*  Handle eccentricity vector case */
  p = a * (1.0 - e * e);
  /*  Special Cases */
  if (muDoubleScalarAbs(e) < 1.0E-6) {
    /*  if circular */
    w *= 0.0;
    if (muDoubleScalarAbs(i) < 1.0E-6) {
      /*  if equatorial */
      /* f0 = acos(r(1, :)./R); % lambda_true */
      Om *= 0.0;
    } else {
      /* f0 = acos(dot(N,r)./(Nnorm.*R)); % u */
    }
  } else if (muDoubleScalarAbs(i) < 1.0E-6) {
    /*  if equatorial */
    /* f0 = acos(r(1, :)./R); % lambda_true */
    Om *= 0.0;
  }
  /*  Storing Variables for Computational Efficiency */
  cosnu = muDoubleScalarCos(f0);
  sinnu = muDoubleScalarSin(f0);
  st.site = &h_emlrtRSI;
  rootmup = mu * (1.0 / p);
  if (rootmup < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  rootmup = muDoubleScalarSqrt(rootmup);
  rPQW_temp_idx_0 = e * cosnu + 1.0;
  rPQW_idx_0 = p * cosnu / rPQW_temp_idx_0;
  rPQW_idx_1 = p * sinnu / rPQW_temp_idx_0;
  vPQW_idx_0 = -rootmup * sinnu;
  vPQW_idx_1 = rootmup * (e + cosnu);
  /*  Computing IJK/PQW */
  rootmup = muDoubleScalarSin(Om);
  cosnu = muDoubleScalarCos(w);
  IJKOverPQW1_tmp_tmp = muDoubleScalarCos(Om);
  rPQW_temp_idx_0 = muDoubleScalarSin(w);
  IJKOverPQW1_tmp = muDoubleScalarCos(i);
  sinnu = IJKOverPQW1_tmp_tmp * cosnu;
  IJKOverPQW1 = sinnu - rootmup * rPQW_temp_idx_0 * IJKOverPQW1_tmp;
  p = rootmup * cosnu;
  IJKOverPQW2 = p + IJKOverPQW1_tmp_tmp * rPQW_temp_idx_0 * IJKOverPQW1_tmp;
  IJKOverPQW3_tmp = muDoubleScalarSin(i);
  IJKOverPQW3 = rPQW_temp_idx_0 * IJKOverPQW3_tmp;
  IJKOverPQW4 = -IJKOverPQW1_tmp_tmp * rPQW_temp_idx_0 - p * IJKOverPQW1_tmp;
  IJKOverPQW5 = -rootmup * rPQW_temp_idx_0 + sinnu * IJKOverPQW1_tmp;
  IJKOverPQW6 = cosnu * IJKOverPQW3_tmp;
  IJKOverPQW7 = rootmup * IJKOverPQW3_tmp;
  p = -IJKOverPQW1_tmp_tmp * IJKOverPQW3_tmp;
  /*  for j = 1:n */
  /*      IJKOverPQW = [IJKOverPQW1(1, j), IJKOverPQW4(1, j), IJKOverPQW7(1, j);
   */
  /*                      IJKOverPQW2(1, j), IJKOverPQW5(1, j), IJKOverPQW8(1,
   * j);  */
  /*                      IJKOverPQW3(1, j), IJKOverPQW6(1, j), IJKOverPQW9(1,
   * j)]; */
  /*      rPQW(:,j) = IJKOverPQW * rPQW(:,j); */
  /*      vPQW(:,j) = IJKOverPQW * vPQW(:,j); */
  /*  end */
  /*  Vectorized rotation of rPQW and vPQW into IJK */
  rPQW_temp_idx_0 = rPQW_idx_0;
  sinnu = vPQW_idx_0;
  cosnu = rPQW_idx_1;
  rootmup = vPQW_idx_1;
  rPQW_idx_0 =
      (IJKOverPQW1 * rPQW_idx_0 + IJKOverPQW4 * rPQW_idx_1) + IJKOverPQW7 * 0.0;
  rPQW_idx_1 =
      (IJKOverPQW2 * rPQW_temp_idx_0 + IJKOverPQW5 * rPQW_idx_1) + p * 0.0;
  vPQW_idx_0 =
      (IJKOverPQW1 * vPQW_idx_0 + IJKOverPQW4 * vPQW_idx_1) + IJKOverPQW7 * 0.0;
  vPQW_idx_1 = (IJKOverPQW2 * sinnu + IJKOverPQW5 * vPQW_idx_1) + p * 0.0;
  y[0] = rPQW_idx_0;
  y[3] = vPQW_idx_0;
  y[1] = rPQW_idx_1;
  y[4] = vPQW_idx_1;
  y[2] = (IJKOverPQW3 * rPQW_temp_idx_0 + IJKOverPQW6 * cosnu) +
         IJKOverPQW1_tmp * 0.0;
  y[5] = (IJKOverPQW3 * sinnu + IJKOverPQW6 * rootmup) + IJKOverPQW1_tmp * 0.0;
}

/* End of code generation (Keplerian2Cartesian.c) */
