/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * lambert.c
 *
 * Code generation for function 'lambert'
 *
 */

/* Include files */
#include "lambert.h"
#include "acos.h"
#include "acosh.h"
#include "asin.h"
#include "asinh.h"
#include "lambert_data.h"
#include "mod.h"
#include "power.h"
#include "rt_nonfinite.h"
#include "sqrt.h"
#include "mwmathutil.h"
#include <math.h>
#include <string.h>

/* Variable Definitions */
static real_T an[25];

static emlrtRSInfo emlrtRSI = {
    201,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo b_emlrtRSI = {
    202,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo c_emlrtRSI = {
    206,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo d_emlrtRSI = {
    208,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo e_emlrtRSI = {
    217,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo f_emlrtRSI = {
    220,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo g_emlrtRSI = {
    224,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo h_emlrtRSI = {
    231,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo i_emlrtRSI = {
    239,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo j_emlrtRSI = {
    240,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo k_emlrtRSI = {
    259,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo l_emlrtRSI = {
    261,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo m_emlrtRSI = {
    264,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo n_emlrtRSI = {
    268,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo o_emlrtRSI = {
    269,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo p_emlrtRSI = {
    290,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo q_emlrtRSI = {
    292,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo r_emlrtRSI = {
    294,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo s_emlrtRSI = {
    296,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo t_emlrtRSI = {
    297,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo u_emlrtRSI = {
    301,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo v_emlrtRSI = {
    303,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo w_emlrtRSI = {
    306,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo x_emlrtRSI = {
    321,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo y_emlrtRSI = {
    339,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo ab_emlrtRSI = {
    343,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo bb_emlrtRSI = {
    345,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo cb_emlrtRSI = {
    347,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo db_emlrtRSI = {
    348,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo eb_emlrtRSI = {
    350,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo fb_emlrtRSI = {
    351,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo gb_emlrtRSI = {
    353,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo hb_emlrtRSI = {
    354,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo ib_emlrtRSI = {
    374,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo jb_emlrtRSI = {
    375,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo kb_emlrtRSI = {
    391,       /* lineNo */
    "lambert", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo
    lb_emlrtRSI =
        {
            44,       /* lineNo */
            "mpower", /* fcnName */
            "/Applications/MATLAB_R2024a.app/toolbox/eml/lib/matlab/matfun/"
            "mpower.m" /* pathName */
};

static emlrtRSInfo mb_emlrtRSI = {
    71,      /* lineNo */
    "power", /* fcnName */
    "/Applications/MATLAB_R2024a.app/toolbox/eml/lib/matlab/ops/power.m" /* pathName
                                                                          */
};

static emlrtRSInfo rb_emlrtRSI = {
    430,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo sb_emlrtRSI = {
    431,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo tb_emlrtRSI = {
    435,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo ub_emlrtRSI = {
    439,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo vb_emlrtRSI = {
    450,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo wb_emlrtRSI = {
    452,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo xb_emlrtRSI = {
    453,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo yb_emlrtRSI = {
    459,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo ac_emlrtRSI = {
    461,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo bc_emlrtRSI = {
    473,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo cc_emlrtRSI = {
    474,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo dc_emlrtRSI = {
    478,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo ec_emlrtRSI = {
    480,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo fc_emlrtRSI = {
    493,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo gc_emlrtRSI = {
    495,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo hc_emlrtRSI = {
    507,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo ic_emlrtRSI = {
    510,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo jc_emlrtRSI = {
    512,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo kc_emlrtRSI = {
    522,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo lc_emlrtRSI = {
    532,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo mc_emlrtRSI = {
    536,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo nc_emlrtRSI = {
    538,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo oc_emlrtRSI = {
    540,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo pc_emlrtRSI = {
    548,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo qc_emlrtRSI = {
    551,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo rc_emlrtRSI = {
    555,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo sc_emlrtRSI = {
    559,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo tc_emlrtRSI = {
    578,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo uc_emlrtRSI = {
    584,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo vc_emlrtRSI = {
    586,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo wc_emlrtRSI = {
    595,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo xc_emlrtRSI = {
    601,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo yc_emlrtRSI = {
    603,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo ad_emlrtRSI = {
    626,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo bd_emlrtRSI = {
    627,                          /* lineNo */
    "lambert_LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo cd_emlrtRSI = {
    671,                  /* lineNo */
    "LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo qd_emlrtRSI = {
    670,                  /* lineNo */
    "LancasterBlanchard", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRSInfo td_emlrtRSI = {
    760,                /* lineNo */
    "minmax_distances", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/lambert.m" /* pathName
                                                                          */
};

static emlrtRTEInfo b_emlrtRTEI = {
    82,         /* lineNo */
    5,          /* colNo */
    "fltpower", /* fName */
    "/Applications/MATLAB_R2024a.app/toolbox/eml/lib/matlab/ops/power.m" /* pName
                                                                          */
};

static emlrtRTEInfo e_emlrtRTEI = {
    14,    /* lineNo */
    9,     /* colNo */
    "log", /* fName */
    "/Applications/MATLAB_R2024a.app/toolbox/eml/lib/matlab/elfun/log.m" /* pName
                                                                          */
};

static const int16_T iv[25] = {0,   2,   6,   12,  20,  30,  42,  56,  72,
                               90,  110, 132, 156, 182, 210, 240, 272, 306,
                               342, 380, 420, 462, 506, 552, 600};

/* Function Declarations */
static real_T LancasterBlanchard(const emlrtStack *sp, real_T q, real_T m);

static real_T b_LancasterBlanchard(const emlrtStack *sp, real_T x, real_T q,
                                   real_T m, real_T *Tp, real_T *Tpp,
                                   real_T *Tppp);

static real_T c_LancasterBlanchard(const emlrtStack *sp, real_T x, real_T q,
                                   real_T m);

static real_T d_LancasterBlanchard(const emlrtStack *sp, real_T x, real_T q,
                                   real_T m, real_T *Tp, real_T *Tpp);

static real_T lambert_LancasterBlanchard(const emlrtStack *sp,
                                         const real_T r1vec[3],
                                         const real_T r2vec[3], real_T tf,
                                         real_T m, real_T muC, real_T V1[3],
                                         real_T V2[3],
                                         real_T extremal_distances[2]);

static void minmax_distances(const emlrtStack *sp, const real_T r1vec[3],
                             real_T r1, const real_T r2vec[3], real_T r2,
                             real_T dth, real_T a, const real_T V1[3],
                             const real_T V2[3], real_T m, real_T muC,
                             real_T extremal_distances[2]);

/* Function Definitions */
static real_T LancasterBlanchard(const emlrtStack *sp, real_T q, real_T m)
{
  emlrtStack st;
  real_T z;
  st.prev = sp;
  st.tls = sp->tls;
  /*  Lancaster & Blanchard's function, and three derivatives thereof */
  /*  protection against idiotic input */
  /*  compute parameter E */
  /*  T(x), T'(x), T''(x) */
  /*  all other cases */
  /*  compute all substitution functions */
  z = q * q;
  st.site = &cd_emlrtRSI;
  if (-z + 1.0 < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  z = muDoubleScalarSqrt(-z + 1.0);
  /*  BUGFIX: (Simon Tardivel) this line is incorrect for E==0 and f+g==0 */
  /*  d  = (E < 0)*(atan2(f, g) + pi*m) + (E > 0)*log( max(0, f + g) ); */
  /*  it should be written out like so: */
  /*  T(x) */
  return -(2.0 *
           ((0.0 - q * z) - (muDoubleScalarAtan2(z - q * 0.0, 0.0 * z - (-q)) +
                             3.1415926535897931 * m)));
  /*   T'(x) */
  /*  T''(x) */
  /*  T'''(x) */
}

static real_T b_LancasterBlanchard(const emlrtStack *sp, real_T x, real_T q,
                                   real_T m, real_T *Tp, real_T *Tpp,
                                   real_T *Tppp)
{
  static const int16_T b_iv[25] = {0,    0,     6,     24,   60,   120,  210,
                                   336,  504,   720,   990,  1320, 1716, 2184,
                                   2730, 3360,  4080,  4896, 5814, 6840, 7980,
                                   9240, 10626, 12144, 13800};
  emlrtStack st;
  real_T E_tmp;
  real_T T;
  int32_T i;
  st.prev = sp;
  st.tls = sp->tls;
  /*  Lancaster & Blanchard's function, and three derivatives thereof */
  /*  protection against idiotic input */
  if (x < -1.0) {
    /*  impossible; negative eccentricity */
    x = muDoubleScalarAbs(x) - 2.0;
  } else if (x == -1.0) {
    /*  impossible; offset x slightly */
    x = -0.99999999999999978;
  }
  /*  compute parameter E */
  E_tmp = x * x;
  /*  T(x), T'(x), T''(x) */
  if (x == 1.0) {
    /*  exactly parabolic; solutions known exactly */
    /*  T(x) */
    T = 1.3333333333333333 * (1.0 - muDoubleScalarPower(q, 3.0));
    /*  T'(x) */
    *Tp = 0.8 * (muDoubleScalarPower(q, 5.0) - 1.0);
    /*  T''(x) */
    *Tpp = *Tp + 1.7142857142857142 * (1.0 - muDoubleScalarPower(q, 7.0));
    /*  T'''(x) */
    *Tppp = 3.0 * (*Tpp - *Tp) +
            2.2222222222222223 * (muDoubleScalarPower(q, 9.0) - 1.0);
  } else if (muDoubleScalarAbs(x - 1.0) < 0.01) {
    real_T b_powers[25];
    real_T dv[25];
    real_T powers[25];
    real_T Tpp_tmp;
    real_T b_Tpp_tmp;
    real_T f;
    real_T g;
    real_T y;
    real_T z;
    real_T z_tmp;
    /*  near-parabolic; compute with series */
    /*  evaluate sigma */
    /*  series approximation to T(x) and its derivatives */
    /*  (used for near-parabolic cases) */
    /*  preload the factors [an] */
    /*  (25 factors is more than enough for 16-digit accuracy) */
    /*  powers of y */
    power(-(E_tmp - 1.0), powers);
    /*  sigma itself */
    /*  dsigma / dx (derivative) */
    /*  d2sigma / dx2 (second derivative) */
    /*  d3sigma / dx3 (third derivative) */
    y = -(E_tmp - 1.0) * q * q;
    /*  series approximation to T(x) and its derivatives */
    /*  (used for near-parabolic cases) */
    /*  preload the factors [an] */
    /*  (25 factors is more than enough for 16-digit accuracy) */
    /*  powers of y */
    power(y, b_powers);
    /*  sigma itself */
    /*  dsigma / dx (derivative) */
    /*  d2sigma / dx2 (second derivative) */
    /*  d3sigma / dx3 (third derivative) */
    /*  T(x) */
    f = 0.0;
    g = 0.0;
    for (i = 0; i < 25; i++) {
      z = an[i];
      f += powers[i] * z;
      g += b_powers[i] * z;
    }
    T = (f + 1.3333333333333333) -
        muDoubleScalarPower(q, 3.0) * (g + 1.3333333333333333);
    /*  T'(x) */
    dv[0] = 1.0;
    for (i = 0; i < 24; i++) {
      dv[i + 1] = ((real_T)i + 2.0) * b_powers[i];
    }
    z = 0.0;
    for (i = 0; i < 25; i++) {
      z += dv[i] * an[i];
    }
    dv[0] = 1.0;
    for (i = 0; i < 24; i++) {
      dv[i + 1] = ((real_T)i + 2.0) * powers[i];
    }
    g = 0.0;
    for (i = 0; i < 25; i++) {
      g += dv[i] * an[i];
    }
    *Tp = 2.0 * x * (muDoubleScalarPower(q, 5.0) * z - g);
    /*  T''(x) */
    z = 1.0 / -(E_tmp - 1.0);
    g = (real_T)iv[0] * z;
    dv[0] = g;
    dv[1] = iv[1];
    for (i = 0; i < 23; i++) {
      dv[i + 2] = (real_T)iv[i + 2] * powers[i];
    }
    b_Tpp_tmp = 0.0;
    for (i = 0; i < 25; i++) {
      b_Tpp_tmp += dv[i] * an[i];
    }
    z_tmp = (real_T)iv[0] * (1.0 / y);
    dv[0] = z_tmp;
    dv[1] = iv[1];
    for (i = 0; i < 23; i++) {
      dv[i + 2] = (real_T)iv[i + 2] * b_powers[i];
    }
    f = 0.0;
    for (i = 0; i < 25; i++) {
      f += dv[i] * an[i];
    }
    Tpp_tmp = *Tp / x;
    *Tpp =
        Tpp_tmp + 4.0 * E_tmp * (b_Tpp_tmp - muDoubleScalarPower(q, 7.0) * f);
    /*  T'''(x) */
    dv[0] = (real_T)b_iv[0] * (1.0 / y / y);
    dv[1] = z_tmp;
    dv[2] = b_iv[2];
    for (i = 0; i < 22; i++) {
      dv[i + 3] = (real_T)b_iv[i + 3] * b_powers[i];
    }
    b_Tpp_tmp = 0.0;
    for (i = 0; i < 25; i++) {
      b_Tpp_tmp += dv[i] * an[i];
    }
    dv[0] = (real_T)b_iv[0] * (z / -(E_tmp - 1.0));
    dv[1] = g;
    dv[2] = b_iv[2];
    for (i = 0; i < 22; i++) {
      dv[i + 3] = (real_T)b_iv[i + 3] * powers[i];
    }
    z = 0.0;
    for (i = 0; i < 25; i++) {
      z += dv[i] * an[i];
    }
    *Tppp = 3.0 * (*Tpp - Tpp_tmp) / x +
            8.0 * x * x * (muDoubleScalarPower(q, 9.0) * b_Tpp_tmp - z);
  } else {
    real_T Tpp_tmp;
    real_T b_Tpp_tmp;
    real_T f;
    real_T g;
    real_T y;
    real_T z;
    real_T z_tmp;
    /*  all other cases */
    /*  compute all substitution functions */
    y = muDoubleScalarAbs(E_tmp - 1.0);
    st.site = &qd_emlrtRSI;
    if (y < 0.0) {
      emlrtErrorWithMessageIdR2018a(
          &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
          "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
    }
    y = muDoubleScalarSqrt(y);
    z_tmp = q * q;
    z = z_tmp * (E_tmp - 1.0) + 1.0;
    st.site = &cd_emlrtRSI;
    if (z < 0.0) {
      emlrtErrorWithMessageIdR2018a(
          &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
          "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
    }
    z = muDoubleScalarSqrt(z);
    f = y * (z - q * x);
    g = x * z - q * (E_tmp - 1.0);
    /*  BUGFIX: (Simon Tardivel) this line is incorrect for E==0 and f+g==0 */
    /*  d  = (E < 0)*(atan2(f, g) + pi*m) + (E > 0)*log( max(0, f + g) ); */
    /*  it should be written out like so: */
    if (E_tmp - 1.0 < 0.0) {
      f = muDoubleScalarAtan2(f, g) + 3.1415926535897931 * m;
    } else if (E_tmp - 1.0 == 0.0) {
      f = 0.0;
    } else {
      f = muDoubleScalarLog(muDoubleScalarMax(0.0, f + g));
    }
    /*  T(x) */
    T = 2.0 * ((x - q * z) - f / y) / (E_tmp - 1.0);
    /*   T'(x) */
    f = muDoubleScalarPower(q, 3.0);
    g = 4.0 * f;
    *Tp = ((4.0 - g * x / z) - 3.0 * x * T) / (E_tmp - 1.0);
    /*  T''(x) */
    Tpp_tmp = z * z;
    b_Tpp_tmp = 1.0 - z_tmp * E_tmp / Tpp_tmp;
    *Tpp =
        ((-4.0 * f / z * b_Tpp_tmp - 3.0 * T) - 3.0 * x * *Tp) / (E_tmp - 1.0);
    /*  T'''(x) */
    *Tppp = ((g / Tpp_tmp * (b_Tpp_tmp + 2.0 * z_tmp * x / Tpp_tmp * (z - x)) -
              8.0 * *Tp) -
             7.0 * x * *Tpp) /
            (E_tmp - 1.0);
  }
  return T;
}

static real_T c_LancasterBlanchard(const emlrtStack *sp, real_T x, real_T q,
                                   real_T m)
{
  emlrtStack st;
  real_T E;
  real_T T;
  int32_T i;
  st.prev = sp;
  st.tls = sp->tls;
  /*  Lancaster & Blanchard's function, and three derivatives thereof */
  /*  protection against idiotic input */
  if (x == -1.0) {
    /*  impossible; offset x slightly */
    x = -0.99999999999999978;
  }
  /*  compute parameter E */
  E = x * x - 1.0;
  /*  T(x), T'(x), T''(x) */
  if (x == 1.0) {
    /*  exactly parabolic; solutions known exactly */
    /*  T(x) */
    T = 1.3333333333333333 * (1.0 - muDoubleScalarPower(q, 3.0));
    /*  T'(x) */
    /*  T''(x) */
    /*  T'''(x) */
  } else if (muDoubleScalarAbs(x - 1.0) < 0.01) {
    real_T b_powers[25];
    real_T powers[25];
    real_T d;
    real_T f;
    /*  near-parabolic; compute with series */
    /*  evaluate sigma */
    /*  series approximation to T(x) and its derivatives */
    /*  (used for near-parabolic cases) */
    /*  preload the factors [an] */
    /*  (25 factors is more than enough for 16-digit accuracy) */
    /*  powers of y */
    power(-E, powers);
    /*  sigma itself */
    /*  dsigma / dx (derivative) */
    /*  d2sigma / dx2 (second derivative) */
    /*  d3sigma / dx3 (third derivative) */
    /*  series approximation to T(x) and its derivatives */
    /*  (used for near-parabolic cases) */
    /*  preload the factors [an] */
    /*  (25 factors is more than enough for 16-digit accuracy) */
    /*  powers of y */
    power(-E * q * q, b_powers);
    /*  sigma itself */
    /*  dsigma / dx (derivative) */
    /*  d2sigma / dx2 (second derivative) */
    /*  d3sigma / dx3 (third derivative) */
    /*  T(x) */
    d = 0.0;
    f = 0.0;
    for (i = 0; i < 25; i++) {
      real_T g;
      g = an[i];
      d += powers[i] * g;
      f += b_powers[i] * g;
    }
    T = (d + 1.3333333333333333) -
        muDoubleScalarPower(q, 3.0) * (f + 1.3333333333333333);
    /*  T'(x) */
    /*  T''(x) */
    /*  T'''(x) */
  } else {
    real_T d;
    real_T f;
    real_T g;
    real_T y;
    real_T z;
    /*  all other cases */
    /*  compute all substitution functions */
    y = muDoubleScalarAbs(E);
    st.site = &qd_emlrtRSI;
    if (y < 0.0) {
      emlrtErrorWithMessageIdR2018a(
          &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
          "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
    }
    y = muDoubleScalarSqrt(y);
    z = q * q * E + 1.0;
    st.site = &cd_emlrtRSI;
    if (z < 0.0) {
      emlrtErrorWithMessageIdR2018a(
          &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
          "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
    }
    z = muDoubleScalarSqrt(z);
    f = y * (z - q * x);
    g = x * z - q * E;
    /*  BUGFIX: (Simon Tardivel) this line is incorrect for E==0 and f+g==0 */
    /*  d  = (E < 0)*(atan2(f, g) + pi*m) + (E > 0)*log( max(0, f + g) ); */
    /*  it should be written out like so: */
    if (E < 0.0) {
      d = muDoubleScalarAtan2(f, g) + 3.1415926535897931 * m;
    } else if (E == 0.0) {
      d = 0.0;
    } else {
      d = muDoubleScalarLog(muDoubleScalarMax(0.0, f + g));
    }
    /*  T(x) */
    T = 2.0 * ((x - q * z) - d / y) / E;
    /*   T'(x) */
    /*  T''(x) */
    /*  T'''(x) */
  }
  return T;
}

static real_T d_LancasterBlanchard(const emlrtStack *sp, real_T x, real_T q,
                                   real_T m, real_T *Tp, real_T *Tpp)
{
  emlrtStack st;
  real_T E_tmp;
  real_T T;
  int32_T i;
  st.prev = sp;
  st.tls = sp->tls;
  /*  Lancaster & Blanchard's function, and three derivatives thereof */
  /*  protection against idiotic input */
  if (x < -1.0) {
    /*  impossible; negative eccentricity */
    x = muDoubleScalarAbs(x) - 2.0;
  } else if (x == -1.0) {
    /*  impossible; offset x slightly */
    x = -0.99999999999999978;
  }
  /*  compute parameter E */
  E_tmp = x * x;
  /*  T(x), T'(x), T''(x) */
  if (x == 1.0) {
    /*  exactly parabolic; solutions known exactly */
    /*  T(x) */
    T = 1.3333333333333333 * (1.0 - muDoubleScalarPower(q, 3.0));
    /*  T'(x) */
    *Tp = 0.8 * (muDoubleScalarPower(q, 5.0) - 1.0);
    /*  T''(x) */
    *Tpp = *Tp + 1.7142857142857142 * (1.0 - muDoubleScalarPower(q, 7.0));
    /*  T'''(x) */
  } else if (muDoubleScalarAbs(x - 1.0) < 0.01) {
    real_T b_powers[25];
    real_T dv[25];
    real_T powers[25];
    real_T f;
    real_T g;
    real_T y;
    real_T z_tmp;
    /*  near-parabolic; compute with series */
    /*  evaluate sigma */
    /*  series approximation to T(x) and its derivatives */
    /*  (used for near-parabolic cases) */
    /*  preload the factors [an] */
    /*  (25 factors is more than enough for 16-digit accuracy) */
    /*  powers of y */
    power(-(E_tmp - 1.0), powers);
    /*  sigma itself */
    /*  dsigma / dx (derivative) */
    /*  d2sigma / dx2 (second derivative) */
    /*  d3sigma / dx3 (third derivative) */
    y = -(E_tmp - 1.0) * q * q;
    /*  series approximation to T(x) and its derivatives */
    /*  (used for near-parabolic cases) */
    /*  preload the factors [an] */
    /*  (25 factors is more than enough for 16-digit accuracy) */
    /*  powers of y */
    power(y, b_powers);
    /*  sigma itself */
    /*  dsigma / dx (derivative) */
    /*  d2sigma / dx2 (second derivative) */
    /*  d3sigma / dx3 (third derivative) */
    /*  T(x) */
    f = 0.0;
    g = 0.0;
    for (i = 0; i < 25; i++) {
      z_tmp = an[i];
      f += powers[i] * z_tmp;
      g += b_powers[i] * z_tmp;
    }
    T = (f + 1.3333333333333333) -
        muDoubleScalarPower(q, 3.0) * (g + 1.3333333333333333);
    /*  T'(x) */
    dv[0] = 1.0;
    for (i = 0; i < 24; i++) {
      dv[i + 1] = ((real_T)i + 2.0) * b_powers[i];
    }
    z_tmp = 0.0;
    for (i = 0; i < 25; i++) {
      z_tmp += dv[i] * an[i];
    }
    dv[0] = 1.0;
    for (i = 0; i < 24; i++) {
      dv[i + 1] = ((real_T)i + 2.0) * powers[i];
    }
    f = 0.0;
    for (i = 0; i < 25; i++) {
      f += dv[i] * an[i];
    }
    *Tp = 2.0 * x * (muDoubleScalarPower(q, 5.0) * z_tmp - f);
    /*  T''(x) */
    dv[0] = (real_T)iv[0] * (1.0 / -(E_tmp - 1.0));
    dv[1] = iv[1];
    for (i = 0; i < 23; i++) {
      dv[i + 2] = (real_T)iv[i + 2] * powers[i];
    }
    z_tmp = 0.0;
    for (i = 0; i < 25; i++) {
      z_tmp += dv[i] * an[i];
    }
    dv[0] = (real_T)iv[0] * (1.0 / y);
    dv[1] = iv[1];
    for (i = 0; i < 23; i++) {
      dv[i + 2] = (real_T)iv[i + 2] * b_powers[i];
    }
    f = 0.0;
    for (i = 0; i < 25; i++) {
      f += dv[i] * an[i];
    }
    *Tpp = *Tp / x + 4.0 * E_tmp * (z_tmp - muDoubleScalarPower(q, 7.0) * f);
    /*  T'''(x) */
  } else {
    real_T f;
    real_T g;
    real_T y;
    real_T z;
    real_T z_tmp;
    /*  all other cases */
    /*  compute all substitution functions */
    y = muDoubleScalarAbs(E_tmp - 1.0);
    st.site = &qd_emlrtRSI;
    if (y < 0.0) {
      emlrtErrorWithMessageIdR2018a(
          &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
          "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
    }
    y = muDoubleScalarSqrt(y);
    z_tmp = q * q;
    z = z_tmp * (E_tmp - 1.0) + 1.0;
    st.site = &cd_emlrtRSI;
    if (z < 0.0) {
      emlrtErrorWithMessageIdR2018a(
          &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
          "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
    }
    z = muDoubleScalarSqrt(z);
    f = y * (z - q * x);
    g = x * z - q * (E_tmp - 1.0);
    /*  BUGFIX: (Simon Tardivel) this line is incorrect for E==0 and f+g==0 */
    /*  d  = (E < 0)*(atan2(f, g) + pi*m) + (E > 0)*log( max(0, f + g) ); */
    /*  it should be written out like so: */
    if (E_tmp - 1.0 < 0.0) {
      f = muDoubleScalarAtan2(f, g) + 3.1415926535897931 * m;
    } else if (E_tmp - 1.0 == 0.0) {
      f = 0.0;
    } else {
      f = muDoubleScalarLog(muDoubleScalarMax(0.0, f + g));
    }
    /*  T(x) */
    T = 2.0 * ((x - q * z) - f / y) / (E_tmp - 1.0);
    /*   T'(x) */
    f = muDoubleScalarPower(q, 3.0);
    *Tp = ((4.0 - 4.0 * f * x / z) - 3.0 * x * T) / (E_tmp - 1.0);
    /*  T''(x) */
    *Tpp = ((-4.0 * f / z * (1.0 - z_tmp * E_tmp / (z * z)) - 3.0 * T) -
            3.0 * x * *Tp) /
           (E_tmp - 1.0);
    /*  T'''(x) */
  }
  return T;
}

static real_T lambert_LancasterBlanchard(const emlrtStack *sp,
                                         const real_T r1vec[3],
                                         const real_T r2vec[3], real_T tf,
                                         real_T m, real_T muC, real_T V1[3],
                                         real_T V2[3],
                                         real_T extremal_distances[2])
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  real_T T;
  real_T T0;
  real_T Td;
  real_T Tp;
  real_T a;
  real_T b_gamma;
  real_T c;
  real_T crsprod_idx_0;
  real_T crsprod_idx_1;
  real_T crsprod_idx_2;
  real_T d;
  real_T dth;
  real_T exitflag;
  real_T mcrsprd;
  real_T phr;
  real_T phr_tmp;
  real_T q;
  real_T r1;
  real_T r1unit_idx_0;
  real_T r1unit_idx_1;
  real_T r1unit_idx_2;
  real_T r2;
  real_T r2unit_idx_0;
  real_T r2unit_idx_1;
  real_T r2unit_idx_2;
  real_T s;
  real_T th1unit_idx_0;
  real_T th1unit_idx_1;
  real_T th1unit_idx_2;
  real_T th2unit_idx_0;
  real_T th2unit_idx_1;
  real_T th2unit_idx_2;
  real_T xM;
  int32_T iterations;
  boolean_T guard1;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  /*  ----------------------------------------------------------------- */
  /*  Lancaster & Blanchard version, with improvements by Gooding */
  /*  Very reliable, moderately fast for both simple and complicated cases */
  /*  ----------------------------------------------------------------- */
  /* #coder */
  /* { */
  /* LAMBERT_LANCASTERBLANCHARD       High-Thrust Lambert-targeter */
  /*  */
  /* lambert_LancasterBlanchard() uses the method developed by */
  /* Lancaster & Blancard, as described in their 1969 paper. Initial values, */
  /* and several details of the procedure, are provided by R.H. Gooding, */
  /* as described in his 1990 paper. */
  /* } */
  /*  Please report bugs and inquiries to: */
  /*  */
  /*  Name       : Rody P.S. Oldenhuis */
  /*  E-mail     : oldenhuis@gmail.com */
  /*  Licence    : 2-clause BSD (see License.txt) */
  /*  If you find this work useful, please consider a donation: */
  /*  https://www.paypal.me/RodyO/3.5 */
  /*  ADJUSTED FOR EML-COMPILATION 29/Sep/2009 */
  /*  manipulate input */
  /*  optimum for numerical noise v.s. actual precision */
  r1 = (r1vec[0] * r1vec[0] + r1vec[1] * r1vec[1]) + r1vec[2] * r1vec[2];
  st.site = &rb_emlrtRSI;
  if (r1 < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  r1 = muDoubleScalarSqrt(r1);
  /*  magnitude of r1vec */
  r2 = (r2vec[0] * r2vec[0] + r2vec[1] * r2vec[1]) + r2vec[2] * r2vec[2];
  st.site = &sb_emlrtRSI;
  if (r2 < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  r2 = muDoubleScalarSqrt(r2);
  /*  magnitude of r2vec */
  /*  unit vector of r1vec */
  /*  unit vector of r2vec */
  crsprod_idx_0 = r1vec[1] * r2vec[2] - r2vec[1] * r1vec[2];
  crsprod_idx_1 = r2vec[0] * r1vec[2] - r1vec[0] * r2vec[2];
  crsprod_idx_2 = r1vec[0] * r2vec[1] - r2vec[0] * r1vec[1];
  /*  cross product of r1vec and r2vec */
  r1unit_idx_0 = r1vec[0] / r1;
  r2unit_idx_0 = r2vec[0] / r2;
  r1unit_idx_1 = r1vec[1] / r1;
  r2unit_idx_1 = r2vec[1] / r2;
  r1unit_idx_2 = r1vec[2] / r1;
  r2unit_idx_2 = r2vec[2] / r2;
  mcrsprd = (crsprod_idx_0 * crsprod_idx_0 + crsprod_idx_1 * crsprod_idx_1) +
            crsprod_idx_2 * crsprod_idx_2;
  st.site = &tb_emlrtRSI;
  if (mcrsprd < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  mcrsprd = muDoubleScalarSqrt(mcrsprd);
  /*  magnitude of that cross product */
  crsprod_idx_0 /= mcrsprd;
  crsprod_idx_1 /= mcrsprd;
  crsprod_idx_2 /= mcrsprd;
  th1unit_idx_0 = crsprod_idx_1 * r1unit_idx_2 - r1unit_idx_1 * crsprod_idx_2;
  th1unit_idx_1 = r1unit_idx_0 * crsprod_idx_2 - crsprod_idx_0 * r1unit_idx_2;
  th1unit_idx_2 = crsprod_idx_0 * r1unit_idx_1 - r1unit_idx_0 * crsprod_idx_1;
  /*  unit vectors in the tangential-directions */
  th2unit_idx_0 = crsprod_idx_1 * r2unit_idx_2 - r2unit_idx_1 * crsprod_idx_2;
  th2unit_idx_1 = r2unit_idx_0 * crsprod_idx_2 - crsprod_idx_0 * r2unit_idx_2;
  th2unit_idx_2 = crsprod_idx_0 * r2unit_idx_1 - r2unit_idx_0 * crsprod_idx_1;
  /*  make 100.4% sure it's in (-1 <= x <= +1) */
  dth = muDoubleScalarMax(-1.0, muDoubleScalarMin(1.0, ((r1vec[0] * r2vec[0] +
                                                         r1vec[1] * r2vec[1]) +
                                                        r1vec[2] * r2vec[2]) /
                                                           r1 / r2));
  st.site = &ub_emlrtRSI;
  if ((dth < -1.0) || (dth > 1.0)) {
    emlrtErrorWithMessageIdR2018a(
        &st, &d_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "acos");
  }
  dth = muDoubleScalarAcos(dth);
  /*  turn angle */
  /*  if the long way was selected, the turn-angle must be negative */
  /*  to take care of the direction of final velocity */
  mcrsprd = muDoubleScalarSign(tf);
  tf = muDoubleScalarAbs(tf);
  if (mcrsprd < 0.0) {
    dth -= 6.2831853071795862;
  }
  /*  left-branch */
  crsprod_idx_1 = muDoubleScalarSign(m);
  m = muDoubleScalarAbs(m);
  /*  define constants */
  st.site = &vb_emlrtRSI;
  b_st.site = &lb_emlrtRSI;
  st.site = &vb_emlrtRSI;
  b_st.site = &lb_emlrtRSI;
  c = (r1 * r1 + r2 * r2) - 2.0 * r1 * r2 * muDoubleScalarCos(dth);
  st.site = &vb_emlrtRSI;
  if (c < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  c = muDoubleScalarSqrt(c);
  s = ((r1 + r2) + c) / 2.0;
  st.site = &wb_emlrtRSI;
  b_st.site = &lb_emlrtRSI;
  d = 8.0 * muC / muDoubleScalarPower(s, 3.0);
  st.site = &wb_emlrtRSI;
  if (d < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  d = muDoubleScalarSqrt(d);
  T = d * tf;
  d = r1 * r2;
  st.site = &xb_emlrtRSI;
  if (d < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  q = muDoubleScalarSqrt(d) / s * muDoubleScalarCos(dth / 2.0);
  /*  general formulae for the initial values (Gooding) */
  /*  ------------------------------------------------- */
  /*  some initial values */
  st.site = &yb_emlrtRSI;
  T0 = LancasterBlanchard(&st, q, m);
  Td = T0 - T;
  st.site = &ac_emlrtRSI;
  b_st.site = &lb_emlrtRSI;
  phr_tmp = q * q;
  phr = b_mod(2.0 * muDoubleScalarAtan2(1.0 - phr_tmp, 2.0 * q));
  /*  initial output is pessimistic */
  V1[0] = rtNaN;
  V2[0] = rtNaN;
  V1[1] = rtNaN;
  V2[1] = rtNaN;
  V1[2] = rtNaN;
  V2[2] = rtNaN;
  extremal_distances[0] = rtNaN;
  extremal_distances[1] = rtNaN;
  /*  single-revolution case */
  guard1 = false;
  if (m == 0.0) {
    if (Td > 0.0) {
      phr = T0 * Td / 4.0 / T;
    } else {
      crsprod_idx_0 = Td / (4.0 - Td);
      crsprod_idx_1 = -Td / (T + T0 / 2.0);
      st.site = &bc_emlrtRSI;
      if (crsprod_idx_1 < 0.0) {
        emlrtErrorWithMessageIdR2018a(
            &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
            "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
      }
      crsprod_idx_1 = muDoubleScalarSqrt(crsprod_idx_1);
      mcrsprd = 2.0 - phr / 3.1415926535897931;
      st.site = &cc_emlrtRSI;
      if (mcrsprd < 0.0) {
        emlrtErrorWithMessageIdR2018a(
            &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
            "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
      }
      mcrsprd = muDoubleScalarSqrt(mcrsprd);
      mcrsprd = crsprod_idx_0 + 1.7 * mcrsprd;
      if (mcrsprd >= 0.0) {
        crsprod_idx_2 = crsprod_idx_0;
      } else {
        st.site = &dc_emlrtRSI;
        b_st.site = &mb_emlrtRSI;
        crsprod_idx_2 = crsprod_idx_0 + muDoubleScalarPower(-mcrsprd, 0.0625) *
                                            (-crsprod_idx_1 - crsprod_idx_0);
      }
      st.site = &ec_emlrtRSI;
      b_st.site = &lb_emlrtRSI;
      st.site = &ec_emlrtRSI;
      if (crsprod_idx_0 + 1.0 < 0.0) {
        emlrtErrorWithMessageIdR2018a(
            &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
            "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
      }
      phr = ((crsprod_idx_2 * (crsprod_idx_0 + 1.0) / 2.0 + 1.0) -
             0.03 * (crsprod_idx_2 * crsprod_idx_2) *
                 muDoubleScalarSqrt(crsprod_idx_0 + 1.0)) *
            crsprod_idx_2;
    }
    /*  this estimate might not give a solution */
    if (phr < -1.0) {
      exitflag = -1.0;
    } else {
      /*  multi-revolution case */
      guard1 = true;
    }
  } else {
    int32_T exitg2;
    /*  determine minimum Tp(x) */
    mcrsprd = 4.0 / (9.42477796076938 * (2.0 * m + 1.0));
    if (phr < 3.1415926535897931) {
      st.site = &fc_emlrtRSI;
      a = phr / 3.1415926535897931;
      b_st.site = &lb_emlrtRSI;
      c_st.site = &mb_emlrtRSI;
      if (a < 0.0) {
        emlrtErrorWithMessageIdR2018a(&c_st, &b_emlrtRTEI,
                                      "Coder:toolbox:power_domainError",
                                      "Coder:toolbox:power_domainError", 0);
      }
      xM = mcrsprd * muDoubleScalarPower(a, 0.125);
    } else if (phr > 3.1415926535897931) {
      st.site = &gc_emlrtRSI;
      a = 2.0 - phr / 3.1415926535897931;
      b_st.site = &lb_emlrtRSI;
      c_st.site = &mb_emlrtRSI;
      if (a < 0.0) {
        emlrtErrorWithMessageIdR2018a(&c_st, &b_emlrtRTEI,
                                      "Coder:toolbox:power_domainError",
                                      "Coder:toolbox:power_domainError", 0);
      }
      xM = mcrsprd * (2.0 - muDoubleScalarPower(a, 0.125));
      /*  EMLMEX requires this one */
    } else {
      xM = 0.0;
    }
    /*  use Halley's method */
    Tp = rtInf;
    iterations = 0;
    do {
      exitg2 = 0;
      if (muDoubleScalarAbs(Tp) > 1.0E-12) {
        /*  iterations */
        iterations++;
        /*  compute first three derivatives */
        st.site = &hc_emlrtRSI;
        b_LancasterBlanchard(&st, xM, q, m, &Tp, &b_gamma, &crsprod_idx_2);
        /*  new value of xM */
        mcrsprd = xM;
        st.site = &ic_emlrtRSI;
        b_st.site = &mb_emlrtRSI;
        xM -= 2.0 * Tp * b_gamma /
              (2.0 * (b_gamma * b_gamma) - Tp * crsprod_idx_2);
        /*  escape clause */
        st.site = &jc_emlrtRSI;
        if (muDoubleScalarRem(iterations, 7.0) != 0.0) {
          xM = (mcrsprd + xM) / 2.0;
        }
        /*  the method might fail. Exit in that case */
        if (*emlrtBreakCheckR2012bFlagVar != 0) {
          emlrtBreakCheckR2012b((emlrtConstCTX)sp);
        }
        if (iterations > 25) {
          exitflag = -2.0;
          exitg2 = 1;
        }
      } else {
        /*  xM should be elliptic (-1 < x < 1) */
        /*  (this should be impossible to go wrong) */
        exitg2 = 2;
      }
    } while (exitg2 == 0);
    if (exitg2 != 1) {
      if ((xM < -1.0) || (xM > 1.0)) {
        exitflag = -1.0;
      } else {
        /*  corresponding time */
        st.site = &kc_emlrtRSI;
        mcrsprd = c_LancasterBlanchard(&st, xM, q, m);
        /*  T should lie above the minimum T */
        if (mcrsprd > T) {
          exitflag = -1.0;
        } else {
          /*  find two initial values for second solution (again with
           * lambda-type patch) */
          /*  --------------------------------------------------------------------------
           */
          /*  some initial values */
          crsprod_idx_2 = T - mcrsprd;
          st.site = &lc_emlrtRSI;
          d_LancasterBlanchard(&st, xM, q, m, &Tp, &b_gamma);
          /*  first estimate (only if m > 0) */
          if (crsprod_idx_1 > 0.0) {
            st.site = &mc_emlrtRSI;
            b_st.site = &lb_emlrtRSI;
            phr = crsprod_idx_2 /
                  (b_gamma / 2.0 + crsprod_idx_2 / ((1.0 - xM) * (1.0 - xM)));
            st.site = &mc_emlrtRSI;
            if (phr < 0.0) {
              emlrtErrorWithMessageIdR2018a(
                  &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
                  "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
            }
            phr = muDoubleScalarSqrt(phr);
            mcrsprd = xM + phr;
            st.site = &nc_emlrtRSI;
            b_st.site = &lb_emlrtRSI;
            mcrsprd = 4.0 * mcrsprd / (crsprod_idx_2 + 4.0) +
                      (1.0 - mcrsprd) * (1.0 - mcrsprd);
            st.site = &oc_emlrtRSI;
            if (mcrsprd < 0.0) {
              emlrtErrorWithMessageIdR2018a(
                  &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
                  "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
            }
            phr = phr * (1.0 - ((m + 1.0) + (dth - 0.5)) / (0.15 * m + 1.0) *
                                   phr *
                                   (mcrsprd / 2.0 +
                                    0.03 * phr * muDoubleScalarSqrt(mcrsprd))) +
                  xM;
            /*  first estimate might not be able to yield possible solution */
            if (phr > 1.0) {
              exitflag = -1.0;
            } else {
              /*  second estimate (only if m > 0) */
              guard1 = true;
            }
          } else {
            if (Td > 0.0) {
              st.site = &pc_emlrtRSI;
              b_st.site = &lb_emlrtRSI;
              crsprod_idx_1 =
                  mcrsprd / (b_gamma / 2.0 -
                             crsprod_idx_2 * (b_gamma / 2.0 / (T0 - mcrsprd) -
                                              1.0 / (xM * xM)));
              st.site = &pc_emlrtRSI;
              if (crsprod_idx_1 < 0.0) {
                emlrtErrorWithMessageIdR2018a(
                    &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
                    "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
              }
              crsprod_idx_1 = muDoubleScalarSqrt(crsprod_idx_1);
              phr = xM - crsprod_idx_1;
            } else {
              crsprod_idx_2 = Td / (4.0 - Td);
              crsprod_idx_1 = 2.0 * (1.0 - phr);
              st.site = &qc_emlrtRSI;
              if (crsprod_idx_1 < 0.0) {
                emlrtErrorWithMessageIdR2018a(
                    &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
                    "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
              }
              crsprod_idx_1 = muDoubleScalarSqrt(crsprod_idx_1);
              mcrsprd = crsprod_idx_2 + 1.7 * crsprod_idx_1;
              if (!(mcrsprd >= 0.0)) {
                st.site = &rc_emlrtRSI;
                b_st.site = &lb_emlrtRSI;
                crsprod_idx_1 = muDoubleScalarPower(-mcrsprd, 0.125);
                st.site = &rc_emlrtRSI;
                if (crsprod_idx_1 < 0.0) {
                  emlrtErrorWithMessageIdR2018a(
                      &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
                      "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
                }
                crsprod_idx_1 = muDoubleScalarSqrt(crsprod_idx_1);
                mcrsprd = -Td / (1.5 * T0 - Td);
                st.site = &rc_emlrtRSI;
                if (mcrsprd < 0.0) {
                  emlrtErrorWithMessageIdR2018a(
                      &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
                      "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
                }
                mcrsprd = muDoubleScalarSqrt(mcrsprd);
                crsprod_idx_2 -= crsprod_idx_1 * (crsprod_idx_2 + mcrsprd);
              }
              mcrsprd = 4.0 / (4.0 - Td);
              st.site = &sc_emlrtRSI;
              if (mcrsprd < 0.0) {
                emlrtErrorWithMessageIdR2018a(
                    &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
                    "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
              }
              phr = crsprod_idx_2 *
                    (((m + 1.0) + 0.24 * (dth - 0.5)) / (0.15 * m + 1.0) *
                         crsprod_idx_2 *
                         (mcrsprd / 2.0 -
                          0.03 * crsprod_idx_2 * muDoubleScalarSqrt(mcrsprd)) +
                     1.0);
            }
            /*  estimate might not give solutions */
            if (phr < -1.0) {
              exitflag = -1.0;
            } else {
              guard1 = true;
            }
          }
        }
      }
    }
  }
  if (guard1) {
    /*  find root of Lancaster & Blancard's function */
    /*  -------------------------------------------- */
    /*  (Halley's method) */
    mcrsprd = rtInf;
    iterations = 0;
    int32_T exitg1;
    do {
      exitg1 = 0;
      if (muDoubleScalarAbs(mcrsprd) > 1.0E-12) {
        /*  iterations */
        iterations++;
        /*  compute function value, and first two derivatives */
        st.site = &tc_emlrtRSI;
        mcrsprd = d_LancasterBlanchard(&st, phr, q, m, &Tp, &b_gamma);
        /*  find the root of the *difference* between the */
        /*  function value [T_x] and the required time [T] */
        mcrsprd -= T;
        /*  new value of x */
        crsprod_idx_2 = phr;
        st.site = &uc_emlrtRSI;
        b_st.site = &lb_emlrtRSI;
        phr -= 2.0 * mcrsprd * Tp / (2.0 * (Tp * Tp) - mcrsprd * b_gamma);
        /*  escape clause */
        st.site = &vc_emlrtRSI;
        if (muDoubleScalarRem(iterations, 7.0) != 0.0) {
          phr = (crsprod_idx_2 + phr) / 2.0;
        }
        /*  Halley's method might fail */
        if (*emlrtBreakCheckR2012bFlagVar != 0) {
          emlrtBreakCheckR2012b((emlrtConstCTX)sp);
        }
        if (iterations > 25) {
          exitflag = -2.0;
          exitg1 = 1;
        }
      } else {
        /*  calculate terminal velocities */
        /*  ----------------------------- */
        /*  constants required for this calculation */
        b_gamma = muC * s / 2.0;
        st.site = &wc_emlrtRSI;
        if (b_gamma < 0.0) {
          emlrtErrorWithMessageIdR2018a(
              &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
              "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
        }
        b_gamma = muDoubleScalarSqrt(b_gamma);
        if (c == 0.0) {
          Tp = 1.0;
          mcrsprd = 0.0;
          crsprod_idx_0 = muDoubleScalarAbs(phr);
        } else {
          st.site = &xc_emlrtRSI;
          b_st.site = &lb_emlrtRSI;
          d /= c * c;
          st.site = &xc_emlrtRSI;
          if (d < 0.0) {
            emlrtErrorWithMessageIdR2018a(
                &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
                "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
          }
          d = muDoubleScalarSqrt(d);
          Tp = 2.0 * d * muDoubleScalarSin(dth / 2.0);
          mcrsprd = (r1 - r2) / c;
          st.site = &yc_emlrtRSI;
          b_st.site = &lb_emlrtRSI;
          st.site = &yc_emlrtRSI;
          b_st.site = &lb_emlrtRSI;
          crsprod_idx_0 = phr_tmp * (phr * phr - 1.0) + 1.0;
          st.site = &yc_emlrtRSI;
          if (crsprod_idx_0 < 0.0) {
            emlrtErrorWithMessageIdR2018a(
                &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
                "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
          }
          crsprod_idx_0 = muDoubleScalarSqrt(crsprod_idx_0);
        }
        /*  radial component */
        xM = q * crsprod_idx_0;
        crsprod_idx_2 = xM - phr;
        xM = mcrsprd * (xM + phr);
        a = b_gamma * (crsprod_idx_2 - xM) / r1;
        crsprod_idx_1 = -b_gamma * (crsprod_idx_2 + xM) / r2;
        /*  tangential component */
        xM = Tp * b_gamma * (crsprod_idx_0 + q * phr);
        crsprod_idx_2 = xM / r1;
        mcrsprd = xM / r2;
        /*  Cartesian velocity */
        V1[0] = crsprod_idx_2 * th1unit_idx_0 + a * r1unit_idx_0;
        V2[0] = mcrsprd * th2unit_idx_0 + crsprod_idx_1 * r2unit_idx_0;
        V1[1] = crsprod_idx_2 * th1unit_idx_1 + a * r1unit_idx_1;
        V2[1] = mcrsprd * th2unit_idx_1 + crsprod_idx_1 * r2unit_idx_1;
        V1[2] = crsprod_idx_2 * th1unit_idx_2 + a * r1unit_idx_2;
        V2[2] = mcrsprd * th2unit_idx_2 + crsprod_idx_1 * r2unit_idx_2;
        /*  exitflag */
        exitflag = 1.0;
        /*  (success) */
        /*  also determine minimum/maximum distance */
        st.site = &ad_emlrtRSI;
        b_st.site = &lb_emlrtRSI;
        /*  semi-major axis */
        st.site = &bd_emlrtRSI;
        minmax_distances(&st, r1vec, r1, r1vec, r2, dth,
                         s / 2.0 / (1.0 - phr * phr), V1, V2, m, muC,
                         extremal_distances);
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }
  return exitflag;
}

static void minmax_distances(const emlrtStack *sp, const real_T r1vec[3],
                             real_T r1, const real_T r2vec[3], real_T r2,
                             real_T dth, real_T a, const real_T V1[3],
                             const real_T V2[3], real_T m, real_T muC,
                             real_T extremal_distances[2])
{
  emlrtStack st;
  real_T absx;
  real_T apocenter;
  real_T d;
  real_T e;
  real_T evec_idx_0;
  real_T evec_idx_1;
  real_T evec_tmp_idx_0;
  real_T evec_tmp_idx_1;
  real_T evec_tmp_idx_2;
  real_T maximum_distance;
  real_T minimum_distance;
  real_T pericenter;
  real_T theta1;
  real_T theta2;
  int32_T exponent;
  st.prev = sp;
  st.tls = sp->tls;
  /*  ----------------------------------------------------------------- */
  /*  Helper functions */
  /*  ----------------------------------------------------------------- */
  /*  compute minimum and maximum distances to the central body */
  /*  default - minimum/maximum of r1,r2 */
  minimum_distance = muDoubleScalarMin(r1, r2);
  maximum_distance = muDoubleScalarMax(r1, r2);
  /*  was the longway used or not? */
  /*  eccentricity vector (use triple product identity) */
  absx = (V1[0] * V1[0] + V1[1] * V1[1]) + V1[2] * V1[2];
  theta1 = (V1[0] * r1vec[0] + V1[1] * r1vec[1]) + V1[2] * r1vec[2];
  /*  eccentricity */
  evec_tmp_idx_2 = r1vec[0] / r1;
  evec_tmp_idx_0 = evec_tmp_idx_2;
  d = (absx * r1vec[0] - theta1 * V1[0]) / muC - evec_tmp_idx_2;
  evec_idx_0 = d;
  theta2 = d * d;
  evec_tmp_idx_2 = r1vec[1] / r1;
  evec_tmp_idx_1 = evec_tmp_idx_2;
  d = (absx * r1vec[1] - theta1 * V1[1]) / muC - evec_tmp_idx_2;
  evec_idx_1 = d;
  theta2 += d * d;
  evec_tmp_idx_2 = r1vec[2] / r1;
  d = (absx * r1vec[2] - theta1 * V1[2]) / muC - evec_tmp_idx_2;
  theta2 += d * d;
  st.site = &td_emlrtRSI;
  if (theta2 < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  e = muDoubleScalarSqrt(theta2);
  /*  apses */
  pericenter = a * (1.0 - e);
  apocenter = rtInf;
  /*  parabolic/hyperbolic case */
  if (e < 1.0) {
    apocenter = a * (e + 1.0);
  }
  /*  elliptic case */
  /*  since we have the eccentricity vector, we know exactly where the */
  /*  pericenter lies. Use this fact, and the given value of [dth], to */
  /*  cross-check if the trajectory goes past it */
  if (m > 0.0) {
    /*  obvious case (always elliptical and both apses are traversed) */
    minimum_distance = pericenter;
    maximum_distance = apocenter;
  } else {
    real_T unnamed_idx_2;
    /*  less obvious case */
    /*  compute theta1&2 ( use (AxB)-(CxD) = (C·B)(D·A) - (C·A)(B·D) )) */
    /*  make 100.4% sure it's in (-1 <= theta12 <= +1) */
    absx = evec_idx_0 / e;
    theta2 = evec_idx_1 / e;
    unnamed_idx_2 = d / e;
    theta1 =
        muDoubleScalarSign(
            r1 * r1 * ((evec_idx_0 * V1[0] + evec_idx_1 * V1[1]) + d * V1[2]) -
            ((r1vec[0] * evec_idx_0 + r1vec[1] * evec_idx_1) + r1vec[2] * d) *
                theta1) *
        muDoubleScalarAcos(muDoubleScalarMax(
            -1.0, muDoubleScalarMin(
                      1.0, (evec_tmp_idx_0 * absx + evec_tmp_idx_1 * theta2) +
                               evec_tmp_idx_2 * unnamed_idx_2)));
    theta2 =
        muDoubleScalarSign(
            r2 * r2 * ((evec_idx_0 * V2[0] + evec_idx_1 * V2[1]) + d * V2[2]) -
            ((r2vec[0] * evec_idx_0 + r2vec[1] * evec_idx_1) + r2vec[2] * d) *
                ((r2vec[0] * V2[0] + r2vec[1] * V2[1]) + r2vec[2] * V2[2])) *
        muDoubleScalarAcos(muDoubleScalarMax(
            -1.0, muDoubleScalarMin(
                      1.0, (r2vec[0] / r2 * absx + r2vec[1] / r2 * theta2) +
                               r2vec[2] / r2 * unnamed_idx_2)));
    /*  points 1&2 are on opposite sides of the symmetry axis -- minimum */
    /*  and maximum distance depends both on the value of [dth], and both */
    /*  [theta1] and [theta2] */
    if (theta1 * theta2 < 0.0) {
      /*  if |th1| + |th2| = turnangle, we know that the pericenter was */
      /*  passed */
      absx = muDoubleScalarAbs(dth);
      if (muDoubleScalarIsInf(absx) || muDoubleScalarIsNaN(absx)) {
        absx = rtNaN;
      } else if (absx < 4.4501477170144028E-308) {
        absx = 4.94065645841247E-324;
      } else {
        frexp(absx, &exponent);
        absx = ldexp(1.0, exponent - 53);
      }
      if (muDoubleScalarAbs(
              (muDoubleScalarAbs(theta1) + muDoubleScalarAbs(theta2)) - dth) <
          5.0 * absx) {
        minimum_distance = pericenter;
        /*  this condition can only be false for elliptic cases, and */
        /*  when it is indeed false, we know that the orbit passed */
        /*  apocenter */
      } else {
        maximum_distance = apocenter;
      }
      /*  points 1&2 are on the same side of the symmetry axis. Only if the */
      /*  long-way was used are the min. and max. distances different from */
      /*  the min. and max. values of the radii (namely, equal to the apses) */
    } else if (muDoubleScalarAbs(dth) > 3.1415926535897931) {
      minimum_distance = pericenter;
      if (e < 1.0) {
        maximum_distance = apocenter;
      }
    }
  }
  /*  output argument */
  extremal_distances[0] = minimum_distance;
  extremal_distances[1] = maximum_distance;
}

void lambert(const emlrtStack *sp, real_T r1vec[3], real_T r2vec[3], real_T tf,
             real_T m, real_T muC, real_T V1[3], real_T V2[3],
             real_T extremal_distances[2], real_T *exitflag)
{
  emlrtStack b_st;
  emlrtStack st;
  real_T crossprd[3];
  real_T aalfa[2];
  real_T bbeta[2];
  real_T dv[2];
  real_T Lambda;
  real_T T;
  real_T V;
  real_T a;
  real_T a_min;
  real_T aa_idx_0;
  real_T b_y1;
  real_T c;
  real_T d;
  real_T d1;
  real_T dth;
  real_T inn1;
  real_T inn2;
  real_T leftbranch;
  real_T logt;
  real_T longway;
  real_T mcr;
  real_T mr2vec;
  real_T r1;
  real_T s;
  real_T x1;
  real_T x2;
  real_T x_tmp;
  real_T x_tmp_tmp;
  real_T xnew;
  real_T xx_idx_0;
  int32_T iterations;
  boolean_T bad;
  boolean_T exitg1;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  /*  LAMBERT            Lambert-targeter for ballistic flights */
  /*                     (Izzo, and Lancaster, Blanchard & Gooding) */
  /*  */
  /*  Usage: */
  /*     [V1, V2, extremal_distances, exitflag] = lambert(r1, r2, tf, m,
   * GM_central) */
  /*  */
  /*  Dimensions: */
  /*              r1, r2 ->  [1x3] */
  /*              V1, V2 ->  [1x3] */
  /*  extremal_distances ->  [1x2] */
  /*               tf, m ->  [1x1] */
  /*          GM_central ->  [1x1] */
  /*  */
  /*  This function solves any Lambert problem *robustly*. It uses two separate
   */
  /*  solvers; the first one tried is a new and unpublished algorithm developed
   */
  /*  by Dr. D. Izzo from the European Space Agency [1]. This version is
   * extremely */
  /*  fast, but especially for larger [m] it still fails quite frequently. In
   * such */
  /*  cases, a MUCH more robust algorithm is started (the one by Lancaster & */
  /*  Blancard [2], with modifications, initial values and other improvements by
   */
  /*  R.Gooding [3]), which is a lot slower partly because of its robustness. */
  /*  */
  /*  INPUT ARGUMENTS: */
  /*  ====================================================================== */
  /*     name        units    description */
  /*  ====================================================================== */
  /*    r1, r1       [km]     position vectors of the two terminal points. */
  /*      tf        [days]    time of flight to solve for */
  /*       m          [-]     specifies the number of complete orbits to
   * complete */
  /*                          (should be an integer) */
  /*  GM_central   [km3/s2]   std. grav. parameter (G�M = mu) of the central
   * body */
  /*  */
  /*  OUTPUT ARGUMENTS: */
  /*  ====================================================================== */
  /*    name             units   description */
  /*  ====================================================================== */
  /*   V1, V2             [km/s]  terminal velocities at the end-points */
  /*   extremal_distances  [km]   minimum(1) and maximum(2) distance of the */
  /*                              spacecraft to the central body. */
  /*   exitflag             [-]   Integer containing information on why the */
  /*                              routine terminated. A value of +1 indicates */
  /*                              success; a normal exit. A value of -1 */
  /*                              indicates that the given problem has no */
  /*                              solution and cannot be solved. A value of -2
   */
  /*                              indicates that both algorithms failed to find
   */
  /*                              a solution. This should never occur since */
  /*                              these problems are well-defined, and at the */
  /*                              very least it can be determined that the */
  /*                              problem has no solution. Nevertheless, it */
  /*                              still occurs sometimes for accidental */
  /*                              erroneous input, so it provides a basic */
  /*                              mechanism to check any application using this
   */
  /*                              algorithm. */
  /*  */
  /*  This routine can be compiled to increase its speed by a factor of about */
  /*  10-15, which is certainly advisable when the complete application requires
   */
  /*  a great number of Lambert problems to be solved. The entire routine is */
  /*  written in embedded MATLAB, so it can be compiled with the emlmex() */
  /*  function (older MATLAB) or codegen() function (MATLAB 2011a and later). */
  /*  */
  /*  To do this using emlmex(), make sure MATLAB's current directory is equal
   */
  /*  to where this file is located. Then, copy-paste and execute the following
   */
  /*  commands to the command window: */
  /*  */
  /*     example_input = {... */
  /*          [0.0, 0.0, 0.0], ...% r1vec */
  /*          [0.0, 0.0, 0.0], ...% r2vec */
  /*           0.0, ...           % tf */
  /*           0.0, ...           % m */
  /*           0.0};              % muC */
  /*     emlmex -eg example_input lambert.m */
  /*  */
  /*  This is of course assuming your compiler is configured correctly. See the
   */
  /*  documentation of emlmex() on how to do that. */
  /*  */
  /*  Using codegen(), the syntax is as follows: */
  /*  */
  /*     example_input = {... */
  /*          [0.0, 0.0, 0.0], ...% r1vec */
  /*          [0.0, 0.0, 0.0], ...% r2vec */
  /*           0.0, ...           % tf */
  /*           0.0, ...           % m */
  /*           0.0};              % muC */
  /*     codegen lambert.m -args example_input */
  /*  */
  /*  Note that in newer MATLAB versions, the code analyzer will complain about
   */
  /*  the pragma "%#eml" after the main function's name, and possibly, issue */
  /*  subsequent warnings related to this issue. To get rid of this problem,
   * simply */
  /*  replace the "%#eml" directive with "%#codegen". */
  /*  */
  /*  References: */
  /*  */
  /*  [1] Izzo, D. ESA Advanced Concepts team. Code used available in MGA.M, on
   */
  /*      http://www.esa.int/gsp/ACT/inf/op/globopt.htm. Last retreived Nov,
   * 2009. */
  /*  [2] Lancaster, E.R. and Blanchard, R.C. "A unified form of Lambert's
   * theorem." */
  /*      NASA technical note TN D-5368,1969. */
  /*  [3] Gooding, R.H. "A procedure for the solution of Lambert's orbital
   * boundary-value */
  /*      problem. Celestial Mechanics and Dynamical Astronomy, 48:145�165,1990.
   */
  /*  */
  /*  See also lambert_low_ExpoSins. */
  /*  Please report bugs and inquiries to: */
  /*  */
  /*  Name       : Rody P.S. Oldenhuis */
  /*  E-mail     : oldenhuis@gmail.com */
  /*  Licence    : 2-clause BSD (see License.txt) */
  /*  If you find this work useful, please consider a donation: */
  /*  https://www.paypal.me/RodyO/3.5 */
  /*  If you want to cite this work in an academic paper, please use */
  /*  the following template: */
  /*  */
  /*  Rody Oldenhuis, orcid.org/0000-0002-3162-3660. "Lambert" <version>, */
  /*  <date you last used it>. MATLAB Robust solver for Lambert's */
  /*  orbital-boundary value problem. */
  /*  https://nl.mathworks.com/matlabcentral/fileexchange/26348 */
  /*  ----------------------------------------------------------------- */
  /*  Izzo's version: */
  /*  Very fast, but not very robust for more complicated cases */
  /*  ----------------------------------------------------------------- */
  /* #coder */
  /*  original documentation: */
  /* { */
  /*  This routine implements a new algorithm that solves Lambert's problem. The
   */
  /*  algorithm has two major characteristics that makes it favorable to other
   */
  /*  existing ones. */
  /*  */
  /*  1) It describes the generic orbit solution of the boundary condition */
  /*  problem through the variable X=log(1+cos(alpha/2)). By doing so the */
  /*  graph of the time of flight become defined in the entire real axis and */
  /*  resembles a straight line. Convergence is granted within few iterations */
  /*  for all the possible geometries (except, of course, when the transfer */
  /*  angle is zero). When multiple revolutions are considered the variable is
   */
  /*  X=tan(cos(alpha/2)*pi/2). */
  /*  */
  /*  2) Once the orbit has been determined in the plane, this routine */
  /*  evaluates the velocity vectors at the two points in a way that is not */
  /*  singular for the transfer angle approaching to pi (Lagrange coefficient */
  /*  based methods are numerically not well suited for this purpose). */
  /*  */
  /*  As a result Lambert's problem is solved (with multiple revolutions */
  /*  being accounted for) with the same computational effort for all */
  /*  possible geometries. The case of near 180 transfers is also solved */
  /*  efficiently. */
  /*  */
  /*   We note here that even when the transfer angle is exactly equal to pi */
  /*  the algorithm does solve the problem in the plane (it finds X), but it */
  /*  is not able to evaluate the plane in which the orbit lies. A solution */
  /*  to this would be to provide the direction of the plane containing the */
  /*  transfer orbit from outside. This has not been implemented in this */
  /*  routine since such a direction would depend on which application the */
  /*  transfer is going to be used in. */
  /*  */
  /*  please report bugs to dario.izzo@esa.int */
  /* } */
  /*  adjusted documentation: */
  /* { */
  /*  By default, the short-way solution is computed. The long way solution */
  /*  may be requested by giving a negative value to the corresponding */
  /*  time-of-flight [tf]. */
  /*  */
  /*  For problems with |m| > 0, there are generally two solutions. By */
  /*  default, the right branch solution will be returned. The left branch */
  /*  may be requested by giving a negative value to the corresponding */
  /*  number of complete revolutions [m]. */
  /* } */
  /*  Authors */
  /*  .-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-.-`-. */
  /*  Name       : Dr. Dario Izzo */
  /*  E-mail     : dario.izzo@esa.int */
  /*  Affiliation: ESA / Advanced Concepts Team (ACT) */
  /*  Made more readible and optimized for speed by Rody P.S. Oldenhuis */
  /*  Code available in MGA.M on   http://www.esa.int/gsp/ACT/inf/op/globopt.htm
   */
  /*  last edited 12/Dec/2009 */
  /*  ADJUSTED FOR EML-COMPILATION 24/Dec/2009 */
  /*  initial values */
  bad = false;
  /*  work with non-dimensional units */
  r1 = (r1vec[0] * r1vec[0] + r1vec[1] * r1vec[1]) + r1vec[2] * r1vec[2];
  st.site = &emlrtRSI;
  if (r1 < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  r1 = muDoubleScalarSqrt(r1);
  r1vec[0] /= r1;
  r1vec[1] /= r1;
  r1vec[2] /= r1;
  V = muC / r1;
  st.site = &b_emlrtRSI;
  if (V < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  V = muDoubleScalarSqrt(V);
  T = r1 / V;
  tf = tf * 86400.0 / T;
  /*  also transform to seconds */
  /*  relevant geometry parameters (non dimensional) */
  d = r2vec[0] / r1;
  r2vec[0] = d;
  inn1 = d * d;
  d = r2vec[1] / r1;
  r2vec[1] = d;
  inn1 += d * d;
  d = r2vec[2] / r1;
  inn1 += d * d;
  st.site = &c_emlrtRSI;
  if (inn1 < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  mr2vec = muDoubleScalarSqrt(inn1);
  /*  make 100% sure it's in (-1 <= dth <= +1) */
  st.site = &d_emlrtRSI;
  dth = muDoubleScalarAcos(muDoubleScalarMax(
      -1.0,
      muDoubleScalarMin(
          1.0, ((r1vec[0] * r2vec[0] + r1vec[1] * r2vec[1]) + r1vec[2] * d) /
                   mr2vec)));
  /*  decide whether to use the left or right branch (for multi-revolution */
  /*  problems), and the long- or short way */
  leftbranch = muDoubleScalarSign(m);
  longway = muDoubleScalarSign(tf);
  m = muDoubleScalarAbs(m);
  tf = muDoubleScalarAbs(tf);
  if (longway < 0.0) {
    dth = 6.2831853071795862 - dth;
  }
  /*  derived quantities */
  st.site = &e_emlrtRSI;
  c = (mr2vec * mr2vec + 1.0) - 2.0 * mr2vec * muDoubleScalarCos(dth);
  st.site = &e_emlrtRSI;
  if (c < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  c = muDoubleScalarSqrt(c);
  /*  non-dimensional chord */
  s = ((mr2vec + 1.0) + c) / 2.0;
  /*  non-dimensional semi-perimeter */
  a_min = s / 2.0;
  /*  minimum energy ellipse semi major axis */
  st.site = &f_emlrtRSI;
  if (mr2vec < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  Lambda = muDoubleScalarSqrt(mr2vec) * muDoubleScalarCos(dth / 2.0) / s;
  /*  lambda parameter (from BATTIN's book) */
  crossprd[0] = r1vec[1] * d - r2vec[1] * r1vec[2];
  crossprd[1] = r2vec[0] * r1vec[2] - r1vec[0] * d;
  crossprd[2] = r1vec[0] * r2vec[1] - r2vec[0] * r1vec[1];
  /* % non-dimensional normal vectors */
  mcr = (crossprd[0] * crossprd[0] + crossprd[1] * crossprd[1]) +
        crossprd[2] * crossprd[2];
  st.site = &g_emlrtRSI;
  if (mcr < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
  }
  mcr = muDoubleScalarSqrt(mcr);
  /*  magnitues thereof */
  /*  unit vector thereof */
  /*  Initial values */
  /*  --------------------------------------------------------- */
  /*  ELMEX requires this variable to be declared OUTSIDE the IF-statement */
  st.site = &h_emlrtRSI;
  if (tf < 0.0) {
    emlrtErrorWithMessageIdR2018a(
        &st, &e_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
        "Coder:toolbox:ElFunDomainError", 3, 4, 3, "log");
  }
  logt = muDoubleScalarLog(tf);
  /*  avoid re-computing the same value */
  /*  single revolution (1 solution) */
  if (m == 0.0) {
    /*  initial values */
    inn1 = -0.5233;
    /*  first initial guess */
    inn2 = 0.5233;
    /*  second initial guess */
    st.site = &i_emlrtRSI;
    x1 = -0.740867916771357;
    /*  transformed first initial guess */
    st.site = &j_emlrtRSI;
    x2 = 0.42087903416051825;
    /*  transformed first second guess */
    /*  multiple revolutions (0, 1 or 2 solutions) */
    /*  the returned soltuion depends on the sign of [m] */
  } else {
    /*  select initial values */
    if (leftbranch < 0.0) {
      inn1 = -0.5234;
      /*  first initial guess, left branch */
      inn2 = -0.2234;
      /*  second initial guess, left branch */
    } else {
      inn1 = 0.7234;
      /*  first initial guess, right branch */
      inn2 = 0.5234;
      /*  second initial guess, right branch */
    }
    x1 = muDoubleScalarTan(inn1 * 3.1415926535897931 / 2.0);
    /*  transformed first initial guess */
    x2 = muDoubleScalarTan(inn2 * 3.1415926535897931 / 2.0);
    /*  transformed first second guess */
  }
  /*  since (inn1, inn2) < 0, initial estimate is always ellipse */
  x_tmp_tmp = s - c;
  x_tmp = x_tmp_tmp / 2.0;
  a = longway * 2.0;
  d1 = a_min / (1.0 - inn1 * inn1);
  aa_idx_0 = d1;
  bbeta[0] = x_tmp / d1;
  d1 = a_min / (1.0 - inn2 * inn2);
  bbeta[1] = x_tmp / d1;
  st.site = &k_emlrtRSI;
  b_sqrt(&st, bbeta);
  st.site = &k_emlrtRSI;
  b_asin(&st, bbeta);
  /*  make 100.4% sure it's in (-1 <= xx <= +1) */
  bbeta[0] *= a;
  aalfa[0] = inn1;
  bbeta[1] *= a;
  aalfa[1] = inn2;
  st.site = &l_emlrtRSI;
  b_acos(&st, aalfa);
  /*  evaluate the time of flight via Lagrange expression */
  inn1 = 2.0 * aalfa[0];
  aalfa[0] = inn1;
  xx_idx_0 = muDoubleScalarSin(inn1);
  dv[0] = aa_idx_0;
  inn1 = 2.0 * aalfa[1];
  dv[1] = d1;
  st.site = &m_emlrtRSI;
  b_sqrt(&st, dv);
  b_y1 = 6.2831853071795862 * m;
  aa_idx_0 =
      aa_idx_0 * dv[0] *
      (((aalfa[0] - xx_idx_0) - (bbeta[0] - muDoubleScalarSin(bbeta[0]))) +
       b_y1);
  inn2 = d1 * dv[1] *
         (((inn1 - muDoubleScalarSin(inn1)) -
           (bbeta[1] - muDoubleScalarSin(bbeta[1]))) +
          b_y1);
  /*  initial estimates for y */
  if (m == 0.0) {
    st.site = &n_emlrtRSI;
    if (aa_idx_0 < 0.0) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
          "Coder:toolbox:ElFunDomainError", 3, 4, 3, "log");
    }
    b_y1 = muDoubleScalarLog(aa_idx_0) - logt;
    st.site = &o_emlrtRSI;
    if (inn2 < 0.0) {
      emlrtErrorWithMessageIdR2018a(
          &st, &e_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
          "Coder:toolbox:ElFunDomainError", 3, 4, 3, "log");
    }
    aa_idx_0 = muDoubleScalarLog(inn2) - logt;
  } else {
    b_y1 = aa_idx_0 - tf;
    aa_idx_0 = inn2 - tf;
  }
  /*  Solve for x */
  /*  --------------------------------------------------------- */
  /*  Newton-Raphson iterations */
  /*  NOTE - the number of iterations will go to infinity in case */
  /*  m > 0  and there is no solution. Start the other routine in */
  /*  that case */
  xx_idx_0 = rtInf;
  iterations = 0;
  xnew = 0.0;
  exitg1 = false;
  while ((!exitg1) && (xx_idx_0 > 1.0E-14)) {
    /*  increment number of iterations */
    iterations++;
    /*  new x */
    xnew = (x1 * aa_idx_0 - b_y1 * x2) / (aa_idx_0 - b_y1);
    /*  copy-pasted code (for performance) */
    if (m == 0.0) {
      inn2 = muDoubleScalarExp(xnew) - 1.0;
    } else {
      inn2 = muDoubleScalarAtan(xnew) * 2.0 / 3.1415926535897931;
    }
    st.site = &p_emlrtRSI;
    a = a_min / (1.0 - inn2 * inn2);
    if (inn2 < 1.0) {
      /*  ellipse */
      st.site = &q_emlrtRSI;
      b_y1 = x_tmp / a;
      b_st.site = &q_emlrtRSI;
      if (b_y1 < 0.0) {
        emlrtErrorWithMessageIdR2018a(
            &b_st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
            "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
      }
      b_y1 = muDoubleScalarSqrt(b_y1);
      if ((b_y1 < -1.0) || (b_y1 > 1.0)) {
        emlrtErrorWithMessageIdR2018a(
            &st, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
            "Coder:toolbox:ElFunDomainError", 3, 4, 4, "asin");
      }
      b_y1 = muDoubleScalarAsin(b_y1);
      inn1 = longway * 2.0 * b_y1;
      /*  make 100.4% sure it's in (-1 <= xx <= +1) */
      st.site = &r_emlrtRSI;
      b_y1 = 2.0 * muDoubleScalarAcos(muDoubleScalarMax(-1.0, inn2));
    } else {
      /*  hyperbola */
      st.site = &s_emlrtRSI;
      b_acosh(&st, &inn2);
      b_y1 = 2.0 * inn2;
      d1 = x_tmp_tmp / (-2.0 * a);
      st.site = &t_emlrtRSI;
      if (d1 < 0.0) {
        emlrtErrorWithMessageIdR2018a(
            &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
            "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
      }
      d1 = muDoubleScalarSqrt(d1);
      b_asinh(&d1);
      inn1 = longway * 2.0 * d1;
    }
    /*  evaluate the time of flight via Lagrange expression */
    if (a > 0.0) {
      st.site = &u_emlrtRSI;
      if (a < 0.0) {
        emlrtErrorWithMessageIdR2018a(
            &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
            "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
      }
      inn1 = a * muDoubleScalarSqrt(a) *
             (((b_y1 - muDoubleScalarSin(b_y1)) -
               (inn1 - muDoubleScalarSin(inn1))) +
              6.2831853071795862 * m);
    } else {
      st.site = &v_emlrtRSI;
      if (-a < 0.0) {
        emlrtErrorWithMessageIdR2018a(
            &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
            "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
      }
      inn1 = -a * muDoubleScalarSqrt(-a) *
             ((muDoubleScalarSinh(b_y1) - b_y1) -
              (muDoubleScalarSinh(inn1) - inn1));
    }
    /*  new value of y */
    if (m == 0.0) {
      st.site = &w_emlrtRSI;
      if (inn1 < 0.0) {
        emlrtErrorWithMessageIdR2018a(
            &st, &e_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
            "Coder:toolbox:ElFunDomainError", 3, 4, 3, "log");
      }
      inn1 = muDoubleScalarLog(inn1) - logt;
    } else {
      inn1 -= tf;
    }
    /*  save previous and current values for the next iterarion */
    /*  (prevents getting stuck between two values) */
    x1 = x2;
    x2 = xnew;
    b_y1 = aa_idx_0;
    aa_idx_0 = inn1;
    /*  update error */
    xx_idx_0 = muDoubleScalarAbs(x1 - xnew);
    /*  escape clause */
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b((emlrtConstCTX)sp);
    }
    if (iterations > 15) {
      bad = true;
      exitg1 = true;
    }
  }
  /*  If the Newton-Raphson scheme failed, try to solve the problem */
  /*  with the other Lambert targeter. */
  if (bad) {
    real_T Vt2[3];
    /*  NOTE: use the original, UN-normalized quantities */
    crossprd[0] = r1vec[0] * r1;
    Vt2[0] = r2vec[0] * r1;
    crossprd[1] = r1vec[1] * r1;
    Vt2[1] = r2vec[1] * r1;
    crossprd[2] = r1vec[2] * r1;
    Vt2[2] = d * r1;
    st.site = &x_emlrtRSI;
    *exitflag = lambert_LancasterBlanchard(&st, crossprd, Vt2, longway * tf * T,
                                           leftbranch * m, muC, V1, V2,
                                           extremal_distances);
  } else {
    real_T Vt2[3];
    /*  convert converged value of x */
    if (m == 0.0) {
      inn2 = muDoubleScalarExp(xnew) - 1.0;
    } else {
      inn2 = muDoubleScalarAtan(xnew) * 2.0 / 3.1415926535897931;
    }
    /*     %{ */
    /*       The solution has been evaluated in terms of log(x+1) or
     * tan(x*pi/2), we */
    /*       now need the conic. As for transfer angles near to pi the Lagrange-
     */
    /*       coefficients technique goes singular (dg approaches a zero/zero
     * that is */
    /*       numerically bad) we here use a different technique for those cases.
     * When */
    /*       the transfer angle is exactly equal to pi, then the ih unit vector
     * is not */
    /*       determined. The remaining equations, though, are still valid. */
    /*     %} */
    /*  Solution for the semi-major axis */
    st.site = &y_emlrtRSI;
    a = a_min / (1.0 - inn2 * inn2);
    /*  Calculate psi */
    if (inn2 < 1.0) {
      /*  ellipse */
      st.site = &ab_emlrtRSI;
      b_y1 = x_tmp / a;
      b_st.site = &ab_emlrtRSI;
      if (b_y1 < 0.0) {
        emlrtErrorWithMessageIdR2018a(
            &b_st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
            "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
      }
      b_y1 = muDoubleScalarSqrt(b_y1);
      if ((b_y1 < -1.0) || (b_y1 > 1.0)) {
        emlrtErrorWithMessageIdR2018a(
            &st, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
            "Coder:toolbox:ElFunDomainError", 3, 4, 4, "asin");
      }
      b_y1 = muDoubleScalarAsin(b_y1);
      /*  make 100.4% sure it's in (-1 <= xx <= +1) */
      st.site = &bb_emlrtRSI;
      st.site = &cb_emlrtRSI;
      b_y1 = muDoubleScalarSin(
          (2.0 * muDoubleScalarAcos(muDoubleScalarMax(-1.0, inn2)) -
           longway * 2.0 * b_y1) /
          2.0);
      xx_idx_0 = 2.0 * a * (b_y1 * b_y1) / s;
      st.site = &db_emlrtRSI;
      if (xx_idx_0 < 0.0) {
        emlrtErrorWithMessageIdR2018a(
            &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
            "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
      }
      inn1 = muDoubleScalarSqrt(xx_idx_0);
    } else {
      /*  hyperbola */
      d1 = (c - s) / 2.0 / a;
      st.site = &eb_emlrtRSI;
      if (d1 < 0.0) {
        emlrtErrorWithMessageIdR2018a(
            &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
            "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
      }
      d1 = muDoubleScalarSqrt(d1);
      b_asinh(&d1);
      inn1 = inn2;
      st.site = &fb_emlrtRSI;
      b_acosh(&st, &inn1);
      st.site = &gb_emlrtRSI;
      b_y1 = muDoubleScalarSinh((2.0 * inn1 - longway * 2.0 * d1) / 2.0);
      xx_idx_0 = -2.0 * a * (b_y1 * b_y1) / s;
      st.site = &hb_emlrtRSI;
      if (xx_idx_0 < 0.0) {
        emlrtErrorWithMessageIdR2018a(
            &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
            "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
      }
      inn1 = muDoubleScalarSqrt(xx_idx_0);
    }
    /*  unit of the normalized normal vector */
    /*  unit vector for normalized [r2vec] */
    crossprd[0] = longway * (crossprd[0] / mcr);
    V2[0] = r2vec[0] / mr2vec;
    crossprd[1] = longway * (crossprd[1] / mcr);
    V2[1] = r2vec[1] / mr2vec;
    crossprd[2] = longway * (crossprd[2] / mcr);
    V2[2] = d / mr2vec;
    /*  cross-products */
    /*  don't use cross() (emlmex() would try to compile it, and this way it */
    /*  also does not create any additional overhead) */
    /*  radial and tangential directions for departure velocity */
    st.site = &ib_emlrtRSI;
    if (a_min < 0.0) {
      emlrtErrorWithMessageIdR2018a(
          &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
          "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
    }
    inn2 = 1.0 / inn1 / muDoubleScalarSqrt(a_min) *
           ((2.0 * Lambda * a_min - Lambda) - inn2 * inn1);
    st.site = &jb_emlrtRSI;
    b_y1 = muDoubleScalarSin(dth / 2.0);
    inn1 = mr2vec / a_min / xx_idx_0 * (b_y1 * b_y1);
    st.site = &jb_emlrtRSI;
    if (inn1 < 0.0) {
      emlrtErrorWithMessageIdR2018a(
          &st, &c_emlrtRTEI, "Coder:toolbox:ElFunDomainError",
          "Coder:toolbox:ElFunDomainError", 3, 4, 4, "sqrt");
    }
    inn1 = muDoubleScalarSqrt(inn1);
    /*  radial and tangential directions for arrival velocity */
    xx_idx_0 = inn1 / mr2vec;
    b_y1 = (inn1 - xx_idx_0) / muDoubleScalarTan(dth / 2.0) - inn2;
    /*  terminal velocities */
    V1[0] = inn1 * (crossprd[1] * r1vec[2] - r1vec[1] * crossprd[2]);
    V1[1] = inn1 * (r1vec[0] * crossprd[2] - crossprd[0] * r1vec[2]);
    V1[2] = inn1 * (crossprd[0] * r1vec[1] - r1vec[0] * crossprd[1]);
    Vt2[0] = xx_idx_0 * (crossprd[1] * V2[2] - V2[1] * crossprd[2]);
    Vt2[1] = xx_idx_0 * (V2[0] * crossprd[2] - crossprd[0] * V2[2]);
    Vt2[2] = xx_idx_0 * (crossprd[0] * V2[1] - V2[0] * crossprd[1]);
    /*  exitflag */
    *exitflag = 1.0;
    /*  (success) */
    /*  also compute minimum distance to central body */
    /*  NOTE: use un-transformed vectors again! */
    V1[0] = (inn2 * r1vec[0] + V1[0]) * V;
    V2[0] = (b_y1 * V2[0] + Vt2[0]) * V;
    crossprd[0] = r1vec[0] * r1;
    Vt2[0] = r2vec[0] * r1;
    V1[1] = (inn2 * r1vec[1] + V1[1]) * V;
    V2[1] = (b_y1 * V2[1] + Vt2[1]) * V;
    crossprd[1] = r1vec[1] * r1;
    Vt2[1] = r2vec[1] * r1;
    V1[2] = (inn2 * r1vec[2] + V1[2]) * V;
    V2[2] = (b_y1 * V2[2] + Vt2[2]) * V;
    crossprd[2] = r1vec[2] * r1;
    Vt2[2] = d * r1;
    st.site = &kb_emlrtRSI;
    minmax_distances(&st, crossprd, r1, Vt2, mr2vec * r1, dth, a * r1, V1, V2,
                     m, muC, extremal_distances);
  }
}

void sigmax_init(void)
{
  static const real_T dv[25] = {0.4,
                                0.2142857142857143,
                                0.0462962962962963,
                                0.006628787878787879,
                                0.00072115384615384609,
                                6.36574074074074E-5,
                                4.7414799253034548E-6,
                                3.0594063283208018E-7,
                                1.74283640925506E-8,
                                8.8924773311095776E-10,
                                4.1101115319865317E-11,
                                1.7367093848414581E-12,
                                6.7597672400414259E-14,
                                2.4391233866140258E-15,
                                8.2034116145380068E-17,
                                2.583771576869575E-18,
                                7.6523313279767163E-20,
                                2.138860629743989E-21,
                                5.6599594511655524E-23,
                                1.4221048338173659E-24,
                                3.4013984832723061E-26,
                                7.7625443047741554E-28,
                                1.693916882090479E-29,
                                3.54129500676686E-31,
                                7.1053361878044024E-33};
  memcpy(&an[0], &dv[0], 25U * sizeof(real_T));
}

/* End of code generation (lambert.c) */
