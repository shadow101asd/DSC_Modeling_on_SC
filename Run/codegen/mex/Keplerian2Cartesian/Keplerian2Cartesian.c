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
#include "Keplerian2Cartesian_emxutil.h"
#include "Keplerian2Cartesian_types.h"
#include "assertCompatibleDims.h"
#include "cos.h"
#include "div.h"
#include "norm.h"
#include "rt_nonfinite.h"
#include "sin.h"
#include "sqrt.h"
#include "mwmathutil.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = {
    25,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo b_emlrtRSI = {
    29,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo c_emlrtRSI = {
    31,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo d_emlrtRSI = {
    38,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo e_emlrtRSI = {
    45,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo f_emlrtRSI = {
    46,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo g_emlrtRSI = {
    47,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo h_emlrtRSI = {
    48,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo i_emlrtRSI = {
    49,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo j_emlrtRSI = {
    52,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo k_emlrtRSI = {
    53,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo l_emlrtRSI = {
    54,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo m_emlrtRSI = {
    55,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo n_emlrtRSI = {
    56,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo o_emlrtRSI = {
    57,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo p_emlrtRSI = {
    58,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo q_emlrtRSI = {
    59,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo r_emlrtRSI = {
    60,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo s_emlrtRSI = {
    83,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo t_emlrtRSI = {
    71,      /* lineNo */
    "power", /* fcnName */
    "/Applications/MATLAB_R2024a.app/toolbox/eml/lib/matlab/ops/power.m" /* pathName
                                                                          */
};

static emlrtRSInfo ab_emlrtRSI = {
    34,               /* lineNo */
    "rdivide_helper", /* fcnName */
    "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/"
    "rdivide_helper.m" /* pathName */
};

static emlrtRSInfo
    bb_emlrtRSI =
        {
            53,    /* lineNo */
            "div", /* fcnName */
            "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/"
            "div.m" /* pathName */
};

static emlrtRSInfo
    cb_emlrtRSI =
        {
            39,    /* lineNo */
            "cat", /* fcnName */
            "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/"
            "cat.m" /* pathName */
};

static emlrtRSInfo
    db_emlrtRSI =
        {
            113,        /* lineNo */
            "cat_impl", /* fcnName */
            "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/"
            "cat.m" /* pathName */
};

static emlrtRSInfo
    eb_emlrtRSI =
        {
            41,    /* lineNo */
            "cat", /* fcnName */
            "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/"
            "cat.m" /* pathName */
};

static emlrtRTEInfo
    emlrtRTEI =
        {
            288,                   /* lineNo */
            27,                    /* colNo */
            "check_non_axis_size", /* fName */
            "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/"
            "cat.m" /* pName */
};

static emlrtECInfo emlrtECI = {
    -1,                    /* nDims */
    80,                    /* lineNo */
    5,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo b_emlrtECI = {
    2,                     /* nDims */
    80,                    /* lineNo */
    17,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo c_emlrtECI = {
    2,                     /* nDims */
    80,                    /* lineNo */
    81,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo d_emlrtECI = {
    2,                     /* nDims */
    80,                    /* lineNo */
    49,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo e_emlrtECI = {
    -1,                    /* nDims */
    79,                    /* lineNo */
    5,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo f_emlrtECI = {
    2,                     /* nDims */
    79,                    /* lineNo */
    17,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo g_emlrtECI = {
    2,                     /* nDims */
    79,                    /* lineNo */
    81,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo h_emlrtECI = {
    2,                     /* nDims */
    79,                    /* lineNo */
    49,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo i_emlrtECI = {
    -1,                    /* nDims */
    78,                    /* lineNo */
    5,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo j_emlrtECI = {
    2,                     /* nDims */
    78,                    /* lineNo */
    17,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo k_emlrtECI = {
    2,                     /* nDims */
    78,                    /* lineNo */
    81,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo l_emlrtECI = {
    2,                     /* nDims */
    78,                    /* lineNo */
    49,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo m_emlrtECI = {
    -1,                    /* nDims */
    76,                    /* lineNo */
    5,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo n_emlrtECI = {
    2,                     /* nDims */
    76,                    /* lineNo */
    17,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo o_emlrtECI = {
    2,                     /* nDims */
    76,                    /* lineNo */
    81,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo p_emlrtECI = {
    2,                     /* nDims */
    76,                    /* lineNo */
    49,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo q_emlrtECI = {
    -1,                    /* nDims */
    75,                    /* lineNo */
    5,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo r_emlrtECI = {
    2,                     /* nDims */
    75,                    /* lineNo */
    17,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo s_emlrtECI = {
    2,                     /* nDims */
    75,                    /* lineNo */
    81,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo t_emlrtECI = {
    2,                     /* nDims */
    75,                    /* lineNo */
    49,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo u_emlrtECI = {
    -1,                    /* nDims */
    74,                    /* lineNo */
    5,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo v_emlrtECI = {
    2,                     /* nDims */
    74,                    /* lineNo */
    17,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo w_emlrtECI = {
    2,                     /* nDims */
    74,                    /* lineNo */
    81,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo x_emlrtECI = {
    2,                     /* nDims */
    74,                    /* lineNo */
    49,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo y_emlrtECI = {
    2,                     /* nDims */
    59,                    /* lineNo */
    19,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo ab_emlrtECI = {
    2,                     /* nDims */
    58,                    /* lineNo */
    19,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo bb_emlrtECI = {
    2,                     /* nDims */
    57,                    /* lineNo */
    19,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo cb_emlrtECI = {
    2,                     /* nDims */
    56,                    /* lineNo */
    19,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo db_emlrtECI = {
    2,                     /* nDims */
    56,                    /* lineNo */
    38,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo eb_emlrtECI = {
    2,                     /* nDims */
    55,                    /* lineNo */
    19,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo fb_emlrtECI = {
    2,                     /* nDims */
    55,                    /* lineNo */
    38,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo gb_emlrtECI = {
    2,                     /* nDims */
    54,                    /* lineNo */
    19,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo hb_emlrtECI = {
    2,                     /* nDims */
    53,                    /* lineNo */
    19,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo ib_emlrtECI = {
    2,                     /* nDims */
    53,                    /* lineNo */
    37,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo jb_emlrtECI = {
    2,                     /* nDims */
    52,                    /* lineNo */
    19,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo kb_emlrtECI = {
    2,                     /* nDims */
    52,                    /* lineNo */
    37,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo lb_emlrtECI = {
    2,                     /* nDims */
    49,                    /* lineNo */
    30,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo mb_emlrtECI = {
    2,                     /* nDims */
    49,                    /* lineNo */
    40,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo nb_emlrtECI = {
    2,                     /* nDims */
    49,                    /* lineNo */
    13,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo ob_emlrtECI = {
    2,                     /* nDims */
    48,                    /* lineNo */
    56,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo pb_emlrtECI = {
    2,                     /* nDims */
    48,                    /* lineNo */
    42,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo qb_emlrtECI = {
    2,                     /* nDims */
    48,                    /* lineNo */
    27,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo rb_emlrtECI = {
    2,                     /* nDims */
    48,                    /* lineNo */
    13,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtECInfo sb_emlrtECI = {
    2,                     /* nDims */
    25,                    /* lineNo */
    9,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo e_emlrtRTEI = {
    25,                    /* lineNo */
    15,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo f_emlrtRTEI = {
    25,                    /* lineNo */
    5,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo g_emlrtRTEI = {
    30,                    /* lineNo */
    9,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo h_emlrtRTEI = {
    40,                    /* lineNo */
    13,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo i_emlrtRTEI = {
    45,                    /* lineNo */
    5,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo j_emlrtRTEI = {
    33,                    /* lineNo */
    13,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo k_emlrtRTEI = {
    46,                    /* lineNo */
    5,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo l_emlrtRTEI = {
    47,                    /* lineNo */
    5,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo m_emlrtRTEI = {
    48,                    /* lineNo */
    13,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo n_emlrtRTEI = {
    48,                    /* lineNo */
    24,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo o_emlrtRTEI = {
    34,               /* lineNo */
    1,                /* colNo */
    "rdivide_helper", /* fName */
    "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/"
    "rdivide_helper.m" /* pName */
};

static emlrtRTEInfo p_emlrtRTEI = {
    48,                    /* lineNo */
    42,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo q_emlrtRTEI = {
    48,                    /* lineNo */
    5,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo r_emlrtRTEI = {
    49,                    /* lineNo */
    13,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo s_emlrtRTEI = {
    49,                    /* lineNo */
    40,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo t_emlrtRTEI = {
    49,                    /* lineNo */
    30,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo u_emlrtRTEI = {
    49,                    /* lineNo */
    5,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo v_emlrtRTEI = {
    52,                    /* lineNo */
    19,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo w_emlrtRTEI = {
    52,                    /* lineNo */
    28,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo x_emlrtRTEI = {
    52,                    /* lineNo */
    37,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo y_emlrtRTEI = {
    52,                    /* lineNo */
    46,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo ab_emlrtRTEI = {
    52,                    /* lineNo */
    54,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo bb_emlrtRTEI = {
    52,                    /* lineNo */
    5,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo cb_emlrtRTEI = {
    53,                    /* lineNo */
    19,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo db_emlrtRTEI = {
    53,                    /* lineNo */
    28,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo eb_emlrtRTEI = {
    53,                    /* lineNo */
    37,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo fb_emlrtRTEI = {
    53,                    /* lineNo */
    46,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo gb_emlrtRTEI = {
    53,                    /* lineNo */
    54,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo hb_emlrtRTEI = {
    53,                    /* lineNo */
    5,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo ib_emlrtRTEI = {
    54,                    /* lineNo */
    19,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo jb_emlrtRTEI = {
    54,                    /* lineNo */
    27,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo kb_emlrtRTEI = {
    54,                    /* lineNo */
    5,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo lb_emlrtRTEI = {
    55,                    /* lineNo */
    20,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo mb_emlrtRTEI = {
    55,                    /* lineNo */
    19,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo nb_emlrtRTEI = {
    55,                    /* lineNo */
    29,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo ob_emlrtRTEI = {
    55,                    /* lineNo */
    38,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo pb_emlrtRTEI = {
    55,                    /* lineNo */
    47,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo qb_emlrtRTEI = {
    55,                    /* lineNo */
    55,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo rb_emlrtRTEI = {
    55,                    /* lineNo */
    5,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo sb_emlrtRTEI = {
    56,                    /* lineNo */
    20,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo tb_emlrtRTEI = {
    56,                    /* lineNo */
    19,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo ub_emlrtRTEI = {
    56,                    /* lineNo */
    29,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo vb_emlrtRTEI = {
    56,                    /* lineNo */
    38,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo wb_emlrtRTEI = {
    56,                    /* lineNo */
    47,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo xb_emlrtRTEI = {
    56,                    /* lineNo */
    55,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo yb_emlrtRTEI = {
    56,                    /* lineNo */
    5,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo ac_emlrtRTEI = {
    57,                    /* lineNo */
    19,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo bc_emlrtRTEI = {
    57,                    /* lineNo */
    27,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo cc_emlrtRTEI = {
    57,                    /* lineNo */
    5,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo dc_emlrtRTEI = {
    58,                    /* lineNo */
    19,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo ec_emlrtRTEI = {
    58,                    /* lineNo */
    28,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo fc_emlrtRTEI = {
    58,                    /* lineNo */
    5,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo gc_emlrtRTEI = {
    59,                    /* lineNo */
    20,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo hc_emlrtRTEI = {
    59,                    /* lineNo */
    19,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo ic_emlrtRTEI = {
    59,                    /* lineNo */
    29,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo jc_emlrtRTEI = {
    59,                    /* lineNo */
    5,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo kc_emlrtRTEI = {
    60,                    /* lineNo */
    5,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo lc_emlrtRTEI = {
    71,                    /* lineNo */
    5,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo mc_emlrtRTEI = {
    72,                    /* lineNo */
    5,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo nc_emlrtRTEI = {
    74,                    /* lineNo */
    17,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo oc_emlrtRTEI = {
    74,                    /* lineNo */
    49,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo pc_emlrtRTEI = {
    74,                    /* lineNo */
    81,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo qc_emlrtRTEI = {
    75,                    /* lineNo */
    17,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo rc_emlrtRTEI = {
    75,                    /* lineNo */
    49,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo sc_emlrtRTEI = {
    75,                    /* lineNo */
    81,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo tc_emlrtRTEI = {
    76,                    /* lineNo */
    17,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo uc_emlrtRTEI = {
    76,                    /* lineNo */
    49,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo vc_emlrtRTEI = {
    76,                    /* lineNo */
    81,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo wc_emlrtRTEI = {
    78,                    /* lineNo */
    17,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo xc_emlrtRTEI = {
    78,                    /* lineNo */
    49,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo yc_emlrtRTEI = {
    78,                    /* lineNo */
    81,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo ad_emlrtRTEI = {
    79,                    /* lineNo */
    17,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo bd_emlrtRTEI = {
    79,                    /* lineNo */
    49,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo cd_emlrtRTEI = {
    79,                    /* lineNo */
    81,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo dd_emlrtRTEI = {
    80,                    /* lineNo */
    17,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo ed_emlrtRTEI = {
    80,                    /* lineNo */
    49,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo fd_emlrtRTEI = {
    80,                    /* lineNo */
    81,                    /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo gd_emlrtRTEI = {
    83,                    /* lineNo */
    5,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRTEInfo jd_emlrtRTEI = {
    25,                    /* lineNo */
    9,                     /* colNo */
    "Keplerian2Cartesian", /* fName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pName */
};

static emlrtRSInfo fb_emlrtRSI = {
    80,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo gb_emlrtRSI = {
    79,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo hb_emlrtRSI = {
    78,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo ib_emlrtRSI = {
    76,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo jb_emlrtRSI = {
    75,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo kb_emlrtRSI = {
    74,                    /* lineNo */
    "Keplerian2Cartesian", /* fcnName */
    "/Users/jpenot/Documents/MATLAB/SC/DSC_Modeling_on_SC/Run/"
    "Keplerian2Cartesian.m" /* pathName */
};

static emlrtRSInfo
    lb_emlrtRSI =
        {
            54,    /* lineNo */
            "div", /* fcnName */
            "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/"
            "div.m" /* pathName */
};

/* Function Declarations */
static void b_plus(const emlrtStack *sp, emxArray_real_T *in1,
                   const emxArray_real_T *in2);

static void b_times(const emlrtStack *sp, emxArray_real_T *in1,
                    const emxArray_real_T *in2, const emxArray_real_T *in3);

static void binary_expand_op(const emlrtStack *sp, emxArray_real_T *in1,
                             const emxArray_real_T *in2);

static void binary_expand_op_1(const emlrtStack *sp, emxArray_real_T *in1,
                               const emxArray_real_T *in2);

static void binary_expand_op_10(const emlrtStack *sp, emxArray_real_T *in1,
                                const emxArray_real_T *in2,
                                const emxArray_real_T *in3);

static void binary_expand_op_11(const emlrtStack *sp, emxArray_real_T *in1,
                                const emxArray_real_T *in2,
                                const emxArray_real_T *in3);

static void binary_expand_op_18(const emlrtStack *sp, emxArray_real_T *in1,
                                const emxArray_real_T *in2,
                                const emxArray_real_T *in3);

static void binary_expand_op_19(const emlrtStack *sp, emxArray_real_T *in1,
                                const emxArray_real_T *in2);

static void binary_expand_op_2(const emlrtStack *sp, emxArray_real_T *in1,
                               const emxArray_real_T *in2);

static void binary_expand_op_9(const emlrtStack *sp, emxArray_real_T *in1,
                               const emxArray_real_T *in2,
                               const emxArray_real_T *in3);

static void minus(const emlrtStack *sp, emxArray_real_T *in1,
                  const emxArray_real_T *in2);

static void plus(const emlrtStack *sp, emxArray_real_T *in1,
                 const emxArray_real_T *in2);

static void times(const emlrtStack *sp, emxArray_real_T *in1,
                  const emxArray_real_T *in2);

/* Function Definitions */
static void b_plus(const emlrtStack *sp, emxArray_real_T *in1,
                   const emxArray_real_T *in2)
{
  emxArray_real_T *b_in2;
  const real_T *in2_data;
  real_T *b_in2_data;
  real_T *in1_data;
  int32_T i;
  int32_T loop_ub;
  int32_T stride_0_1;
  int32_T stride_1_1;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real_T(sp, &b_in2, &s_emlrtRTEI);
  i = b_in2->size[0] * b_in2->size[1];
  b_in2->size[0] = 1;
  if (in1->size[1] == 1) {
    loop_ub = in2->size[1];
  } else {
    loop_ub = in1->size[1];
  }
  b_in2->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, b_in2, i, &s_emlrtRTEI);
  b_in2_data = b_in2->data;
  stride_0_1 = (in2->size[1] != 1);
  stride_1_1 = (in1->size[1] != 1);
  for (i = 0; i < loop_ub; i++) {
    b_in2_data[i] = in2_data[i * stride_0_1] + in1_data[i * stride_1_1];
  }
  i = in1->size[0] * in1->size[1];
  in1->size[0] = 1;
  in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, in1, i, &s_emlrtRTEI);
  in1_data = in1->data;
  for (i = 0; i < loop_ub; i++) {
    in1_data[i] = b_in2_data[i];
  }
  emxFree_real_T(sp, &b_in2);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void b_times(const emlrtStack *sp, emxArray_real_T *in1,
                    const emxArray_real_T *in2, const emxArray_real_T *in3)
{
  const real_T *in2_data;
  const real_T *in3_data;
  real_T *in1_data;
  int32_T i;
  int32_T loop_ub;
  int32_T stride_0_1;
  int32_T stride_1_1;
  in3_data = in3->data;
  in2_data = in2->data;
  i = in1->size[0] * in1->size[1];
  in1->size[0] = 1;
  emxEnsureCapacity_real_T(sp, in1, i, &m_emlrtRTEI);
  if (in3->size[1] == 1) {
    loop_ub = in2->size[1];
  } else {
    loop_ub = in3->size[1];
  }
  i = in1->size[0] * in1->size[1];
  in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, in1, i, &m_emlrtRTEI);
  in1_data = in1->data;
  stride_0_1 = (in2->size[1] != 1);
  stride_1_1 = (in3->size[1] != 1);
  for (i = 0; i < loop_ub; i++) {
    in1_data[i] = in2_data[i * stride_0_1] * in3_data[i * stride_1_1];
  }
}

static void binary_expand_op(const emlrtStack *sp, emxArray_real_T *in1,
                             const emxArray_real_T *in2)
{
  emxArray_real_T *b_in1;
  const real_T *in2_data;
  real_T *b_in1_data;
  real_T *in1_data;
  int32_T i;
  int32_T loop_ub;
  int32_T stride_0_1;
  int32_T stride_1_1;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real_T(sp, &b_in1, &fd_emlrtRTEI);
  i = b_in1->size[0] * b_in1->size[1];
  b_in1->size[0] = 1;
  if (in2->size[1] == 1) {
    loop_ub = in1->size[1];
  } else {
    loop_ub = in2->size[1];
  }
  b_in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, b_in1, i, &fd_emlrtRTEI);
  b_in1_data = b_in1->data;
  stride_0_1 = (in1->size[1] != 1);
  stride_1_1 = (in2->size[1] != 1);
  for (i = 0; i < loop_ub; i++) {
    b_in1_data[i] =
        in1_data[i * stride_0_1] * in2_data[3 * (i * stride_1_1) + 2];
  }
  i = in1->size[0] * in1->size[1];
  in1->size[0] = 1;
  in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, in1, i, &fd_emlrtRTEI);
  in1_data = in1->data;
  for (i = 0; i < loop_ub; i++) {
    in1_data[i] = b_in1_data[i];
  }
  emxFree_real_T(sp, &b_in1);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void binary_expand_op_1(const emlrtStack *sp, emxArray_real_T *in1,
                               const emxArray_real_T *in2)
{
  emxArray_real_T *b_in1;
  const real_T *in2_data;
  real_T *b_in1_data;
  real_T *in1_data;
  int32_T i;
  int32_T loop_ub;
  int32_T stride_0_1;
  int32_T stride_1_1;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real_T(sp, &b_in1, &ed_emlrtRTEI);
  i = b_in1->size[0] * b_in1->size[1];
  b_in1->size[0] = 1;
  if (in2->size[1] == 1) {
    loop_ub = in1->size[1];
  } else {
    loop_ub = in2->size[1];
  }
  b_in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, b_in1, i, &ed_emlrtRTEI);
  b_in1_data = b_in1->data;
  stride_0_1 = (in1->size[1] != 1);
  stride_1_1 = (in2->size[1] != 1);
  for (i = 0; i < loop_ub; i++) {
    b_in1_data[i] =
        in1_data[i * stride_0_1] * in2_data[3 * (i * stride_1_1) + 1];
  }
  i = in1->size[0] * in1->size[1];
  in1->size[0] = 1;
  in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, in1, i, &ed_emlrtRTEI);
  in1_data = in1->data;
  for (i = 0; i < loop_ub; i++) {
    in1_data[i] = b_in1_data[i];
  }
  emxFree_real_T(sp, &b_in1);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void binary_expand_op_10(const emlrtStack *sp, emxArray_real_T *in1,
                                const emxArray_real_T *in2,
                                const emxArray_real_T *in3)
{
  const real_T *in2_data;
  const real_T *in3_data;
  real_T *in1_data;
  int32_T i;
  int32_T loop_ub;
  int32_T stride_0_1;
  int32_T stride_1_1;
  in3_data = in3->data;
  in2_data = in2->data;
  i = in1->size[0] * in1->size[1];
  in1->size[0] = 1;
  emxEnsureCapacity_real_T(sp, in1, i, &uc_emlrtRTEI);
  if (in3->size[1] == 1) {
    loop_ub = in2->size[1];
  } else {
    loop_ub = in3->size[1];
  }
  i = in1->size[0] * in1->size[1];
  in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, in1, i, &uc_emlrtRTEI);
  in1_data = in1->data;
  stride_0_1 = (in2->size[1] != 1);
  stride_1_1 = (in3->size[1] != 1);
  for (i = 0; i < loop_ub; i++) {
    in1_data[i] = in2_data[i * stride_0_1] * in3_data[3 * (i * stride_1_1) + 1];
  }
}

static void binary_expand_op_11(const emlrtStack *sp, emxArray_real_T *in1,
                                const emxArray_real_T *in2,
                                const emxArray_real_T *in3)
{
  const real_T *in2_data;
  const real_T *in3_data;
  real_T *in1_data;
  int32_T i;
  int32_T loop_ub;
  int32_T stride_0_1;
  int32_T stride_1_1;
  in3_data = in3->data;
  in2_data = in2->data;
  i = in1->size[0] * in1->size[1];
  in1->size[0] = 1;
  emxEnsureCapacity_real_T(sp, in1, i, &tc_emlrtRTEI);
  if (in3->size[1] == 1) {
    loop_ub = in2->size[1];
  } else {
    loop_ub = in3->size[1];
  }
  i = in1->size[0] * in1->size[1];
  in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, in1, i, &tc_emlrtRTEI);
  in1_data = in1->data;
  stride_0_1 = (in2->size[1] != 1);
  stride_1_1 = (in3->size[1] != 1);
  for (i = 0; i < loop_ub; i++) {
    in1_data[i] = in2_data[i * stride_0_1] * in3_data[3 * (i * stride_1_1)];
  }
}

static void binary_expand_op_18(const emlrtStack *sp, emxArray_real_T *in1,
                                const emxArray_real_T *in2,
                                const emxArray_real_T *in3)
{
  const real_T *in2_data;
  const real_T *in3_data;
  real_T *in1_data;
  int32_T i;
  int32_T loop_ub;
  int32_T stride_0_1;
  int32_T stride_1_1;
  in3_data = in3->data;
  in2_data = in2->data;
  i = in1->size[0] * in1->size[1];
  in1->size[0] = 1;
  emxEnsureCapacity_real_T(sp, in1, i, &n_emlrtRTEI);
  if (in3->size[1] == 1) {
    loop_ub = in2->size[1];
  } else {
    loop_ub = in3->size[1];
  }
  i = in1->size[0] * in1->size[1];
  in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, in1, i, &n_emlrtRTEI);
  in1_data = in1->data;
  stride_0_1 = (in2->size[1] != 1);
  stride_1_1 = (in3->size[1] != 1);
  for (i = 0; i < loop_ub; i++) {
    in1_data[i] = in2_data[i * stride_0_1] * in3_data[i * stride_1_1] + 1.0;
  }
}

static void binary_expand_op_19(const emlrtStack *sp, emxArray_real_T *in1,
                                const emxArray_real_T *in2)
{
  emxArray_real_T *b_in2;
  const real_T *in2_data;
  real_T *b_in2_data;
  real_T *in1_data;
  int32_T i;
  int32_T loop_ub;
  int32_T stride_0_1;
  int32_T stride_1_1;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real_T(sp, &b_in2, &jd_emlrtRTEI);
  i = b_in2->size[0] * b_in2->size[1];
  b_in2->size[0] = 1;
  if (in1->size[1] == 1) {
    loop_ub = in2->size[1];
  } else {
    loop_ub = in1->size[1];
  }
  b_in2->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, b_in2, i, &jd_emlrtRTEI);
  b_in2_data = b_in2->data;
  stride_0_1 = (in2->size[1] != 1);
  stride_1_1 = (in1->size[1] != 1);
  for (i = 0; i < loop_ub; i++) {
    b_in2_data[i] = in2_data[i * stride_0_1] * (1.0 - in1_data[i * stride_1_1]);
  }
  i = in1->size[0] * in1->size[1];
  in1->size[0] = 1;
  in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, in1, i, &jd_emlrtRTEI);
  in1_data = in1->data;
  for (i = 0; i < loop_ub; i++) {
    in1_data[i] = b_in2_data[i];
  }
  emxFree_real_T(sp, &b_in2);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void binary_expand_op_2(const emlrtStack *sp, emxArray_real_T *in1,
                               const emxArray_real_T *in2)
{
  emxArray_real_T *b_in1;
  const real_T *in2_data;
  real_T *b_in1_data;
  real_T *in1_data;
  int32_T i;
  int32_T loop_ub;
  int32_T stride_0_1;
  int32_T stride_1_1;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real_T(sp, &b_in1, &dd_emlrtRTEI);
  i = b_in1->size[0] * b_in1->size[1];
  b_in1->size[0] = 1;
  if (in2->size[1] == 1) {
    loop_ub = in1->size[1];
  } else {
    loop_ub = in2->size[1];
  }
  b_in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, b_in1, i, &dd_emlrtRTEI);
  b_in1_data = b_in1->data;
  stride_0_1 = (in1->size[1] != 1);
  stride_1_1 = (in2->size[1] != 1);
  for (i = 0; i < loop_ub; i++) {
    b_in1_data[i] = in1_data[i * stride_0_1] * in2_data[3 * (i * stride_1_1)];
  }
  i = in1->size[0] * in1->size[1];
  in1->size[0] = 1;
  in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, in1, i, &dd_emlrtRTEI);
  in1_data = in1->data;
  for (i = 0; i < loop_ub; i++) {
    in1_data[i] = b_in1_data[i];
  }
  emxFree_real_T(sp, &b_in1);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void binary_expand_op_9(const emlrtStack *sp, emxArray_real_T *in1,
                               const emxArray_real_T *in2,
                               const emxArray_real_T *in3)
{
  const real_T *in2_data;
  const real_T *in3_data;
  real_T *in1_data;
  int32_T i;
  int32_T loop_ub;
  int32_T stride_0_1;
  int32_T stride_1_1;
  in3_data = in3->data;
  in2_data = in2->data;
  i = in1->size[0] * in1->size[1];
  in1->size[0] = 1;
  emxEnsureCapacity_real_T(sp, in1, i, &vc_emlrtRTEI);
  if (in3->size[1] == 1) {
    loop_ub = in2->size[1];
  } else {
    loop_ub = in3->size[1];
  }
  i = in1->size[0] * in1->size[1];
  in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, in1, i, &vc_emlrtRTEI);
  in1_data = in1->data;
  stride_0_1 = (in2->size[1] != 1);
  stride_1_1 = (in3->size[1] != 1);
  for (i = 0; i < loop_ub; i++) {
    in1_data[i] = in2_data[i * stride_0_1] * in3_data[3 * (i * stride_1_1) + 2];
  }
}

static void minus(const emlrtStack *sp, emxArray_real_T *in1,
                  const emxArray_real_T *in2)
{
  emxArray_real_T *b_in1;
  const real_T *in2_data;
  real_T *b_in1_data;
  real_T *in1_data;
  int32_T i;
  int32_T loop_ub;
  int32_T stride_0_1;
  int32_T stride_1_1;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real_T(sp, &b_in1, &mb_emlrtRTEI);
  i = b_in1->size[0] * b_in1->size[1];
  b_in1->size[0] = 1;
  if (in2->size[1] == 1) {
    loop_ub = in1->size[1];
  } else {
    loop_ub = in2->size[1];
  }
  b_in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, b_in1, i, &mb_emlrtRTEI);
  b_in1_data = b_in1->data;
  stride_0_1 = (in1->size[1] != 1);
  stride_1_1 = (in2->size[1] != 1);
  for (i = 0; i < loop_ub; i++) {
    b_in1_data[i] = in1_data[i * stride_0_1] - in2_data[i * stride_1_1];
  }
  i = in1->size[0] * in1->size[1];
  in1->size[0] = 1;
  in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, in1, i, &mb_emlrtRTEI);
  in1_data = in1->data;
  for (i = 0; i < loop_ub; i++) {
    in1_data[i] = b_in1_data[i];
  }
  emxFree_real_T(sp, &b_in1);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void plus(const emlrtStack *sp, emxArray_real_T *in1,
                 const emxArray_real_T *in2)
{
  emxArray_real_T *b_in1;
  const real_T *in2_data;
  real_T *b_in1_data;
  real_T *in1_data;
  int32_T i;
  int32_T loop_ub;
  int32_T stride_0_1;
  int32_T stride_1_1;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real_T(sp, &b_in1, &dd_emlrtRTEI);
  i = b_in1->size[0] * b_in1->size[1];
  b_in1->size[0] = 1;
  if (in2->size[1] == 1) {
    loop_ub = in1->size[1];
  } else {
    loop_ub = in2->size[1];
  }
  b_in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, b_in1, i, &dd_emlrtRTEI);
  b_in1_data = b_in1->data;
  stride_0_1 = (in1->size[1] != 1);
  stride_1_1 = (in2->size[1] != 1);
  for (i = 0; i < loop_ub; i++) {
    b_in1_data[i] = in1_data[i * stride_0_1] + in2_data[i * stride_1_1];
  }
  i = in1->size[0] * in1->size[1];
  in1->size[0] = 1;
  in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, in1, i, &dd_emlrtRTEI);
  in1_data = in1->data;
  for (i = 0; i < loop_ub; i++) {
    in1_data[i] = b_in1_data[i];
  }
  emxFree_real_T(sp, &b_in1);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

static void times(const emlrtStack *sp, emxArray_real_T *in1,
                  const emxArray_real_T *in2)
{
  emxArray_real_T *b_in1;
  const real_T *in2_data;
  real_T *b_in1_data;
  real_T *in1_data;
  int32_T i;
  int32_T loop_ub;
  int32_T stride_0_1;
  int32_T stride_1_1;
  in2_data = in2->data;
  in1_data = in1->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
  emxInit_real_T(sp, &b_in1, &hc_emlrtRTEI);
  i = b_in1->size[0] * b_in1->size[1];
  b_in1->size[0] = 1;
  if (in2->size[1] == 1) {
    loop_ub = in1->size[1];
  } else {
    loop_ub = in2->size[1];
  }
  b_in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, b_in1, i, &hc_emlrtRTEI);
  b_in1_data = b_in1->data;
  stride_0_1 = (in1->size[1] != 1);
  stride_1_1 = (in2->size[1] != 1);
  for (i = 0; i < loop_ub; i++) {
    b_in1_data[i] = in1_data[i * stride_0_1] * in2_data[i * stride_1_1];
  }
  i = in1->size[0] * in1->size[1];
  in1->size[0] = 1;
  in1->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, in1, i, &hc_emlrtRTEI);
  in1_data = in1->data;
  for (i = 0; i < loop_ub; i++) {
    in1_data[i] = b_in1_data[i];
  }
  emxFree_real_T(sp, &b_in1);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

void Keplerian2Cartesian(const emlrtStack *sp, const emxArray_real_T *a,
                         const emxArray_real_T *e, const emxArray_real_T *i,
                         emxArray_real_T *Om, emxArray_real_T *w,
                         const emxArray_real_T *f0, real_T mu,
                         emxArray_real_T *y)
{
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  emxArray_real_T *IJKOverPQW1;
  emxArray_real_T *IJKOverPQW2;
  emxArray_real_T *IJKOverPQW3;
  emxArray_real_T *IJKOverPQW4;
  emxArray_real_T *IJKOverPQW5;
  emxArray_real_T *IJKOverPQW7;
  emxArray_real_T *IJKOverPQW8;
  emxArray_real_T *cosnu;
  emxArray_real_T *p;
  emxArray_real_T *rPQW;
  emxArray_real_T *rPQW_temp;
  emxArray_real_T *rootmup;
  emxArray_real_T *sinnu;
  emxArray_real_T *vPQW;
  emxArray_real_T *vPQW_temp;
  const real_T *a_data;
  const real_T *e_data;
  const real_T *f0_data;
  const real_T *i_data;
  real_T varargin_1;
  real_T *IJKOverPQW1_data;
  real_T *IJKOverPQW2_data;
  real_T *IJKOverPQW3_data;
  real_T *IJKOverPQW4_data;
  real_T *IJKOverPQW5_data;
  real_T *IJKOverPQW7_data;
  real_T *IJKOverPQW8_data;
  real_T *Om_data;
  real_T *cosnu_data;
  real_T *p_data;
  real_T *rPQW_data;
  real_T *rootmup_data;
  real_T *sinnu_data;
  real_T *vPQW_data;
  real_T *w_data;
  int32_T input_sizes[2];
  int32_T b_i;
  int32_T b_loop_ub;
  int32_T c_loop_ub;
  int32_T d_loop_ub;
  int32_T e_loop_ub;
  int32_T f_loop_ub;
  int32_T loop_ub;
  int32_T loop_ub_tmp;
  int8_T input_sizes_idx_0;
  int8_T sizes_idx_0;
  boolean_T empty_non_axis_sizes;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  f0_data = f0->data;
  w_data = w->data;
  Om_data = Om->data;
  i_data = i->data;
  e_data = e->data;
  a_data = a->data;
  emlrtHeapReferenceStackEnterFcnR2012b((emlrtConstCTX)sp);
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
  st.site = &emlrtRSI;
  b_st.site = &t_emlrtRSI;
  emxInit_real_T(&b_st, &p, &f_emlrtRTEI);
  b_i = p->size[0] * p->size[1];
  p->size[0] = 1;
  loop_ub = e->size[1];
  p->size[1] = e->size[1];
  emxEnsureCapacity_real_T(&b_st, p, b_i, &e_emlrtRTEI);
  p_data = p->data;
  for (b_i = 0; b_i < loop_ub; b_i++) {
    varargin_1 = e_data[b_i];
    p_data[b_i] = varargin_1 * varargin_1;
  }
  b_loop_ub = a->size[1];
  if ((a->size[1] != e->size[1]) && ((a->size[1] != 1) && (e->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(a->size[1], e->size[1], &sb_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if (a->size[1] == p->size[1]) {
    c_loop_ub = a->size[1] - 1;
    b_i = p->size[0] * p->size[1];
    p->size[0] = 1;
    p->size[1] = a->size[1];
    emxEnsureCapacity_real_T(sp, p, b_i, &f_emlrtRTEI);
    p_data = p->data;
    for (b_i = 0; b_i <= c_loop_ub; b_i++) {
      p_data[b_i] = a_data[b_i] * (1.0 - p_data[b_i]);
    }
  } else {
    st.site = &emlrtRSI;
    binary_expand_op_19(&st, p, a);
    p_data = p->data;
  }
  /*  Special Cases */
  st.site = &b_emlrtRSI;
  varargin_1 = 1.0E-6 * muDoubleScalarSqrt(a->size[1]);
  if (b_norm(e) < varargin_1) {
    /*  if circular */
    b_i = w->size[0] * w->size[1];
    w->size[0] = 1;
    emxEnsureCapacity_real_T(sp, w, b_i, &g_emlrtRTEI);
    w_data = w->data;
    c_loop_ub = w->size[1] - 1;
    for (b_i = 0; b_i <= c_loop_ub; b_i++) {
      w_data[b_i] *= 0.0;
    }
    st.site = &c_emlrtRSI;
    if (b_norm(i) < varargin_1) {
      /*  if equatorial */
      /* f0 = acos(r(1, :)./R); % lambda_true */
      b_i = Om->size[0] * Om->size[1];
      Om->size[0] = 1;
      emxEnsureCapacity_real_T(sp, Om, b_i, &j_emlrtRTEI);
      Om_data = Om->data;
      c_loop_ub = Om->size[1] - 1;
      for (b_i = 0; b_i <= c_loop_ub; b_i++) {
        Om_data[b_i] *= 0.0;
      }
    } else {
      /* f0 = acos(dot(N,r)./(Nnorm.*R)); % u */
    }
  } else {
    st.site = &d_emlrtRSI;
    if (b_norm(i) < varargin_1) {
      /*  if equatorial */
      /* f0 = acos(r(1, :)./R); % lambda_true */
      b_i = Om->size[0] * Om->size[1];
      Om->size[0] = 1;
      emxEnsureCapacity_real_T(sp, Om, b_i, &h_emlrtRTEI);
      Om_data = Om->data;
      c_loop_ub = Om->size[1] - 1;
      for (b_i = 0; b_i <= c_loop_ub; b_i++) {
        Om_data[b_i] *= 0.0;
      }
    }
  }
  /*  Storing Variables for Computational Efficiency */
  emxInit_real_T(sp, &cosnu, &i_emlrtRTEI);
  b_i = cosnu->size[0] * cosnu->size[1];
  cosnu->size[0] = 1;
  c_loop_ub = f0->size[1];
  cosnu->size[1] = f0->size[1];
  emxEnsureCapacity_real_T(sp, cosnu, b_i, &i_emlrtRTEI);
  cosnu_data = cosnu->data;
  for (b_i = 0; b_i < c_loop_ub; b_i++) {
    cosnu_data[b_i] = f0_data[b_i];
  }
  st.site = &e_emlrtRSI;
  b_cos(&st, cosnu);
  cosnu_data = cosnu->data;
  emxInit_real_T(sp, &sinnu, &k_emlrtRTEI);
  b_i = sinnu->size[0] * sinnu->size[1];
  sinnu->size[0] = 1;
  sinnu->size[1] = f0->size[1];
  emxEnsureCapacity_real_T(sp, sinnu, b_i, &k_emlrtRTEI);
  sinnu_data = sinnu->data;
  for (b_i = 0; b_i < c_loop_ub; b_i++) {
    sinnu_data[b_i] = f0_data[b_i];
  }
  st.site = &f_emlrtRSI;
  b_sin(&st, sinnu);
  sinnu_data = sinnu->data;
  st.site = &g_emlrtRSI;
  b_st.site = &t_emlrtRSI;
  emxInit_real_T(sp, &rootmup, &l_emlrtRTEI);
  b_i = rootmup->size[0] * rootmup->size[1];
  rootmup->size[0] = 1;
  c_loop_ub = p->size[1];
  rootmup->size[1] = p->size[1];
  emxEnsureCapacity_real_T(sp, rootmup, b_i, &l_emlrtRTEI);
  rootmup_data = rootmup->data;
  for (b_i = 0; b_i < c_loop_ub; b_i++) {
    varargin_1 = p_data[b_i];
    rootmup_data[b_i] = mu * (1.0 / varargin_1);
  }
  st.site = &g_emlrtRSI;
  b_sqrt(&st, rootmup);
  rootmup_data = rootmup->data;
  if ((p->size[1] != cosnu->size[1]) &&
      ((p->size[1] != 1) && (cosnu->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(p->size[1], cosnu->size[1], &rb_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if ((e->size[1] != cosnu->size[1]) &&
      ((e->size[1] != 1) && (cosnu->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(e->size[1], cosnu->size[1], &qb_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if ((p->size[1] != sinnu->size[1]) &&
      ((p->size[1] != 1) && (sinnu->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(p->size[1], sinnu->size[1], &pb_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if ((e->size[1] != cosnu->size[1]) &&
      ((e->size[1] != 1) && (cosnu->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(e->size[1], cosnu->size[1], &ob_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  emxInit_real_T(sp, &IJKOverPQW7, &fc_emlrtRTEI);
  st.site = &h_emlrtRSI;
  if (p->size[1] == cosnu->size[1]) {
    b_i = IJKOverPQW7->size[0] * IJKOverPQW7->size[1];
    IJKOverPQW7->size[0] = 1;
    IJKOverPQW7->size[1] = p->size[1];
    emxEnsureCapacity_real_T(&st, IJKOverPQW7, b_i, &m_emlrtRTEI);
    IJKOverPQW7_data = IJKOverPQW7->data;
    for (b_i = 0; b_i < c_loop_ub; b_i++) {
      IJKOverPQW7_data[b_i] = p_data[b_i] * cosnu_data[b_i];
    }
  } else {
    b_st.site = &h_emlrtRSI;
    b_times(&b_st, IJKOverPQW7, p, cosnu);
  }
  emxInit_real_T(&st, &IJKOverPQW8, &jc_emlrtRTEI);
  if (e->size[1] == cosnu->size[1]) {
    b_i = IJKOverPQW8->size[0] * IJKOverPQW8->size[1];
    IJKOverPQW8->size[0] = 1;
    IJKOverPQW8->size[1] = e->size[1];
    emxEnsureCapacity_real_T(&st, IJKOverPQW8, b_i, &n_emlrtRTEI);
    IJKOverPQW8_data = IJKOverPQW8->data;
    for (b_i = 0; b_i < loop_ub; b_i++) {
      IJKOverPQW8_data[b_i] = e_data[b_i] * cosnu_data[b_i] + 1.0;
    }
  } else {
    b_st.site = &h_emlrtRSI;
    binary_expand_op_18(&b_st, IJKOverPQW8, e, cosnu);
    IJKOverPQW8_data = IJKOverPQW8->data;
  }
  b_st.site = &ab_emlrtRSI;
  c_st.site = &bb_emlrtRSI;
  assertCompatibleDims(&c_st, IJKOverPQW7, IJKOverPQW8);
  if (IJKOverPQW7->size[1] == IJKOverPQW8->size[1]) {
    loop_ub = IJKOverPQW7->size[1] - 1;
    b_i = IJKOverPQW7->size[0] * IJKOverPQW7->size[1];
    IJKOverPQW7->size[0] = 1;
    emxEnsureCapacity_real_T(&b_st, IJKOverPQW7, b_i, &o_emlrtRTEI);
    IJKOverPQW7_data = IJKOverPQW7->data;
    for (b_i = 0; b_i <= loop_ub; b_i++) {
      IJKOverPQW7_data[b_i] /= IJKOverPQW8_data[b_i];
    }
  } else {
    c_st.site = &lb_emlrtRSI;
    rdivide(&c_st, IJKOverPQW7, IJKOverPQW8);
    IJKOverPQW7_data = IJKOverPQW7->data;
  }
  st.site = &h_emlrtRSI;
  if (p->size[1] == sinnu->size[1]) {
    loop_ub = p->size[1] - 1;
    b_i = p->size[0] * p->size[1];
    p->size[0] = 1;
    emxEnsureCapacity_real_T(&st, p, b_i, &p_emlrtRTEI);
    p_data = p->data;
    for (b_i = 0; b_i <= loop_ub; b_i++) {
      p_data[b_i] *= sinnu_data[b_i];
    }
  } else {
    b_st.site = &h_emlrtRSI;
    times(&b_st, p, sinnu);
  }
  b_st.site = &ab_emlrtRSI;
  c_st.site = &bb_emlrtRSI;
  assertCompatibleDims(&c_st, p, IJKOverPQW8);
  if (p->size[1] == IJKOverPQW8->size[1]) {
    loop_ub = p->size[1] - 1;
    b_i = p->size[0] * p->size[1];
    p->size[0] = 1;
    emxEnsureCapacity_real_T(&b_st, p, b_i, &o_emlrtRTEI);
    p_data = p->data;
    for (b_i = 0; b_i <= loop_ub; b_i++) {
      p_data[b_i] /= IJKOverPQW8_data[b_i];
    }
  } else {
    c_st.site = &lb_emlrtRSI;
    rdivide(&c_st, p, IJKOverPQW8);
    p_data = p->data;
  }
  st.site = &h_emlrtRSI;
  b_st.site = &cb_emlrtRSI;
  c_st.site = &db_emlrtRSI;
  if (p->size[1] != IJKOverPQW7->size[1]) {
    emlrtErrorWithMessageIdR2018a(&c_st, &emlrtRTEI,
                                  "MATLAB:catenate:matrixDimensionMismatch",
                                  "MATLAB:catenate:matrixDimensionMismatch", 0);
  }
  if (a->size[1] != IJKOverPQW7->size[1]) {
    emlrtErrorWithMessageIdR2018a(&c_st, &emlrtRTEI,
                                  "MATLAB:catenate:matrixDimensionMismatch",
                                  "MATLAB:catenate:matrixDimensionMismatch", 0);
  }
  emxInit_real_T(&b_st, &rPQW, &q_emlrtRTEI);
  b_i = rPQW->size[0] * rPQW->size[1];
  rPQW->size[0] = 3;
  loop_ub = IJKOverPQW7->size[1];
  rPQW->size[1] = IJKOverPQW7->size[1];
  emxEnsureCapacity_real_T(&b_st, rPQW, b_i, &q_emlrtRTEI);
  rPQW_data = rPQW->data;
  for (b_i = 0; b_i < loop_ub; b_i++) {
    rPQW_data[3 * b_i] = IJKOverPQW7_data[b_i];
  }
  c_loop_ub = p->size[1];
  for (b_i = 0; b_i < c_loop_ub; b_i++) {
    rPQW_data[3 * b_i + 1] = p_data[b_i];
  }
  for (b_i = 0; b_i < b_loop_ub; b_i++) {
    rPQW_data[3 * b_i + 2] = 0.0;
  }
  b_i = p->size[0] * p->size[1];
  p->size[0] = 1;
  c_loop_ub = rootmup->size[1];
  p->size[1] = rootmup->size[1];
  emxEnsureCapacity_real_T(sp, p, b_i, &r_emlrtRTEI);
  p_data = p->data;
  for (b_i = 0; b_i < c_loop_ub; b_i++) {
    p_data[b_i] = -rootmup_data[b_i];
  }
  if ((rootmup->size[1] != sinnu->size[1]) &&
      ((rootmup->size[1] != 1) && (sinnu->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(rootmup->size[1], sinnu->size[1], &nb_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if ((e->size[1] != cosnu->size[1]) &&
      ((e->size[1] != 1) && (cosnu->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(e->size[1], cosnu->size[1], &mb_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if (e->size[1] == cosnu->size[1]) {
    c_loop_ub = e->size[1] - 1;
    b_i = cosnu->size[0] * cosnu->size[1];
    cosnu->size[0] = 1;
    cosnu->size[1] = e->size[1];
    emxEnsureCapacity_real_T(sp, cosnu, b_i, &s_emlrtRTEI);
    cosnu_data = cosnu->data;
    for (b_i = 0; b_i <= c_loop_ub; b_i++) {
      cosnu_data[b_i] += e_data[b_i];
    }
  } else {
    st.site = &i_emlrtRSI;
    b_plus(&st, cosnu, e);
    cosnu_data = cosnu->data;
  }
  if ((rootmup->size[1] != cosnu->size[1]) &&
      ((rootmup->size[1] != 1) && (cosnu->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(rootmup->size[1], cosnu->size[1], &lb_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  st.site = &i_emlrtRSI;
  if (p->size[1] == sinnu->size[1]) {
    c_loop_ub = p->size[1] - 1;
    b_i = p->size[0] * p->size[1];
    p->size[0] = 1;
    emxEnsureCapacity_real_T(&st, p, b_i, &r_emlrtRTEI);
    p_data = p->data;
    for (b_i = 0; b_i <= c_loop_ub; b_i++) {
      p_data[b_i] *= sinnu_data[b_i];
    }
  } else {
    b_st.site = &i_emlrtRSI;
    times(&b_st, p, sinnu);
    p_data = p->data;
  }
  if (rootmup->size[1] == cosnu->size[1]) {
    c_loop_ub = rootmup->size[1] - 1;
    b_i = rootmup->size[0] * rootmup->size[1];
    rootmup->size[0] = 1;
    emxEnsureCapacity_real_T(&st, rootmup, b_i, &t_emlrtRTEI);
    rootmup_data = rootmup->data;
    for (b_i = 0; b_i <= c_loop_ub; b_i++) {
      rootmup_data[b_i] *= cosnu_data[b_i];
    }
  } else {
    b_st.site = &i_emlrtRSI;
    times(&b_st, rootmup, cosnu);
    rootmup_data = rootmup->data;
  }
  b_st.site = &cb_emlrtRSI;
  c_st.site = &db_emlrtRSI;
  if (rootmup->size[1] != p->size[1]) {
    emlrtErrorWithMessageIdR2018a(&c_st, &emlrtRTEI,
                                  "MATLAB:catenate:matrixDimensionMismatch",
                                  "MATLAB:catenate:matrixDimensionMismatch", 0);
  }
  if (a->size[1] != p->size[1]) {
    emlrtErrorWithMessageIdR2018a(&c_st, &emlrtRTEI,
                                  "MATLAB:catenate:matrixDimensionMismatch",
                                  "MATLAB:catenate:matrixDimensionMismatch", 0);
  }
  emxInit_real_T(&b_st, &vPQW, &u_emlrtRTEI);
  b_i = vPQW->size[0] * vPQW->size[1];
  vPQW->size[0] = 3;
  c_loop_ub = p->size[1];
  vPQW->size[1] = p->size[1];
  emxEnsureCapacity_real_T(&b_st, vPQW, b_i, &u_emlrtRTEI);
  vPQW_data = vPQW->data;
  for (b_i = 0; b_i < c_loop_ub; b_i++) {
    vPQW_data[3 * b_i] = p_data[b_i];
  }
  d_loop_ub = rootmup->size[1];
  for (b_i = 0; b_i < d_loop_ub; b_i++) {
    vPQW_data[3 * b_i + 1] = rootmup_data[b_i];
  }
  for (b_i = 0; b_i < b_loop_ub; b_i++) {
    vPQW_data[3 * b_i + 2] = 0.0;
  }
  /*  Computing IJK/PQW */
  emxInit_real_T(sp, &IJKOverPQW1, &bb_emlrtRTEI);
  b_i = IJKOverPQW1->size[0] * IJKOverPQW1->size[1];
  IJKOverPQW1->size[0] = 1;
  b_loop_ub = Om->size[1];
  IJKOverPQW1->size[1] = Om->size[1];
  emxEnsureCapacity_real_T(sp, IJKOverPQW1, b_i, &v_emlrtRTEI);
  IJKOverPQW1_data = IJKOverPQW1->data;
  for (b_i = 0; b_i < b_loop_ub; b_i++) {
    IJKOverPQW1_data[b_i] = Om_data[b_i];
  }
  st.site = &j_emlrtRSI;
  b_cos(&st, IJKOverPQW1);
  b_i = sinnu->size[0] * sinnu->size[1];
  sinnu->size[0] = 1;
  d_loop_ub = w->size[1];
  sinnu->size[1] = w->size[1];
  emxEnsureCapacity_real_T(sp, sinnu, b_i, &w_emlrtRTEI);
  sinnu_data = sinnu->data;
  for (b_i = 0; b_i < d_loop_ub; b_i++) {
    sinnu_data[b_i] = w_data[b_i];
  }
  st.site = &j_emlrtRSI;
  b_cos(&st, sinnu);
  sinnu_data = sinnu->data;
  if ((IJKOverPQW1->size[1] != sinnu->size[1]) &&
      ((IJKOverPQW1->size[1] != 1) && (sinnu->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW1->size[1], sinnu->size[1],
                                &jb_emlrtECI, (emlrtConstCTX)sp);
  }
  b_i = rootmup->size[0] * rootmup->size[1];
  rootmup->size[0] = 1;
  rootmup->size[1] = Om->size[1];
  emxEnsureCapacity_real_T(sp, rootmup, b_i, &x_emlrtRTEI);
  rootmup_data = rootmup->data;
  for (b_i = 0; b_i < b_loop_ub; b_i++) {
    rootmup_data[b_i] = Om_data[b_i];
  }
  st.site = &j_emlrtRSI;
  b_sin(&st, rootmup);
  b_i = p->size[0] * p->size[1];
  p->size[0] = 1;
  p->size[1] = w->size[1];
  emxEnsureCapacity_real_T(sp, p, b_i, &y_emlrtRTEI);
  p_data = p->data;
  for (b_i = 0; b_i < d_loop_ub; b_i++) {
    p_data[b_i] = w_data[b_i];
  }
  st.site = &j_emlrtRSI;
  b_sin(&st, p);
  p_data = p->data;
  if ((rootmup->size[1] != p->size[1]) &&
      ((rootmup->size[1] != 1) && (p->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(rootmup->size[1], p->size[1], &kb_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if (rootmup->size[1] == p->size[1]) {
    e_loop_ub = rootmup->size[1] - 1;
    b_i = rootmup->size[0] * rootmup->size[1];
    rootmup->size[0] = 1;
    emxEnsureCapacity_real_T(sp, rootmup, b_i, &x_emlrtRTEI);
    rootmup_data = rootmup->data;
    for (b_i = 0; b_i <= e_loop_ub; b_i++) {
      rootmup_data[b_i] *= p_data[b_i];
    }
  } else {
    st.site = &j_emlrtRSI;
    times(&st, rootmup, p);
  }
  b_i = p->size[0] * p->size[1];
  p->size[0] = 1;
  e_loop_ub = i->size[1];
  p->size[1] = i->size[1];
  emxEnsureCapacity_real_T(sp, p, b_i, &ab_emlrtRTEI);
  p_data = p->data;
  for (b_i = 0; b_i < e_loop_ub; b_i++) {
    p_data[b_i] = i_data[b_i];
  }
  st.site = &j_emlrtRSI;
  b_cos(&st, p);
  p_data = p->data;
  if ((rootmup->size[1] != p->size[1]) &&
      ((rootmup->size[1] != 1) && (p->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(rootmup->size[1], p->size[1], &kb_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if (IJKOverPQW1->size[1] == sinnu->size[1]) {
    f_loop_ub = IJKOverPQW1->size[1] - 1;
    b_i = IJKOverPQW1->size[0] * IJKOverPQW1->size[1];
    IJKOverPQW1->size[0] = 1;
    emxEnsureCapacity_real_T(sp, IJKOverPQW1, b_i, &v_emlrtRTEI);
    IJKOverPQW1_data = IJKOverPQW1->data;
    for (b_i = 0; b_i <= f_loop_ub; b_i++) {
      IJKOverPQW1_data[b_i] *= sinnu_data[b_i];
    }
  } else {
    st.site = &j_emlrtRSI;
    times(&st, IJKOverPQW1, sinnu);
  }
  if (rootmup->size[1] == p->size[1]) {
    f_loop_ub = rootmup->size[1] - 1;
    b_i = rootmup->size[0] * rootmup->size[1];
    rootmup->size[0] = 1;
    emxEnsureCapacity_real_T(sp, rootmup, b_i, &x_emlrtRTEI);
    rootmup_data = rootmup->data;
    for (b_i = 0; b_i <= f_loop_ub; b_i++) {
      rootmup_data[b_i] *= p_data[b_i];
    }
  } else {
    st.site = &j_emlrtRSI;
    times(&st, rootmup, p);
    rootmup_data = rootmup->data;
  }
  if ((IJKOverPQW1->size[1] != rootmup->size[1]) &&
      ((IJKOverPQW1->size[1] != 1) && (rootmup->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW1->size[1], rootmup->size[1],
                                &jb_emlrtECI, (emlrtConstCTX)sp);
  }
  if (IJKOverPQW1->size[1] == rootmup->size[1]) {
    f_loop_ub = IJKOverPQW1->size[1] - 1;
    b_i = IJKOverPQW1->size[0] * IJKOverPQW1->size[1];
    IJKOverPQW1->size[0] = 1;
    emxEnsureCapacity_real_T(sp, IJKOverPQW1, b_i, &bb_emlrtRTEI);
    IJKOverPQW1_data = IJKOverPQW1->data;
    for (b_i = 0; b_i <= f_loop_ub; b_i++) {
      IJKOverPQW1_data[b_i] -= rootmup_data[b_i];
    }
  } else {
    st.site = &j_emlrtRSI;
    minus(&st, IJKOverPQW1, rootmup);
    IJKOverPQW1_data = IJKOverPQW1->data;
  }
  emxInit_real_T(sp, &IJKOverPQW2, &hb_emlrtRTEI);
  b_i = IJKOverPQW2->size[0] * IJKOverPQW2->size[1];
  IJKOverPQW2->size[0] = 1;
  IJKOverPQW2->size[1] = Om->size[1];
  emxEnsureCapacity_real_T(sp, IJKOverPQW2, b_i, &cb_emlrtRTEI);
  IJKOverPQW2_data = IJKOverPQW2->data;
  for (b_i = 0; b_i < b_loop_ub; b_i++) {
    IJKOverPQW2_data[b_i] = Om_data[b_i];
  }
  st.site = &k_emlrtRSI;
  b_sin(&st, IJKOverPQW2);
  b_i = sinnu->size[0] * sinnu->size[1];
  sinnu->size[0] = 1;
  sinnu->size[1] = w->size[1];
  emxEnsureCapacity_real_T(sp, sinnu, b_i, &db_emlrtRTEI);
  sinnu_data = sinnu->data;
  for (b_i = 0; b_i < d_loop_ub; b_i++) {
    sinnu_data[b_i] = w_data[b_i];
  }
  st.site = &k_emlrtRSI;
  b_cos(&st, sinnu);
  sinnu_data = sinnu->data;
  if ((IJKOverPQW2->size[1] != sinnu->size[1]) &&
      ((IJKOverPQW2->size[1] != 1) && (sinnu->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW2->size[1], sinnu->size[1],
                                &hb_emlrtECI, (emlrtConstCTX)sp);
  }
  b_i = rootmup->size[0] * rootmup->size[1];
  rootmup->size[0] = 1;
  rootmup->size[1] = Om->size[1];
  emxEnsureCapacity_real_T(sp, rootmup, b_i, &eb_emlrtRTEI);
  rootmup_data = rootmup->data;
  for (b_i = 0; b_i < b_loop_ub; b_i++) {
    rootmup_data[b_i] = Om_data[b_i];
  }
  st.site = &k_emlrtRSI;
  b_cos(&st, rootmup);
  b_i = p->size[0] * p->size[1];
  p->size[0] = 1;
  p->size[1] = w->size[1];
  emxEnsureCapacity_real_T(sp, p, b_i, &fb_emlrtRTEI);
  p_data = p->data;
  for (b_i = 0; b_i < d_loop_ub; b_i++) {
    p_data[b_i] = w_data[b_i];
  }
  st.site = &k_emlrtRSI;
  b_sin(&st, p);
  p_data = p->data;
  if ((rootmup->size[1] != p->size[1]) &&
      ((rootmup->size[1] != 1) && (p->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(rootmup->size[1], p->size[1], &ib_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if (rootmup->size[1] == p->size[1]) {
    f_loop_ub = rootmup->size[1] - 1;
    b_i = rootmup->size[0] * rootmup->size[1];
    rootmup->size[0] = 1;
    emxEnsureCapacity_real_T(sp, rootmup, b_i, &eb_emlrtRTEI);
    rootmup_data = rootmup->data;
    for (b_i = 0; b_i <= f_loop_ub; b_i++) {
      rootmup_data[b_i] *= p_data[b_i];
    }
  } else {
    st.site = &k_emlrtRSI;
    times(&st, rootmup, p);
  }
  b_i = p->size[0] * p->size[1];
  p->size[0] = 1;
  p->size[1] = i->size[1];
  emxEnsureCapacity_real_T(sp, p, b_i, &gb_emlrtRTEI);
  p_data = p->data;
  for (b_i = 0; b_i < e_loop_ub; b_i++) {
    p_data[b_i] = i_data[b_i];
  }
  st.site = &k_emlrtRSI;
  b_cos(&st, p);
  p_data = p->data;
  if ((rootmup->size[1] != p->size[1]) &&
      ((rootmup->size[1] != 1) && (p->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(rootmup->size[1], p->size[1], &ib_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if (IJKOverPQW2->size[1] == sinnu->size[1]) {
    f_loop_ub = IJKOverPQW2->size[1] - 1;
    b_i = IJKOverPQW2->size[0] * IJKOverPQW2->size[1];
    IJKOverPQW2->size[0] = 1;
    emxEnsureCapacity_real_T(sp, IJKOverPQW2, b_i, &cb_emlrtRTEI);
    IJKOverPQW2_data = IJKOverPQW2->data;
    for (b_i = 0; b_i <= f_loop_ub; b_i++) {
      IJKOverPQW2_data[b_i] *= sinnu_data[b_i];
    }
  } else {
    st.site = &k_emlrtRSI;
    times(&st, IJKOverPQW2, sinnu);
  }
  if (rootmup->size[1] == p->size[1]) {
    f_loop_ub = rootmup->size[1] - 1;
    b_i = rootmup->size[0] * rootmup->size[1];
    rootmup->size[0] = 1;
    emxEnsureCapacity_real_T(sp, rootmup, b_i, &eb_emlrtRTEI);
    rootmup_data = rootmup->data;
    for (b_i = 0; b_i <= f_loop_ub; b_i++) {
      rootmup_data[b_i] *= p_data[b_i];
    }
  } else {
    st.site = &k_emlrtRSI;
    times(&st, rootmup, p);
    rootmup_data = rootmup->data;
  }
  if ((IJKOverPQW2->size[1] != rootmup->size[1]) &&
      ((IJKOverPQW2->size[1] != 1) && (rootmup->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW2->size[1], rootmup->size[1],
                                &hb_emlrtECI, (emlrtConstCTX)sp);
  }
  if (IJKOverPQW2->size[1] == rootmup->size[1]) {
    f_loop_ub = IJKOverPQW2->size[1] - 1;
    b_i = IJKOverPQW2->size[0] * IJKOverPQW2->size[1];
    IJKOverPQW2->size[0] = 1;
    emxEnsureCapacity_real_T(sp, IJKOverPQW2, b_i, &hb_emlrtRTEI);
    IJKOverPQW2_data = IJKOverPQW2->data;
    for (b_i = 0; b_i <= f_loop_ub; b_i++) {
      IJKOverPQW2_data[b_i] += rootmup_data[b_i];
    }
  } else {
    st.site = &k_emlrtRSI;
    plus(&st, IJKOverPQW2, rootmup);
    IJKOverPQW2_data = IJKOverPQW2->data;
  }
  emxInit_real_T(sp, &IJKOverPQW3, &kb_emlrtRTEI);
  b_i = IJKOverPQW3->size[0] * IJKOverPQW3->size[1];
  IJKOverPQW3->size[0] = 1;
  IJKOverPQW3->size[1] = w->size[1];
  emxEnsureCapacity_real_T(sp, IJKOverPQW3, b_i, &ib_emlrtRTEI);
  IJKOverPQW3_data = IJKOverPQW3->data;
  for (b_i = 0; b_i < d_loop_ub; b_i++) {
    IJKOverPQW3_data[b_i] = w_data[b_i];
  }
  st.site = &l_emlrtRSI;
  b_sin(&st, IJKOverPQW3);
  b_i = sinnu->size[0] * sinnu->size[1];
  sinnu->size[0] = 1;
  sinnu->size[1] = i->size[1];
  emxEnsureCapacity_real_T(sp, sinnu, b_i, &jb_emlrtRTEI);
  sinnu_data = sinnu->data;
  for (b_i = 0; b_i < e_loop_ub; b_i++) {
    sinnu_data[b_i] = i_data[b_i];
  }
  st.site = &l_emlrtRSI;
  b_sin(&st, sinnu);
  sinnu_data = sinnu->data;
  if ((IJKOverPQW3->size[1] != sinnu->size[1]) &&
      ((IJKOverPQW3->size[1] != 1) && (sinnu->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW3->size[1], sinnu->size[1],
                                &gb_emlrtECI, (emlrtConstCTX)sp);
  }
  if (IJKOverPQW3->size[1] == sinnu->size[1]) {
    f_loop_ub = IJKOverPQW3->size[1] - 1;
    b_i = IJKOverPQW3->size[0] * IJKOverPQW3->size[1];
    IJKOverPQW3->size[0] = 1;
    emxEnsureCapacity_real_T(sp, IJKOverPQW3, b_i, &kb_emlrtRTEI);
    IJKOverPQW3_data = IJKOverPQW3->data;
    for (b_i = 0; b_i <= f_loop_ub; b_i++) {
      IJKOverPQW3_data[b_i] *= sinnu_data[b_i];
    }
  } else {
    st.site = &l_emlrtRSI;
    times(&st, IJKOverPQW3, sinnu);
    IJKOverPQW3_data = IJKOverPQW3->data;
  }
  emxInit_real_T(sp, &IJKOverPQW4, &rb_emlrtRTEI);
  b_i = IJKOverPQW4->size[0] * IJKOverPQW4->size[1];
  IJKOverPQW4->size[0] = 1;
  IJKOverPQW4->size[1] = Om->size[1];
  emxEnsureCapacity_real_T(sp, IJKOverPQW4, b_i, &lb_emlrtRTEI);
  IJKOverPQW4_data = IJKOverPQW4->data;
  for (b_i = 0; b_i < b_loop_ub; b_i++) {
    IJKOverPQW4_data[b_i] = Om_data[b_i];
  }
  st.site = &m_emlrtRSI;
  b_cos(&st, IJKOverPQW4);
  b_i = IJKOverPQW4->size[0] * IJKOverPQW4->size[1];
  IJKOverPQW4->size[0] = 1;
  emxEnsureCapacity_real_T(sp, IJKOverPQW4, b_i, &mb_emlrtRTEI);
  IJKOverPQW4_data = IJKOverPQW4->data;
  loop_ub_tmp = IJKOverPQW4->size[1] - 1;
  for (b_i = 0; b_i <= loop_ub_tmp; b_i++) {
    IJKOverPQW4_data[b_i] = -IJKOverPQW4_data[b_i];
  }
  b_i = sinnu->size[0] * sinnu->size[1];
  sinnu->size[0] = 1;
  sinnu->size[1] = w->size[1];
  emxEnsureCapacity_real_T(sp, sinnu, b_i, &nb_emlrtRTEI);
  sinnu_data = sinnu->data;
  for (b_i = 0; b_i < d_loop_ub; b_i++) {
    sinnu_data[b_i] = w_data[b_i];
  }
  st.site = &m_emlrtRSI;
  b_sin(&st, sinnu);
  sinnu_data = sinnu->data;
  if ((IJKOverPQW4->size[1] != sinnu->size[1]) &&
      ((IJKOverPQW4->size[1] != 1) && (sinnu->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW4->size[1], sinnu->size[1],
                                &eb_emlrtECI, (emlrtConstCTX)sp);
  }
  b_i = rootmup->size[0] * rootmup->size[1];
  rootmup->size[0] = 1;
  rootmup->size[1] = Om->size[1];
  emxEnsureCapacity_real_T(sp, rootmup, b_i, &ob_emlrtRTEI);
  rootmup_data = rootmup->data;
  for (b_i = 0; b_i < b_loop_ub; b_i++) {
    rootmup_data[b_i] = Om_data[b_i];
  }
  st.site = &m_emlrtRSI;
  b_sin(&st, rootmup);
  b_i = p->size[0] * p->size[1];
  p->size[0] = 1;
  p->size[1] = w->size[1];
  emxEnsureCapacity_real_T(sp, p, b_i, &pb_emlrtRTEI);
  p_data = p->data;
  for (b_i = 0; b_i < d_loop_ub; b_i++) {
    p_data[b_i] = w_data[b_i];
  }
  st.site = &m_emlrtRSI;
  b_cos(&st, p);
  p_data = p->data;
  if ((rootmup->size[1] != p->size[1]) &&
      ((rootmup->size[1] != 1) && (p->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(rootmup->size[1], p->size[1], &fb_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if (rootmup->size[1] == p->size[1]) {
    f_loop_ub = rootmup->size[1] - 1;
    b_i = rootmup->size[0] * rootmup->size[1];
    rootmup->size[0] = 1;
    emxEnsureCapacity_real_T(sp, rootmup, b_i, &ob_emlrtRTEI);
    rootmup_data = rootmup->data;
    for (b_i = 0; b_i <= f_loop_ub; b_i++) {
      rootmup_data[b_i] *= p_data[b_i];
    }
  } else {
    st.site = &m_emlrtRSI;
    times(&st, rootmup, p);
  }
  b_i = p->size[0] * p->size[1];
  p->size[0] = 1;
  p->size[1] = i->size[1];
  emxEnsureCapacity_real_T(sp, p, b_i, &qb_emlrtRTEI);
  p_data = p->data;
  for (b_i = 0; b_i < e_loop_ub; b_i++) {
    p_data[b_i] = i_data[b_i];
  }
  st.site = &m_emlrtRSI;
  b_cos(&st, p);
  p_data = p->data;
  if ((rootmup->size[1] != p->size[1]) &&
      ((rootmup->size[1] != 1) && (p->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(rootmup->size[1], p->size[1], &fb_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if (IJKOverPQW4->size[1] == sinnu->size[1]) {
    b_i = IJKOverPQW4->size[0] * IJKOverPQW4->size[1];
    IJKOverPQW4->size[0] = 1;
    emxEnsureCapacity_real_T(sp, IJKOverPQW4, b_i, &mb_emlrtRTEI);
    IJKOverPQW4_data = IJKOverPQW4->data;
    for (b_i = 0; b_i <= loop_ub_tmp; b_i++) {
      IJKOverPQW4_data[b_i] *= sinnu_data[b_i];
    }
  } else {
    st.site = &m_emlrtRSI;
    times(&st, IJKOverPQW4, sinnu);
  }
  if (rootmup->size[1] == p->size[1]) {
    f_loop_ub = rootmup->size[1] - 1;
    b_i = rootmup->size[0] * rootmup->size[1];
    rootmup->size[0] = 1;
    emxEnsureCapacity_real_T(sp, rootmup, b_i, &ob_emlrtRTEI);
    rootmup_data = rootmup->data;
    for (b_i = 0; b_i <= f_loop_ub; b_i++) {
      rootmup_data[b_i] *= p_data[b_i];
    }
  } else {
    st.site = &m_emlrtRSI;
    times(&st, rootmup, p);
    rootmup_data = rootmup->data;
  }
  if ((IJKOverPQW4->size[1] != rootmup->size[1]) &&
      ((IJKOverPQW4->size[1] != 1) && (rootmup->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW4->size[1], rootmup->size[1],
                                &eb_emlrtECI, (emlrtConstCTX)sp);
  }
  if (IJKOverPQW4->size[1] == rootmup->size[1]) {
    f_loop_ub = IJKOverPQW4->size[1] - 1;
    b_i = IJKOverPQW4->size[0] * IJKOverPQW4->size[1];
    IJKOverPQW4->size[0] = 1;
    emxEnsureCapacity_real_T(sp, IJKOverPQW4, b_i, &rb_emlrtRTEI);
    IJKOverPQW4_data = IJKOverPQW4->data;
    for (b_i = 0; b_i <= f_loop_ub; b_i++) {
      IJKOverPQW4_data[b_i] -= rootmup_data[b_i];
    }
  } else {
    st.site = &m_emlrtRSI;
    minus(&st, IJKOverPQW4, rootmup);
    IJKOverPQW4_data = IJKOverPQW4->data;
  }
  emxInit_real_T(sp, &IJKOverPQW5, &yb_emlrtRTEI);
  b_i = IJKOverPQW5->size[0] * IJKOverPQW5->size[1];
  IJKOverPQW5->size[0] = 1;
  IJKOverPQW5->size[1] = Om->size[1];
  emxEnsureCapacity_real_T(sp, IJKOverPQW5, b_i, &sb_emlrtRTEI);
  IJKOverPQW5_data = IJKOverPQW5->data;
  for (b_i = 0; b_i < b_loop_ub; b_i++) {
    IJKOverPQW5_data[b_i] = Om_data[b_i];
  }
  st.site = &n_emlrtRSI;
  b_sin(&st, IJKOverPQW5);
  b_i = IJKOverPQW5->size[0] * IJKOverPQW5->size[1];
  IJKOverPQW5->size[0] = 1;
  emxEnsureCapacity_real_T(sp, IJKOverPQW5, b_i, &tb_emlrtRTEI);
  IJKOverPQW5_data = IJKOverPQW5->data;
  loop_ub_tmp = IJKOverPQW5->size[1] - 1;
  for (b_i = 0; b_i <= loop_ub_tmp; b_i++) {
    IJKOverPQW5_data[b_i] = -IJKOverPQW5_data[b_i];
  }
  b_i = sinnu->size[0] * sinnu->size[1];
  sinnu->size[0] = 1;
  sinnu->size[1] = w->size[1];
  emxEnsureCapacity_real_T(sp, sinnu, b_i, &ub_emlrtRTEI);
  sinnu_data = sinnu->data;
  for (b_i = 0; b_i < d_loop_ub; b_i++) {
    sinnu_data[b_i] = w_data[b_i];
  }
  st.site = &n_emlrtRSI;
  b_sin(&st, sinnu);
  sinnu_data = sinnu->data;
  if ((IJKOverPQW5->size[1] != sinnu->size[1]) &&
      ((IJKOverPQW5->size[1] != 1) && (sinnu->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW5->size[1], sinnu->size[1],
                                &cb_emlrtECI, (emlrtConstCTX)sp);
  }
  b_i = rootmup->size[0] * rootmup->size[1];
  rootmup->size[0] = 1;
  rootmup->size[1] = Om->size[1];
  emxEnsureCapacity_real_T(sp, rootmup, b_i, &vb_emlrtRTEI);
  rootmup_data = rootmup->data;
  for (b_i = 0; b_i < b_loop_ub; b_i++) {
    rootmup_data[b_i] = Om_data[b_i];
  }
  st.site = &n_emlrtRSI;
  b_cos(&st, rootmup);
  b_i = p->size[0] * p->size[1];
  p->size[0] = 1;
  p->size[1] = w->size[1];
  emxEnsureCapacity_real_T(sp, p, b_i, &wb_emlrtRTEI);
  p_data = p->data;
  for (b_i = 0; b_i < d_loop_ub; b_i++) {
    p_data[b_i] = w_data[b_i];
  }
  st.site = &n_emlrtRSI;
  b_cos(&st, p);
  p_data = p->data;
  if ((rootmup->size[1] != p->size[1]) &&
      ((rootmup->size[1] != 1) && (p->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(rootmup->size[1], p->size[1], &db_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if (rootmup->size[1] == p->size[1]) {
    f_loop_ub = rootmup->size[1] - 1;
    b_i = rootmup->size[0] * rootmup->size[1];
    rootmup->size[0] = 1;
    emxEnsureCapacity_real_T(sp, rootmup, b_i, &vb_emlrtRTEI);
    rootmup_data = rootmup->data;
    for (b_i = 0; b_i <= f_loop_ub; b_i++) {
      rootmup_data[b_i] *= p_data[b_i];
    }
  } else {
    st.site = &n_emlrtRSI;
    times(&st, rootmup, p);
  }
  b_i = p->size[0] * p->size[1];
  p->size[0] = 1;
  p->size[1] = i->size[1];
  emxEnsureCapacity_real_T(sp, p, b_i, &xb_emlrtRTEI);
  p_data = p->data;
  for (b_i = 0; b_i < e_loop_ub; b_i++) {
    p_data[b_i] = i_data[b_i];
  }
  st.site = &n_emlrtRSI;
  b_cos(&st, p);
  p_data = p->data;
  if ((rootmup->size[1] != p->size[1]) &&
      ((rootmup->size[1] != 1) && (p->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(rootmup->size[1], p->size[1], &db_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if (IJKOverPQW5->size[1] == sinnu->size[1]) {
    b_i = IJKOverPQW5->size[0] * IJKOverPQW5->size[1];
    IJKOverPQW5->size[0] = 1;
    emxEnsureCapacity_real_T(sp, IJKOverPQW5, b_i, &tb_emlrtRTEI);
    IJKOverPQW5_data = IJKOverPQW5->data;
    for (b_i = 0; b_i <= loop_ub_tmp; b_i++) {
      IJKOverPQW5_data[b_i] *= sinnu_data[b_i];
    }
  } else {
    st.site = &n_emlrtRSI;
    times(&st, IJKOverPQW5, sinnu);
  }
  if (rootmup->size[1] == p->size[1]) {
    f_loop_ub = rootmup->size[1] - 1;
    b_i = rootmup->size[0] * rootmup->size[1];
    rootmup->size[0] = 1;
    emxEnsureCapacity_real_T(sp, rootmup, b_i, &vb_emlrtRTEI);
    rootmup_data = rootmup->data;
    for (b_i = 0; b_i <= f_loop_ub; b_i++) {
      rootmup_data[b_i] *= p_data[b_i];
    }
  } else {
    st.site = &n_emlrtRSI;
    times(&st, rootmup, p);
    rootmup_data = rootmup->data;
  }
  if ((IJKOverPQW5->size[1] != rootmup->size[1]) &&
      ((IJKOverPQW5->size[1] != 1) && (rootmup->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW5->size[1], rootmup->size[1],
                                &cb_emlrtECI, (emlrtConstCTX)sp);
  }
  if (IJKOverPQW5->size[1] == rootmup->size[1]) {
    f_loop_ub = IJKOverPQW5->size[1] - 1;
    b_i = IJKOverPQW5->size[0] * IJKOverPQW5->size[1];
    IJKOverPQW5->size[0] = 1;
    emxEnsureCapacity_real_T(sp, IJKOverPQW5, b_i, &yb_emlrtRTEI);
    IJKOverPQW5_data = IJKOverPQW5->data;
    for (b_i = 0; b_i <= f_loop_ub; b_i++) {
      IJKOverPQW5_data[b_i] += rootmup_data[b_i];
    }
  } else {
    st.site = &n_emlrtRSI;
    plus(&st, IJKOverPQW5, rootmup);
    IJKOverPQW5_data = IJKOverPQW5->data;
  }
  b_i = p->size[0] * p->size[1];
  p->size[0] = 1;
  p->size[1] = w->size[1];
  emxEnsureCapacity_real_T(sp, p, b_i, &ac_emlrtRTEI);
  p_data = p->data;
  for (b_i = 0; b_i < d_loop_ub; b_i++) {
    p_data[b_i] = w_data[b_i];
  }
  st.site = &o_emlrtRSI;
  b_cos(&st, p);
  b_i = sinnu->size[0] * sinnu->size[1];
  sinnu->size[0] = 1;
  sinnu->size[1] = i->size[1];
  emxEnsureCapacity_real_T(sp, sinnu, b_i, &bc_emlrtRTEI);
  sinnu_data = sinnu->data;
  for (b_i = 0; b_i < e_loop_ub; b_i++) {
    sinnu_data[b_i] = i_data[b_i];
  }
  st.site = &o_emlrtRSI;
  b_sin(&st, sinnu);
  sinnu_data = sinnu->data;
  if ((p->size[1] != sinnu->size[1]) &&
      ((p->size[1] != 1) && (sinnu->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(p->size[1], sinnu->size[1], &bb_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if (p->size[1] == sinnu->size[1]) {
    d_loop_ub = p->size[1] - 1;
    b_i = p->size[0] * p->size[1];
    p->size[0] = 1;
    emxEnsureCapacity_real_T(sp, p, b_i, &cc_emlrtRTEI);
    p_data = p->data;
    for (b_i = 0; b_i <= d_loop_ub; b_i++) {
      p_data[b_i] *= sinnu_data[b_i];
    }
  } else {
    st.site = &o_emlrtRSI;
    times(&st, p, sinnu);
    p_data = p->data;
  }
  b_i = IJKOverPQW7->size[0] * IJKOverPQW7->size[1];
  IJKOverPQW7->size[0] = 1;
  IJKOverPQW7->size[1] = Om->size[1];
  emxEnsureCapacity_real_T(sp, IJKOverPQW7, b_i, &dc_emlrtRTEI);
  IJKOverPQW7_data = IJKOverPQW7->data;
  for (b_i = 0; b_i < b_loop_ub; b_i++) {
    IJKOverPQW7_data[b_i] = Om_data[b_i];
  }
  st.site = &p_emlrtRSI;
  b_sin(&st, IJKOverPQW7);
  b_i = sinnu->size[0] * sinnu->size[1];
  sinnu->size[0] = 1;
  sinnu->size[1] = i->size[1];
  emxEnsureCapacity_real_T(sp, sinnu, b_i, &ec_emlrtRTEI);
  sinnu_data = sinnu->data;
  for (b_i = 0; b_i < e_loop_ub; b_i++) {
    sinnu_data[b_i] = i_data[b_i];
  }
  st.site = &p_emlrtRSI;
  b_sin(&st, sinnu);
  sinnu_data = sinnu->data;
  if ((IJKOverPQW7->size[1] != sinnu->size[1]) &&
      ((IJKOverPQW7->size[1] != 1) && (sinnu->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW7->size[1], sinnu->size[1],
                                &ab_emlrtECI, (emlrtConstCTX)sp);
  }
  if (IJKOverPQW7->size[1] == sinnu->size[1]) {
    d_loop_ub = IJKOverPQW7->size[1] - 1;
    b_i = IJKOverPQW7->size[0] * IJKOverPQW7->size[1];
    IJKOverPQW7->size[0] = 1;
    emxEnsureCapacity_real_T(sp, IJKOverPQW7, b_i, &fc_emlrtRTEI);
    IJKOverPQW7_data = IJKOverPQW7->data;
    for (b_i = 0; b_i <= d_loop_ub; b_i++) {
      IJKOverPQW7_data[b_i] *= sinnu_data[b_i];
    }
  } else {
    st.site = &p_emlrtRSI;
    times(&st, IJKOverPQW7, sinnu);
    IJKOverPQW7_data = IJKOverPQW7->data;
  }
  b_i = IJKOverPQW8->size[0] * IJKOverPQW8->size[1];
  IJKOverPQW8->size[0] = 1;
  IJKOverPQW8->size[1] = Om->size[1];
  emxEnsureCapacity_real_T(sp, IJKOverPQW8, b_i, &gc_emlrtRTEI);
  IJKOverPQW8_data = IJKOverPQW8->data;
  for (b_i = 0; b_i < b_loop_ub; b_i++) {
    IJKOverPQW8_data[b_i] = Om_data[b_i];
  }
  st.site = &q_emlrtRSI;
  b_cos(&st, IJKOverPQW8);
  b_i = IJKOverPQW8->size[0] * IJKOverPQW8->size[1];
  IJKOverPQW8->size[0] = 1;
  emxEnsureCapacity_real_T(sp, IJKOverPQW8, b_i, &hc_emlrtRTEI);
  IJKOverPQW8_data = IJKOverPQW8->data;
  loop_ub_tmp = IJKOverPQW8->size[1] - 1;
  for (b_i = 0; b_i <= loop_ub_tmp; b_i++) {
    IJKOverPQW8_data[b_i] = -IJKOverPQW8_data[b_i];
  }
  b_i = sinnu->size[0] * sinnu->size[1];
  sinnu->size[0] = 1;
  sinnu->size[1] = i->size[1];
  emxEnsureCapacity_real_T(sp, sinnu, b_i, &ic_emlrtRTEI);
  sinnu_data = sinnu->data;
  for (b_i = 0; b_i < e_loop_ub; b_i++) {
    sinnu_data[b_i] = i_data[b_i];
  }
  st.site = &q_emlrtRSI;
  b_sin(&st, sinnu);
  sinnu_data = sinnu->data;
  if ((IJKOverPQW8->size[1] != sinnu->size[1]) &&
      ((IJKOverPQW8->size[1] != 1) && (sinnu->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW8->size[1], sinnu->size[1],
                                &y_emlrtECI, (emlrtConstCTX)sp);
  }
  if (IJKOverPQW8->size[1] == sinnu->size[1]) {
    b_i = IJKOverPQW8->size[0] * IJKOverPQW8->size[1];
    IJKOverPQW8->size[0] = 1;
    emxEnsureCapacity_real_T(sp, IJKOverPQW8, b_i, &jc_emlrtRTEI);
    IJKOverPQW8_data = IJKOverPQW8->data;
    for (b_i = 0; b_i <= loop_ub_tmp; b_i++) {
      IJKOverPQW8_data[b_i] *= sinnu_data[b_i];
    }
  } else {
    st.site = &q_emlrtRSI;
    times(&st, IJKOverPQW8, sinnu);
    IJKOverPQW8_data = IJKOverPQW8->data;
  }
  b_i = cosnu->size[0] * cosnu->size[1];
  cosnu->size[0] = 1;
  cosnu->size[1] = i->size[1];
  emxEnsureCapacity_real_T(sp, cosnu, b_i, &kc_emlrtRTEI);
  cosnu_data = cosnu->data;
  for (b_i = 0; b_i < e_loop_ub; b_i++) {
    cosnu_data[b_i] = i_data[b_i];
  }
  st.site = &r_emlrtRSI;
  b_cos(&st, cosnu);
  cosnu_data = cosnu->data;
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
  emxInit_real_T(sp, &rPQW_temp, &lc_emlrtRTEI);
  b_i = rPQW_temp->size[0] * rPQW_temp->size[1];
  rPQW_temp->size[0] = 3;
  rPQW_temp->size[1] = loop_ub;
  emxEnsureCapacity_real_T(sp, rPQW_temp, b_i, &lc_emlrtRTEI);
  w_data = rPQW_temp->data;
  b_loop_ub = 3 * rPQW->size[1];
  for (b_i = 0; b_i < b_loop_ub; b_i++) {
    w_data[b_i] = rPQW_data[b_i];
  }
  emxInit_real_T(sp, &vPQW_temp, &mc_emlrtRTEI);
  b_i = vPQW_temp->size[0] * vPQW_temp->size[1];
  vPQW_temp->size[0] = 3;
  vPQW_temp->size[1] = c_loop_ub;
  emxEnsureCapacity_real_T(sp, vPQW_temp, b_i, &mc_emlrtRTEI);
  Om_data = vPQW_temp->data;
  b_loop_ub = 3 * vPQW->size[1];
  for (b_i = 0; b_i < b_loop_ub; b_i++) {
    Om_data[b_i] = vPQW_data[b_i];
  }
  b_loop_ub = IJKOverPQW1->size[1];
  if ((IJKOverPQW1->size[1] != loop_ub) &&
      ((IJKOverPQW1->size[1] != 1) && (loop_ub != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW1->size[1], loop_ub, &v_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  d_loop_ub = IJKOverPQW4->size[1];
  if ((IJKOverPQW4->size[1] != loop_ub) &&
      ((IJKOverPQW4->size[1] != 1) && (loop_ub != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW4->size[1], loop_ub, &x_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if (IJKOverPQW1->size[1] == rPQW->size[1]) {
    b_i = sinnu->size[0] * sinnu->size[1];
    sinnu->size[0] = 1;
    sinnu->size[1] = IJKOverPQW1->size[1];
    emxEnsureCapacity_real_T(sp, sinnu, b_i, &nc_emlrtRTEI);
    sinnu_data = sinnu->data;
    for (b_i = 0; b_i < b_loop_ub; b_i++) {
      sinnu_data[b_i] = IJKOverPQW1_data[b_i] * rPQW_data[3 * b_i];
    }
  } else {
    st.site = &kb_emlrtRSI;
    binary_expand_op_11(&st, sinnu, IJKOverPQW1, rPQW);
  }
  if (IJKOverPQW4->size[1] == rPQW->size[1]) {
    b_i = rootmup->size[0] * rootmup->size[1];
    rootmup->size[0] = 1;
    rootmup->size[1] = IJKOverPQW4->size[1];
    emxEnsureCapacity_real_T(sp, rootmup, b_i, &oc_emlrtRTEI);
    rootmup_data = rootmup->data;
    for (b_i = 0; b_i < d_loop_ub; b_i++) {
      rootmup_data[b_i] = IJKOverPQW4_data[b_i] * rPQW_data[3 * b_i + 1];
    }
  } else {
    st.site = &kb_emlrtRSI;
    binary_expand_op_10(&st, rootmup, IJKOverPQW4, rPQW);
    rootmup_data = rootmup->data;
  }
  if ((sinnu->size[1] != rootmup->size[1]) &&
      ((sinnu->size[1] != 1) && (rootmup->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(sinnu->size[1], rootmup->size[1], &v_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  b_loop_ub = IJKOverPQW7->size[1];
  if ((IJKOverPQW7->size[1] != rPQW->size[1]) &&
      ((IJKOverPQW7->size[1] != 1) && (rPQW->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW7->size[1], rPQW->size[1],
                                &w_emlrtECI, (emlrtConstCTX)sp);
  }
  if (sinnu->size[1] == rootmup->size[1]) {
    d_loop_ub = sinnu->size[1] - 1;
    b_i = sinnu->size[0] * sinnu->size[1];
    sinnu->size[0] = 1;
    emxEnsureCapacity_real_T(sp, sinnu, b_i, &nc_emlrtRTEI);
    sinnu_data = sinnu->data;
    for (b_i = 0; b_i <= d_loop_ub; b_i++) {
      sinnu_data[b_i] += rootmup_data[b_i];
    }
  } else {
    st.site = &kb_emlrtRSI;
    plus(&st, sinnu, rootmup);
  }
  if (IJKOverPQW7->size[1] == rPQW->size[1]) {
    b_i = rootmup->size[0] * rootmup->size[1];
    rootmup->size[0] = 1;
    rootmup->size[1] = IJKOverPQW7->size[1];
    emxEnsureCapacity_real_T(sp, rootmup, b_i, &pc_emlrtRTEI);
    rootmup_data = rootmup->data;
    for (b_i = 0; b_i < b_loop_ub; b_i++) {
      rootmup_data[b_i] = IJKOverPQW7_data[b_i] * rPQW_data[3 * b_i + 2];
    }
  } else {
    st.site = &kb_emlrtRSI;
    binary_expand_op_9(&st, rootmup, IJKOverPQW7, rPQW);
    rootmup_data = rootmup->data;
  }
  if ((sinnu->size[1] != rootmup->size[1]) &&
      ((sinnu->size[1] != 1) && (rootmup->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(sinnu->size[1], rootmup->size[1], &v_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if (sinnu->size[1] == rootmup->size[1]) {
    b_loop_ub = sinnu->size[1] - 1;
    b_i = sinnu->size[0] * sinnu->size[1];
    sinnu->size[0] = 1;
    emxEnsureCapacity_real_T(sp, sinnu, b_i, &nc_emlrtRTEI);
    sinnu_data = sinnu->data;
    for (b_i = 0; b_i <= b_loop_ub; b_i++) {
      sinnu_data[b_i] += rootmup_data[b_i];
    }
  } else {
    st.site = &kb_emlrtRSI;
    plus(&st, sinnu, rootmup);
    sinnu_data = sinnu->data;
  }
  input_sizes[0] = 1;
  f_loop_ub = rPQW->size[1];
  input_sizes[1] = rPQW->size[1];
  emlrtSubAssignSizeCheckR2012b(&input_sizes[0], 2, &sinnu->size[0], 2,
                                &u_emlrtECI, (emlrtCTX)sp);
  for (b_i = 0; b_i < f_loop_ub; b_i++) {
    rPQW_data[3 * b_i] = sinnu_data[b_i];
  }
  b_loop_ub = IJKOverPQW2->size[1];
  if ((IJKOverPQW2->size[1] != loop_ub) &&
      ((IJKOverPQW2->size[1] != 1) && (loop_ub != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW2->size[1], loop_ub, &r_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  d_loop_ub = IJKOverPQW5->size[1];
  if ((IJKOverPQW5->size[1] != loop_ub) &&
      ((IJKOverPQW5->size[1] != 1) && (loop_ub != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW5->size[1], loop_ub, &t_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if (IJKOverPQW2->size[1] == rPQW_temp->size[1]) {
    b_i = sinnu->size[0] * sinnu->size[1];
    sinnu->size[0] = 1;
    sinnu->size[1] = IJKOverPQW2->size[1];
    emxEnsureCapacity_real_T(sp, sinnu, b_i, &qc_emlrtRTEI);
    sinnu_data = sinnu->data;
    for (b_i = 0; b_i < b_loop_ub; b_i++) {
      sinnu_data[b_i] = IJKOverPQW2_data[b_i] * w_data[3 * b_i];
    }
  } else {
    st.site = &jb_emlrtRSI;
    binary_expand_op_11(&st, sinnu, IJKOverPQW2, rPQW_temp);
  }
  if (IJKOverPQW5->size[1] == rPQW_temp->size[1]) {
    b_i = rootmup->size[0] * rootmup->size[1];
    rootmup->size[0] = 1;
    rootmup->size[1] = IJKOverPQW5->size[1];
    emxEnsureCapacity_real_T(sp, rootmup, b_i, &rc_emlrtRTEI);
    rootmup_data = rootmup->data;
    for (b_i = 0; b_i < d_loop_ub; b_i++) {
      rootmup_data[b_i] = IJKOverPQW5_data[b_i] * w_data[3 * b_i + 1];
    }
  } else {
    st.site = &jb_emlrtRSI;
    binary_expand_op_10(&st, rootmup, IJKOverPQW5, rPQW_temp);
    rootmup_data = rootmup->data;
  }
  if ((sinnu->size[1] != rootmup->size[1]) &&
      ((sinnu->size[1] != 1) && (rootmup->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(sinnu->size[1], rootmup->size[1], &r_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  loop_ub = IJKOverPQW8->size[1];
  if ((IJKOverPQW8->size[1] != rPQW_temp->size[1]) &&
      ((IJKOverPQW8->size[1] != 1) && (rPQW_temp->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW8->size[1], rPQW_temp->size[1],
                                &s_emlrtECI, (emlrtConstCTX)sp);
  }
  if (sinnu->size[1] == rootmup->size[1]) {
    b_loop_ub = sinnu->size[1] - 1;
    b_i = sinnu->size[0] * sinnu->size[1];
    sinnu->size[0] = 1;
    emxEnsureCapacity_real_T(sp, sinnu, b_i, &qc_emlrtRTEI);
    sinnu_data = sinnu->data;
    for (b_i = 0; b_i <= b_loop_ub; b_i++) {
      sinnu_data[b_i] += rootmup_data[b_i];
    }
  } else {
    st.site = &jb_emlrtRSI;
    plus(&st, sinnu, rootmup);
  }
  if (IJKOverPQW8->size[1] == rPQW_temp->size[1]) {
    b_i = rootmup->size[0] * rootmup->size[1];
    rootmup->size[0] = 1;
    rootmup->size[1] = IJKOverPQW8->size[1];
    emxEnsureCapacity_real_T(sp, rootmup, b_i, &sc_emlrtRTEI);
    rootmup_data = rootmup->data;
    for (b_i = 0; b_i < loop_ub; b_i++) {
      rootmup_data[b_i] = IJKOverPQW8_data[b_i] * w_data[3 * b_i + 2];
    }
  } else {
    st.site = &jb_emlrtRSI;
    binary_expand_op_9(&st, rootmup, IJKOverPQW8, rPQW_temp);
    rootmup_data = rootmup->data;
  }
  if ((sinnu->size[1] != rootmup->size[1]) &&
      ((sinnu->size[1] != 1) && (rootmup->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(sinnu->size[1], rootmup->size[1], &r_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if (sinnu->size[1] == rootmup->size[1]) {
    loop_ub = sinnu->size[1] - 1;
    b_i = sinnu->size[0] * sinnu->size[1];
    sinnu->size[0] = 1;
    emxEnsureCapacity_real_T(sp, sinnu, b_i, &qc_emlrtRTEI);
    sinnu_data = sinnu->data;
    for (b_i = 0; b_i <= loop_ub; b_i++) {
      sinnu_data[b_i] += rootmup_data[b_i];
    }
  } else {
    st.site = &jb_emlrtRSI;
    plus(&st, sinnu, rootmup);
    sinnu_data = sinnu->data;
  }
  input_sizes[0] = 1;
  input_sizes[1] = rPQW->size[1];
  emlrtSubAssignSizeCheckR2012b(&input_sizes[0], 2, &sinnu->size[0], 2,
                                &q_emlrtECI, (emlrtCTX)sp);
  for (b_i = 0; b_i < f_loop_ub; b_i++) {
    rPQW_data[3 * b_i + 1] = sinnu_data[b_i];
  }
  loop_ub = IJKOverPQW3->size[1];
  if ((IJKOverPQW3->size[1] != rPQW_temp->size[1]) &&
      ((IJKOverPQW3->size[1] != 1) && (rPQW_temp->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW3->size[1], rPQW_temp->size[1],
                                &n_emlrtECI, (emlrtConstCTX)sp);
  }
  b_loop_ub = p->size[1];
  if ((p->size[1] != rPQW_temp->size[1]) &&
      ((p->size[1] != 1) && (rPQW_temp->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(p->size[1], rPQW_temp->size[1], &p_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if (IJKOverPQW3->size[1] == rPQW_temp->size[1]) {
    b_i = sinnu->size[0] * sinnu->size[1];
    sinnu->size[0] = 1;
    sinnu->size[1] = IJKOverPQW3->size[1];
    emxEnsureCapacity_real_T(sp, sinnu, b_i, &tc_emlrtRTEI);
    sinnu_data = sinnu->data;
    for (b_i = 0; b_i < loop_ub; b_i++) {
      sinnu_data[b_i] = IJKOverPQW3_data[b_i] * w_data[3 * b_i];
    }
  } else {
    st.site = &ib_emlrtRSI;
    binary_expand_op_11(&st, sinnu, IJKOverPQW3, rPQW_temp);
  }
  if (p->size[1] == rPQW_temp->size[1]) {
    b_i = rootmup->size[0] * rootmup->size[1];
    rootmup->size[0] = 1;
    rootmup->size[1] = p->size[1];
    emxEnsureCapacity_real_T(sp, rootmup, b_i, &uc_emlrtRTEI);
    rootmup_data = rootmup->data;
    for (b_i = 0; b_i < b_loop_ub; b_i++) {
      rootmup_data[b_i] = p_data[b_i] * w_data[3 * b_i + 1];
    }
  } else {
    st.site = &ib_emlrtRSI;
    binary_expand_op_10(&st, rootmup, p, rPQW_temp);
    rootmup_data = rootmup->data;
  }
  if ((sinnu->size[1] != rootmup->size[1]) &&
      ((sinnu->size[1] != 1) && (rootmup->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(sinnu->size[1], rootmup->size[1], &n_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  loop_ub = cosnu->size[1];
  if ((cosnu->size[1] != rPQW_temp->size[1]) &&
      ((cosnu->size[1] != 1) && (rPQW_temp->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(cosnu->size[1], rPQW_temp->size[1], &o_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if (sinnu->size[1] == rootmup->size[1]) {
    b_loop_ub = sinnu->size[1] - 1;
    b_i = sinnu->size[0] * sinnu->size[1];
    sinnu->size[0] = 1;
    emxEnsureCapacity_real_T(sp, sinnu, b_i, &tc_emlrtRTEI);
    sinnu_data = sinnu->data;
    for (b_i = 0; b_i <= b_loop_ub; b_i++) {
      sinnu_data[b_i] += rootmup_data[b_i];
    }
  } else {
    st.site = &ib_emlrtRSI;
    plus(&st, sinnu, rootmup);
  }
  if (cosnu->size[1] == rPQW_temp->size[1]) {
    b_i = rootmup->size[0] * rootmup->size[1];
    rootmup->size[0] = 1;
    rootmup->size[1] = cosnu->size[1];
    emxEnsureCapacity_real_T(sp, rootmup, b_i, &vc_emlrtRTEI);
    rootmup_data = rootmup->data;
    for (b_i = 0; b_i < loop_ub; b_i++) {
      rootmup_data[b_i] = cosnu_data[b_i] * w_data[3 * b_i + 2];
    }
  } else {
    st.site = &ib_emlrtRSI;
    binary_expand_op_9(&st, rootmup, cosnu, rPQW_temp);
    rootmup_data = rootmup->data;
  }
  emxFree_real_T(sp, &rPQW_temp);
  if ((sinnu->size[1] != rootmup->size[1]) &&
      ((sinnu->size[1] != 1) && (rootmup->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(sinnu->size[1], rootmup->size[1], &n_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if (sinnu->size[1] == rootmup->size[1]) {
    loop_ub = sinnu->size[1] - 1;
    b_i = sinnu->size[0] * sinnu->size[1];
    sinnu->size[0] = 1;
    emxEnsureCapacity_real_T(sp, sinnu, b_i, &tc_emlrtRTEI);
    sinnu_data = sinnu->data;
    for (b_i = 0; b_i <= loop_ub; b_i++) {
      sinnu_data[b_i] += rootmup_data[b_i];
    }
  } else {
    st.site = &ib_emlrtRSI;
    plus(&st, sinnu, rootmup);
    sinnu_data = sinnu->data;
  }
  emxFree_real_T(sp, &rootmup);
  input_sizes[0] = 1;
  input_sizes[1] = rPQW->size[1];
  emlrtSubAssignSizeCheckR2012b(&input_sizes[0], 2, &sinnu->size[0], 2,
                                &m_emlrtECI, (emlrtCTX)sp);
  for (b_i = 0; b_i < f_loop_ub; b_i++) {
    rPQW_data[3 * b_i + 2] = sinnu_data[b_i];
  }
  emxFree_real_T(sp, &sinnu);
  loop_ub = IJKOverPQW1->size[1];
  if ((IJKOverPQW1->size[1] != c_loop_ub) &&
      ((IJKOverPQW1->size[1] != 1) && (c_loop_ub != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW1->size[1], c_loop_ub, &j_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  b_loop_ub = IJKOverPQW4->size[1];
  if ((IJKOverPQW4->size[1] != c_loop_ub) &&
      ((IJKOverPQW4->size[1] != 1) && (c_loop_ub != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW4->size[1], c_loop_ub, &l_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if (IJKOverPQW1->size[1] == vPQW->size[1]) {
    b_i = IJKOverPQW1->size[0] * IJKOverPQW1->size[1];
    IJKOverPQW1->size[0] = 1;
    emxEnsureCapacity_real_T(sp, IJKOverPQW1, b_i, &wc_emlrtRTEI);
    IJKOverPQW1_data = IJKOverPQW1->data;
    for (b_i = 0; b_i < loop_ub; b_i++) {
      IJKOverPQW1_data[b_i] *= vPQW_data[3 * b_i];
    }
  } else {
    st.site = &hb_emlrtRSI;
    binary_expand_op_2(&st, IJKOverPQW1, vPQW);
  }
  if (IJKOverPQW4->size[1] == vPQW->size[1]) {
    b_i = IJKOverPQW4->size[0] * IJKOverPQW4->size[1];
    IJKOverPQW4->size[0] = 1;
    emxEnsureCapacity_real_T(sp, IJKOverPQW4, b_i, &xc_emlrtRTEI);
    IJKOverPQW4_data = IJKOverPQW4->data;
    for (b_i = 0; b_i < b_loop_ub; b_i++) {
      IJKOverPQW4_data[b_i] *= vPQW_data[3 * b_i + 1];
    }
  } else {
    st.site = &hb_emlrtRSI;
    binary_expand_op_1(&st, IJKOverPQW4, vPQW);
    IJKOverPQW4_data = IJKOverPQW4->data;
  }
  if ((IJKOverPQW1->size[1] != IJKOverPQW4->size[1]) &&
      ((IJKOverPQW1->size[1] != 1) && (IJKOverPQW4->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW1->size[1], IJKOverPQW4->size[1],
                                &j_emlrtECI, (emlrtConstCTX)sp);
  }
  loop_ub = IJKOverPQW7->size[1];
  if ((IJKOverPQW7->size[1] != vPQW->size[1]) &&
      ((IJKOverPQW7->size[1] != 1) && (vPQW->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW7->size[1], vPQW->size[1],
                                &k_emlrtECI, (emlrtConstCTX)sp);
  }
  if (IJKOverPQW1->size[1] == IJKOverPQW4->size[1]) {
    b_loop_ub = IJKOverPQW1->size[1] - 1;
    b_i = IJKOverPQW1->size[0] * IJKOverPQW1->size[1];
    IJKOverPQW1->size[0] = 1;
    emxEnsureCapacity_real_T(sp, IJKOverPQW1, b_i, &wc_emlrtRTEI);
    IJKOverPQW1_data = IJKOverPQW1->data;
    for (b_i = 0; b_i <= b_loop_ub; b_i++) {
      IJKOverPQW1_data[b_i] += IJKOverPQW4_data[b_i];
    }
  } else {
    st.site = &hb_emlrtRSI;
    plus(&st, IJKOverPQW1, IJKOverPQW4);
  }
  emxFree_real_T(sp, &IJKOverPQW4);
  if (IJKOverPQW7->size[1] == vPQW->size[1]) {
    b_i = IJKOverPQW7->size[0] * IJKOverPQW7->size[1];
    IJKOverPQW7->size[0] = 1;
    emxEnsureCapacity_real_T(sp, IJKOverPQW7, b_i, &yc_emlrtRTEI);
    IJKOverPQW7_data = IJKOverPQW7->data;
    for (b_i = 0; b_i < loop_ub; b_i++) {
      IJKOverPQW7_data[b_i] *= vPQW_data[3 * b_i + 2];
    }
  } else {
    st.site = &hb_emlrtRSI;
    binary_expand_op(&st, IJKOverPQW7, vPQW);
    IJKOverPQW7_data = IJKOverPQW7->data;
  }
  if ((IJKOverPQW1->size[1] != IJKOverPQW7->size[1]) &&
      ((IJKOverPQW1->size[1] != 1) && (IJKOverPQW7->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW1->size[1], IJKOverPQW7->size[1],
                                &j_emlrtECI, (emlrtConstCTX)sp);
  }
  if (IJKOverPQW1->size[1] == IJKOverPQW7->size[1]) {
    loop_ub = IJKOverPQW1->size[1] - 1;
    b_i = IJKOverPQW1->size[0] * IJKOverPQW1->size[1];
    IJKOverPQW1->size[0] = 1;
    emxEnsureCapacity_real_T(sp, IJKOverPQW1, b_i, &wc_emlrtRTEI);
    IJKOverPQW1_data = IJKOverPQW1->data;
    for (b_i = 0; b_i <= loop_ub; b_i++) {
      IJKOverPQW1_data[b_i] += IJKOverPQW7_data[b_i];
    }
  } else {
    st.site = &hb_emlrtRSI;
    plus(&st, IJKOverPQW1, IJKOverPQW7);
    IJKOverPQW1_data = IJKOverPQW1->data;
  }
  emxFree_real_T(sp, &IJKOverPQW7);
  input_sizes[0] = 1;
  f_loop_ub = vPQW->size[1];
  input_sizes[1] = vPQW->size[1];
  emlrtSubAssignSizeCheckR2012b(&input_sizes[0], 2, &IJKOverPQW1->size[0], 2,
                                &i_emlrtECI, (emlrtCTX)sp);
  for (b_i = 0; b_i < f_loop_ub; b_i++) {
    vPQW_data[3 * b_i] = IJKOverPQW1_data[b_i];
  }
  emxFree_real_T(sp, &IJKOverPQW1);
  loop_ub = IJKOverPQW2->size[1];
  if ((IJKOverPQW2->size[1] != c_loop_ub) &&
      ((IJKOverPQW2->size[1] != 1) && (c_loop_ub != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW2->size[1], c_loop_ub, &f_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  b_loop_ub = IJKOverPQW5->size[1];
  if ((IJKOverPQW5->size[1] != c_loop_ub) &&
      ((IJKOverPQW5->size[1] != 1) && (c_loop_ub != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW5->size[1], c_loop_ub, &h_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if (IJKOverPQW2->size[1] == vPQW_temp->size[1]) {
    b_i = IJKOverPQW2->size[0] * IJKOverPQW2->size[1];
    IJKOverPQW2->size[0] = 1;
    emxEnsureCapacity_real_T(sp, IJKOverPQW2, b_i, &ad_emlrtRTEI);
    IJKOverPQW2_data = IJKOverPQW2->data;
    for (b_i = 0; b_i < loop_ub; b_i++) {
      IJKOverPQW2_data[b_i] *= Om_data[3 * b_i];
    }
  } else {
    st.site = &gb_emlrtRSI;
    binary_expand_op_2(&st, IJKOverPQW2, vPQW_temp);
  }
  if (IJKOverPQW5->size[1] == vPQW_temp->size[1]) {
    b_i = IJKOverPQW5->size[0] * IJKOverPQW5->size[1];
    IJKOverPQW5->size[0] = 1;
    emxEnsureCapacity_real_T(sp, IJKOverPQW5, b_i, &bd_emlrtRTEI);
    IJKOverPQW5_data = IJKOverPQW5->data;
    for (b_i = 0; b_i < b_loop_ub; b_i++) {
      IJKOverPQW5_data[b_i] *= Om_data[3 * b_i + 1];
    }
  } else {
    st.site = &gb_emlrtRSI;
    binary_expand_op_1(&st, IJKOverPQW5, vPQW_temp);
    IJKOverPQW5_data = IJKOverPQW5->data;
  }
  if ((IJKOverPQW2->size[1] != IJKOverPQW5->size[1]) &&
      ((IJKOverPQW2->size[1] != 1) && (IJKOverPQW5->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW2->size[1], IJKOverPQW5->size[1],
                                &f_emlrtECI, (emlrtConstCTX)sp);
  }
  loop_ub = IJKOverPQW8->size[1];
  if ((IJKOverPQW8->size[1] != vPQW_temp->size[1]) &&
      ((IJKOverPQW8->size[1] != 1) && (vPQW_temp->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW8->size[1], vPQW_temp->size[1],
                                &g_emlrtECI, (emlrtConstCTX)sp);
  }
  if (IJKOverPQW2->size[1] == IJKOverPQW5->size[1]) {
    b_loop_ub = IJKOverPQW2->size[1] - 1;
    b_i = IJKOverPQW2->size[0] * IJKOverPQW2->size[1];
    IJKOverPQW2->size[0] = 1;
    emxEnsureCapacity_real_T(sp, IJKOverPQW2, b_i, &ad_emlrtRTEI);
    IJKOverPQW2_data = IJKOverPQW2->data;
    for (b_i = 0; b_i <= b_loop_ub; b_i++) {
      IJKOverPQW2_data[b_i] += IJKOverPQW5_data[b_i];
    }
  } else {
    st.site = &gb_emlrtRSI;
    plus(&st, IJKOverPQW2, IJKOverPQW5);
  }
  emxFree_real_T(sp, &IJKOverPQW5);
  if (IJKOverPQW8->size[1] == vPQW_temp->size[1]) {
    b_i = IJKOverPQW8->size[0] * IJKOverPQW8->size[1];
    IJKOverPQW8->size[0] = 1;
    emxEnsureCapacity_real_T(sp, IJKOverPQW8, b_i, &cd_emlrtRTEI);
    IJKOverPQW8_data = IJKOverPQW8->data;
    for (b_i = 0; b_i < loop_ub; b_i++) {
      IJKOverPQW8_data[b_i] *= Om_data[3 * b_i + 2];
    }
  } else {
    st.site = &gb_emlrtRSI;
    binary_expand_op(&st, IJKOverPQW8, vPQW_temp);
    IJKOverPQW8_data = IJKOverPQW8->data;
  }
  if ((IJKOverPQW2->size[1] != IJKOverPQW8->size[1]) &&
      ((IJKOverPQW2->size[1] != 1) && (IJKOverPQW8->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW2->size[1], IJKOverPQW8->size[1],
                                &f_emlrtECI, (emlrtConstCTX)sp);
  }
  if (IJKOverPQW2->size[1] == IJKOverPQW8->size[1]) {
    loop_ub = IJKOverPQW2->size[1] - 1;
    b_i = IJKOverPQW2->size[0] * IJKOverPQW2->size[1];
    IJKOverPQW2->size[0] = 1;
    emxEnsureCapacity_real_T(sp, IJKOverPQW2, b_i, &ad_emlrtRTEI);
    IJKOverPQW2_data = IJKOverPQW2->data;
    for (b_i = 0; b_i <= loop_ub; b_i++) {
      IJKOverPQW2_data[b_i] += IJKOverPQW8_data[b_i];
    }
  } else {
    st.site = &gb_emlrtRSI;
    plus(&st, IJKOverPQW2, IJKOverPQW8);
    IJKOverPQW2_data = IJKOverPQW2->data;
  }
  emxFree_real_T(sp, &IJKOverPQW8);
  input_sizes[0] = 1;
  input_sizes[1] = vPQW->size[1];
  emlrtSubAssignSizeCheckR2012b(&input_sizes[0], 2, &IJKOverPQW2->size[0], 2,
                                &e_emlrtECI, (emlrtCTX)sp);
  for (b_i = 0; b_i < f_loop_ub; b_i++) {
    vPQW_data[3 * b_i + 1] = IJKOverPQW2_data[b_i];
  }
  emxFree_real_T(sp, &IJKOverPQW2);
  loop_ub = IJKOverPQW3->size[1];
  if ((IJKOverPQW3->size[1] != vPQW_temp->size[1]) &&
      ((IJKOverPQW3->size[1] != 1) && (vPQW_temp->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW3->size[1], vPQW_temp->size[1],
                                &b_emlrtECI, (emlrtConstCTX)sp);
  }
  b_loop_ub = p->size[1];
  if ((p->size[1] != vPQW_temp->size[1]) &&
      ((p->size[1] != 1) && (vPQW_temp->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(p->size[1], vPQW_temp->size[1], &d_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if (IJKOverPQW3->size[1] == vPQW_temp->size[1]) {
    b_i = IJKOverPQW3->size[0] * IJKOverPQW3->size[1];
    IJKOverPQW3->size[0] = 1;
    emxEnsureCapacity_real_T(sp, IJKOverPQW3, b_i, &dd_emlrtRTEI);
    IJKOverPQW3_data = IJKOverPQW3->data;
    for (b_i = 0; b_i < loop_ub; b_i++) {
      IJKOverPQW3_data[b_i] *= Om_data[3 * b_i];
    }
  } else {
    st.site = &fb_emlrtRSI;
    binary_expand_op_2(&st, IJKOverPQW3, vPQW_temp);
  }
  if (p->size[1] == vPQW_temp->size[1]) {
    b_i = p->size[0] * p->size[1];
    p->size[0] = 1;
    emxEnsureCapacity_real_T(sp, p, b_i, &ed_emlrtRTEI);
    p_data = p->data;
    for (b_i = 0; b_i < b_loop_ub; b_i++) {
      p_data[b_i] *= Om_data[3 * b_i + 1];
    }
  } else {
    st.site = &fb_emlrtRSI;
    binary_expand_op_1(&st, p, vPQW_temp);
    p_data = p->data;
  }
  if ((IJKOverPQW3->size[1] != p->size[1]) &&
      ((IJKOverPQW3->size[1] != 1) && (p->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW3->size[1], p->size[1], &b_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  loop_ub = cosnu->size[1];
  if ((cosnu->size[1] != vPQW_temp->size[1]) &&
      ((cosnu->size[1] != 1) && (vPQW_temp->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(cosnu->size[1], vPQW_temp->size[1], &c_emlrtECI,
                                (emlrtConstCTX)sp);
  }
  if (IJKOverPQW3->size[1] == p->size[1]) {
    b_loop_ub = IJKOverPQW3->size[1] - 1;
    b_i = IJKOverPQW3->size[0] * IJKOverPQW3->size[1];
    IJKOverPQW3->size[0] = 1;
    emxEnsureCapacity_real_T(sp, IJKOverPQW3, b_i, &dd_emlrtRTEI);
    IJKOverPQW3_data = IJKOverPQW3->data;
    for (b_i = 0; b_i <= b_loop_ub; b_i++) {
      IJKOverPQW3_data[b_i] += p_data[b_i];
    }
  } else {
    st.site = &fb_emlrtRSI;
    plus(&st, IJKOverPQW3, p);
  }
  emxFree_real_T(sp, &p);
  if (cosnu->size[1] == vPQW_temp->size[1]) {
    b_i = cosnu->size[0] * cosnu->size[1];
    cosnu->size[0] = 1;
    emxEnsureCapacity_real_T(sp, cosnu, b_i, &fd_emlrtRTEI);
    cosnu_data = cosnu->data;
    for (b_i = 0; b_i < loop_ub; b_i++) {
      cosnu_data[b_i] *= Om_data[3 * b_i + 2];
    }
  } else {
    st.site = &fb_emlrtRSI;
    binary_expand_op(&st, cosnu, vPQW_temp);
    cosnu_data = cosnu->data;
  }
  emxFree_real_T(sp, &vPQW_temp);
  if ((IJKOverPQW3->size[1] != cosnu->size[1]) &&
      ((IJKOverPQW3->size[1] != 1) && (cosnu->size[1] != 1))) {
    emlrtDimSizeImpxCheckR2021b(IJKOverPQW3->size[1], cosnu->size[1],
                                &b_emlrtECI, (emlrtConstCTX)sp);
  }
  if (IJKOverPQW3->size[1] == cosnu->size[1]) {
    loop_ub = IJKOverPQW3->size[1] - 1;
    b_i = IJKOverPQW3->size[0] * IJKOverPQW3->size[1];
    IJKOverPQW3->size[0] = 1;
    emxEnsureCapacity_real_T(sp, IJKOverPQW3, b_i, &dd_emlrtRTEI);
    IJKOverPQW3_data = IJKOverPQW3->data;
    for (b_i = 0; b_i <= loop_ub; b_i++) {
      IJKOverPQW3_data[b_i] += cosnu_data[b_i];
    }
  } else {
    st.site = &fb_emlrtRSI;
    plus(&st, IJKOverPQW3, cosnu);
    IJKOverPQW3_data = IJKOverPQW3->data;
  }
  emxFree_real_T(sp, &cosnu);
  input_sizes[0] = 1;
  input_sizes[1] = vPQW->size[1];
  emlrtSubAssignSizeCheckR2012b(&input_sizes[0], 2, &IJKOverPQW3->size[0], 2,
                                &emlrtECI, (emlrtCTX)sp);
  for (b_i = 0; b_i < f_loop_ub; b_i++) {
    vPQW_data[3 * b_i + 2] = IJKOverPQW3_data[b_i];
  }
  emxFree_real_T(sp, &IJKOverPQW3);
  st.site = &s_emlrtRSI;
  b_st.site = &eb_emlrtRSI;
  if (rPQW->size[1] != 0) {
    f_loop_ub = rPQW->size[1];
  } else if (vPQW->size[1] != 0) {
    f_loop_ub = vPQW->size[1];
  } else {
    f_loop_ub = 0;
  }
  c_st.site = &db_emlrtRSI;
  if ((rPQW->size[1] != f_loop_ub) && (rPQW->size[1] != 0)) {
    emlrtErrorWithMessageIdR2018a(&c_st, &emlrtRTEI,
                                  "MATLAB:catenate:matrixDimensionMismatch",
                                  "MATLAB:catenate:matrixDimensionMismatch", 0);
  }
  if ((vPQW->size[1] != f_loop_ub) && (vPQW->size[1] != 0)) {
    emlrtErrorWithMessageIdR2018a(&c_st, &emlrtRTEI,
                                  "MATLAB:catenate:matrixDimensionMismatch",
                                  "MATLAB:catenate:matrixDimensionMismatch", 0);
  }
  empty_non_axis_sizes = (f_loop_ub == 0);
  if (empty_non_axis_sizes || (rPQW->size[1] != 0)) {
    input_sizes_idx_0 = 3;
  } else {
    input_sizes_idx_0 = 0;
  }
  if (empty_non_axis_sizes || (vPQW->size[1] != 0)) {
    sizes_idx_0 = 3;
  } else {
    sizes_idx_0 = 0;
  }
  e_loop_ub = sizes_idx_0;
  b_i = y->size[0] * y->size[1];
  y->size[0] = input_sizes_idx_0 + sizes_idx_0;
  y->size[1] = f_loop_ub;
  emxEnsureCapacity_real_T(&b_st, y, b_i, &gd_emlrtRTEI);
  w_data = y->data;
  for (b_i = 0; b_i < f_loop_ub; b_i++) {
    loop_ub = input_sizes_idx_0;
    for (loop_ub_tmp = 0; loop_ub_tmp < loop_ub; loop_ub_tmp++) {
      w_data[loop_ub_tmp + y->size[0] * b_i] =
          rPQW_data[loop_ub_tmp + input_sizes_idx_0 * b_i];
    }
    for (loop_ub_tmp = 0; loop_ub_tmp < e_loop_ub; loop_ub_tmp++) {
      w_data[(loop_ub_tmp + input_sizes_idx_0) + y->size[0] * b_i] =
          vPQW_data[loop_ub_tmp + sizes_idx_0 * b_i];
    }
  }
  emxFree_real_T(&b_st, &vPQW);
  emxFree_real_T(&b_st, &rPQW);
  emlrtHeapReferenceStackLeaveFcnR2012b((emlrtConstCTX)sp);
}

/* End of code generation (Keplerian2Cartesian.c) */
