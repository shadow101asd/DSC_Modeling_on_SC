/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * _coder_propagateFromKeplerians_info.c
 *
 * Code generation for function 'propagateFromKeplerians'
 *
 */

/* Include files */
#include "_coder_propagateFromKeplerians_info.h"
#include "emlrt.h"
#include "tmwtypes.h"

/* Function Declarations */
static const mxArray *c_emlrtMexFcnResolvedFunctionsI(void);

/* Function Definitions */
static const mxArray *c_emlrtMexFcnResolvedFunctionsI(void)
{
  const mxArray *nameCaptureInfo;
  const char_T *data[6] = {
      "789ced56cb6ed340149dd082d800a99090f889ba3c4ae9d2755ad4a686e2b8422d42e960"
      "df246e3c0ff9913a1fd07559b2e31760c707f001fc05ea0276ec8913"
      "db892d46ae6431088bbbb9be733c73eec33e1ad4d8d51b08a1db6866ef5766fe56123713"
      "7f0de5ad883712bf548853bb8e9673fb52fc22f116a30144c12ca098",
      "40b6d366c4a19806e69803f2c067ee08ec29d2735c301d029dc5e0791c919d05280b6228"
      "7ed606600d3b2141dec09f67e82e06593f3e0bea5dbe623f0e05fd68"
      "16f0d7db6f94431f3c5f39e54059a0b4981512a081afe8aab9af6e291d4d6975b4aece6c"
      "701ddaef32da9d2c192155b8c738eee300763c46dac05df01c4cfd55",
      "b25007af58c79d923a523ca41c5bc32c898cffa422ff0d21ff0cb159f8d68579bd1f2bf2"
      "e942be3c5e656e855ee5e67522c82fed7359fe453f7fffe6d47ff879"
      "d94312f9d0bbaf4b52f912fb5b7c91e0bcab7e7ff7047ccd02de3f8bc69bbcb56e18fa68"
      "1f6ff4d4f6f163f26c9ec741094f591e4810cb3abf2efa5bf57bb85b",
      "52478a67f40f35ec05e0c71a8c50f97ff0a774b8eafc5e0af9f27895f9fda667d3d9c526"
      "4b3fbe7c3f1bc8e41b5e7eba2f932fb5baeb313e8d86bda7c6934d5b"
      "dba38f760fc89a3adadafeafc775d3e395923a523ce4f62411d30b41a58c60779cacffab"
      "f7e21742be3c5ee95e5cecd96472b27463f5875c1d3efa76fe40265f",
      "6a75d761e7d83d3237f6e015d622c332db6bf6405baf810eff02b3c135a7",
      ""};
  nameCaptureInfo = NULL;
  emlrtNameCaptureMxArrayR2016a(&data[0], 4528U, &nameCaptureInfo);
  return nameCaptureInfo;
}

mxArray *emlrtMexFcnProperties(void)
{
  mxArray *xEntryPoints;
  mxArray *xInputs;
  mxArray *xResult;
  const char_T *propFieldName[9] = {"Version",
                                    "ResolvedFunctions",
                                    "Checksum",
                                    "EntryPoints",
                                    "CoverageInfo",
                                    "IsPolymorphic",
                                    "PropertyList",
                                    "UUID",
                                    "ClassEntryPointIsHandle"};
  const char_T *epFieldName[8] = {
      "Name",     "NumberOfInputs", "NumberOfOutputs", "ConstantInputs",
      "FullPath", "TimeStamp",      "Constructor",     "Visible"};
  xEntryPoints =
      emlrtCreateStructMatrix(1, 1, 8, (const char_T **)&epFieldName[0]);
  xInputs = emlrtCreateLogicalMatrix(1, 3);
  emlrtSetField(xEntryPoints, 0, "Name",
                emlrtMxCreateString("propagateFromKeplerians"));
  emlrtSetField(xEntryPoints, 0, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(3.0));
  emlrtSetField(xEntryPoints, 0, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(1.0));
  emlrtSetField(xEntryPoints, 0, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 0, "FullPath",
      emlrtMxCreateString("/Users/jpenot/Documents/MATLAB/SC/"
                          "DSC_Modeling_on_SC/Run/propagateFromKeplerians.m"));
  emlrtSetField(xEntryPoints, 0, "TimeStamp",
                emlrtMxCreateDoubleScalar(739814.87836805556));
  emlrtSetField(xEntryPoints, 0, "Constructor",
                emlrtMxCreateLogicalScalar(false));
  emlrtSetField(xEntryPoints, 0, "Visible", emlrtMxCreateLogicalScalar(true));
  xResult =
      emlrtCreateStructMatrix(1, 1, 9, (const char_T **)&propFieldName[0]);
  emlrtSetField(xResult, 0, "Version",
                emlrtMxCreateString("24.1.0.2837808 (R2024a) Update 7"));
  emlrtSetField(xResult, 0, "ResolvedFunctions",
                (mxArray *)c_emlrtMexFcnResolvedFunctionsI());
  emlrtSetField(xResult, 0, "Checksum",
                emlrtMxCreateString("kefC8xvjx09Cs3CzqndwfG"));
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

/* End of code generation (_coder_propagateFromKeplerians_info.c) */
