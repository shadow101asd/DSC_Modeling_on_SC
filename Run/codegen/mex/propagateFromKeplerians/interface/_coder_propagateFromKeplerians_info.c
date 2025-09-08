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
      "789ced56cb6ed340149dd082d814522121f5277029a5344bd749abb431b471ba28154aa7"
      "f14de2c4f3905f4956fd0458b263cf0a76ddd30fe037bae8b27be2c4"
      "76628b912b198c887a37d7778e67ce7dd847830a55b580107a82a6f67975ea5782b818f8"
      "07286e49bc10f8a5441cda43b41cdb17e21f03df62d481a1330d2826",
      "10edd4193128a64e63c401596033d3037d82b40d131a06016d3e78eb4764770e8a021ff2"
      "9f952eb4fa9a4b90d5b567199af341d48f4b41bdcb77ecc7b1a01fc5"
      "047e5af9201ddb60d9528f03658e54662d9700756c49951b357947d214a9ac294d95e960"
      "1ab4d364b4395eaabb54e216e3b8831dd8b51839006e8265606abf20",
      "7375f08c753c4da923c45dca71ab1f2511f19f65e47f24e49f223a73cf4d98d5fb2d239f"
      "2ae48be359e696e8556c5e6782fcc23ea7e59ff4b3f71f4ffc97dbeb"
      "36ca910f7dfab9942b5f60ff8a6f2838efaedfdf73015f31817706c35189975fd7ebaa57"
      "c36fdaf2c1fb4db237cbe33085272d0f2488f33a7f51f437ebf7f02c",
      "a58e108fe837146c3960fb1a8c50fa7ff0b77438ebfc8e847c713ccbfc7ed3b3c9ec7ccb"
      "4b3fae6e06dd3cf9fad7dfd7f2e40b6dd1f518f786fdf6767daba42b"
      "fbf455f590accbde4ee55e8f174d8f5753ea087197ebe3441a960b3265049ba360fd7fbd"
      "17bf13f2c5f14cf7e264cfc693cb4b37a4ca45ae3afcb5fae3fe5e8c",
      "febc0e57555d33e4935aa9e4d92f654f1b9c6c9f6f2ec0bdf81765f9349f",
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
                emlrtMxCreateString("ADfbSCEaPu8BZenB3K3kJE"));
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

/* End of code generation (_coder_propagateFromKeplerians_info.c) */
