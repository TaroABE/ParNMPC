//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// _coder_Simu_Matlab_mex.cpp
//
// Code generation for function 'Simu_Matlab'
//

// Include files
#include "_coder_Simu_Matlab_mex.h"
#include "_coder_Simu_Matlab_api.h"

// Function Definitions
void mexFunction(int32_T nlhs, mxArray *[], int32_T nrhs, const mxArray *[])
{
  mexAtExit(&Simu_Matlab_atexit);
  // Module initialization.
  Simu_Matlab_initialize();
  // Dispatch the entry-point.
  unsafe_Simu_Matlab_mexFunction(nlhs, nrhs);
  // Module termination.
  Simu_Matlab_terminate();
}

emlrtCTX mexFunctionCreateRootTLS()
{
  emlrtCreateRootTLSR2022a(&emlrtRootTLSGlobal, &emlrtContextGlobal, nullptr, 1,
                           nullptr, "UTF-8", true);
  return emlrtRootTLSGlobal;
}

void unsafe_Simu_Matlab_mexFunction(int32_T nlhs, int32_T nrhs)
{
  emlrtStack st{
      nullptr, // site
      nullptr, // tls
      nullptr  // prev
  };
  st.tls = emlrtRootTLSGlobal;
  // Check for proper number of arguments.
  if (nrhs != 0) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 0, 4,
                        11, "Simu_Matlab");
  }
  if (nlhs > 0) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 11,
                        "Simu_Matlab");
  }
  // Call the function.
  Simu_Matlab_api();
}

// End of code generation (_coder_Simu_Matlab_mex.cpp)
