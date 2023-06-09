//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// Simu_Matlab.h
//
// Code generation for function 'Simu_Matlab'
//

#ifndef SIMU_MATLAB_H
#define SIMU_MATLAB_H

// Include files
#include "rtwtypes.h"
#include "omp.h"
#include <cstddef>
#include <cstdlib>

// Variable Declarations
extern omp_nest_lock_t Simu_Matlab_nestLockGlobal;

// Function Declarations
extern void Simu_Matlab();

extern void Simu_Matlab_initialize();

extern void Simu_Matlab_terminate();

#endif
// End of code generation (Simu_Matlab.h)
