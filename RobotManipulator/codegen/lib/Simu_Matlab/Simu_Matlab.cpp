//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// Simu_Matlab.cpp
//
// Code generation for function 'Simu_Matlab'
//

// Include files
#include "Simu_Matlab.h"
#include "coder_bounded_array.h"
#include "iiwa14.h"
#include "omp.h"
#include "stdio.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>

// Variable Definitions
static boolean_T xInit_not_empty;

static std::FILE *eml_openfiles[20];

static boolean_T eml_autoflush[20];

omp_nest_lock_t Simu_Matlab_nestLockGlobal;

static boolean_T isInitialized_Simu_Matlab{false};

// Function Declarations
static double NMPC_Solve(const double x0[14], const double p[168],
                         double solution_lambda[336], double solution_u[192],
                         double solution_x[336], double solution_z[696],
                         double solution_LAMBDA[4704], double &output_cost,
                         double &output_KKTError,
                         double &output_timeElapsed_searchDirection,
                         double &output_timeElapsed_lineSearch,
                         double &output_timeElapsed_KKTErrorCheck,
                         double &output_timeElapsed_total, double &output_rho,
                         double &output_iterTotal, double &output_exitflag);

static double NMPC_Solve_SearchDirection(
    const double x0[14], const double p_data[], double rho,
    double lambda_data[], int lambda_size[3], int mu_size[3], double u_data[],
    int u_size[3], double x_data[], int x_size[3], double z_data[],
    int z_size[3], double LAMBDA_data[], double &KKTError_C,
    double &KKTError_Hu, double &KKTError_costateEquation, double &costL,
    double &timeElapsed);

static void NMPC_Solve_init();

static void SIM_Plant_RK4(const double u[7], const double x[14],
                          double xNext[14]);

static void coarse_update_func(
    double lambda_i_data[], double u_i_data[], double x_i_data[],
    const double z_i_data[], const double p_i_data[], double xPrev_i_data[],
    double lambdaNext_i_data[], double LAMBDA_i_data[], double rho, double i,
    double p_muu_F_i_data[], int p_muu_F_i_size[3],
    double p_muu_Lambda_i_data[], int p_muu_Lambda_i_size[3],
    double p_lambda_Lambda_i_data[], int p_lambda_Lambda_i_size[3],
    double p_x_Lambda_i_data[], int p_x_Lambda_i_size[3], double p_x_F_i_data[],
    int p_x_F_i_size[3], double KKTxEquation_i_data[],
    int KKTxEquation_i_size[2], double KKTC_i_data[], int KKTC_i_size[2],
    double KKTHu_i_data[], int KKTHu_i_size[2],
    double KKTlambdaEquation_i_data[], int KKTlambdaEquation_i_size[2],
    double L_i_data[], int L_i_size[2], double LB_i_data[], int LB_i_size[2]);

namespace coder {
namespace internal {
static int cfclose(double fid);

static signed char cfopen();

static std::FILE *fileManager(double varargin_1, boolean_T &a);

static signed char filedata();

namespace reflapack {
static int xzgetrf(double A[196], int ipiv[14]);

}
} // namespace internal
static void inv(const double x[64], double y[64]);

static void mldivide(double B[14]);

} // namespace coder
static void filedata_init();

static double fraction_to_boundary_parallel_func(
    const double u_i_data[], const double x_i_data[], double z_i_data[],
    const double u_k_i_data[], const double x_k_i_data[],
    const double z_k_i_data[], double rho, double &stepSizeMaxG_i);

// Function Definitions
static double NMPC_Solve(const double x0[14], const double p[168],
                         double solution_lambda[336], double solution_u[192],
                         double solution_x[336], double solution_z[696],
                         double solution_LAMBDA[4704], double &output_cost,
                         double &output_KKTError,
                         double &output_timeElapsed_searchDirection,
                         double &output_timeElapsed_lineSearch,
                         double &output_timeElapsed_KKTErrorCheck,
                         double &output_timeElapsed_total, double &output_rho,
                         double &output_iterTotal, double &output_exitflag)
{
  static coder::bounded_array<double, 4704U, 4U> LAMBDAInit;
  static coder::bounded_array<double, 336U, 3U> lambdaInit;
  static coder::bounded_array<double, 336U, 3U> xInit;
  static coder::bounded_array<double, 192U, 3U> uInit;
  static const signed char iv[192]{
      0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1};
  coder::bounded_array<double, 696U, 3U> zInit;
  double KKTError_Hu;
  double KKTError_costateEquation;
  double lambda_L1Norm;
  double output_iterInit;
  double tSearchDirection;
  double tStart;
  int lambdaSplit_size[3];
  int muSplit_size[3];
  int uSplit_size[3];
  int xSplit_size[3];
  int zSplit_size[3];
  int b_iter;
  int iter;
  int mode;
  boolean_T exitg1;
  tStart = omp_get_wtime();
  //  only for the very first problem
  //  Reshape
  if (!xInit_not_empty) {
    lambdaInit.size[0] = 14;
    lambdaInit.size[1] = 6;
    lambdaInit.size[2] = 4;
    std::memset(&lambdaInit.data[0], 0, 336U * sizeof(double));
    uInit.size[0] = 8;
    uInit.size[1] = 6;
    uInit.size[2] = 4;
    for (mode = 0; mode < 192; mode++) {
      uInit.data[mode] = iv[mode];
    }
    xInit.size[0] = 14;
    xInit.size[1] = 6;
    xInit.size[2] = 4;
    std::memset(&xInit.data[0], 0, 336U * sizeof(double));
    xInit_not_empty = true;
    LAMBDAInit.size[0] = 14;
    LAMBDAInit.size[1] = 14;
    LAMBDAInit.size[2] = 6;
    LAMBDAInit.size[3] = 4;
    std::memset(&LAMBDAInit.data[0], 0, 4704U * sizeof(double));
    //  assert
  }
  zInit.size[0] = 29;
  zInit.size[1] = 6;
  zInit.size[2] = 4;
  for (mode = 0; mode < 696; mode++) {
    zInit.data[mode] = 1.0;
  }
  //  Init
  lambdaSplit_size[0] = 14;
  lambdaSplit_size[1] = 6;
  lambdaSplit_size[2] = 4;
  std::copy(&lambdaInit.data[0], &lambdaInit.data[336], &solution_lambda[0]);
  muSplit_size[0] = 0;
  muSplit_size[1] = 6;
  muSplit_size[2] = 4;
  uSplit_size[0] = 8;
  uSplit_size[1] = 6;
  uSplit_size[2] = 4;
  std::copy(&uInit.data[0], &uInit.data[192], &solution_u[0]);
  xSplit_size[0] = 14;
  xSplit_size[1] = 6;
  xSplit_size[2] = 4;
  std::copy(&xInit.data[0], &xInit.data[336], &solution_x[0]);
  zSplit_size[0] = 29;
  zSplit_size[1] = 6;
  zSplit_size[2] = 4;
  std::copy(&zInit.data[0], &zInit.data[696], &solution_z[0]);
  std::copy(&LAMBDAInit.data[0], &LAMBDAInit.data[4704], &solution_LAMBDA[0]);
  output_iterInit = 0.0;
  output_timeElapsed_searchDirection = 0.0;
  output_timeElapsed_lineSearch = 0.0;
  output_timeElapsed_KKTErrorCheck = 0.0;
  output_exitflag = 0.0;
  output_cost = 0.0;
  output_KKTError = 0.0;
  mode = 1;
  //  initialize a filter [LAll;xEq+C;flag]
  //  Iteration
  iter = 1;
  b_iter = 0;
  exitg1 = false;
  while ((!exitg1) && (b_iter < 10)) {
    double KKTError_stateEquation;
    iter = b_iter + 1;
    //  backup
    //         %% Search direction
    KKTError_stateEquation = NMPC_Solve_SearchDirection(
        x0, p, 0.001, solution_lambda, lambdaSplit_size, muSplit_size,
        solution_u, uSplit_size, solution_x, xSplit_size, solution_z,
        zSplit_size, solution_LAMBDA, lambda_L1Norm, KKTError_Hu,
        KKTError_costateEquation, output_cost, tSearchDirection);
    output_timeElapsed_searchDirection += tSearchDirection;
    //         %% Line search
    //         %% KKT error
    //  avoid chattering of the num of iter
    lambda_L1Norm = 0.0;
    for (int k{0}; k < 336; k++) {
      lambda_L1Norm += std::abs(solution_lambda[k]);
    }
    lambda_L1Norm = std::fmax(100.0, lambda_L1Norm / 14.0 / 24.0) / 100.0;
    tSearchDirection = KKTError_Hu / lambda_L1Norm / 5.0;
    lambda_L1Norm = KKTError_costateEquation / lambda_L1Norm / 5.0;
    output_KKTError = KKTError_stateEquation * 5.0;
    if (KKTError_stateEquation * 5.0 < 0.0) {
      output_KKTError = 0.0;
    }
    if (output_KKTError < tSearchDirection) {
      output_KKTError = tSearchDirection;
    }
    if (output_KKTError < lambda_L1Norm) {
      output_KKTError = lambda_L1Norm;
    }
    //         %% print
    //         %% barrier parameter update
    if (mode == 1) {
      //  init
      if ((output_KKTError < 0.001) || (b_iter + 1 >= 10)) {
        mode = 2;
        lambdaInit.size[0] = 14;
        lambdaInit.size[1] = 6;
        lambdaInit.size[2] = 4;
        std::copy(&solution_lambda[0], &solution_lambda[336],
                  &lambdaInit.data[0]);
        uInit.size[0] = 8;
        uInit.size[1] = 6;
        uInit.size[2] = 4;
        std::copy(&solution_u[0], &solution_u[192], &uInit.data[0]);
        xInit.size[0] = 14;
        xInit.size[1] = 6;
        xInit.size[2] = 4;
        std::copy(&solution_x[0], &solution_x[336], &xInit.data[0]);
        LAMBDAInit.size[0] = 14;
        LAMBDAInit.size[1] = 14;
        LAMBDAInit.size[2] = 6;
        LAMBDAInit.size[3] = 4;
        std::copy(&solution_LAMBDA[0], &solution_LAMBDA[4704],
                  &LAMBDAInit.data[0]);
        output_iterInit = static_cast<double>(b_iter) + 1.0;
        if (output_KKTError < 0.001) {
          output_exitflag = 1.0;
          exitg1 = true;
        } else {
          b_iter++;
        }
      } else {
        b_iter++;
      }

      //  decay
    } else if (output_KKTError < 0.001) {
      output_exitflag = 1.0;
      exitg1 = true;
    } else {
      b_iter++;
    }
  }
  //     %% Reshape
  output_rho = 0.001;
  output_iterTotal = iter;
  lambda_L1Norm = omp_get_wtime();
  output_timeElapsed_total = lambda_L1Norm - tStart;
  return output_iterInit;
}

static double NMPC_Solve_SearchDirection(
    const double x0[14], const double p_data[], double rho,
    double lambda_data[], int lambda_size[3], int mu_size[3], double u_data[],
    int u_size[3], double x_data[], int x_size[3], double z_data[],
    int z_size[3], double LAMBDA_data[], double &KKTError_C,
    double &KKTError_Hu, double &KKTError_costateEquation, double &costL,
    double &timeElapsed)
{
  double p_lambda_Lambda_data[4704];
  double p_x_F_data[4704];
  double p_x_Lambda_data[4704];
  double p_muu_F_data[2688];
  double p_muu_Lambda_data[2688];
  double LAMBDA_i_data[1176];
  double p_lambda_Lambda_i_data[1176];
  double p_x_F_i_data[1176];
  double p_x_Lambda_i_data[1176];
  double z_k_data[696];
  double p_muu_F_i_data[672];
  double p_muu_Lambda_i_data[672];
  double dlambda_data[336];
  double lambdaNext_data[336];
  double lambda_k_data[336];
  double xPrev_data[336];
  double x_k_data[336];
  double u_k_data[192];
  double z_i_data[174];
  double z_k_i_data[174];
  double lambda_i_data[84];
  double xPrev_i_data[84];
  double x_i_data[84];
  double x_k_i_data[84];
  double u_i_data[48];
  double u_k_i_data[48];
  double p_i_data[42];
  double KKTHu_data[24];
  double KKTlambdaEquation_data[24];
  double KKTxEquation_data[24];
  double L_data[24];
  double dlambda_j_i[14];
  double lambda_next[14];
  double dmu_u_j_i[8];
  double KKTC_i_data[6];
  double KKTHu_i_data[6];
  double KKTlambdaEquation_i_data[6];
  double KKTxEquation_i_data[6];
  double LB_i_data[6];
  double L_i_data[6];
  double stepSizeG_data[4];
  double stepSizeZ_data[4];
  double KKTError_stateEquation;
  double stepSizeG_i;
  double stepSizeZ_i;
  double tEnd;
  double tStart;
  int p_lambda_Lambda_i_size[3];
  int p_muu_F_i_size[3];
  int p_muu_Lambda_i_size[3];
  int p_x_F_i_size[3];
  int p_x_Lambda_i_size[3];
  int KKTC_i_size[2];
  int KKTHu_i_size[2];
  int KKTlambdaEquation_i_size[2];
  int KKTxEquation_i_size[2];
  int LB_i_size[2];
  int L_i_size[2];
  int b_xPrev_data_tmp;
  int i1;
  int i2;
  int i3;
  int j;
  int lambda_i_data_tmp;
  int xPrev_data_tmp;
  tStart = omp_get_wtime();
  //  global variables
  //  parallel seg
  //  Code generation
  //  backup
  std::copy(&lambda_data[0], &lambda_data[336], &lambda_k_data[0]);
  std::copy(&u_data[0], &u_data[192], &u_k_data[0]);
  std::copy(&x_data[0], &x_data[336], &x_k_data[0]);
  std::copy(&z_data[0], &z_data[696], &z_k_data[0]);
  //  line search parameters
  //  local variables
  std::memset(&lambdaNext_data[0], 0, 336U * sizeof(double));
  std::memset(&xPrev_data[0], 0, 336U * sizeof(double));
  std::memset(&dlambda_data[0], 0, 336U * sizeof(double));
  //  coupling variable for each segment
  std::copy(&x0[0], &x0[14], &xPrev_data[0]);
  for (int i{0}; i < 3; i++) {
    for (int b_i{0}; b_i < 14; b_i++) {
      xPrev_data_tmp = b_i + 84 * (i + 1);
      b_xPrev_data_tmp = (b_i + 84 * i) + 70;
      xPrev_data[xPrev_data_tmp] = x_data[b_xPrev_data_tmp];
      lambdaNext_data[b_xPrev_data_tmp] = lambda_data[xPrev_data_tmp];
      for (i1 = 0; i1 < 14; i1++) {
        b_xPrev_data_tmp = i1 + 14 * b_i;
        LAMBDA_data[(b_xPrev_data_tmp + 1176 * i) + 980] =
            LAMBDA_data[b_xPrev_data_tmp + 1176 * (i + 1)];
      }
    }
  }
  std::memset(&LAMBDA_data[4508], 0, 196U * sizeof(double));
  //     %% V(:,index_inside_sig_j,which_sigment_i)
#pragma omp parallel for num_threads(                                          \
    4 > omp_get_max_threads()                                                  \
        ? omp_get_max_threads()                                                \
        : 4) private(LAMBDA_i_data, x_k_i_data, xPrev_i_data, x_i_data,        \
                     u_i_data, lambda_i_data, LB_i_data, LB_i_size, L_i_data,  \
                     L_i_size, KKTlambdaEquation_i_data,                       \
                     KKTlambdaEquation_i_size, KKTHu_i_data, KKTHu_i_size,     \
                     KKTC_i_data, KKTC_i_size, KKTxEquation_i_data,            \
                     KKTxEquation_i_size, p_x_F_i_data, p_x_F_i_size,          \
                     p_x_Lambda_i_data, p_x_Lambda_i_size,                     \
                     p_lambda_Lambda_i_data, p_lambda_Lambda_i_size,           \
                     p_muu_Lambda_i_data, p_muu_Lambda_i_size, p_muu_F_i_data, \
                     p_muu_F_i_size, p_i_data, z_i_data, i2, i3,               \
                     lambda_i_data_tmp, j)

  for (int c_i = 0; c_i < 4; c_i++) {
    //      for i=1:1:DoP
    for (i2 = 0; i2 < 6; i2++) {
      std::copy(&lambda_data[c_i * 84 + i2 * 14],
                &lambda_data[static_cast<int>(
                    static_cast<unsigned int>(c_i * 84 + i2 * 14) + 14U)],
                &lambda_i_data[i2 * 14]);
      std::copy(&u_data[c_i * 48 + i2 * 8],
                &u_data[static_cast<int>(
                    static_cast<unsigned int>(c_i * 48 + i2 * 8) + 8U)],
                &u_i_data[i2 * 8]);
      std::copy(&x_data[c_i * 84 + i2 * 14],
                &x_data[static_cast<int>(
                    static_cast<unsigned int>(c_i * 84 + i2 * 14) + 14U)],
                &x_i_data[i2 * 14]);
      std::copy(&z_data[c_i * 174 + i2 * 29],
                &z_data[static_cast<int>(
                    static_cast<unsigned int>(c_i * 174 + i2 * 29) + 29U)],
                &z_i_data[i2 * 29]);
      for (i3 = 0; i3 < 7; i3++) {
        lambda_i_data_tmp = i3 + 7 * i2;
        p_i_data[lambda_i_data_tmp] = p_data[lambda_i_data_tmp + 42 * c_i];
      }
      for (i3 = 0; i3 < 14; i3++) {
        std::copy(
            &LAMBDA_data[(c_i * 1176 + i2 * 196) + i3 * 14],
            &LAMBDA_data[static_cast<int>(
                static_cast<unsigned int>((c_i * 1176 + i2 * 196) + i3 * 14) +
                14U)],
            &LAMBDA_i_data[i2 * 196 + i3 * 14]);
        lambda_i_data_tmp = i3 + 14 * i2;
        j = lambda_i_data_tmp + 84 * c_i;
        xPrev_i_data[lambda_i_data_tmp] = xPrev_data[j];
        x_k_i_data[lambda_i_data_tmp] = lambdaNext_data[j];
      }
    }
    coarse_update_func(
        lambda_i_data, u_i_data, x_i_data, z_i_data, p_i_data, xPrev_i_data,
        x_k_i_data, LAMBDA_i_data, rho, static_cast<double>(c_i) + 1.0,
        p_muu_F_i_data, p_muu_F_i_size, p_muu_Lambda_i_data,
        p_muu_Lambda_i_size, p_lambda_Lambda_i_data, p_lambda_Lambda_i_size,
        p_x_Lambda_i_data, p_x_Lambda_i_size, p_x_F_i_data, p_x_F_i_size,
        KKTxEquation_i_data, KKTxEquation_i_size, KKTC_i_data, KKTC_i_size,
        KKTHu_i_data, KKTHu_i_size, KKTlambdaEquation_i_data,
        KKTlambdaEquation_i_size, L_i_data, L_i_size, LB_i_data, LB_i_size);
    //  Recover
    //
    for (i2 = 0; i2 < 6; i2++) {
      std::copy(&lambda_i_data[i2 * 14],
                &lambda_i_data[static_cast<int>(
                    static_cast<unsigned int>(i2 * 14) + 14U)],
                &lambda_data[c_i * 84 + i2 * 14]);
      std::copy(
          &u_i_data[i2 * 8],
          &u_i_data[static_cast<int>(static_cast<unsigned int>(i2 * 8) + 8U)],
          &u_data[c_i * 48 + i2 * 8]);
      for (i3 = 0; i3 < 14; i3++) {
        lambda_i_data_tmp = i3 + 14 * i2;
        j = lambda_i_data_tmp + 84 * c_i;
        x_data[j] = x_i_data[lambda_i_data_tmp];
        xPrev_data[j] = xPrev_i_data[lambda_i_data_tmp];
        std::copy(&LAMBDA_i_data[i2 * 196 + i3 * 14],
                  &LAMBDA_i_data[static_cast<int>(
                      static_cast<unsigned int>(i2 * 196 + i3 * 14) + 14U)],
                  &LAMBDA_data[(c_i * 1176 + i2 * 196) + i3 * 14]);
        lambdaNext_data[j] = x_k_i_data[lambda_i_data_tmp];
        std::copy(&p_muu_F_i_data[i2 * 112 + i3 * 8],
                  &p_muu_F_i_data[static_cast<int>(
                      static_cast<unsigned int>(i2 * 112 + i3 * 8) + 8U)],
                  &p_muu_F_data[(c_i * 672 + i2 * 112) + i3 * 8]);
        std::copy(&p_muu_Lambda_i_data[i2 * 112 + i3 * 8],
                  &p_muu_Lambda_i_data[static_cast<int>(
                      static_cast<unsigned int>(i2 * 112 + i3 * 8) + 8U)],
                  &p_muu_Lambda_data[(c_i * 672 + i2 * 112) + i3 * 8]);
        std::copy(&p_lambda_Lambda_i_data[i2 * 196 + i3 * 14],
                  &p_lambda_Lambda_i_data[static_cast<int>(
                      static_cast<unsigned int>(i2 * 196 + i3 * 14) + 14U)],
                  &p_lambda_Lambda_data[(c_i * 1176 + i2 * 196) + i3 * 14]);
        std::copy(&p_x_Lambda_i_data[i2 * 196 + i3 * 14],
                  &p_x_Lambda_i_data[static_cast<int>(
                      static_cast<unsigned int>(i2 * 196 + i3 * 14) + 14U)],
                  &p_x_Lambda_data[(c_i * 1176 + i2 * 196) + i3 * 14]);
        std::copy(&p_x_F_i_data[i2 * 196 + i3 * 14],
                  &p_x_F_i_data[static_cast<int>(
                      static_cast<unsigned int>(i2 * 196 + i3 * 14) + 14U)],
                  &p_x_F_data[(c_i * 1176 + i2 * 196) + i3 * 14]);
      }
      lambda_i_data_tmp = i2 + 6 * c_i;
      KKTxEquation_data[lambda_i_data_tmp] = KKTxEquation_i_data[i2];
      KKTHu_data[lambda_i_data_tmp] = KKTHu_i_data[i2];
      KKTlambdaEquation_data[lambda_i_data_tmp] = KKTlambdaEquation_i_data[i2];
      L_data[lambda_i_data_tmp] = L_i_data[i2];
    }
  }
  //     %%
  KKTError_stateEquation = KKTxEquation_data[0];
  KKTError_C = 0.0;
  KKTError_Hu = KKTHu_data[0];
  KKTError_costateEquation = KKTlambdaEquation_data[0];
  costL = L_data[0];
  for (xPrev_data_tmp = 0; xPrev_data_tmp < 23; xPrev_data_tmp++) {
    tEnd = KKTxEquation_data[xPrev_data_tmp + 1];
    if (KKTError_stateEquation < tEnd) {
      KKTError_stateEquation = tEnd;
    }
    tEnd = KKTHu_data[xPrev_data_tmp + 1];
    if (KKTError_Hu < tEnd) {
      KKTError_Hu = tEnd;
    }
    tEnd = KKTlambdaEquation_data[xPrev_data_tmp + 1];
    if (KKTError_costateEquation < tEnd) {
      KKTError_costateEquation = tEnd;
    }
    costL += L_data[xPrev_data_tmp + 1];
  }
  //     %% Step 2: Backward correction due to the approximation of lambda
  for (int i{0}; i < 3; i++) {
    for (b_xPrev_data_tmp = 0; b_xPrev_data_tmp < 6; b_xPrev_data_tmp++) {
      if (6 - b_xPrev_data_tmp == 6) {
        std::copy(&lambda_data[i * -84 + 252], &lambda_data[i * -84 + 266],
                  &lambda_next[0]);
      } else {
        std::copy(&lambda_data[(i * -84 + b_xPrev_data_tmp * -14) + 252],
                  &lambda_data[(i * -84 + b_xPrev_data_tmp * -14) + 266],
                  &lambda_next[0]);
      }
      for (int b_i{0}; b_i < 14; b_i++) {
        xPrev_data_tmp = (b_i + 14 * (5 - b_xPrev_data_tmp)) + 84 * (2 - i);
        dlambda_data[xPrev_data_tmp] =
            lambda_next[b_i] - lambdaNext_data[xPrev_data_tmp];
      }
      for (int b_i{0}; b_i < 14; b_i++) {
        tEnd = 0.0;
        for (i1 = 0; i1 < 14; i1++) {
          tEnd +=
              p_lambda_Lambda_data[((b_i + 14 * i1) +
                                    196 * (5 - b_xPrev_data_tmp)) +
                                   1176 * (2 - i)] *
              dlambda_data[(i1 + 14 * (5 - b_xPrev_data_tmp)) + 84 * (2 - i)];
        }
        i1 = (b_i + 14 * (5 - b_xPrev_data_tmp)) + 84 * (2 - i);
        lambda_data[i1] -= tEnd;
      }
    }
  }
#pragma omp parallel for num_threads(                                          \
    4 > omp_get_max_threads() ? omp_get_max_threads()                          \
                              : 4) private(dmu_u_j_i, dlambda_j_i, x_i_data,   \
                                           u_i_data, i2, j, stepSizeG_i, i3)

  for (int c_i = 0; c_i < 4; c_i++) {
    //      for i=1:1:DoP
    for (i2 = 0; i2 < 6; i2++) {
      std::copy(&u_data[c_i * 48 + i2 * 8],
                &u_data[static_cast<int>(
                    static_cast<unsigned int>(c_i * 48 + i2 * 8) + 8U)],
                &u_i_data[i2 * 8]);
      std::copy(&x_data[c_i * 84 + i2 * 14],
                &x_data[static_cast<int>(
                    static_cast<unsigned int>(c_i * 84 + i2 * 14) + 14U)],
                &x_i_data[i2 * 14]);
    }
    for (j = 0; j < 6; j++) {
      for (i2 = 0; i2 < 8; i2++) {
        stepSizeG_i = 0.0;
        for (i3 = 0; i3 < 14; i3++) {
          stepSizeG_i +=
              p_muu_Lambda_data[((i2 + 8 * i3) + 112 * (5 - j)) + 672 * c_i] *
              dlambda_data[(i3 + 14 * (5 - j)) + 84 * c_i];
        }
        dmu_u_j_i[i2] = stepSizeG_i;
      }
      for (i2 = 0; i2 < 14; i2++) {
        stepSizeG_i = 0.0;
        for (i3 = 0; i3 < 14; i3++) {
          stepSizeG_i +=
              p_x_Lambda_data[((i2 + 14 * i3) + 196 * (5 - j)) + 1176 * c_i] *
              dlambda_data[(i3 + 14 * (5 - j)) + 84 * c_i];
        }
        dlambda_j_i[i2] = stepSizeG_i;
      }
      for (i2 = 0; i2 < 8; i2++) {
        i3 = i2 + 8 * (5 - j);
        stepSizeG_i = u_i_data[i3] - dmu_u_j_i[i2];
        dmu_u_j_i[i2] = stepSizeG_i;
        u_i_data[i3] = stepSizeG_i;
      }
      for (i2 = 0; i2 < 14; i2++) {
        i3 = i2 + 14 * (5 - j);
        stepSizeG_i = x_i_data[i3] - dlambda_j_i[i2];
        dlambda_j_i[i2] = stepSizeG_i;
        x_i_data[i3] = stepSizeG_i;
      }
    }
    for (i2 = 0; i2 < 6; i2++) {
      std::copy(
          &u_i_data[i2 * 8],
          &u_i_data[static_cast<int>(static_cast<unsigned int>(i2 * 8) + 8U)],
          &u_data[c_i * 48 + i2 * 8]);
      std::copy(
          &x_i_data[i2 * 14],
          &x_i_data[static_cast<int>(static_cast<unsigned int>(i2 * 14) + 14U)],
          &x_data[c_i * 84 + i2 * 14]);
    }
  }
  //     %% Step 3: Forward correction due to the approximation of x
  for (int i{0}; i < 4; i++) {
    for (b_xPrev_data_tmp = 0; b_xPrev_data_tmp < 6; b_xPrev_data_tmp++) {
      if (b_xPrev_data_tmp + 1 == 1) {
        if (i + 1 == 1) {
          std::copy(&x0[0], &x0[14], &lambda_next[0]);
        } else {
          std::copy(&x_data[i * 84 + -14], &x_data[i * 84], &lambda_next[0]);
        }
      } else {
        std::copy(&x_data[(i * 84 + b_xPrev_data_tmp * 14) + -14],
                  &x_data[i * 84 + b_xPrev_data_tmp * 14], &lambda_next[0]);
      }
      for (int b_i{0}; b_i < 14; b_i++) {
        xPrev_data_tmp = (b_i + 14 * b_xPrev_data_tmp) + 84 * i;
        lambdaNext_data[xPrev_data_tmp] =
            lambda_next[b_i] - xPrev_data[xPrev_data_tmp];
      }
      for (int b_i{0}; b_i < 14; b_i++) {
        tEnd = 0.0;
        for (i1 = 0; i1 < 14; i1++) {
          tEnd += p_x_F_data[((b_i + 14 * i1) + 196 * b_xPrev_data_tmp) +
                             1176 * i] *
                  lambdaNext_data[(i1 + 14 * b_xPrev_data_tmp) + 84 * i];
        }
        i1 = (b_i + 14 * b_xPrev_data_tmp) + 84 * i;
        x_data[i1] -= tEnd;
      }
    }
  }
#pragma omp parallel for num_threads(                                          \
    4 > omp_get_max_threads()                                                  \
        ? omp_get_max_threads()                                                \
        : 4) private(z_i_data, stepSizeG_i, stepSizeZ_i, dmu_u_j_i,            \
                     z_k_i_data, x_k_i_data, u_k_i_data, x_i_data, u_i_data,   \
                     lambda_i_data, i2, j, i3, lambda_i_data_tmp)

  for (int c_i = 0; c_i < 4; c_i++) {
    //      for i=1:1:DoP
    //  line search variables
    for (i2 = 0; i2 < 6; i2++) {
      std::copy(&lambda_data[c_i * 84 + i2 * 14],
                &lambda_data[static_cast<int>(
                    static_cast<unsigned int>(c_i * 84 + i2 * 14) + 14U)],
                &lambda_i_data[i2 * 14]);
      std::copy(&u_data[c_i * 48 + i2 * 8],
                &u_data[static_cast<int>(
                    static_cast<unsigned int>(c_i * 48 + i2 * 8) + 8U)],
                &u_i_data[i2 * 8]);
      std::copy(&x_data[c_i * 84 + i2 * 14],
                &x_data[static_cast<int>(
                    static_cast<unsigned int>(c_i * 84 + i2 * 14) + 14U)],
                &x_i_data[i2 * 14]);
      std::copy(&z_data[c_i * 174 + i2 * 29],
                &z_data[static_cast<int>(
                    static_cast<unsigned int>(c_i * 174 + i2 * 29) + 29U)],
                &z_i_data[i2 * 29]);
      std::copy(&u_k_data[c_i * 48 + i2 * 8],
                &u_k_data[static_cast<int>(
                    static_cast<unsigned int>(c_i * 48 + i2 * 8) + 8U)],
                &u_k_i_data[i2 * 8]);
      std::copy(&x_k_data[c_i * 84 + i2 * 14],
                &x_k_data[static_cast<int>(
                    static_cast<unsigned int>(c_i * 84 + i2 * 14) + 14U)],
                &x_k_i_data[i2 * 14]);
      std::copy(&z_k_data[c_i * 174 + i2 * 29],
                &z_k_data[static_cast<int>(
                    static_cast<unsigned int>(c_i * 174 + i2 * 29) + 29U)],
                &z_k_i_data[i2 * 29]);
    }
    for (j = 0; j < 6; j++) {
      for (i2 = 0; i2 < 8; i2++) {
        stepSizeG_i = 0.0;
        for (i3 = 0; i3 < 14; i3++) {
          stepSizeG_i +=
              p_muu_F_data[((i2 + 8 * i3) + 112 * (5 - j)) + 672 * c_i] *
              lambdaNext_data[(i3 + 14 * (5 - j)) + 84 * c_i];
        }
        dmu_u_j_i[i2] = stepSizeG_i;
      }
      for (i2 = 0; i2 < 14; i2++) {
        stepSizeG_i = 0.0;
        for (i3 = 0; i3 < 14; i3++) {
          stepSizeG_i +=
              LAMBDA_data[((i2 + 14 * i3) + 196 * (5 - j)) + 1176 * c_i] *
              lambdaNext_data[(i3 + 14 * (5 - j)) + 84 * c_i];
        }
        lambda_i_data_tmp = i2 + 14 * (5 - j);
        lambda_i_data[lambda_i_data_tmp] -= stepSizeG_i;
      }
      for (i2 = 0; i2 < 8; i2++) {
        i3 = i2 + 8 * (5 - j);
        stepSizeG_i = u_i_data[i3] - dmu_u_j_i[i2];
        dmu_u_j_i[i2] = stepSizeG_i;
        u_i_data[i3] = stepSizeG_i;
      }
    }
    for (i2 = 0; i2 < 6; i2++) {
      std::copy(&lambda_i_data[i2 * 14],
                &lambda_i_data[static_cast<int>(
                    static_cast<unsigned int>(i2 * 14) + 14U)],
                &lambda_data[c_i * 84 + i2 * 14]);
      std::copy(
          &u_i_data[i2 * 8],
          &u_i_data[static_cast<int>(static_cast<unsigned int>(i2 * 8) + 8U)],
          &u_data[c_i * 48 + i2 * 8]);
    }
    stepSizeZ_i = fraction_to_boundary_parallel_func(
        u_i_data, x_i_data, z_i_data, u_k_i_data, x_k_i_data, z_k_i_data, rho,
        stepSizeG_i);
    //  Recover
    for (i2 = 0; i2 < 6; i2++) {
      std::copy(&lambda_i_data[i2 * 14],
                &lambda_i_data[static_cast<int>(
                    static_cast<unsigned int>(i2 * 14) + 14U)],
                &lambda_data[c_i * 84 + i2 * 14]);
      std::copy(
          &u_i_data[i2 * 8],
          &u_i_data[static_cast<int>(static_cast<unsigned int>(i2 * 8) + 8U)],
          &u_data[c_i * 48 + i2 * 8]);
      std::copy(
          &z_i_data[i2 * 29],
          &z_i_data[static_cast<int>(static_cast<unsigned int>(i2 * 29) + 29U)],
          &z_data[c_i * 174 + i2 * 29]);
    }
    stepSizeZ_data[c_i] = stepSizeZ_i;
    stepSizeG_data[c_i] = stepSizeG_i;
  }
  //     %% Line Search to Guarantee Primal Stability
  tEnd = stepSizeG_data[0];
  if (stepSizeG_data[0] > stepSizeG_data[1]) {
    tEnd = stepSizeG_data[1];
  }
  if (tEnd > stepSizeG_data[2]) {
    tEnd = stepSizeG_data[2];
  }
  if (tEnd > stepSizeG_data[3]) {
    tEnd = stepSizeG_data[3];
  }
  if (tEnd != 1.0) {
    lambda_size[0] = 14;
    lambda_size[1] = 6;
    lambda_size[2] = 4;
    mu_size[0] = 0;
    mu_size[1] = 6;
    mu_size[2] = 4;
    u_size[0] = 8;
    u_size[1] = 6;
    u_size[2] = 4;
    x_size[0] = 14;
    x_size[1] = 6;
    x_size[2] = 4;
    for (int b_i{0}; b_i < 4; b_i++) {
      for (i1 = 0; i1 < 6; i1++) {
        for (b_xPrev_data_tmp = 0; b_xPrev_data_tmp < 14; b_xPrev_data_tmp++) {
          xPrev_data_tmp = (b_xPrev_data_tmp + 14 * i1) + 84 * b_i;
          lambda_data[xPrev_data_tmp] =
              (1.0 - tEnd) * lambda_k_data[xPrev_data_tmp] +
              tEnd * lambda_data[xPrev_data_tmp];
        }
        for (b_xPrev_data_tmp = 0; b_xPrev_data_tmp < 8; b_xPrev_data_tmp++) {
          xPrev_data_tmp = (b_xPrev_data_tmp + 8 * i1) + 48 * b_i;
          u_data[xPrev_data_tmp] = (1.0 - tEnd) * u_k_data[xPrev_data_tmp] +
                                   tEnd * u_data[xPrev_data_tmp];
        }
        for (b_xPrev_data_tmp = 0; b_xPrev_data_tmp < 14; b_xPrev_data_tmp++) {
          xPrev_data_tmp = (b_xPrev_data_tmp + 14 * i1) + 84 * b_i;
          x_data[xPrev_data_tmp] = (1.0 - tEnd) * x_k_data[xPrev_data_tmp] +
                                   tEnd * x_data[xPrev_data_tmp];
        }
      }
    }
  }
  tEnd = stepSizeZ_data[0];
  if (stepSizeZ_data[0] > stepSizeZ_data[1]) {
    tEnd = stepSizeZ_data[1];
  }
  if (tEnd > stepSizeZ_data[2]) {
    tEnd = stepSizeZ_data[2];
  }
  if (tEnd > stepSizeZ_data[3]) {
    tEnd = stepSizeZ_data[3];
  }
  if (tEnd != 1.0) {
    z_size[0] = 29;
    z_size[1] = 6;
    z_size[2] = 4;
    for (int b_i{0}; b_i < 696; b_i++) {
      z_data[b_i] = (1.0 - tEnd) * z_k_data[b_i] + tEnd * z_data[b_i];
    }
  }
  //     %%
  tEnd = omp_get_wtime();
  timeElapsed = tEnd - tStart;
  return KKTError_stateEquation;
}

static void NMPC_Solve_init()
{
  xInit_not_empty = false;
}

static void SIM_Plant_RK4(const double u[7], const double x[14],
                          double xNext[14])
{
  double b_x[14];
  double k2[14];
  double k3[14];
  double q[7];
  double qd[7];
  double qdd[7];
  double tau[7];
  double d;
  //  dx = f(u,x,p)
  //  Specify your own f(u,x,p) function for code generation
  for (int i{0}; i < 7; i++) {
    q[i] = x[i];
    qd[i] = x[i + 7];
    qdd[i] = 0.0;
    tau[i] = u[i];
  }
  sim_qdd_cal(&q[0], &qd[0], &qdd[0], &tau[0]);
  for (int i{0}; i < 7; i++) {
    xNext[i] = qd[i];
    xNext[i + 7] = qdd[i];
  }
  coder::mldivide(xNext);
  for (int i{0}; i < 14; i++) {
    d = 0.001 * xNext[i];
    xNext[i] = d;
    b_x[i] = x[i] + d / 2.0;
  }
  //  dx = f(u,x,p)
  //  Specify your own f(u,x,p) function for code generation
  for (int i{0}; i < 7; i++) {
    q[i] = b_x[i];
    qd[i] = b_x[i + 7];
    qdd[i] = 0.0;
    tau[i] = u[i];
  }
  sim_qdd_cal(&q[0], &qd[0], &qdd[0], &tau[0]);
  for (int i{0}; i < 7; i++) {
    k2[i] = qd[i];
    k2[i + 7] = qdd[i];
  }
  coder::mldivide(k2);
  for (int i{0}; i < 14; i++) {
    d = 0.001 * k2[i];
    k2[i] = d;
    b_x[i] = x[i] + d / 2.0;
  }
  //  dx = f(u,x,p)
  //  Specify your own f(u,x,p) function for code generation
  for (int i{0}; i < 7; i++) {
    q[i] = b_x[i];
    qd[i] = b_x[i + 7];
    qdd[i] = 0.0;
    tau[i] = u[i];
  }
  sim_qdd_cal(&q[0], &qd[0], &qdd[0], &tau[0]);
  for (int i{0}; i < 7; i++) {
    k3[i] = qd[i];
    k3[i + 7] = qdd[i];
  }
  coder::mldivide(k3);
  for (int i{0}; i < 14; i++) {
    d = 0.001 * k3[i];
    k3[i] = d;
    b_x[i] = x[i] + d;
  }
  //  dx = f(u,x,p)
  //  Specify your own f(u,x,p) function for code generation
  for (int i{0}; i < 7; i++) {
    q[i] = b_x[i];
    qd[i] = b_x[i + 7];
    qdd[i] = 0.0;
    tau[i] = u[i];
  }
  sim_qdd_cal(&q[0], &qd[0], &qdd[0], &tau[0]);
  for (int i{0}; i < 7; i++) {
    b_x[i] = qd[i];
    b_x[i + 7] = qdd[i];
  }
  coder::mldivide(b_x);
  for (int i{0}; i < 14; i++) {
    xNext[i] =
        x[i] +
        (((xNext[i] + 2.0 * k2[i]) + 2.0 * k3[i]) + 0.001 * b_x[i]) / 6.0;
  }
}

static void coarse_update_func(
    double lambda_i_data[], double u_i_data[], double x_i_data[],
    const double z_i_data[], const double p_i_data[], double xPrev_i_data[],
    double lambdaNext_i_data[], double LAMBDA_i_data[], double rho, double i,
    double p_muu_F_i_data[], int p_muu_F_i_size[3],
    double p_muu_Lambda_i_data[], int p_muu_Lambda_i_size[3],
    double p_lambda_Lambda_i_data[], int p_lambda_Lambda_i_size[3],
    double p_x_Lambda_i_data[], int p_x_Lambda_i_size[3], double p_x_F_i_data[],
    int p_x_F_i_size[3], double KKTxEquation_i_data[],
    int KKTxEquation_i_size[2], double KKTC_i_data[], int KKTC_i_size[2],
    double KKTHu_i_data[], int KKTHu_i_size[2],
    double KKTlambdaEquation_i_data[], int KKTlambdaEquation_i_size[2],
    double L_i_data[], int L_i_size[2], double LB_i_data[], int LB_i_size[2])
{
  static const signed char iv[98]{
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1};
  double Ix[196];
  double dv[196];
  double invFx_j_i[196];
  double Aux_invFx_j_i[112];
  double FuT_invFxT_j_i[112];
  double p_lambda_muu_j_i[112];
  double p_muu_Lambda_j_i[112];
  double p_x_muu_j_i[112];
  double z_i[64];
  double dq[49];
  double dqd[49];
  double dtau[49];
  double b_Ix[14];
  //  regularization
  p_muu_F_i_size[0] = 8;
  p_muu_F_i_size[1] = 14;
  p_muu_F_i_size[2] = 6;
  p_muu_Lambda_i_size[0] = 8;
  p_muu_Lambda_i_size[1] = 14;
  p_muu_Lambda_i_size[2] = 6;
  p_lambda_Lambda_i_size[0] = 14;
  p_lambda_Lambda_i_size[1] = 14;
  p_lambda_Lambda_i_size[2] = 6;
  p_x_Lambda_i_size[0] = 14;
  p_x_Lambda_i_size[1] = 14;
  p_x_Lambda_i_size[2] = 6;
  p_x_F_i_size[0] = 14;
  p_x_F_i_size[1] = 14;
  p_x_F_i_size[2] = 6;
  L_i_size[0] = 1;
  L_i_size[1] = 6;
  LB_i_size[0] = 1;
  LB_i_size[1] = 6;
  KKTxEquation_i_size[0] = 1;
  KKTxEquation_i_size[1] = 6;
  KKTC_i_size[0] = 1;
  KKTC_i_size[1] = 6;
  KKTHu_i_size[0] = 1;
  KKTHu_i_size[1] = 6;
  KKTlambdaEquation_i_size[0] = 1;
  KKTlambdaEquation_i_size[1] = 6;
  dv[0] = 1.0;
  std::memset(&dv[1], 0, 14U * sizeof(double));
  dv[15] = 1.0;
  std::memset(&dv[16], 0, 14U * sizeof(double));
  dv[30] = 1.0;
  std::memset(&dv[31], 0, 14U * sizeof(double));
  dv[45] = 1.0;
  std::memset(&dv[46], 0, 14U * sizeof(double));
  dv[60] = 1.0;
  std::memset(&dv[61], 0, 14U * sizeof(double));
  dv[75] = 1.0;
  std::memset(&dv[76], 0, 14U * sizeof(double));
  dv[90] = 1.0;
  std::memset(&dv[91], 0, 13U * sizeof(double));
  dv[104] = 0.0;
  std::memset(&dv[106], 0, 13U * sizeof(double));
  dv[119] = 0.0;
  std::memset(&dv[121], 0, 13U * sizeof(double));
  dv[134] = 0.0;
  std::memset(&dv[136], 0, 13U * sizeof(double));
  dv[149] = 0.0;
  std::memset(&dv[151], 0, 13U * sizeof(double));
  dv[164] = 0.0;
  std::memset(&dv[166], 0, 13U * sizeof(double));
  dv[179] = 0.0;
  std::memset(&dv[181], 0, 14U * sizeof(double));
  z_i[1] = 0.0;
  z_i[2] = 0.0;
  z_i[3] = 0.0;
  z_i[4] = 0.0;
  z_i[5] = 0.0;
  z_i[6] = 0.0;
  z_i[7] = 0.0;
  z_i[8] = 0.0;
  z_i[10] = 0.0;
  z_i[11] = 0.0;
  z_i[12] = 0.0;
  z_i[13] = 0.0;
  z_i[14] = 0.0;
  z_i[15] = 0.0;
  z_i[16] = 0.0;
  z_i[17] = 0.0;
  z_i[19] = 0.0;
  z_i[20] = 0.0;
  z_i[21] = 0.0;
  z_i[22] = 0.0;
  z_i[23] = 0.0;
  z_i[24] = 0.0;
  z_i[25] = 0.0;
  z_i[26] = 0.0;
  z_i[28] = 0.0;
  z_i[29] = 0.0;
  z_i[30] = 0.0;
  z_i[31] = 0.0;
  z_i[32] = 0.0;
  z_i[33] = 0.0;
  z_i[34] = 0.0;
  z_i[35] = 0.0;
  z_i[37] = 0.0;
  z_i[38] = 0.0;
  z_i[39] = 0.0;
  z_i[40] = 0.0;
  z_i[41] = 0.0;
  z_i[42] = 0.0;
  z_i[43] = 0.0;
  z_i[44] = 0.0;
  z_i[46] = 0.0;
  z_i[47] = 0.0;
  z_i[48] = 0.0;
  z_i[49] = 0.0;
  z_i[50] = 0.0;
  z_i[51] = 0.0;
  z_i[52] = 0.0;
  z_i[53] = 0.0;
  std::memset(&z_i[55], 0, 8U * sizeof(double));
  for (int j{0}; j < 6; j++) {
    double LAMBDAUncrt_j_i[196];
    double LAMBDA_i[196];
    double c_Ix[196];
    double FuT_LAMBDAUncrt_m_Aux_invFx_j_i[112];
    double reshapes_f2[112];
    double Inv_dKKT23_mu_u_j_i[64];
    double b_z_i[64];
    double lambdaEq_j_i[14];
    double xEq_j_i[14];
    double HAlluT_j_i[8];
    double b_HAlluT_j_i[8];
    double dmu_u_j_i[8];
    double b_qd[7];
    double q[7];
    double qd[7];
    double qdd[7];
    double tau[7];
    double d;
    double d1;
    double d10;
    double d11;
    double d12;
    double d13;
    double d14;
    double d15;
    double d2;
    double d3;
    double d4;
    double d5;
    double d6;
    double d7;
    double d8;
    double d9;
    double t10;
    double t12;
    double t14;
    double t17;
    double t18;
    double t19;
    double t2;
    double t20;
    double t21;
    double t22;
    double t23;
    double t24;
    double t25;
    double t26;
    double t27;
    double t28;
    double t29;
    double t30;
    double t4;
    double t41_tmp;
    double t43_tmp;
    double t45_tmp;
    double t47_tmp;
    double t49_tmp;
    double t51_tmp;
    double t6;
    double t8;
    int ipiv[14];
    int b_i;
    int i1;
    int i2;
    int i3;
    int invFx_j_i_tmp;
    int kAcol;
    int pipk;
    int t31_tmp;
    signed char p[14];
    //  Function and Jacobian
    // OCP_GEN_L_Lu_Lx
    //     [L,Lu,Lx] = OCP_GEN_L_Lu_Lx(IN1,IN2,IN3)
    //     This function was generated by the Symbolic Math Toolbox version 9.3.
    //     2023/06/09 12:24:34
    b_i = 14 * (5 - j);
    d = x_i_data[b_i];
    pipk = 7 * (5 - j);
    t17 = p_i_data[pipk];
    t2 = t17 - d;
    t29 = x_i_data[b_i + 1];
    t18 = p_i_data[pipk + 1];
    t4 = t18 - t29;
    t19 = x_i_data[b_i + 2];
    t20 = p_i_data[pipk + 2];
    t6 = t20 - t19;
    t21 = x_i_data[b_i + 3];
    t22 = p_i_data[pipk + 3];
    t8 = t22 - t21;
    t23 = x_i_data[b_i + 4];
    t24 = p_i_data[pipk + 4];
    t10 = t24 - t23;
    t25 = x_i_data[b_i + 5];
    t26 = p_i_data[pipk + 5];
    t12 = t26 - t25;
    t27 = x_i_data[b_i + 6];
    t28 = p_i_data[pipk + 6];
    t14 = t28 - t27;
    d1 = x_i_data[b_i + 7];
    d2 = x_i_data[b_i + 8];
    d3 = x_i_data[b_i + 9];
    d4 = x_i_data[b_i + 10];
    d5 = x_i_data[b_i + 11];
    d6 = x_i_data[b_i + 12];
    d7 = x_i_data[b_i + 13];
    i1 = 8 * (5 - j);
    d8 = u_i_data[i1];
    d9 = u_i_data[i1 + 1];
    d10 = u_i_data[i1 + 2];
    d11 = u_i_data[i1 + 3];
    d12 = u_i_data[i1 + 4];
    d13 = u_i_data[i1 + 5];
    d14 = u_i_data[i1 + 6];
    d15 = u_i_data[i1 + 7];
    L_i_data[5 - j] =
        ((((((((((((((((((d8 * d8 / 2000.0 + d9 * d9 / 2000.0) +
                         d10 * d10 / 2000.0) +
                        d11 * d11 / 2000.0) +
                       d12 * d12 / 2000.0) +
                      d13 * d13 / 2000.0) +
                     d14 * d14 / 2000.0) +
                    d15 * d15 * 1000.0) +
                   d1 * d1 / 20.0) +
                  d2 * d2 / 20.0) +
                 d3 * d3 / 20.0) +
                d4 * d4 / 20.0) +
               d5 * d5 / 20.0) +
              d6 * d6 / 20.0) +
             d7 * d7 / 20.0) +
            t2 * (t17 / 2.0 - d / 2.0)) +
           t4 * (t18 / 2.0 - t29 / 2.0)) +
          t6 * (t20 / 2.0 - t19 / 2.0)) +
         t8 * (t22 / 2.0 - t21 / 2.0)) +
        ((t10 * (t24 / 2.0 - t23 / 2.0) + t12 * (t26 / 2.0 - t25 / 2.0)) +
         t14 * (t28 / 2.0 - t27 / 2.0));
    // OCP_GEN_LB_LBu_LBx
    //     [LB,LBu,LBx] = OCP_GEN_LB_LBu_LBx(IN1,IN2,IN3)
    //     This function was generated by the Symbolic Math Toolbox version 9.3.
    //     2023/06/09 12:24:34
    t17 = (d15 + 1.5707963267948966) + d1;
    t18 = (d15 + 1.5707963267948966) + d2;
    t19 = (d15 + 1.5707963267948966) + d3;
    t20 = (d15 + 1.5707963267948966) + d4;
    t21 = (d15 + 1.5707963267948966) + d5;
    t22 = (d15 + 1.5707963267948966) + d6;
    t23 = (d15 + 1.5707963267948966) + d7;
    t24 = (-d1 + 1.5707963267948966) + d15;
    t25 = (-d2 + 1.5707963267948966) + d15;
    t26 = (-d3 + 1.5707963267948966) + d15;
    t27 = (-d4 + 1.5707963267948966) + d15;
    t28 = (-d5 + 1.5707963267948966) + d15;
    t29 = (-d6 + 1.5707963267948966) + d15;
    t30 = (-d7 + 1.5707963267948966) + d15;
    LB_i_data[5 - j] =
        (((((((((((((((((((((((((((((d15 * 0.0015 + 0.0021991148575128553) -
                                    std::log(-d8 + 10.0)) -
                                   std::log(-d9 + 10.0)) -
                                  std::log(-d10 + 10.0)) -
                                 std::log(-d11 + 10.0)) -
                                std::log(-d12 + 10.0)) -
                               std::log(-d13 + 10.0)) -
                              std::log(-d14 + 10.0)) -
                             std::log(d8 + 10.0)) -
                            std::log(d9 + 10.0)) -
                           std::log(d10 + 10.0)) -
                          std::log(d11 + 10.0)) -
                         std::log(d12 + 10.0)) -
                        std::log(d13 + 10.0)) -
                       std::log(d14 + 10.0)) -
                      std::log(t17)) -
                     std::log(t18)) -
                    std::log(t19)) -
                   std::log(t20)) -
                  std::log(t21)) -
                 std::log(t22)) -
                std::log(t23)) -
               std::log(t24)) -
              std::log(t25)) -
             std::log(t26)) -
            std::log(t27)) -
           std::log(t28)) -
          std::log(t29)) -
         std::log(t30)) +
        (-std::log(d15) + 0.014);
    t17 = 1.0 / t17;
    t18 = 1.0 / t18;
    t19 = 1.0 / t19;
    t20 = 1.0 / t20;
    t21 = 1.0 / t21;
    t22 = 1.0 / t22;
    t23 = 1.0 / t23;
    t24 = 1.0 / t24;
    t41_tmp = 1.0 / t25;
    t43_tmp = 1.0 / t26;
    t45_tmp = 1.0 / t27;
    t47_tmp = 1.0 / t28;
    t49_tmp = 1.0 / t29;
    t51_tmp = 1.0 / t30;
    std::memset(&Ix[0], 0, 196U * sizeof(double));
    for (int k{0}; k < 14; k++) {
      Ix[k + 14 * k] = 1.0;
    }
    //  M disabled
    //  'Euler'
    //  dx = f(u,x,p)
    //  Specify your own f(u,x,p) function for code generation
    for (int c_i{0}; c_i < 7; c_i++) {
      pipk = c_i + b_i;
      q[c_i] = x_i_data[pipk];
      qd[c_i] = x_i_data[pipk + 7];
      qdd[c_i] = 0.0;
      tau[c_i] = u_i_data[c_i + i1];
    }
    qdd_cal(&q[0], &qd[0], &qdd[0], &tau[0], i);
    std::memset(&p_lambda_muu_j_i[0], 0, 112U * sizeof(double));
    for (i2 = 0; i2 < 7; i2++) {
      pipk = i2 + b_i;
      q[i2] = x_i_data[pipk];
      b_qd[i2] = x_i_data[pipk + 7];
      tau[i2] = u_i_data[i2 + i1];
    }
    std::memset(&dq[0], 0, 49U * sizeof(double));
    std::memset(&dqd[0], 0, 49U * sizeof(double));
    std::memset(&dtau[0], 0, 49U * sizeof(double));
    derivatives_cal(&q[0], &b_qd[0], &tau[0], &dq[0], &dqd[0], &dtau[0], i);
    for (i2 = 0; i2 < 7; i2++) {
      for (i3 = 0; i3 < 7; i3++) {
        p_lambda_muu_j_i[(i3 + 14 * i2) + 7] = dtau[i3 + 7 * i2];
      }
    }
    for (i2 = 0; i2 < 112; i2++) {
      p_lambda_muu_j_i[i2] *= 0.041667;
    }
    for (i2 = 0; i2 < 14; i2++) {
      for (i3 = 0; i3 < 7; i3++) {
        pipk = i3 + 14 * i2;
        LAMBDAUncrt_j_i[pipk] =
            static_cast<double>(iv[i3 + 7 * i2]) * 0.041667 - Ix[pipk];
      }
    }
    for (i2 = 0; i2 < 7; i2++) {
      for (i3 = 0; i3 < 7; i3++) {
        pipk = (i3 + 14 * i2) + 7;
        kAcol = i3 + 7 * i2;
        LAMBDAUncrt_j_i[pipk] = dq[kAcol] * 0.041667 - Ix[pipk];
        pipk = (i3 + 14 * (i2 + 7)) + 7;
        LAMBDAUncrt_j_i[pipk] = dqd[kAcol] * 0.041667 - Ix[pipk];
      }
    }
    //  KKT
    if (6 - j > 1) {
      std::copy(&x_i_data[j * -14 + 56], &x_i_data[j * -14 + 70],
                &xPrev_i_data[b_i]);
    }
    if (6 - j < 6) {
      std::copy(&lambda_i_data[j * -14 + 84], &lambda_i_data[j * -14 + 98],
                &lambdaNext_i_data[b_i]);
    }
    for (i2 = 0; i2 < 7; i2++) {
      xEq_j_i[i2] = qd[i2] * 0.041667;
      xEq_j_i[i2 + 7] = qdd[i2] * 0.041667;
    }
    for (i2 = 0; i2 < 14; i2++) {
      pipk = i2 + b_i;
      xEq_j_i[i2] = xPrev_i_data[pipk] + (xEq_j_i[i2] - x_i_data[pipk]);
      for (i3 = 0; i3 < 8; i3++) {
        FuT_invFxT_j_i[i3 + (i2 << 3)] = p_lambda_muu_j_i[i2 + 14 * i3];
      }
    }
    std::copy(&FuT_invFxT_j_i[0], &FuT_invFxT_j_i[112], &p_x_muu_j_i[0]);
    lambdaEq_j_i[0] = -t2;
    lambdaEq_j_i[1] = -t4;
    lambdaEq_j_i[2] = -t6;
    lambdaEq_j_i[3] = -t8;
    lambdaEq_j_i[4] = -t10;
    lambdaEq_j_i[5] = -t12;
    lambdaEq_j_i[6] = -t14;
    lambdaEq_j_i[7] = d1 / 10.0;
    lambdaEq_j_i[8] = d2 / 10.0;
    lambdaEq_j_i[9] = d3 / 10.0;
    lambdaEq_j_i[10] = d4 / 10.0;
    lambdaEq_j_i[11] = d5 / 10.0;
    lambdaEq_j_i[12] = d6 / 10.0;
    lambdaEq_j_i[13] = d7 / 10.0;
    b_Ix[0] = 0.0;
    b_Ix[1] = 0.0;
    b_Ix[2] = 0.0;
    b_Ix[3] = 0.0;
    b_Ix[4] = 0.0;
    b_Ix[5] = 0.0;
    b_Ix[6] = 0.0;
    b_Ix[7] = rho * (-t17 + t24);
    b_Ix[8] = rho * (-t18 + t41_tmp);
    b_Ix[9] = rho * (-t19 + t43_tmp);
    b_Ix[10] = rho * (-t20 + t45_tmp);
    b_Ix[11] = rho * (-t21 + t47_tmp);
    b_Ix[12] = rho * (-t22 + t49_tmp);
    b_Ix[13] = rho * (-t23 + t51_tmp);
    for (i2 = 0; i2 < 14; i2++) {
      d = 0.0;
      for (i3 = 0; i3 < 14; i3++) {
        d += LAMBDAUncrt_j_i[i3 + 14 * i2] * lambda_i_data[i3 + b_i];
      }
      lambdaEq_j_i[i2] =
          (lambdaNext_i_data[i2 + b_i] + (lambdaEq_j_i[i2] + d)) + b_Ix[i2];
    }
    HAlluT_j_i[0] = d8 / 1000.0;
    HAlluT_j_i[1] = d9 / 1000.0;
    HAlluT_j_i[2] = d10 / 1000.0;
    HAlluT_j_i[3] = d11 / 1000.0;
    HAlluT_j_i[4] = d12 / 1000.0;
    HAlluT_j_i[5] = d13 / 1000.0;
    HAlluT_j_i[6] = d14 / 1000.0;
    HAlluT_j_i[7] = d15 * 2000.0;
    dmu_u_j_i[0] = rho * (-1.0 / (d8 - 10.0) - 1.0 / (d8 + 10.0));
    dmu_u_j_i[1] = rho * (-1.0 / (d9 - 10.0) - 1.0 / (d9 + 10.0));
    dmu_u_j_i[2] = rho * (-1.0 / (d10 - 10.0) - 1.0 / (d10 + 10.0));
    dmu_u_j_i[3] = rho * (-1.0 / (d11 - 10.0) - 1.0 / (d11 + 10.0));
    dmu_u_j_i[4] = rho * (-1.0 / (d12 - 10.0) - 1.0 / (d12 + 10.0));
    dmu_u_j_i[5] = rho * (-1.0 / (d13 - 10.0) - 1.0 / (d13 + 10.0));
    dmu_u_j_i[6] = rho * (-1.0 / (d14 - 10.0) - 1.0 / (d14 + 10.0));
    dmu_u_j_i[7] =
        rho *
        (((((((((((((((-t17 - t24) - t18) - t41_tmp) - t19) - t43_tmp) - t20) -
                 t45_tmp) -
                t21) -
               t47_tmp) -
              t22) -
             t49_tmp) -
            t23) -
           t51_tmp) -
          1.0 / d15) +
         0.0015);
    for (i2 = 0; i2 < 8; i2++) {
      d = 0.0;
      for (i3 = 0; i3 < 14; i3++) {
        d += p_x_muu_j_i[i2 + (i3 << 3)] * lambda_i_data[i3 + b_i];
      }
      HAlluT_j_i[i2] = (HAlluT_j_i[i2] + d) + dmu_u_j_i[i2];
    }
    //  Hessian
    //  --- AuuCondensed_j_i = Auu_j_i + Gu_j_i.'*(z_j_i./G_j_i.*Gu_j_i);
    //  --- AuxCondensed_j_i = Aux_j_i + Gu_j_i.'*(z_j_i./G_j_i.*Gx_j_i);
    //  --- AxxCondensed_j_i = Axx_j_i + Gx_j_i.'*(z_j_i./G_j_i.*Gx_j_i);
    // OCP_GEN_Auu_Aux_Axx_Condensed
    //     [AuuCondensed,AuxCondensed,AxxCondensed] =
    //     OCP_GEN_Auu_Aux_Axx_Condensed(IN1,IN2,IN3,IN4,IN5,IN6) This function
    //     was generated by the Symbolic Math Toolbox version 9.3. 2023/06/09
    //     12:24:35
    t31_tmp = 29 * (5 - j);
    t2 = t17 * z_i_data[t31_tmp + 15];
    t30 = t18 * z_i_data[t31_tmp + 16];
    t29 = t19 * z_i_data[t31_tmp + 17];
    t28 = t20 * z_i_data[t31_tmp + 18];
    t27 = t21 * z_i_data[t31_tmp + 19];
    t26 = t22 * z_i_data[t31_tmp + 20];
    t25 = t23 * z_i_data[t31_tmp + 21];
    t24 *= z_i_data[t31_tmp + 22];
    t23 = t41_tmp * z_i_data[t31_tmp + 23];
    t22 = t43_tmp * z_i_data[t31_tmp + 24];
    t21 = t45_tmp * z_i_data[t31_tmp + 25];
    t20 = t47_tmp * z_i_data[t31_tmp + 26];
    t19 = t49_tmp * z_i_data[t31_tmp + 27];
    t18 = t51_tmp * z_i_data[t31_tmp + 28];
    //  descent regularization
    //  descent regularization
    //  nonsingular regularization
    //  Intermediate Variables
    if (6 - j < 6) {
      for (i2 = 0; i2 < 196; i2++) {
        i3 = i2 + j * -196;
        LAMBDA_i_data[i3 + 980] = LAMBDA_i_data[i3 + 1176];
      }
    }
    std::memset(&invFx_j_i[0], 0, 196U * sizeof(double));
    coder::internal::reflapack::xzgetrf(LAMBDAUncrt_j_i, ipiv);
    for (i2 = 0; i2 < 14; i2++) {
      p[i2] = static_cast<signed char>(i2 + 1);
    }
    for (int k{0}; k < 13; k++) {
      i2 = ipiv[k];
      if (i2 > k + 1) {
        pipk = p[i2 - 1];
        p[i2 - 1] = p[k];
        p[k] = static_cast<signed char>(pipk);
      }
    }
    for (int k{0}; k < 14; k++) {
      invFx_j_i_tmp = 14 * (p[k] - 1);
      invFx_j_i[k + invFx_j_i_tmp] = 1.0;
      for (int b_j{k + 1}; b_j < 15; b_j++) {
        i2 = (b_j + invFx_j_i_tmp) - 1;
        if (invFx_j_i[i2] != 0.0) {
          i3 = b_j + 1;
          for (int c_i{i3}; c_i < 15; c_i++) {
            pipk = (c_i + invFx_j_i_tmp) - 1;
            invFx_j_i[pipk] -=
                invFx_j_i[i2] * LAMBDAUncrt_j_i[(c_i + 14 * (b_j - 1)) - 1];
          }
        }
      }
    }
    for (int b_j{0}; b_j < 14; b_j++) {
      pipk = 14 * b_j;
      for (int k{13}; k >= 0; k--) {
        kAcol = 14 * k;
        i2 = k + pipk;
        d = invFx_j_i[i2];
        if (d != 0.0) {
          invFx_j_i[i2] = d / LAMBDAUncrt_j_i[k + kAcol];
          for (int c_i{0}; c_i < k; c_i++) {
            invFx_j_i_tmp = c_i + pipk;
            invFx_j_i[invFx_j_i_tmp] -=
                invFx_j_i[i2] * LAMBDAUncrt_j_i[c_i + kAcol];
          }
        }
      }
    }
    dv[105] = (t2 + t24) + 0.1;
    dv[120] = (t30 + t23) + 0.1;
    dv[135] = (t29 + t22) + 0.1;
    dv[150] = (t28 + t21) + 0.1;
    dv[165] = (t27 + t20) + 0.1;
    dv[180] = (t26 + t19) + 0.1;
    dv[195] = (t25 + t18) + 0.1;
    for (i2 = 0; i2 < 14; i2++) {
      for (i3 = 0; i3 < 14; i3++) {
        kAcol = i3 + 14 * i2;
        Ix[kAcol] = invFx_j_i[i2 + 14 * i3];
        LAMBDA_i[kAcol] = dv[kAcol] - LAMBDA_i_data[kAcol + 196 * (5 - j)];
      }
    }
    for (i2 = 0; i2 < 14; i2++) {
      for (i3 = 0; i3 < 14; i3++) {
        d = 0.0;
        for (invFx_j_i_tmp = 0; invFx_j_i_tmp < 14; invFx_j_i_tmp++) {
          d += Ix[i2 + 14 * invFx_j_i_tmp] * LAMBDA_i[invFx_j_i_tmp + 14 * i3];
        }
        c_Ix[i2 + 14 * i3] = d;
      }
      for (i3 = 0; i3 < 14; i3++) {
        d = 0.0;
        for (invFx_j_i_tmp = 0; invFx_j_i_tmp < 14; invFx_j_i_tmp++) {
          d += c_Ix[i2 + 14 * invFx_j_i_tmp] *
               invFx_j_i[invFx_j_i_tmp + 14 * i3];
        }
        LAMBDAUncrt_j_i[i2 + 14 * i3] = d;
      }
    }
    std::memset(&FuT_invFxT_j_i[0], 0, 63U * sizeof(double));
    FuT_invFxT_j_i[63] = t2 - t24;
    FuT_invFxT_j_i[64] = 0.0;
    FuT_invFxT_j_i[65] = 0.0;
    FuT_invFxT_j_i[66] = 0.0;
    FuT_invFxT_j_i[67] = 0.0;
    FuT_invFxT_j_i[68] = 0.0;
    FuT_invFxT_j_i[69] = 0.0;
    FuT_invFxT_j_i[70] = 0.0;
    FuT_invFxT_j_i[71] = t30 - t23;
    FuT_invFxT_j_i[72] = 0.0;
    FuT_invFxT_j_i[73] = 0.0;
    FuT_invFxT_j_i[74] = 0.0;
    FuT_invFxT_j_i[75] = 0.0;
    FuT_invFxT_j_i[76] = 0.0;
    FuT_invFxT_j_i[77] = 0.0;
    FuT_invFxT_j_i[78] = 0.0;
    FuT_invFxT_j_i[79] = t29 - t22;
    FuT_invFxT_j_i[80] = 0.0;
    FuT_invFxT_j_i[81] = 0.0;
    FuT_invFxT_j_i[82] = 0.0;
    FuT_invFxT_j_i[83] = 0.0;
    FuT_invFxT_j_i[84] = 0.0;
    FuT_invFxT_j_i[85] = 0.0;
    FuT_invFxT_j_i[86] = 0.0;
    FuT_invFxT_j_i[87] = t28 - t21;
    FuT_invFxT_j_i[88] = 0.0;
    FuT_invFxT_j_i[89] = 0.0;
    FuT_invFxT_j_i[90] = 0.0;
    FuT_invFxT_j_i[91] = 0.0;
    FuT_invFxT_j_i[92] = 0.0;
    FuT_invFxT_j_i[93] = 0.0;
    FuT_invFxT_j_i[94] = 0.0;
    FuT_invFxT_j_i[95] = t27 - t20;
    FuT_invFxT_j_i[96] = 0.0;
    FuT_invFxT_j_i[97] = 0.0;
    FuT_invFxT_j_i[98] = 0.0;
    FuT_invFxT_j_i[99] = 0.0;
    FuT_invFxT_j_i[100] = 0.0;
    FuT_invFxT_j_i[101] = 0.0;
    FuT_invFxT_j_i[102] = 0.0;
    FuT_invFxT_j_i[103] = t26 - t19;
    FuT_invFxT_j_i[104] = 0.0;
    FuT_invFxT_j_i[105] = 0.0;
    FuT_invFxT_j_i[106] = 0.0;
    FuT_invFxT_j_i[107] = 0.0;
    FuT_invFxT_j_i[108] = 0.0;
    FuT_invFxT_j_i[109] = 0.0;
    FuT_invFxT_j_i[110] = 0.0;
    FuT_invFxT_j_i[111] = t25 - t18;
    for (i2 = 0; i2 < 8; i2++) {
      for (i3 = 0; i3 < 14; i3++) {
        d = 0.0;
        for (invFx_j_i_tmp = 0; invFx_j_i_tmp < 14; invFx_j_i_tmp++) {
          d += FuT_invFxT_j_i[i2 + (invFx_j_i_tmp << 3)] *
               invFx_j_i[invFx_j_i_tmp + 14 * i3];
        }
        Aux_invFx_j_i[i2 + (i3 << 3)] = d;
      }
    }
    for (int b_j{0}; b_j < 14; b_j++) {
      kAcol = b_j << 3;
      for (int c_i{0}; c_i < 8; c_i++) {
        pipk = c_i * 14;
        t17 = 0.0;
        for (int k{0}; k < 14; k++) {
          t17 += p_lambda_muu_j_i[pipk + k] * invFx_j_i[k * 14 + b_j];
        }
        FuT_invFxT_j_i[kAcol + c_i] = t17;
      }
    }
    for (int b_j{0}; b_j < 8; b_j++) {
      for (i2 = 0; i2 < 14; i2++) {
        d = 0.0;
        for (i3 = 0; i3 < 14; i3++) {
          d += p_x_muu_j_i[b_j + (i3 << 3)] * LAMBDAUncrt_j_i[i3 + 14 * i2];
        }
        kAcol = b_j + (i2 << 3);
        FuT_LAMBDAUncrt_m_Aux_invFx_j_i[kAcol] = d - Aux_invFx_j_i[kAcol];
      }
      kAcol = b_j << 3;
      for (int c_i{0}; c_i < 8; c_i++) {
        pipk = c_i * 14;
        t17 = 0.0;
        for (int k{0}; k < 14; k++) {
          t17 += p_lambda_muu_j_i[pipk + k] * Aux_invFx_j_i[(k << 3) + b_j];
        }
        Inv_dKKT23_mu_u_j_i[kAcol + c_i] = t17;
      }
    }
    z_i[0] = (z_i_data[t31_tmp] / (d8 + 10.0) -
              z_i_data[t31_tmp + 7] / (d8 - 10.0)) +
             0.001;
    z_i[9] = (z_i_data[t31_tmp + 1] / (d9 + 10.0) -
              z_i_data[t31_tmp + 8] / (d9 - 10.0)) +
             0.001;
    z_i[18] = (z_i_data[t31_tmp + 2] / (d10 + 10.0) -
               z_i_data[t31_tmp + 9] / (d10 - 10.0)) +
              0.001;
    z_i[27] = (z_i_data[t31_tmp + 3] / (d11 + 10.0) -
               z_i_data[t31_tmp + 10] / (d11 - 10.0)) +
              0.001;
    z_i[36] = (z_i_data[t31_tmp + 4] / (d12 + 10.0) -
               z_i_data[t31_tmp + 11] / (d12 - 10.0)) +
              0.001;
    z_i[45] = (z_i_data[t31_tmp + 5] / (d13 + 10.0) -
               z_i_data[t31_tmp + 12] / (d13 - 10.0)) +
              0.001;
    z_i[54] = (z_i_data[t31_tmp + 6] / (d14 + 10.0) -
               z_i_data[t31_tmp + 13] / (d14 - 10.0)) +
              0.001;
    z_i[63] =
        ((((((((((((((t2 + t30) + t29) + t28) + t27) + t26) + t25) + t24) +
               t23) +
              t22) +
             t21) +
            t20) +
           t19) +
          t18) +
         z_i_data[t31_tmp + 14] / d15) +
        2000.0;
    //  Sensitivities
    for (i2 = 0; i2 < 8; i2++) {
      for (i3 = 0; i3 < 8; i3++) {
        d = 0.0;
        for (invFx_j_i_tmp = 0; invFx_j_i_tmp < 14; invFx_j_i_tmp++) {
          d += FuT_LAMBDAUncrt_m_Aux_invFx_j_i[i2 + (invFx_j_i_tmp << 3)] *
               p_lambda_muu_j_i[invFx_j_i_tmp + 14 * i3];
        }
        pipk = i2 + (i3 << 3);
        b_z_i[pipk] = z_i[pipk] + (d - Inv_dKKT23_mu_u_j_i[pipk]);
      }
      for (i3 = 0; i3 < 14; i3++) {
        reshapes_f2[i3 + 14 * i2] =
            FuT_LAMBDAUncrt_m_Aux_invFx_j_i[i2 + (i3 << 3)];
      }
    }
    coder::inv(b_z_i, Inv_dKKT23_mu_u_j_i);
    for (i2 = 0; i2 < 8; i2++) {
      for (i3 = 0; i3 < 14; i3++) {
        pipk = i2 + (i3 << 3);
        kAcol = i3 + 14 * i2;
        p_lambda_muu_j_i[kAcol] = FuT_LAMBDAUncrt_m_Aux_invFx_j_i[pipk];
        p_x_muu_j_i[kAcol] = -FuT_invFxT_j_i[pipk];
        d = 0.0;
        for (invFx_j_i_tmp = 0; invFx_j_i_tmp < 8; invFx_j_i_tmp++) {
          d += Inv_dKKT23_mu_u_j_i[i2 + (invFx_j_i_tmp << 3)] *
               reshapes_f2[i3 + 14 * invFx_j_i_tmp];
        }
        Aux_invFx_j_i[pipk] = d;
      }
    }
    for (i2 = 0; i2 < 8; i2++) {
      for (i3 = 0; i3 < 14; i3++) {
        d = 0.0;
        for (invFx_j_i_tmp = 0; invFx_j_i_tmp < 8; invFx_j_i_tmp++) {
          d += Inv_dKKT23_mu_u_j_i[i2 + (invFx_j_i_tmp << 3)] *
               p_x_muu_j_i[i3 + 14 * invFx_j_i_tmp];
        }
        p_muu_Lambda_j_i[i2 + (i3 << 3)] = d;
      }
    }
    //  --- LAMBDA = p_lambda_F
    for (i2 = 0; i2 < 14; i2++) {
      for (i3 = 0; i3 < 14; i3++) {
        d = 0.0;
        t29 = 0.0;
        t19 = 0.0;
        t21 = 0.0;
        for (invFx_j_i_tmp = 0; invFx_j_i_tmp < 8; invFx_j_i_tmp++) {
          kAcol = i2 + 14 * invFx_j_i_tmp;
          t23 = p_lambda_muu_j_i[kAcol];
          pipk = invFx_j_i_tmp + (i3 << 3);
          t25 = p_muu_Lambda_j_i[pipk];
          d += t23 * t25;
          t27 = p_x_muu_j_i[kAcol];
          t21 += t27 * t25;
          t25 = Aux_invFx_j_i[pipk];
          t29 += t27 * t25;
          t19 += t23 * t25;
        }
        kAcol = i2 + 14 * i3;
        pipk = kAcol + 196 * (5 - j);
        p_x_Lambda_i_data[pipk] = t21;
        p_lambda_Lambda_i_data[pipk] = d + Ix[kAcol];
        p_x_F_i_data[pipk] = t29 + invFx_j_i[kAcol];
        LAMBDA_i_data[pipk] = t19 - LAMBDAUncrt_j_i[kAcol];
      }
    }
    //  Coarse Iteration
    for (i2 = 0; i2 < 8; i2++) {
      d = 0.0;
      t29 = 0.0;
      for (i3 = 0; i3 < 14; i3++) {
        invFx_j_i_tmp = i2 + (i3 << 3);
        d += FuT_invFxT_j_i[invFx_j_i_tmp] * lambdaEq_j_i[i3];
        t29 += FuT_LAMBDAUncrt_m_Aux_invFx_j_i[invFx_j_i_tmp] * xEq_j_i[i3];
      }
      b_HAlluT_j_i[i2] = (HAlluT_j_i[i2] - d) + t29;
    }
    for (i2 = 0; i2 < 8; i2++) {
      d = 0.0;
      for (i3 = 0; i3 < 8; i3++) {
        d += Inv_dKKT23_mu_u_j_i[i2 + (i3 << 3)] * b_HAlluT_j_i[i3];
      }
      dmu_u_j_i[i2] = d;
    }
    for (i2 = 0; i2 < 196; i2++) {
      LAMBDAUncrt_j_i[i2] = -LAMBDAUncrt_j_i[i2];
    }
    for (i2 = 0; i2 < 14; i2++) {
      d = 0.0;
      t29 = 0.0;
      for (i3 = 0; i3 < 14; i3++) {
        invFx_j_i_tmp = i2 + 14 * i3;
        d += LAMBDAUncrt_j_i[invFx_j_i_tmp] * xEq_j_i[i3];
        t29 += Ix[invFx_j_i_tmp] * lambdaEq_j_i[i3];
      }
      t19 = 0.0;
      for (i3 = 0; i3 < 8; i3++) {
        t19 += reshapes_f2[i2 + 14 * i3] * dmu_u_j_i[i3];
      }
      i3 = i2 + b_i;
      lambda_i_data[i3] -= (d + t29) + t19;
    }
    for (i2 = 0; i2 < 8; i2++) {
      i3 = i2 + i1;
      u_i_data[i3] -= dmu_u_j_i[i2];
    }
    //  Recover
    //
    t18 = 0.0;
    for (int k{0}; k < 14; k++) {
      d = 0.0;
      for (i1 = 0; i1 < 8; i1++) {
        d += p_x_muu_j_i[k + 14 * i1] * dmu_u_j_i[i1];
      }
      t29 = 0.0;
      for (i1 = 0; i1 < 14; i1++) {
        t29 += invFx_j_i[k + 14 * i1] * xEq_j_i[i1];
      }
      i1 = k + b_i;
      x_i_data[i1] -= d + t29;
      std::copy(&Aux_invFx_j_i[k * 8],
                &Aux_invFx_j_i[static_cast<int>(
                    static_cast<unsigned int>(k * 8) + 8U)],
                &p_muu_F_i_data[(j * -112 + k * 8) + 560]);
      std::copy(&p_muu_Lambda_j_i[k * 8],
                &p_muu_Lambda_j_i[static_cast<int>(
                    static_cast<unsigned int>(k * 8) + 8U)],
                &p_muu_Lambda_i_data[(j * -112 + k * 8) + 560]);
      t17 = std::abs(xEq_j_i[k]);
      if (t17 > t18) {
        t18 = t17;
      }
    }
    KKTxEquation_i_data[5 - j] = t18;
    KKTC_i_data[5 - j] = 0.0;
    t18 = 0.0;
    for (int k{0}; k < 8; k++) {
      t17 = std::abs(HAlluT_j_i[k]);
      if (t17 > t18) {
        t18 = t17;
      }
    }
    KKTHu_i_data[5 - j] = t18;
    t18 = 0.0;
    for (int k{0}; k < 14; k++) {
      t17 = std::abs(lambdaEq_j_i[k]);
      if (t17 > t18) {
        t18 = t17;
      }
    }
    KKTlambdaEquation_i_data[5 - j] = t18;
  }
}

namespace coder {
namespace internal {
static int cfclose(double fid)
{
  std::FILE *f;
  int st;
  signed char b_fileid;
  signed char fileid;
  st = -1;
  fileid = static_cast<signed char>(fid);
  if ((static_cast<signed char>(fid) < 0) ||
      (fid != static_cast<signed char>(fid))) {
    fileid = -1;
  }
  b_fileid = fileid;
  if (fileid < 0) {
    b_fileid = -1;
  }
  if (b_fileid >= 3) {
    f = eml_openfiles[b_fileid - 3];
  } else if (b_fileid == 0) {
    f = stdin;
  } else if (b_fileid == 1) {
    f = stdout;
  } else if (b_fileid == 2) {
    f = stderr;
  } else {
    f = nullptr;
  }
  if ((f != nullptr) && (fileid >= 3)) {
    int cst;
    cst = std::fclose(f);
    if (cst == 0) {
      st = 0;
      eml_openfiles[fileid - 3] = nullptr;
      eml_autoflush[fileid - 3] = true;
    }
  }
  return st;
}

static signed char cfopen()
{
  std::FILE *filestar;
  signed char fileid;
  signed char j;
  fileid = -1;
  j = filedata();
  if (j >= 1) {
    filestar = std::fopen("GEN_log_rec.txt", "wb");
    if (filestar != nullptr) {
      int i;
      eml_openfiles[j - 1] = filestar;
      eml_autoflush[j - 1] = true;
      i = j + 2;
      if (j + 2 > 127) {
        i = 127;
      }
      fileid = static_cast<signed char>(i);
    }
  }
  return fileid;
}

static std::FILE *fileManager(double varargin_1, boolean_T &a)
{
  std::FILE *f;
  signed char fileid;
  fileid = static_cast<signed char>(varargin_1);
  if ((static_cast<signed char>(varargin_1) < 0) ||
      (varargin_1 != static_cast<signed char>(varargin_1))) {
    fileid = -1;
  }
  if (fileid >= 3) {
    a = eml_autoflush[fileid - 3];
    f = eml_openfiles[fileid - 3];
  } else if (fileid == 0) {
    f = stdin;
    a = true;
  } else if (fileid == 1) {
    f = stdout;
    a = true;
  } else if (fileid == 2) {
    f = stderr;
    a = true;
  } else {
    f = nullptr;
    a = true;
  }
  return f;
}

static signed char filedata()
{
  int k;
  signed char f;
  boolean_T exitg1;
  f = 0;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 20)) {
    if (eml_openfiles[k] == nullptr) {
      f = static_cast<signed char>(k + 1);
      exitg1 = true;
    } else {
      k++;
    }
  }
  return f;
}

namespace reflapack {
static int xzgetrf(double A[196], int ipiv[14])
{
  int i;
  int info;
  for (i = 0; i < 14; i++) {
    ipiv[i] = i + 1;
  }
  info = 0;
  for (int j{0}; j < 13; j++) {
    double smax;
    int a;
    int b_tmp;
    int jA;
    int jp1j;
    int mmj_tmp;
    mmj_tmp = 12 - j;
    b_tmp = j * 15;
    jp1j = b_tmp + 2;
    jA = 14 - j;
    a = 0;
    smax = std::abs(A[b_tmp]);
    for (int k{2}; k <= jA; k++) {
      double s;
      s = std::abs(A[(b_tmp + k) - 1]);
      if (s > smax) {
        a = k - 1;
        smax = s;
      }
    }
    if (A[b_tmp + a] != 0.0) {
      if (a != 0) {
        jA = j + a;
        ipiv[j] = jA + 1;
        for (int k{0}; k < 14; k++) {
          a = j + k * 14;
          smax = A[a];
          i = jA + k * 14;
          A[a] = A[i];
          A[i] = smax;
        }
      }
      i = (b_tmp - j) + 14;
      for (jA = jp1j; jA <= i; jA++) {
        A[jA - 1] /= A[b_tmp];
      }
    } else {
      info = j + 1;
    }
    jA = b_tmp;
    for (jp1j = 0; jp1j <= mmj_tmp; jp1j++) {
      smax = A[(b_tmp + jp1j * 14) + 14];
      if (smax != 0.0) {
        i = jA + 16;
        a = (jA - j) + 28;
        for (int k{i}; k <= a; k++) {
          A[k - 1] += A[((b_tmp + k) - jA) - 15] * -smax;
        }
      }
      jA += 14;
    }
  }
  if ((info == 0) && (A[195] == 0.0)) {
    info = 14;
  }
  return info;
}

} // namespace reflapack
} // namespace internal
static void inv(const double x[64], double y[64])
{
  double b_x[64];
  double smax;
  int i;
  int jA;
  int jp1j;
  int kAcol;
  int temp_tmp;
  signed char ipiv[8];
  signed char p[8];
  for (i = 0; i < 64; i++) {
    y[i] = 0.0;
    b_x[i] = x[i];
  }
  for (i = 0; i < 8; i++) {
    ipiv[i] = static_cast<signed char>(i + 1);
  }
  for (int j{0}; j < 7; j++) {
    int b_tmp;
    int mmj_tmp;
    mmj_tmp = 6 - j;
    b_tmp = j * 9;
    jp1j = b_tmp + 2;
    jA = 8 - j;
    kAcol = 0;
    smax = std::abs(b_x[b_tmp]);
    for (int k{2}; k <= jA; k++) {
      double s;
      s = std::abs(b_x[(b_tmp + k) - 1]);
      if (s > smax) {
        kAcol = k - 1;
        smax = s;
      }
    }
    if (b_x[b_tmp + kAcol] != 0.0) {
      if (kAcol != 0) {
        kAcol += j;
        ipiv[j] = static_cast<signed char>(kAcol + 1);
        for (int k{0}; k < 8; k++) {
          jA = k << 3;
          temp_tmp = j + jA;
          smax = b_x[temp_tmp];
          jA += kAcol;
          b_x[temp_tmp] = b_x[jA];
          b_x[jA] = smax;
        }
      }
      i = (b_tmp - j) + 8;
      for (int b_i{jp1j}; b_i <= i; b_i++) {
        b_x[b_i - 1] /= b_x[b_tmp];
      }
    }
    jA = b_tmp;
    for (kAcol = 0; kAcol <= mmj_tmp; kAcol++) {
      smax = b_x[(b_tmp + (kAcol << 3)) + 8];
      if (smax != 0.0) {
        i = jA + 10;
        jp1j = (jA - j) + 16;
        for (temp_tmp = i; temp_tmp <= jp1j; temp_tmp++) {
          b_x[temp_tmp - 1] += b_x[((b_tmp + temp_tmp) - jA) - 9] * -smax;
        }
      }
      jA += 8;
    }
  }
  for (i = 0; i < 8; i++) {
    p[i] = static_cast<signed char>(i + 1);
  }
  for (int k{0}; k < 7; k++) {
    signed char i1;
    i1 = ipiv[k];
    if (i1 > k + 1) {
      jA = p[i1 - 1];
      p[i1 - 1] = p[k];
      p[k] = static_cast<signed char>(jA);
    }
  }
  for (int k{0}; k < 8; k++) {
    temp_tmp = (p[k] - 1) << 3;
    y[k + temp_tmp] = 1.0;
    for (int j{k + 1}; j < 9; j++) {
      i = (j + temp_tmp) - 1;
      if (y[i] != 0.0) {
        jp1j = j + 1;
        for (int b_i{jp1j}; b_i < 9; b_i++) {
          jA = (b_i + temp_tmp) - 1;
          y[jA] -= y[i] * b_x[(b_i + ((j - 1) << 3)) - 1];
        }
      }
    }
  }
  for (int j{0}; j < 8; j++) {
    jA = j << 3;
    for (int k{7}; k >= 0; k--) {
      kAcol = k << 3;
      i = k + jA;
      smax = y[i];
      if (smax != 0.0) {
        y[i] = smax / b_x[k + kAcol];
        for (int b_i{0}; b_i < k; b_i++) {
          temp_tmp = b_i + jA;
          y[temp_tmp] -= y[i] * b_x[b_i + kAcol];
        }
      }
    }
  }
}

static void mldivide(double B[14])
{
  static const signed char b_A[196]{
      1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};
  double A[196];
  double temp;
  int ipiv[14];
  int i;
  int kAcol;
  for (i = 0; i < 196; i++) {
    A[i] = b_A[i];
  }
  internal::reflapack::xzgetrf(A, ipiv);
  for (int b_i{0}; b_i < 13; b_i++) {
    i = ipiv[b_i];
    if (i != b_i + 1) {
      temp = B[b_i];
      B[b_i] = B[i - 1];
      B[i - 1] = temp;
    }
  }
  for (int k{0}; k < 14; k++) {
    kAcol = 14 * k;
    if (B[k] != 0.0) {
      i = k + 2;
      for (int b_i{i}; b_i < 15; b_i++) {
        B[b_i - 1] -= B[k] * A[(b_i + kAcol) - 1];
      }
    }
  }
  for (int k{13}; k >= 0; k--) {
    kAcol = 14 * k;
    temp = B[k];
    if (temp != 0.0) {
      temp /= A[k + kAcol];
      B[k] = temp;
      for (int b_i{0}; b_i < k; b_i++) {
        B[b_i] -= B[k] * A[b_i + kAcol];
      }
    }
  }
}

} // namespace coder
static void filedata_init()
{
  for (int i{0}; i < 20; i++) {
    eml_autoflush[i] = false;
  }
  for (int i{0}; i < 20; i++) {
    eml_openfiles[i] = nullptr;
  }
}

static double fraction_to_boundary_parallel_func(
    const double u_i_data[], const double x_i_data[], double z_i_data[],
    const double u_k_i_data[], const double x_k_i_data[],
    const double z_k_i_data[], double rho, double &stepSizeMaxG_i)
{
  double b_u_i_data[174];
  double stepSizeG_i_data[174];
  double stepSizeZ_i_data[174];
  double d;
  double d1;
  double d2;
  double d3;
  double d4;
  double d5;
  double d6;
  double d7;
  double stepSizeMaxZ_i;
  int i;
  int i1;
  int i10;
  int i11;
  int i12;
  int i13;
  int i2;
  int i3;
  int i4;
  int i5;
  int i6;
  int i7;
  int i8;
  int i9;
  int j;
  //  line search variables
  //  update z
  for (j = 0; j < 6; j++) {
    double u_k_i[29];
    double z_i[29];
    // OCP_GEN_G
    //     G = OCP_GEN_G(IN1,IN2,IN3)
    //     This function was generated by the Symbolic Math Toolbox version 9.3.
    //     2023/06/09 12:24:34
    // OCP_GEN_G
    //     G = OCP_GEN_G(IN1,IN2,IN3)
    //     This function was generated by the Symbolic Math Toolbox version 9.3.
    //     2023/06/09 12:24:34
    d = u_i_data[8 * j];
    z_i[0] = z_i_data[29 * j] * (d + 10.0) - rho;
    i = 8 * j + 1;
    d1 = u_i_data[i];
    z_i[1] = z_i_data[29 * j + 1] * (d1 + 10.0) - rho;
    i1 = 8 * j + 2;
    d2 = u_i_data[i1];
    z_i[2] = z_i_data[29 * j + 2] * (d2 + 10.0) - rho;
    i2 = 8 * j + 3;
    d3 = u_i_data[i2];
    z_i[3] = z_i_data[29 * j + 3] * (d3 + 10.0) - rho;
    i3 = 8 * j + 4;
    d4 = u_i_data[i3];
    z_i[4] = z_i_data[29 * j + 4] * (d4 + 10.0) - rho;
    i4 = 8 * j + 5;
    d5 = u_i_data[i4];
    z_i[5] = z_i_data[29 * j + 5] * (d5 + 10.0) - rho;
    i5 = 8 * j + 6;
    d6 = u_i_data[i5];
    z_i[6] = z_i_data[29 * j + 6] * (d6 + 10.0) - rho;
    z_i[7] = z_i_data[29 * j + 7] * (-d + 10.0) - rho;
    z_i[8] = z_i_data[29 * j + 8] * (-d1 + 10.0) - rho;
    z_i[9] = z_i_data[29 * j + 9] * (-d2 + 10.0) - rho;
    z_i[10] = z_i_data[29 * j + 10] * (-d3 + 10.0) - rho;
    z_i[11] = z_i_data[29 * j + 11] * (-d4 + 10.0) - rho;
    z_i[12] = z_i_data[29 * j + 12] * (-d5 + 10.0) - rho;
    z_i[13] = z_i_data[29 * j + 13] * (-d6 + 10.0) - rho;
    i6 = 8 * j + 7;
    d = u_i_data[i6];
    z_i[14] = d * z_i_data[29 * j + 14] - rho;
    i7 = 14 * j + 7;
    d1 = x_i_data[i7];
    z_i[15] = z_i_data[29 * j + 15] * ((d + 1.5707963267948966) + d1) - rho;
    i8 = 14 * j + 8;
    d2 = x_i_data[i8];
    z_i[16] = z_i_data[29 * j + 16] * ((d + 1.5707963267948966) + d2) - rho;
    i9 = 14 * j + 9;
    d3 = x_i_data[i9];
    z_i[17] = z_i_data[29 * j + 17] * ((d + 1.5707963267948966) + d3) - rho;
    i10 = 14 * j + 10;
    d4 = x_i_data[i10];
    z_i[18] = z_i_data[29 * j + 18] * ((d + 1.5707963267948966) + d4) - rho;
    i11 = 14 * j + 11;
    d5 = x_i_data[i11];
    z_i[19] = z_i_data[29 * j + 19] * ((d + 1.5707963267948966) + d5) - rho;
    i12 = 14 * j + 12;
    d6 = x_i_data[i12];
    z_i[20] = z_i_data[29 * j + 20] * ((d + 1.5707963267948966) + d6) - rho;
    i13 = 14 * j + 13;
    d7 = x_i_data[i13];
    z_i[21] = z_i_data[29 * j + 21] * ((d + 1.5707963267948966) + d7) - rho;
    z_i[22] = z_i_data[29 * j + 22] * ((d + 1.5707963267948966) - d1) - rho;
    z_i[23] = z_i_data[29 * j + 23] * ((d + 1.5707963267948966) - d2) - rho;
    z_i[24] = z_i_data[29 * j + 24] * ((d + 1.5707963267948966) - d3) - rho;
    z_i[25] = z_i_data[29 * j + 25] * ((d + 1.5707963267948966) - d4) - rho;
    z_i[26] = z_i_data[29 * j + 26] * ((d + 1.5707963267948966) - d5) - rho;
    z_i[27] = z_i_data[29 * j + 27] * ((d + 1.5707963267948966) - d6) - rho;
    z_i[28] = z_i_data[29 * j + 28] * ((d + 1.5707963267948966) - d7) - rho;
    d = u_k_i_data[8 * j];
    u_k_i[0] = d + 10.0;
    d1 = u_k_i_data[i];
    u_k_i[1] = d1 + 10.0;
    d2 = u_k_i_data[i1];
    u_k_i[2] = d2 + 10.0;
    d3 = u_k_i_data[i2];
    u_k_i[3] = d3 + 10.0;
    d4 = u_k_i_data[i3];
    u_k_i[4] = d4 + 10.0;
    d5 = u_k_i_data[i4];
    u_k_i[5] = d5 + 10.0;
    d6 = u_k_i_data[i5];
    u_k_i[6] = d6 + 10.0;
    u_k_i[7] = -d + 10.0;
    u_k_i[8] = -d1 + 10.0;
    u_k_i[9] = -d2 + 10.0;
    u_k_i[10] = -d3 + 10.0;
    u_k_i[11] = -d4 + 10.0;
    u_k_i[12] = -d5 + 10.0;
    u_k_i[13] = -d6 + 10.0;
    d = u_k_i_data[i6];
    u_k_i[14] = d;
    d1 = x_k_i_data[i7];
    u_k_i[15] = (d + 1.5707963267948966) + d1;
    d2 = x_k_i_data[i8];
    u_k_i[16] = (d + 1.5707963267948966) + d2;
    d3 = x_k_i_data[i9];
    u_k_i[17] = (d + 1.5707963267948966) + d3;
    d4 = x_k_i_data[i10];
    u_k_i[18] = (d + 1.5707963267948966) + d4;
    d5 = x_k_i_data[i11];
    u_k_i[19] = (d + 1.5707963267948966) + d5;
    d6 = x_k_i_data[i12];
    u_k_i[20] = (d + 1.5707963267948966) + d6;
    d7 = x_k_i_data[i13];
    u_k_i[21] = (d + 1.5707963267948966) + d7;
    u_k_i[22] = (d + 1.5707963267948966) - d1;
    u_k_i[23] = (d + 1.5707963267948966) - d2;
    u_k_i[24] = (d + 1.5707963267948966) - d3;
    u_k_i[25] = (d + 1.5707963267948966) - d4;
    u_k_i[26] = (d + 1.5707963267948966) - d5;
    u_k_i[27] = (d + 1.5707963267948966) - d6;
    u_k_i[28] = (d + 1.5707963267948966) - d7;
    for (i = 0; i < 29; i++) {
      i1 = i + 29 * j;
      d = z_i_data[i1] - z_i[i] / u_k_i[i];
      z_i[i] = d;
      z_i_data[i1] = d;
    }
  }
  //         %% Line Search for feasibility
  //  z
  for (j = 0; j < 174; j++) {
    d = z_k_i_data[j];
    d = -0.95 * (d / (z_i_data[j] - d));
    stepSizeZ_i_data[j] = d;
    if ((d > 1.0) || (d < 0.0)) {
      stepSizeZ_i_data[j] = 1.0;
    }
  }
  stepSizeMaxZ_i = stepSizeZ_i_data[0];
  for (j = 0; j < 173; j++) {
    d = stepSizeZ_i_data[j + 1];
    if (stepSizeMaxZ_i > d) {
      stepSizeMaxZ_i = d;
    }
  }
  //  G
  // OCP_GEN_G
  //     G = OCP_GEN_G(IN1,IN2,IN3)
  //     This function was generated by the Symbolic Math Toolbox version 9.3.
  //     2023/06/09 12:24:34
  //  GMin
  // OCP_GEN_G
  //     G = OCP_GEN_G(IN1,IN2,IN3)
  //     This function was generated by the Symbolic Math Toolbox version 9.3.
  //     2023/06/09 12:24:34
  for (i = 0; i < 6; i++) {
    int G_k_i_data_tmp;
    int ab_G_k_i_data_tmp;
    int b_G_k_i_data_tmp;
    int bb_G_k_i_data_tmp;
    int c_G_k_i_data_tmp;
    int d_G_k_i_data_tmp;
    int e_G_k_i_data_tmp;
    int f_G_k_i_data_tmp;
    int g_G_k_i_data_tmp;
    int h_G_k_i_data_tmp;
    int i14;
    int i_G_k_i_data_tmp;
    int j_G_k_i_data_tmp;
    int k_G_k_i_data_tmp;
    int l_G_k_i_data_tmp;
    int m_G_k_i_data_tmp;
    int n_G_k_i_data_tmp;
    int o_G_k_i_data_tmp;
    int p_G_k_i_data_tmp;
    int q_G_k_i_data_tmp;
    int r_G_k_i_data_tmp;
    int s_G_k_i_data_tmp;
    int t_G_k_i_data_tmp;
    int u_G_k_i_data_tmp;
    int v_G_k_i_data_tmp;
    int w_G_k_i_data_tmp;
    int x_G_k_i_data_tmp;
    int y_G_k_i_data_tmp;
    d = u_k_i_data[8 * i];
    stepSizeZ_i_data[29 * i] = d + 10.0;
    i1 = 8 * i + 1;
    d1 = u_k_i_data[i1];
    j = 29 * i + 1;
    stepSizeZ_i_data[j] = d1 + 10.0;
    i2 = 8 * i + 2;
    d2 = u_k_i_data[i2];
    G_k_i_data_tmp = 29 * i + 2;
    stepSizeZ_i_data[G_k_i_data_tmp] = d2 + 10.0;
    i3 = 8 * i + 3;
    d3 = u_k_i_data[i3];
    b_G_k_i_data_tmp = 29 * i + 3;
    stepSizeZ_i_data[b_G_k_i_data_tmp] = d3 + 10.0;
    i4 = 8 * i + 4;
    d4 = u_k_i_data[i4];
    c_G_k_i_data_tmp = 29 * i + 4;
    stepSizeZ_i_data[c_G_k_i_data_tmp] = d4 + 10.0;
    i5 = 8 * i + 5;
    d5 = u_k_i_data[i5];
    d_G_k_i_data_tmp = 29 * i + 5;
    stepSizeZ_i_data[d_G_k_i_data_tmp] = d5 + 10.0;
    i6 = 8 * i + 6;
    d6 = u_k_i_data[i6];
    e_G_k_i_data_tmp = 29 * i + 6;
    stepSizeZ_i_data[e_G_k_i_data_tmp] = d6 + 10.0;
    f_G_k_i_data_tmp = 29 * i + 7;
    stepSizeZ_i_data[f_G_k_i_data_tmp] = -d + 10.0;
    g_G_k_i_data_tmp = 29 * i + 8;
    stepSizeZ_i_data[g_G_k_i_data_tmp] = -d1 + 10.0;
    h_G_k_i_data_tmp = 29 * i + 9;
    stepSizeZ_i_data[h_G_k_i_data_tmp] = -d2 + 10.0;
    i_G_k_i_data_tmp = 29 * i + 10;
    stepSizeZ_i_data[i_G_k_i_data_tmp] = -d3 + 10.0;
    j_G_k_i_data_tmp = 29 * i + 11;
    stepSizeZ_i_data[j_G_k_i_data_tmp] = -d4 + 10.0;
    k_G_k_i_data_tmp = 29 * i + 12;
    stepSizeZ_i_data[k_G_k_i_data_tmp] = -d5 + 10.0;
    l_G_k_i_data_tmp = 29 * i + 13;
    stepSizeZ_i_data[l_G_k_i_data_tmp] = -d6 + 10.0;
    i7 = 8 * i + 7;
    d = u_k_i_data[i7];
    m_G_k_i_data_tmp = 29 * i + 14;
    stepSizeZ_i_data[m_G_k_i_data_tmp] = d;
    i8 = 14 * i + 7;
    d1 = x_k_i_data[i8];
    n_G_k_i_data_tmp = 29 * i + 15;
    stepSizeZ_i_data[n_G_k_i_data_tmp] = (d + 1.5707963267948966) + d1;
    i9 = 14 * i + 8;
    d2 = x_k_i_data[i9];
    o_G_k_i_data_tmp = 29 * i + 16;
    stepSizeZ_i_data[o_G_k_i_data_tmp] = (d + 1.5707963267948966) + d2;
    i10 = 14 * i + 9;
    d3 = x_k_i_data[i10];
    p_G_k_i_data_tmp = 29 * i + 17;
    stepSizeZ_i_data[p_G_k_i_data_tmp] = (d + 1.5707963267948966) + d3;
    i11 = 14 * i + 10;
    d4 = x_k_i_data[i11];
    q_G_k_i_data_tmp = 29 * i + 18;
    stepSizeZ_i_data[q_G_k_i_data_tmp] = (d + 1.5707963267948966) + d4;
    i12 = 14 * i + 11;
    d5 = x_k_i_data[i12];
    r_G_k_i_data_tmp = 29 * i + 19;
    stepSizeZ_i_data[r_G_k_i_data_tmp] = (d + 1.5707963267948966) + d5;
    i13 = 14 * i + 12;
    d6 = x_k_i_data[i13];
    s_G_k_i_data_tmp = 29 * i + 20;
    stepSizeZ_i_data[s_G_k_i_data_tmp] = (d + 1.5707963267948966) + d6;
    i14 = 14 * i + 13;
    d7 = x_k_i_data[i14];
    t_G_k_i_data_tmp = 29 * i + 21;
    stepSizeZ_i_data[t_G_k_i_data_tmp] = (d + 1.5707963267948966) + d7;
    u_G_k_i_data_tmp = 29 * i + 22;
    stepSizeZ_i_data[u_G_k_i_data_tmp] = (d + 1.5707963267948966) - d1;
    v_G_k_i_data_tmp = 29 * i + 23;
    stepSizeZ_i_data[v_G_k_i_data_tmp] = (d + 1.5707963267948966) - d2;
    w_G_k_i_data_tmp = 29 * i + 24;
    stepSizeZ_i_data[w_G_k_i_data_tmp] = (d + 1.5707963267948966) - d3;
    x_G_k_i_data_tmp = 29 * i + 25;
    stepSizeZ_i_data[x_G_k_i_data_tmp] = (d + 1.5707963267948966) - d4;
    y_G_k_i_data_tmp = 29 * i + 26;
    stepSizeZ_i_data[y_G_k_i_data_tmp] = (d + 1.5707963267948966) - d5;
    ab_G_k_i_data_tmp = 29 * i + 27;
    stepSizeZ_i_data[ab_G_k_i_data_tmp] = (d + 1.5707963267948966) - d6;
    bb_G_k_i_data_tmp = 29 * i + 28;
    stepSizeZ_i_data[bb_G_k_i_data_tmp] = (d + 1.5707963267948966) - d7;
    d = u_i_data[8 * i];
    b_u_i_data[29 * i] = d + 10.0;
    d1 = u_i_data[i1];
    b_u_i_data[j] = d1 + 10.0;
    d2 = u_i_data[i2];
    b_u_i_data[G_k_i_data_tmp] = d2 + 10.0;
    d3 = u_i_data[i3];
    b_u_i_data[b_G_k_i_data_tmp] = d3 + 10.0;
    d4 = u_i_data[i4];
    b_u_i_data[c_G_k_i_data_tmp] = d4 + 10.0;
    d5 = u_i_data[i5];
    b_u_i_data[d_G_k_i_data_tmp] = d5 + 10.0;
    d6 = u_i_data[i6];
    b_u_i_data[e_G_k_i_data_tmp] = d6 + 10.0;
    b_u_i_data[f_G_k_i_data_tmp] = -d + 10.0;
    b_u_i_data[g_G_k_i_data_tmp] = -d1 + 10.0;
    b_u_i_data[h_G_k_i_data_tmp] = -d2 + 10.0;
    b_u_i_data[i_G_k_i_data_tmp] = -d3 + 10.0;
    b_u_i_data[j_G_k_i_data_tmp] = -d4 + 10.0;
    b_u_i_data[k_G_k_i_data_tmp] = -d5 + 10.0;
    b_u_i_data[l_G_k_i_data_tmp] = -d6 + 10.0;
    d = u_i_data[i7];
    b_u_i_data[m_G_k_i_data_tmp] = d;
    d1 = x_i_data[i8];
    b_u_i_data[n_G_k_i_data_tmp] = (d + 1.5707963267948966) + d1;
    d2 = x_i_data[i9];
    b_u_i_data[o_G_k_i_data_tmp] = (d + 1.5707963267948966) + d2;
    d3 = x_i_data[i10];
    b_u_i_data[p_G_k_i_data_tmp] = (d + 1.5707963267948966) + d3;
    d4 = x_i_data[i11];
    b_u_i_data[q_G_k_i_data_tmp] = (d + 1.5707963267948966) + d4;
    d5 = x_i_data[i12];
    b_u_i_data[r_G_k_i_data_tmp] = (d + 1.5707963267948966) + d5;
    d6 = x_i_data[i13];
    b_u_i_data[s_G_k_i_data_tmp] = (d + 1.5707963267948966) + d6;
    d7 = x_i_data[i14];
    b_u_i_data[t_G_k_i_data_tmp] = (d + 1.5707963267948966) + d7;
    b_u_i_data[u_G_k_i_data_tmp] = (d + 1.5707963267948966) - d1;
    b_u_i_data[v_G_k_i_data_tmp] = (d + 1.5707963267948966) - d2;
    b_u_i_data[w_G_k_i_data_tmp] = (d + 1.5707963267948966) - d3;
    b_u_i_data[x_G_k_i_data_tmp] = (d + 1.5707963267948966) - d4;
    b_u_i_data[y_G_k_i_data_tmp] = (d + 1.5707963267948966) - d5;
    b_u_i_data[ab_G_k_i_data_tmp] = (d + 1.5707963267948966) - d6;
    b_u_i_data[bb_G_k_i_data_tmp] = (d + 1.5707963267948966) - d7;
  }
  for (j = 0; j < 174; j++) {
    d = stepSizeZ_i_data[j];
    d = -0.95 * (d / (b_u_i_data[j] - d));
    stepSizeG_i_data[j] = d;
    if ((d > 1.0) || (d < 0.0)) {
      stepSizeG_i_data[j] = 1.0;
    }
  }
  stepSizeMaxG_i = stepSizeG_i_data[0];
  for (j = 0; j < 173; j++) {
    d = stepSizeG_i_data[j + 1];
    if (stepSizeMaxG_i > d) {
      stepSizeMaxG_i = d;
    }
  }
  return stepSizeMaxZ_i;
}

void Simu_Matlab()
{
  static const double b[7]{1.5707963267948966, 0.0, 1.5707963267948966, 0.0,
                           1.5707963267948966, 0.0, 1.5707963267948966};
  static const double b_b[7]{0.0, 1.5707963267948966, 0.0, 1.5707963267948966,
                             0.0, 1.5707963267948966, 0.0};
  static double rec_x[112014];
  static double rec_u[64000];
  static double rec_cpuTime[8000];
  static double rec_cpuTimeSearchDirection[8000];
  static double rec_error[8000];
  static double rec_numIter[8000];
  static double d_expl_temp[4704];
  std::FILE *filestar;
  double x0[14];
  double x0Measured[14];
  double e_expl_temp;
  double f_expl_temp;
  double g_expl_temp;
  double h_expl_temp;
  double output_KKTError;
  double output_iterTotal;
  double output_timeElapsed_searchDirection;
  double output_timeElapsed_total;
  double tTotal;
  char varargin_1[3];
  signed char fileid;
  boolean_T autoflush;
  if (!isInitialized_Simu_Matlab) {
    Simu_Matlab_initialize();
  }
  //  For code generation
  //  sampling interval
  //  Load data
  std::memset(&x0[0], 0, 14U * sizeof(double));
  //  define record variables
  for (int j{0}; j < 112014; j++) {
    rec_x[j] = 1.0;
  }
  for (int j{0}; j < 14; j++) {
    rec_x[8001 * j] = 0.0;
  }
  //  Simulation
  //  Create options for NMPC_Solve
  //     %% options for the initial rho
  //  initial rho
  //  max number of iterations for the initial rho problem
  //  KKT tolerence for the initial rho problem
  //     %% barrier parameter decaying rate
  //     %% options for the target rho
  //  target/end rho
  //  max number of iterations for the initial rho problem
  //  KKT tolerence for the end rho problem
  //     %% line search parameters
  //  enable or disable line search
  //  line search method
  //     %% degree of parallism
  //  1: in serial, otherwise in parallel
  //     %% display
  //     %% check KKT error after iteration
  //  whether to check the KKT error after iteration
  //  degree of parallism: 1 = in serial, otherwise in parallel
  //  do not check the KKT error after iteration
  //  init
  tTotal = 0.0;
  for (int step{0}; step < 8000; step++) {
    double c_expl_temp[696];
    double b_expl_temp[336];
    double expl_temp[336];
    double solution_u[192];
    double p[168];
    double unnamed_idx_0_tmp;
    // simulation steps
    //  Set the reference angle q_ref with a rate limit
    //  0<t<4: q_ref = [0,   pi/2,0,    pi/2, 0,    pi/2, 0]
    //  4<t<8: q_ref = [pi/2,0,   pi/2, 0,    pi/2, 0,    pi/2]
    if (step + 1 < 4000) {
      double ref[7];
      unnamed_idx_0_tmp = (static_cast<double>(step) + 1.0) * 0.008;
      if (unnamed_idx_0_tmp > 1.0) {
        unnamed_idx_0_tmp = 1.0;
      }
      for (int iRef{0}; iRef < 7; iRef++) {
        ref[iRef] = unnamed_idx_0_tmp * b_b[iRef];
        for (int j{0}; j < 24; j++) {
          p[iRef + 7 * j] = ref[iRef];
        }
      }
    } else {
      double ref[7];
      unnamed_idx_0_tmp = ((static_cast<double>(step) + 1.0) - 4000.0) * 0.008;
      if (unnamed_idx_0_tmp > 1.0) {
        unnamed_idx_0_tmp = 1.0;
      }
      for (int iRef{0}; iRef < 7; iRef++) {
        ref[iRef] =
            unnamed_idx_0_tmp * b[iRef] + (1.0 - unnamed_idx_0_tmp) * b_b[iRef];
        for (int j{0}; j < 24; j++) {
          p[iRef + 7 * j] = ref[iRef];
        }
      }
    }
    //  Solve the optimal control problem
    std::copy(&x0[0], &x0[14], &x0Measured[0]);
    NMPC_Solve(x0, p, expl_temp, solution_u, b_expl_temp, c_expl_temp,
               d_expl_temp, unnamed_idx_0_tmp, output_KKTError,
               output_timeElapsed_searchDirection, e_expl_temp, f_expl_temp,
               output_timeElapsed_total, g_expl_temp, output_iterTotal,
               h_expl_temp);
    tTotal += output_timeElapsed_total;
    //  Obtain the first optimal control input
    //  System simulation by the 4th-order Explicit Runge-Kutta Method
    SIM_Plant_RK4(&(*(double(*)[8]) & solution_u[0])[0], x0Measured, x0);
    //  Record data
    for (int j{0}; j < 14; j++) {
      rec_x[(step + 8001 * j) + 1] = x0Measured[j];
    }
    for (int j{0}; j < 8; j++) {
      rec_u[step + 8000 * j] = solution_u[j];
    }
    rec_error[step] = output_KKTError;
    rec_cpuTime[step] = output_timeElapsed_total * 1.0E+6;
    rec_cpuTimeSearchDirection[step] =
        output_timeElapsed_searchDirection * 1.0E+6;
    rec_numIter[step] = output_iterTotal;
  }
  const char *fmt1;
  //  Log to file
  //  Code generation
  //  show Time Elapsed for RTI
  fmt1 = "Time Elapsed for NMPC_Solve: %f s\r\n";
  printf(fmt1, tTotal);
  //  Log to file
  fileid = coder::internal::cfopen();
  //  printf header
  for (int j{0}; j < 14; j++) {
    varargin_1[0] = 'x';
    varargin_1[1] = static_cast<char>((static_cast<double>(j) + 1.0) + 48.0);
    varargin_1[2] = '\x00';
    filestar =
        coder::internal::fileManager(static_cast<double>(fileid), autoflush);
    if (!(filestar == nullptr)) {
      std::fprintf(filestar, "%s\t", &varargin_1[0]);
      if (autoflush) {
        std::fflush(filestar);
      }
    }
  }
  for (int j{0}; j < 8; j++) {
    varargin_1[0] = 'u';
    varargin_1[1] = static_cast<char>((static_cast<double>(j) + 1.0) + 48.0);
    varargin_1[2] = '\x00';
    filestar =
        coder::internal::fileManager(static_cast<double>(fileid), autoflush);
    if (!(filestar == nullptr)) {
      std::fprintf(filestar, "%s\t", &varargin_1[0]);
      if (autoflush) {
        std::fflush(filestar);
      }
    }
  }
  filestar =
      coder::internal::fileManager(static_cast<double>(fileid), autoflush);
  if (!(filestar == nullptr)) {
    std::fprintf(filestar, "%s\t", "error");
    if (autoflush) {
      std::fflush(filestar);
    }
  }
  filestar =
      coder::internal::fileManager(static_cast<double>(fileid), autoflush);
  if (!(filestar == nullptr)) {
    std::fprintf(filestar, "%s\t", "numIter");
    if (autoflush) {
      std::fflush(filestar);
    }
  }
  filestar =
      coder::internal::fileManager(static_cast<double>(fileid), autoflush);
  if (!(filestar == nullptr)) {
    std::fprintf(filestar, "%s\t", "cpuTime");
    if (autoflush) {
      std::fflush(filestar);
    }
  }
  filestar =
      coder::internal::fileManager(static_cast<double>(fileid), autoflush);
  if (!(filestar == nullptr)) {
    std::fprintf(filestar, "%s\t", "cpuTimeSearchDirection");
    if (autoflush) {
      std::fflush(filestar);
    }
  }
  filestar =
      coder::internal::fileManager(static_cast<double>(fileid), autoflush);
  if (!(filestar == nullptr)) {
    std::fprintf(filestar, "%s\t", "cpuTimeLineSearch");
    if (autoflush) {
      std::fflush(filestar);
    }
  }
  filestar =
      coder::internal::fileManager(static_cast<double>(fileid), autoflush);
  if (!(filestar == nullptr)) {
    std::fprintf(filestar, "%s\n", "cpuTimeKKTError");
    if (autoflush) {
      std::fflush(filestar);
    }
  }
  //  printf data
  for (int step{0}; step < 8000; step++) {
    for (int j{0}; j < 14; j++) {
      filestar =
          coder::internal::fileManager(static_cast<double>(fileid), autoflush);
      if (!(filestar == nullptr)) {
        std::fprintf(filestar, "%f\t", rec_x[step + 8001 * j]);
        if (autoflush) {
          std::fflush(filestar);
        }
      }
    }
    for (int j{0}; j < 8; j++) {
      filestar =
          coder::internal::fileManager(static_cast<double>(fileid), autoflush);
      if (!(filestar == nullptr)) {
        std::fprintf(filestar, "%f\t", rec_u[step + 8000 * j]);
        if (autoflush) {
          std::fflush(filestar);
        }
      }
    }
    filestar =
        coder::internal::fileManager(static_cast<double>(fileid), autoflush);
    if (!(filestar == nullptr)) {
      std::fprintf(filestar, "%f\t", rec_error[step]);
      if (autoflush) {
        std::fflush(filestar);
      }
    }
    filestar =
        coder::internal::fileManager(static_cast<double>(fileid), autoflush);
    if (!(filestar == nullptr)) {
      std::fprintf(filestar, "%f\t", rec_numIter[step]);
      if (autoflush) {
        std::fflush(filestar);
      }
    }
    filestar =
        coder::internal::fileManager(static_cast<double>(fileid), autoflush);
    if (!(filestar == nullptr)) {
      std::fprintf(filestar, "%f\t", rec_cpuTime[step]);
      if (autoflush) {
        std::fflush(filestar);
      }
    }
    filestar =
        coder::internal::fileManager(static_cast<double>(fileid), autoflush);
    if (!(filestar == nullptr)) {
      std::fprintf(filestar, "%f\t", rec_cpuTimeSearchDirection[step]);
      if (autoflush) {
        std::fflush(filestar);
      }
    }
    filestar =
        coder::internal::fileManager(static_cast<double>(fileid), autoflush);
    if (!(filestar == nullptr)) {
      std::fprintf(filestar, "%f\t", 0.0);
      if (autoflush) {
        std::fflush(filestar);
      }
    }
    filestar =
        coder::internal::fileManager(static_cast<double>(fileid), autoflush);
    if (!(filestar == nullptr)) {
      std::fprintf(filestar, "%f\n", 0.0);
      if (autoflush) {
        std::fflush(filestar);
      }
    }
  }
  coder::internal::cfclose(static_cast<double>(fileid));
}

void Simu_Matlab_initialize()
{
  omp_init_nest_lock(&Simu_Matlab_nestLockGlobal);
  NMPC_Solve_init();
  filedata_init();
  // user code (Initialize function Body)
  {
    iiwa14_init();
  }

  isInitialized_Simu_Matlab = true;
}

void Simu_Matlab_terminate()
{
  omp_destroy_nest_lock(&Simu_Matlab_nestLockGlobal);
  isInitialized_Simu_Matlab = false;
}

// End of code generation (Simu_Matlab.cpp)
