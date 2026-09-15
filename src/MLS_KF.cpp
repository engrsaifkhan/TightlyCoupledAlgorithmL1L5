#include "MLS_KF.h"
#include "ECEF_to_NED.h"
#include "Skew_Symmetric.h"
#include "Gravity_ECEF.h"
#include <math.h>
#include <string.h>
#include <cmath>
#include <stddef.h>
#include <Arduino.h>

#include "All_Configs.h"
#include "Globals.h"

static constexpr size_t delta_z_DIM = 2 * MAX_GNSS_ROWS_PER_EPOCH + 1;
static constexpr size_t MLS_U_AS_E_T_SIZE  = 3 * MAX_GNSS_ROWS_PER_EPOCH;
static constexpr size_t MLS_PRED_MEAS_SIZE = 2 * MAX_GNSS_ROWS_PER_EPOCH;
static constexpr size_t MLS_HROWS_MAX      = 2 * MAX_GNSS_ROWS_PER_EPOCH + 1;
static constexpr size_t MLS_H_MATRIX_SIZE  = MLS_HROWS_MAX * 17;

static double KG_PHt[17 * MLS_HROWS_MAX];
static double KG_S[MLS_HROWS_MAX * MLS_HROWS_MAX];
static double KG_rhs[MLS_HROWS_MAX];
static int KG_piv[MLS_HROWS_MAX];
static double MLS_K_matrix[17 * MLS_HROWS_MAX];
// MLS_KF_Epoch temporary PSRAM workspaces
static double MLS_Q_prime_matrix[P_DIM] = {0.0};                    //  P_DIM = 17*17
static double MLS_P_matrix_propagated[P_DIM] = {0.0};               //  P_DIM = 17*17
static double MLS_u_as_e_T[MLS_U_AS_E_T_SIZE] = {0.0};
static double MLS_pred_meas[MLS_PRED_MEAS_SIZE] = {0.0};
static double MLS_H_matrix[MLS_H_MATRIX_SIZE] = {0.0};
static double MLS_delta_z[delta_z_DIM] = {0.0};

#define DEBUG_KALMAN_GAIN 0
// Column-major index: A(row, col) = A[row + ld * col]
static inline size_t CM(const size_t row, const size_t col, const size_t ld) {
  return row + ld * col;
}

static inline void Mat3MulCM(const double A[9],
                             const double B[9],
                             double C[9])
{
  for (int col = 0; col < 3; ++col) {
    for (int row = 0; row < 3; ++row) {
      C[row + 3 * col] =
          A[row + 3 * 0] * B[0 + 3 * col] +
          A[row + 3 * 1] * B[1 + 3 * col] +
          A[row + 3 * 2] * B[2 + 3 * col];
    }
  }
}


static inline void Mat3VecCM(const double A[9],
                             const double x[3],
                             double y[3])
{
  y[0] = A[0] * x[0] + A[3] * x[1] + A[6] * x[2];
  y[1] = A[1] * x[0] + A[4] * x[1] + A[7] * x[2];
  y[2] = A[2] * x[0] + A[5] * x[1] + A[8] * x[2];
}

void PrintKalmanGainMatrix(const double* K_matrix,
                           size_t state_dim,
                           size_t meas_dim) {
 #if DEBUG_KALMAN_GAIN

  Serial.println(F("========== Kalman Gain Matrix Debug =========="));

  Serial.print(F("state_dim = "));
  Serial.println(state_dim);

  Serial.print(F("meas_dim = "));
  Serial.println(meas_dim);

  for (size_t r = 0; r < state_dim; r++) {
    Serial.print(F("K row "));
    Serial.print(r);
    Serial.print(F(": "));

    for (size_t c = 0; c < meas_dim; c++) {
      double val = K_matrix[r + state_dim * c];

      Serial.print(val, 9);

      if (c + 1 < meas_dim) {
        Serial.print(F(", "));
      }
    }

    Serial.println();
  }

  Serial.println(F("=============================================="));

 #endif
}
static bool Cholesky_Lower_InPlace(double *A, const size_t n, const size_t ld) {
  for (size_t j = 0; j < n; ++j) {
    double d = A[CM(j, j, ld)];

    for (size_t k = 0; k < j; ++k) {
      const double Ljk = A[CM(j, k, ld)];
      d -= Ljk * Ljk;
    }

    if (!(d > 0.0) || !isfinite(d)) {
      return false;
    }

    const double Ljj = sqrt(d);
    A[CM(j, j, ld)] = Ljj;

    for (size_t i = j + 1; i < n; ++i) {
      double s = A[CM(i, j, ld)];

      for (size_t k = 0; k < j; ++k) {
        s -= A[CM(i, k, ld)] * A[CM(j, k, ld)];
      }

      A[CM(i, j, ld)] = s / Ljj;
    }

    for (size_t k = j + 1; k < n; ++k) {
      A[CM(j, k, ld)] = 0.0;
    }
  }

  return true;
}

static void Cholesky_Lower_Solve_InPlace(const double *L, const size_t n, const size_t ld, double *b) {
  // Solve L * y = b
  for (size_t i = 0; i < n; ++i) {
    double s = b[i];

    for (size_t k = 0; k < i; ++k) {
      s -= L[CM(i, k, ld)] * b[k];
    }

    b[i] = s / L[CM(i, i, ld)];
  }

  // Solve L' * x = y
  for (size_t ii = n; ii-- > 0;) {
    double s = b[ii];

    for (size_t k = ii + 1; k < n; ++k) {
      s -= L[CM(k, ii, ld)] * b[k];
    }

    b[ii] = s / L[CM(ii, ii, ld)];
  }
}

static bool LU_Factor_InPlace(double *A, const size_t n, const size_t ld, int *piv) {
  for (size_t k = 0; k < n; ++k) {
    size_t pivot = k;
    double max_abs = fabs(A[CM(k, k, ld)]);

    for (size_t i = k + 1; i < n; ++i) {
      const double v = fabs(A[CM(i, k, ld)]);

      if (v > max_abs) {
        max_abs = v;
        pivot = i;
      }
    }

    if (!(max_abs > 0.0) || !isfinite(max_abs)) {
      return false;
    }

    piv[k] = static_cast<int>(pivot);

    if (pivot != k) {
      for (size_t col = 0; col < n; ++col) {
        const size_t idx1 = CM(k, col, ld);
        const size_t idx2 = CM(pivot, col, ld);

        const double tmp = A[idx1];
        A[idx1] = A[idx2];
        A[idx2] = tmp;
      }
    }

    for (size_t i = k + 1; i < n; ++i) {
      A[CM(i, k, ld)] /= A[CM(k, k, ld)];

      const double Lik = A[CM(i, k, ld)];

      for (size_t col = k + 1; col < n; ++col) {
        A[CM(i, col, ld)] -= Lik * A[CM(k, col, ld)];
      }
    }
  }

  return true;
}

static bool LU_Solve_InPlace(const double *LU, const size_t n, const size_t ld, const int *piv, double *b) {
  for (size_t k = 0; k < n; ++k) {
    const size_t p = static_cast<size_t>(piv[k]);

    if (p != k) {
      const double tmp = b[k];
      b[k] = b[p];
      b[p] = tmp;
    }
  }

  for (size_t i = 0; i < n; ++i) {
    double s = b[i];

    for (size_t k = 0; k < i; ++k) {
      s -= LU[CM(i, k, ld)] * b[k];
    }

    b[i] = s;
  }

  for (size_t ii = n; ii-- > 0;) {
    double s = b[ii];

    for (size_t k = ii + 1; k < n; ++k) {
      s -= LU[CM(ii, k, ld)] * b[k];
    }

    const double diag = LU[CM(ii, ii, ld)];

    if (!(fabs(diag) > 0.0) || !isfinite(diag)) {
      return false;
    }

    b[ii] = s / diag;
  }

  return true;
}

static void Build_PHt(const double *P_matrix_propagated, const double *H_matrix, const size_t H_rows) {
  // PHt = P_matrix_propagated * H_matrix'
  // PHt size = 17 x H_rows

  for (size_t col = 0; col < H_rows; ++col) {
    for (size_t row = 0; row < 17; ++row) {
      double sum = 0.0;

      for (size_t k = 0; k < 17; ++k) {
        // H_matrix'(k, col) = H_matrix(col, k)
        sum += P_matrix_propagated[CM(row, k, 17)] *  H_matrix[CM(col, k, H_rows)];
      }

      KG_PHt[CM(row, col, 17)] = sum;
    }
  }
}

static void Build_S_Matrix(const double *H_matrix, const double *R_matrix, const size_t H_rows, const size_t R_ld) {
  for (size_t col = 0; col < H_rows; ++col) {
    for (size_t row = 0; row < H_rows; ++row) {
      double sum = R_matrix[CM(row, col, R_ld)];

      for (size_t k = 0; k < 17; ++k) {
        sum += H_matrix[CM(row, k, H_rows)] *
               KG_PHt[CM(k, col, 17)];
      }

      // KG_S is still stored using MLS_HROWS_MAX as leading dimension
      KG_S[CM(row, col, MLS_HROWS_MAX)] = sum;
    }
  }
}

static void Symmetrize_S(const size_t H_rows) {
    for (size_t col = 0; col < H_rows; ++col) {
      for (size_t row = col + 1; row < H_rows; ++row) {
        const double a = 0.5 * (
          KG_S[CM(row, col, MLS_HROWS_MAX)] +
          KG_S[CM(col, row, MLS_HROWS_MAX)]
        );

        KG_S[CM(row, col, MLS_HROWS_MAX)] = a;
        KG_S[CM(col, row, MLS_HROWS_MAX)] = a;
      }
    }
}

static bool Kalman_Gain_Solver_PSram(const double *P_matrix_propagated, const double *H_matrix, const double *R_matrix, const size_t H_rows, const size_t R_ld, double *K_out) {
  if (H_rows == 0 || H_rows > MLS_HROWS_MAX) {
    Serial.println("Kalman Gain Error: invalid H_rows");
    return false;
  }
  if (R_ld < H_rows) {
    Serial.println("Kalman Gain Error: invalid R leading dimension");
    return false;
  }

  // ------------------------------------------------------------
  // 1. PHt = P * H'
  // ------------------------------------------------------------
  Build_PHt(P_matrix_propagated, H_matrix, H_rows);

  // ------------------------------------------------------------
  // 2. S = H * PHt + R
  // ------------------------------------------------------------
  Build_S_Matrix(H_matrix, R_matrix, H_rows, R_ld);

  // ------------------------------------------------------------
  // 3. Try Cholesky first
  // ------------------------------------------------------------
  Symmetrize_S(H_rows);

  // Dimensionality Check
  for (size_t i = 0; i < H_rows; ++i) {
    const double sii = KG_S[CM(i, i, MLS_HROWS_MAX)];
    if (!(sii > 0.0) || !isfinite(sii)) {
      Serial.print("Kalman Gain Error: bad S diagonal at ");
      Serial.println(i);
      return false;
    }
  }

  if (Cholesky_Lower_InPlace(KG_S, H_rows, MLS_HROWS_MAX)) {
    for (size_t state = 0; state < 17; ++state) {
      for (size_t r = 0; r < H_rows; ++r) {
        KG_rhs[r] = KG_PHt[CM(state, r, 17)];
      }
      Cholesky_Lower_Solve_InPlace(KG_S, H_rows, MLS_HROWS_MAX, KG_rhs);
      for (size_t col = 0; col < H_rows; ++col) {
        K_out[CM(state, col, 17)] = KG_rhs[col];
      }
    }
    return true;
  }

  // ------------------------------------------------------------
  // 4. Reject update if S is not SPD
  // ------------------------------------------------------------
  Serial.println("Kalman Gain Error: Cholesky failed. GNSS update rejected.");
  return false;
}



bool MLS_KF_predict(
  const double tor_i,
  double C_interpt[9],
  double v_interpt[3],
  double r_interpt[3],
  double meas_f_ib_b[3],
  TC_KFConfig &TC_KF_Config,
  double meas_omega_ib_b[3],
  MLSConfig &MLS_Config,
  double P_matrix[289]) {

  // Serial.println("Kalman filter prediction epoch");
  const double tor  = tor_i;
  const double tor2 = tor * tor;
  const double tor3 = tor2 * tor;

  const double half_tor2 = 0.5 * tor2;
  const double one_third_tor3 = (1.0 / 3.0) * tor3;
  const double one_sixth_tor3 = (1.0 / 6.0) * tor3;
  ///---------------------------------Initializing Earth Parameters------------------/////
  const double omega_ie = 7.292115E-5;        //  Earth Rotation Rate in rad/s
  const double R_0 = 6378137.0;               //  WGS84 Equatorial Radius in Meters
  const double e = 0.0818191908425;           //  Eccentricity
  double est_L_b_old;
  double est_lambda_b_old;
  double h_b;
  double est_v_eb_n_old[3] = {0.0, 0.0, 0.0};
  double est_C_b_n_old[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

 

  if (!ECEF_to_NED(C_interpt, v_interpt, r_interpt, est_L_b_old, est_lambda_b_old, h_b, est_v_eb_n_old, est_C_b_n_old)) {
    Serial.println("MLS_KF_predict Error: ECEF to NED Failed");
    return false;
  }

  ///-----------------------------------Skew Symmetric Matrix of Earth Rate-----------------------------////
  double omega_vec[3] = {0.0, 0.0, omega_ie};
  double Omega_ie[9];
  Skew_Symmetric(omega_vec, Omega_ie);

  double Cv[3];
  double S[9];

  double OC[9];     // Omega_ie * C_interpt
  double SO[9];     // S * Omega_ie
  double OS[9];     // Omega_ie * S
  double SC[9];     // S * C_interpt
  double SOC[9];    // S * Omega_ie * C_interpt
  double OSC[9];    // Omega_ie * S * C_interpt

  Mat3VecCM(C_interpt, meas_f_ib_b, Cv);
  Skew_Symmetric(Cv, S);

  Mat3MulCM(Omega_ie, C_interpt, OC);
  Mat3MulCM(S, Omega_ie, SO);
  Mat3MulCM(Omega_ie, S, OS);
  Mat3MulCM(S, C_interpt, SC);
  Mat3MulCM(SO, C_interpt, SOC);
  Mat3MulCM(OS, C_interpt, OSC);
  ///----------------------------------Earth Radius update Calculation From 2.137--------------------------------/////
  const double sinL = sin(est_L_b_old);
  const double cosL = cos(est_L_b_old);
  const double e2 = e * e;
  const double one_minus_e2 = 1.0 - e2;

  const double denom = sqrt(1.0 - e2 * sinL * sinL);

  double geocentric_radius =
      (R_0 / denom) *
      sqrt(cosL * cosL +
          one_minus_e2 * one_minus_e2 * sinL * sinL);
  ////------------------------------------------SYSTEM PROPAGATION PHASE-------------------------------------------------------------------//////////
  ///-------------------------------------1. Determine ECEF State Transition Matrix (17x17) Using (I.1) (Third-Order Approx)-------------------------//////////////
  
  // MATALB = Phi_matrix = eye(17);
  static double Phi_matrix[289] = {0.0};
  memset(Phi_matrix, 0, 289 * sizeof(double));


  for (int i = 0; i < 17; i++) {
    Phi_matrix[i + 17 * i] = 1.0;
  }

  //  MATLAB = Phi_matrix(1:3,1:3) = Phi_matrix(1:3,1:3) - Omega_ie * tor_i;
  for (int col = 0; col < 3; ++col) {
    for (int row = 0; row < 3; ++row) {
      Phi_matrix[row + 17 * col] -= Omega_ie[row + 3 * col] * tor_i;
    }
  }

  // Phi(1:3,13:15) = C*tor - 0.5*Omega*C*tor^2
  for (int col = 0; col < 3; ++col) {
    for (int row = 0; row < 3; ++row) {
      Phi_matrix[row + 17 * (col + 12)] =
          C_interpt[row + 3 * col] * tor -
          OC[row + 3 * col] * half_tor2;
    }
  }

  // Phi(4:6,1:3)
  for (int col = 0; col < 3; ++col) {
    for (int row = 0; row < 3; ++row) {
      Phi_matrix[(row + 3) + 17 * col] =
          -S[row + 3 * col] * tor +
          0.5 * SO[row + 3 * col] * tor2 +
          OS[row + 3 * col] * tor2;
    }
  }

  // Phi(4:6,4:6) = I - 2*Omega_ie*tor
  for (int col = 0; col < 3; ++col) {
    for (int row = 0; row < 3; ++row) {
      Phi_matrix[(row + 3) + 17 * (col + 3)] -=
          2.0 * Omega_ie[row + 3 * col] * tor;
    }
  }
  
  // MATLAB = Phi_matrix(4:6,7:9) = -tor_s * 2 * Gravity_ECEF(est_r_eb_e_old) / geocentric_radius * est_r_eb_e_old' / sqrt (est_r_eb_e_old' * est_r_eb_e_old);
  double gravity[3];
  Gravity_ECEF(r_interpt, gravity);
  double r_norm = sqrt(r_interpt[0] * r_interpt[0] +
                       r_interpt[1] * r_interpt[1] +
                       r_interpt[2] * r_interpt[2]
                      );

  
  //=================SAFETY CHECK===========================//
  if (!(r_norm > 0.0) || !isfinite(r_norm)) {
    Serial.println("MLS_KF_predict Error: invalid ECEF position norm");
    return false;
  }

  if (!(geocentric_radius > 0.0) || !isfinite(geocentric_radius)) {
    Serial.println("MLS_KF_predict Error: invalid geocentric radius");
    return false;
  }
  //=================SAFETY CHECK===========================//

  double factor = (-tor_i * 2.0) / geocentric_radius;
  for (int col = 0; col < 3; ++col) {
    for (int row = 0; row < 3; ++row) {
      Phi_matrix[(row + 3) + 17 * (col + 6)] = factor * gravity[row] * r_interpt[col] / r_norm;
    }
  }

  // Phi(4:6,10:12) = C*tor - Omega*C*tor^2
  for (int col = 0; col < 3; ++col) {
    for (int row = 0; row < 3; ++row) {
      Phi_matrix[(row + 3) + 17 * (col + 9)] =
          C_interpt[row + 3 * col] * tor -
          OC[row + 3 * col] * tor2;
    }
  }

  // Phi(4:6,13:15)
  for (int col = 0; col < 3; ++col) {
    for (int row = 0; row < 3; ++row) {
      Phi_matrix[(row + 3) + 17 * (col + 12)] =
          -SC[row + 3 * col] * half_tor2 +
          SOC[row + 3 * col] * one_sixth_tor3 +
          OSC[row + 3 * col] * one_third_tor3;
    }
  }

  // Phi(7:9,1:3)
  for (int col = 0; col < 3; ++col) {
    for (int row = 0; row < 3; ++row) {
      Phi_matrix[(row + 6) + 17 * col] =
          -S[row + 3 * col] * half_tor2 +
          SO[row + 3 * col] * one_sixth_tor3 +
          OS[row + 3 * col] * one_third_tor3;
    }
  }


  // Phi(7:9,4:6)
  for (int col = 0; col < 3; ++col) {
    for (int row = 0; row < 3; ++row) {
      const double I = (row == col) ? 1.0 : 0.0;

      Phi_matrix[(row + 6) + 17 * (col + 3)] =
          I * tor - Omega_ie[row + 3 * col] * tor2;
    }
  }

  // Phi(7:9,10:12)
  for (int col = 0; col < 3; ++col) {
    for (int row = 0; row < 3; ++row) {
      Phi_matrix[(row + 6) + 17 * (col + 9)] =
          0.5 * C_interpt[row + 3 * col] * tor2 -
          OC[row + 3 * col] * one_third_tor3;
    }
  }


  // Phi(7:9,13:15)
  for (int col = 0; col < 3; ++col) {
    for (int row = 0; row < 3; ++row) {
      Phi_matrix[(row + 6) + 17 * (col + 12)] =
          -SC[row + 3 * col] * one_sixth_tor3;
    }
  }

  // Phi(16,17) = tor
  Phi_matrix[15 + 17 * 16] = tor;

  ///-------------------------------------2. Determine Approximate System Noise Covariance Matrix (17x17) Using (14.82)-------------------------//////////////
  double Q_prime_matrix[17] = {0.0};

  const double q_gyro       = TC_KF_Config.gyro_noise_PSD  * tor;
  const double q_accel      = TC_KF_Config.accel_noise_PSD * tor;
  const double q_accel_bias = TC_KF_Config.accel_bias_PSD * tor;
  const double q_gyro_bias  = TC_KF_Config.gyro_bias_PSD  * tor;
  const double q_clk_phase  = TC_KF_Config.clock_phase_PSD * tor;
  const double q_clk_freq   = TC_KF_Config.clock_freq_PSD  * tor;

  for (int i = 0; i < 3; i++) {
    Q_prime_matrix[i] = q_gyro;
  }

  for (int i = 3; i < 6; i++) {
    Q_prime_matrix[i] = q_accel;
  }

  for (int i = 9; i < 12; i++) {
    Q_prime_matrix[i] = q_accel_bias;
  }

  for (int i = 12; i < 15; i++) {
    Q_prime_matrix[i] = q_gyro_bias;
  }

  Q_prime_matrix[15] = q_clk_phase;
  Q_prime_matrix[16] = q_clk_freq;
 
  ///-----------------------------------3. Propagate State Estimation Error Covariance Matrix (17x17)------------------------------------///////////
  static double P_matrix_propagated[289] = {0.0};
  static double A[289] = {0.0};
  static double B[289] = {0.0};
  


   // MATLAB = P_matrix_propagated = Phi_matrix * (P_matrix_old + 0.5 * Q_prime_matrix) * Phi_matrix' + 0.5 * Q_prime_matrix;
  memcpy(A, P_matrix, 289 * sizeof(double));

  for (int i = 0; i < 17; i++) {
    A[i + 17 * i] += 0.5 * Q_prime_matrix[i];
  }
  for (int col = 0; col < 17; col++) {
    for (int row = 0; row < 17; row++) {
      double sum = 0.0;
      for (int k = 0; k < 17; k++) {
        sum += Phi_matrix[row + 17 * k] * A[k + 17 * col];
      }
      B[row + 17 * col] = sum;
    }
  }
  for (int col = 0; col < 17; col++) {
    for (int row = 0; row < 17; row++) {
      double sum = 0.0;
      for (int k = 0; k < 17; k++) {
        sum += B[row + 17 * k] * Phi_matrix[col + 17 * k];
      }
      P_matrix_propagated[row + 17 * col] = sum ;
    }
  }
  for (int i = 0; i < 17; i++) {
    P_matrix_propagated[i + 17 * i] += 0.5 * Q_prime_matrix[i];
  }
  // Write propagated covariance back to P_matrix for IMU-rate prediction.
  for (int col = 0; col < 17; col++) {
    for (int row = 0; row < 17; row++) {
      P_matrix[row + 17 * col] = P_matrix_propagated[row + 17 * col];
    }
  }

  // Force covariance symmetry after propagation.
  for (int col = 0; col < 17; col++) {
    for (int row = col + 1; row < 17; row++) {
      const double a = 0.5 * (P_matrix[row + 17 * col] + P_matrix[col + 17 * row]);
      P_matrix[row + 17 * col] = a;
      P_matrix[col + 17 * row] = a;
    }
  }

  (void)v_interpt;
  (void)meas_omega_ib_b;
  (void)MLS_Config;

  return true;
}


bool MLS_KF_update(
  GNSSRow gnssRows[],
  size_t gnssCount,
  double C_interpt[9],
  double v_interpt[3],
  double r_interpt[3],
  double meas_f_ib_b[3],
  TC_KFConfig &TC_KF_Config,
  const double L_ba_b[3],
  double meas_omega_ib_b[3],
  MLSConfig &MLS_Config,
  double P_matrix[289],
  double *R_matrix,
  double *R_matrix_NextLoop,
  size_t R_dim,
  GNSSConfig &gnssConfig,
  double est_C_b_e[9],
  double est_v_eb_e[3],
  double est_r_eb_e[3],
  double est_IMU_Bias[6]) {
    

  ///---------------------------------Initializing Earth Parameters------------------/////
  const double omega_ie = 7.292115E-5;        // Earth rotation rate in rad/s
  const double c = 299792458.0;               // Speed of light in vacuum

  if (gnssRows == nullptr || gnssCount == 0 || gnssCount > MAX_GNSS_ROWS_PER_EPOCH) {
    Serial.println("MLS_KF_update Error: invalid GNSS measurement count");
    return false;
  }

  // Clear update workspaces
  memset(KG_PHt, 0, sizeof(KG_PHt));
  memset(KG_S, 0, sizeof(KG_S));
  memset(KG_rhs, 0, sizeof(KG_rhs));
  memset(KG_piv, 0, sizeof(KG_piv));
  memset(MLS_K_matrix, 0, sizeof(MLS_K_matrix));
  memset(MLS_u_as_e_T, 0, sizeof(MLS_u_as_e_T));
  memset(MLS_pred_meas, 0, sizeof(MLS_pred_meas));
  memset(MLS_H_matrix, 0, sizeof(MLS_H_matrix));
  memset(MLS_delta_z, 0, sizeof(MLS_delta_z));

  // In the split architecture, P_matrix entering MLS_KF_update() is already the
  // predicted/prior covariance. Keep a copy as P_matrix_propagated because the
  // original MLS residual logic uses this variable name.
  double *P_matrix_propagated = MLS_P_matrix_propagated;
  for (int i = 0; i < 289; i++) {
    P_matrix_propagated[i] = P_matrix[i];
  }

  // Error state prior is zero in an error-state Kalman filter.
  double x_est_propagated[17] = {0.0};

  ////---------------------------------------------------------------------------------MEASUREMENT UPDATE PHASE-----------------------------------------------------------/////////////
  
  ///--------------------------------------------Receiver Position & Velocity In ECEF Without Lever Arm Compensation--------------------------------//////////////
  double *u_as_e_T = MLS_u_as_e_T;
  double *pred_meas = MLS_pred_meas;
  memset(u_as_e_T, 0, MLS_U_AS_E_T_SIZE * sizeof(double));
  memset(pred_meas, 0, MLS_PRED_MEAS_SIZE * sizeof(double));

  //  MATLAB = est_r_ea_e_old = est_r_eb_e_old + est_C_b_e_old * L_ba_b;
  double est_r_ea_e_old[3] = {0.0};
  est_r_ea_e_old[0] = r_interpt[0] + (C_interpt[0] * L_ba_b[0] + C_interpt[3] * L_ba_b[1] + C_interpt[6] * L_ba_b[2]);
  est_r_ea_e_old[1] = r_interpt[1] + (C_interpt[1] * L_ba_b[0] + C_interpt[4] * L_ba_b[1] + C_interpt[7] * L_ba_b[2]);
  est_r_ea_e_old[2] = r_interpt[2] + (C_interpt[2] * L_ba_b[0] + C_interpt[5] * L_ba_b[1] + C_interpt[8] * L_ba_b[2]);

  //  MATLAB = est_v_ea_e_old = est_v_eb_e_old + est_C_b_e_old * (Skew_symmetric(meas_omega_ib_b)*L_ba_b);
  double est_v_ea_e_old[3] = {0.0};

  double omega_ie_b[3];

  omega_ie_b[0] = C_interpt[2] * omega_ie;
  omega_ie_b[1] = C_interpt[5] * omega_ie;
  omega_ie_b[2] = C_interpt[8] * omega_ie;

  double omega_eb_b[3];

  omega_eb_b[0] = meas_omega_ib_b[0] - omega_ie_b[0];
  omega_eb_b[1] = meas_omega_ib_b[1] - omega_ie_b[1];
  omega_eb_b[2] = meas_omega_ib_b[2] - omega_ie_b[2];
  
  double temp_skew_symmetric[9] = {0.0};
  Skew_Symmetric(omega_eb_b, temp_skew_symmetric);
  est_v_ea_e_old[0]= v_interpt[0] + (C_interpt[0] * (temp_skew_symmetric[0] * L_ba_b[0] + temp_skew_symmetric[3] * L_ba_b[1] + temp_skew_symmetric[6] * L_ba_b[2]) +
                     C_interpt[3] * (temp_skew_symmetric[1] * L_ba_b[0] + temp_skew_symmetric[4] * L_ba_b[1] + temp_skew_symmetric[7] * L_ba_b[2]) +
                     C_interpt[6] * (temp_skew_symmetric[2] * L_ba_b[0] + temp_skew_symmetric[5] * L_ba_b[1] + temp_skew_symmetric[8] * L_ba_b[2]));
  
  est_v_ea_e_old[1]= v_interpt[1] + (C_interpt[1] * (temp_skew_symmetric[0] * L_ba_b[0] + temp_skew_symmetric[3] * L_ba_b[1] + temp_skew_symmetric[6] * L_ba_b[2]) +
                     C_interpt[4] * (temp_skew_symmetric[1] * L_ba_b[0] + temp_skew_symmetric[4] * L_ba_b[1] + temp_skew_symmetric[7] * L_ba_b[2]) +
                     C_interpt[7] * (temp_skew_symmetric[2] * L_ba_b[0] + temp_skew_symmetric[5] * L_ba_b[1] + temp_skew_symmetric[8] * L_ba_b[2]));
  
  est_v_ea_e_old[2]= v_interpt[2] + (C_interpt[2] * (temp_skew_symmetric[0] * L_ba_b[0] + temp_skew_symmetric[3] * L_ba_b[1] + temp_skew_symmetric[6] * L_ba_b[2]) +
                     C_interpt[5] * (temp_skew_symmetric[1] * L_ba_b[0] + temp_skew_symmetric[4] * L_ba_b[1] + temp_skew_symmetric[7] * L_ba_b[2]) +
                     C_interpt[8] * (temp_skew_symmetric[2] * L_ba_b[0] + temp_skew_symmetric[5] * L_ba_b[1] + temp_skew_symmetric[8] * L_ba_b[2]));

  
  ///------------------------------------------------------Loop Measurements--------------------------------------------------------/////////////////
  for (size_t i = 0; i < gnssCount; i++) {

    ///------------------------Predict Pseudo-Range Using (9.165)----------------------------//////////
    // MATLAB = delta_r = GNSS_measurements(j,3:5)' - est_r_ea_e_old;

    double delta_r[3] = {0.0};

    delta_r[0] = gnssRows[i].sat_x - est_r_ea_e_old[0];
    delta_r[1] = gnssRows[i].sat_y - est_r_ea_e_old[1];
    delta_r[2] = gnssRows[i].sat_z - est_r_ea_e_old[2];

    // MATLAB = range = sqrt(delta_r' * delta_r);

    double range = sqrt(
        (delta_r[0] * delta_r[0]) +
        (delta_r[1] * delta_r[1]) +
        (delta_r[2] * delta_r[2])
    );

    if (range <= 0.0) {
        pred_meas[i] = 0.0;
        pred_meas[i + gnssCount] = 0.0;

        u_as_e_T[i + gnssCount * 0] = 0.0;
        u_as_e_T[i + gnssCount * 1] = 0.0;
        u_as_e_T[i + gnssCount * 2] = 0.0;

        continue;
    }

    ////-------------------------Predict Line of Sight--------------------------////////////

    // MATLAB = u_as_e_T(j,1:3) = delta_r' / range;

    u_as_e_T[i + gnssCount * 0] = delta_r[0] / range;  // MATLAB u_as_e_T(j+1,1)
    u_as_e_T[i + gnssCount * 1] = delta_r[1] / range;  // MATLAB u_as_e_T(j+1,2)
    u_as_e_T[i + gnssCount * 2] = delta_r[2] / range;  // MATLAB u_as_e_T(j+1,3)

    ///------------------------Sagnac Correction for Pseudo-Range----------------------------//////////
    // MATLAB:
    // sag_range = (omega_ie / c) * ...
    //     (sat_x * rec_y - sat_y * rec_x);

    double sag_range = (omega_ie / c) *
        (
            gnssRows[i].sat_x * est_r_ea_e_old[1] -
            gnssRows[i].sat_y * est_r_ea_e_old[0]
        );

    // MATLAB:
    // pred_meas(j,1) = range + sag_range + est_clock_old(1,1);

    pred_meas[i] = range + sag_range + gnssConfig.est_clock[0];

    ///------------------------Sagnac Correction for Pseudo-Range Rate----------------------------//////////
    // MATLAB:
    // sag_cor = (omega_ie/c) * ...
    //     (sat_vx * rec_y + sat_x * rec_vy ...
    //    - sat_vy * rec_x - sat_y * rec_vx);

    double sag_cor = (omega_ie / c) *
        (
            gnssRows[i].sat_vx * est_r_ea_e_old[1] +
            gnssRows[i].sat_x  * est_v_ea_e_old[1] -
            gnssRows[i].sat_vy * est_r_ea_e_old[0] -
            gnssRows[i].sat_y  * est_v_ea_e_old[0]
        );

    ///-----------------------Predict Pseudo-Range Rate Using (8.45) With Sagnac Effect Compensation---------------------------------///////////////
    // MATLAB:
    // range_rate = u_as_e_T(j,1:3) * ...
    //     (GNSS_measurements(j,6:8)' - est_v_ea_e_old) + sag_cor;

    double range_rate =
        u_as_e_T[i + gnssCount * 0] * (gnssRows[i].sat_vx - est_v_ea_e_old[0]) +
        u_as_e_T[i + gnssCount * 1] * (gnssRows[i].sat_vy - est_v_ea_e_old[1]) +
        u_as_e_T[i + gnssCount * 2] * (gnssRows[i].sat_vz - est_v_ea_e_old[2]) +
        sag_cor;

    // MATLAB:
    // pred_meas(j,2) = range_rate + est_clock_old(1,2);

    pred_meas[i + gnssCount] = range_rate + gnssConfig.est_clock[1];
  }


  ///----------------------------------5. Set-UP Measurement Matrix Using (14.126)----------------------------------------//////////////////
  double *H_matrix = MLS_H_matrix;
  size_t H_rows = 2 * gnssCount;            // No of Rows in Matrix
  if (TC_KF_Config.StationaryFlag == 1) {
    H_rows += 1;                          // Increase Row with 1
  }

  if (R_dim != H_rows) {
    Serial.println("MLS_KF_update Error: R_dim must equal H_rows");
    return false;
  }

  if (H_rows > MLS_HROWS_MAX) {
    Serial.println("MLS_KF_update Error: H_rows exceeds MLS_HROWS_MAX");
    return false;
  }

  // MATLAB equivalent: H_matrix = zeros(H_rows,17)
  memset(H_matrix, 0, H_rows * 17 * sizeof(double));
 
  // MATLAB = H_matrix(1:no_meas,7:9) = u_as_e_T(1:no_meas,1:3);
  for (size_t i = 0 ; i < gnssCount ; i++) {
    H_matrix[i + H_rows * 6] = u_as_e_T[i + gnssCount * 0];
    H_matrix[i + H_rows * 7] = u_as_e_T[i + gnssCount * 1];
    H_matrix[i + H_rows * 8] = u_as_e_T[i + gnssCount * 2];

    // MATLAB = H_matrix(1:no_meas,16) = ones(no_meas,1);
    H_matrix[i + H_rows * 15] = 1; 
  }

  // MATLAB = H_matrix((no_meas + 1):(2 * no_meas), 4:6) = u_as_e_T(1:no_meas, 1:3);
  for (size_t row = 0; row < gnssCount; ++row) {
    H_matrix[(gnssCount + row) + H_rows * 3] = u_as_e_T[row + gnssCount * 0];
    H_matrix[(gnssCount + row) + H_rows * 4] = u_as_e_T[row + gnssCount * 1];
    H_matrix[(gnssCount + row) + H_rows * 5] = u_as_e_T[row + gnssCount * 2];
  }

  // MATLAB = H_matrix((no_meas + 1):(2 * no_meas),17) = ones(no_meas,1);
  for (size_t row = 0; row < gnssCount; ++row) {
    H_matrix[(gnssCount + row) + H_rows * 16] = 1;
  }


  // MATLAB = if TC_KF_config.StationarityFlag == 1
          //     H_matrix((2 * no_meas)+1,13:15) = -1;
          // end
  if (TC_KF_Config.StationaryFlag == 1) {
    size_t row = 2 * gnssCount;
   
    H_matrix[row + H_rows * 14] = -1.0;
  } 
  
  ///-------------------------------------------------------8. Formulate Measurement Innovations Using (14.119)------------------------------------------////////////
  
  // MATLAB = delta_z=zeros(no_meas*2,1);
           // if TC_KF_config.StationarityFlag == 1
           //     delta_z=zeros(no_meas*2+1,1);
           // end
  double *delta_z = MLS_delta_z;
  memset(delta_z, 0, delta_z_DIM * sizeof(double));

  size_t delta_z_len = 2 * gnssCount;
  if (TC_KF_Config.StationaryFlag == 1) {
    delta_z_len += 1;
  }

  //  MATLAB = delta_z(1:no_meas,1) = GNSS_measurements(1:no_meas,1) - pred_meas(1:no_meas,1);
  //  MATLAB = delta_z((no_meas + 1):(2 * no_meas),1) = GNSS_measurements(1:no_meas,2) - pred_meas(1:no_meas,2);
  for (size_t i = 0; i < gnssCount; i++) {
    
    // MATLAB = delta_z(1:no_meas,1) = GNSS_measurements(1:no_meas,1) - pred_meas(1:no_meas,1);
    delta_z[i] = gnssRows[i].pr_corr - pred_meas[i];
    
    // MATLAB = delta_z((no_meas + 1):(2 * no_meas),1) = GNSS_measurements(1:no_meas,2) - pred_meas(1:no_meas,2);
    delta_z[gnssCount + i] = gnssRows[i].prrate_corr - pred_meas[gnssCount + i];
  }

  // MATLAB =  if TC_KF_config.StationarityFlag == 1
          //   delta_z(2*no_meas + 1,1)=-meas_omega_ib_b(3);
          // end
  if (TC_KF_Config.StationaryFlag == 1) {
    delta_z[2 * gnssCount] = -meas_omega_ib_b[2];
  }
  
  Serial.println(F("========== GNSS Innovation delta_z Debug =========="));

  Serial.print("GNSS_sec = ");
  Serial.print(packet.gnss_obssec);
  Serial.print(F("gnssCount = "));
  Serial.println(gnssCount);

  for (size_t i = 0; i < gnssCount; i++) {
  Serial.print(F("SV "));
  Serial.print(i);

  Serial.print(F(", PR_meas = "));
  Serial.print(gnssRows[i].pr_corr, 6);

  Serial.print(F(", PR_pred = "));
  Serial.print(pred_meas[i], 6);

  Serial.print(F(", PR_delta_z = "));
  Serial.print(delta_z[i], 6);

  Serial.print(F(", PRR_meas = "));
  Serial.print(gnssRows[i].prrate_corr, 6);

  Serial.print(F(", PRR_pred = "));
  Serial.print(pred_meas[gnssCount + i], 6);

  Serial.print(F(", PRR_delta_z = "));
  Serial.println(delta_z[gnssCount + i], 6);
  }
  ///-------------------------------------------------------7. Calculate Kalman Gain Using (3.21)------------------------------------------////////////
  // MATLAB = K_matrix = P_matrix_propagated * H_matrix' / (H_matrix * P_matrix_propagated * H_matrix' + R_matrix);
  double *K_matrix = MLS_K_matrix;
  memset(K_matrix, 0, 17 * MLS_HROWS_MAX * sizeof(double));
  if (!Kalman_Gain_Solver_PSram(
        P_matrix_propagated,   // 17 x 17
        H_matrix,              // HROWS_MAX x 17
        R_matrix,              // R_dim x R_dim
        H_rows,                // active measurement rows
        R_dim,                 // leading dimension of R_matrix
        K_matrix               // output: 17 x H_rows
      )) {
    Serial.println("Kalman Gain Solver Failed");
    return false;
  }
  
  // PrintKalmanGainMatrix(K_matrix, STATE_DIM, R_dim);
  // Serial.println("Kalman filter update epoch");
  ///-------------------------------------------------------9. Update State Estimates Using (3.24)---------------------------------------///////////////

  static double x_est_new[17] = {0.0};

  for (size_t row = 0; row < 17; ++row) {
    double sum = 0.0;
    for (size_t col = 0; col < H_rows; ++col) {
      sum += K_matrix[CM(row, col, 17)] * delta_z[col];
    }
    x_est_new[row] = x_est_propagated[row] + sum;
  }

  ///-----------------------------------10. Update State Estimation Error Covariance Matrix Using Joseph Form(3.58)--------------------------------///////////////
  // MATLAB:
  // P_matrix = (I - K_matrix * H_matrix) * P_matrix_propagated * ...
  //            (I - K_matrix * H_matrix)' + K_matrix * R_matrix * K_matrix';
  //
  // Result is written directly into P_matrix.

  static double A_joseph[289] = {0.0};      // A = I - K*H
  static double AP_joseph[289] = {0.0};     // A * P_propagated
  static double APA_joseph[289] = {0.0};    // A * P_propagated * A'

  static double KR_joseph[17 * MLS_HROWS_MAX] = {0.0};  // K * R
  static double KRK_joseph[289] = {0.0};                // K * R * K'

  // ------------------------------------------------------------
  // 1. A_joseph = I - K_matrix * H_matrix
  // ------------------------------------------------------------
  for (size_t col = 0; col < 17; ++col) {
    for (size_t row = 0; row < 17; ++row) {
      double sum = 0.0;
      for (size_t k = 0; k < H_rows; ++k) {
        sum += K_matrix[CM(row, k, 17)] * H_matrix[CM(k, col, H_rows)];
      }
      A_joseph[CM(row, col, 17)] = ((row == col) ? 1.0 : 0.0) - sum;
    }
  }

  // ------------------------------------------------------------
  // 2. AP_joseph = A_joseph * P_matrix_propagated
  // ------------------------------------------------------------
  for (size_t col = 0; col < 17; ++col) {
    for (size_t row = 0; row < 17; ++row) {
      double sum = 0.0;
      for (size_t k = 0; k < 17; ++k) {
        sum += A_joseph[CM(row, k, 17)] *  P_matrix_propagated[CM(k, col, 17)];
      }
      AP_joseph[CM(row, col, 17)] = sum;
    }
  }

  // ------------------------------------------------------------
  // 3. APA_joseph = AP_joseph * A_joseph'
  // ------------------------------------------------------------
  for (size_t col = 0; col < 17; ++col) {
    for (size_t row = 0; row < 17; ++row) {

      double sum = 0.0;

      for (size_t k = 0; k < 17; ++k) {
        // A_joseph'(k,col) = A_joseph(col,k)
        sum += AP_joseph[CM(row, k, 17)] * A_joseph[CM(col, k, 17)];
      }
      APA_joseph[CM(row, col, 17)] = sum;
    }
  }

  // ------------------------------------------------------------
  // 4. KR_joseph = K_matrix * R_matrix
  // ------------------------------------------------------------
  for (size_t col = 0; col < H_rows; ++col) {
    for (size_t row = 0; row < 17; ++row) {
      double sum = 0.0;
      for (size_t k = 0; k < H_rows; ++k) {
        sum += K_matrix[CM(row, k, 17)] * R_matrix[CM(k, col, R_dim)];
      }
      KR_joseph[CM(row, col, 17)] = sum;
    }
  }

  // ------------------------------------------------------------
  // 5. KRK_joseph = KR_joseph * K_matrix'
  // ------------------------------------------------------------
  for (size_t col = 0; col < 17; ++col) {
    for (size_t row = 0; row < 17; ++row) {

      double sum = 0.0;

      for (size_t k = 0; k < H_rows; ++k) {
        // K_matrix'(k,col) = K_matrix(col,k)
        sum += KR_joseph[CM(row, k, 17)] *
              K_matrix[CM(col, k, 17)];
      }

      KRK_joseph[CM(row, col, 17)] = sum;
    }
  }

  // ------------------------------------------------------------
  // 6. Final covariance update directly into P_matrix
  // ------------------------------------------------------------
  for (size_t col = 0; col < 17; ++col) {
    for (size_t row = 0; row < 17; ++row) {
      P_matrix[CM(row, col, 17)] =
          APA_joseph[CM(row, col, 17)] +
          KRK_joseph[CM(row, col, 17)];
    }
  }

  ///-------------------------------------------------------6. Update Measurement Noise Covariance Matrix (R) Using M Gain Residual-------------------------------------------------------//////////////////////

  // MATLAB:
  //
  // Resi = (eye(length(delta_z)) - H_matrix*K_matrix) * delta_z;
  //
  // UnitWeightVar = ...
  //   (delta_z' / (R_matrix + H_matrix*P_matrix_propagated*H_matrix')) ...
  //   * delta_z / length(delta_z);
  //
  // CovarianceOfResi = R_matrix - H_matrix*P_matrix_new*H_matrix';
  //
  // StdResi = Resi / sqrt(UnitWeightVar) ./ sqrt(diag(CovarianceOfResi));
  //
  // MGain = eye(length(delta_z));
  // index_outlier = find(abs(StdResi) > MLS.K1);
  // index_buffer  = find(abs(StdResi) > MLS.K0 & abs(StdResi) <= MLS.K1);
  //
  // R_matrix_NextLoop = MGain * R_matrix;
  //
  // In this C++ version:
  // P_matrix already contains P_matrix_new from Joseph update.
  // R_matrix is updated in-place and becomes R_matrix_NextLoop.

  const size_t m = H_rows;
  const size_t H_ld = H_rows;

  // Your current R_matrix initialization uses:
  // R_Matrix[col * R_dim + row]
  // so R_ld = R_dim is correct.
  const size_t R_ld = R_dim;

  if (m == 0 || m > MLS_HROWS_MAX) {
    Serial.println("MLS residual update error: invalid m");
    return false;
  }

  if (R_dim != m) {
    Serial.println("MLS residual update error: R_dim != H_rows");
    return false;
  }

  static double HK_mls[MLS_HROWS_MAX * MLS_HROWS_MAX] = {0.0};  // H*K, size m x m
  static double HP_mls[MLS_HROWS_MAX * 17] = {0.0};             // H*P, size m x 17

  static double Resi[MLS_HROWS_MAX] = {0.0};
  static double StdResi[MLS_HROWS_MAX] = {0.0};
  static double MGainDiag[MLS_HROWS_MAX] = {0.0};

  //
  // ------------------------------------------------------------
  // 1. HK_mls = H_matrix * K_matrix
  // MATLAB operation order:
  // Resi = (eye(length(delta_z)) - H_matrix*K_matrix) * delta_z;
  // ------------------------------------------------------------
  // H_matrix : m x 17
  // K_matrix : 17 x m
  // HK_mls   : m x m
  //

  for (size_t col = 0; col < m; ++col) {
    for (size_t row = 0; row < m; ++row) {
      double sum = 0.0;
      for (size_t k = 0; k < 17; ++k) {
        sum += H_matrix[CM(row, k, H_ld)] * K_matrix[CM(k, col, 17)];
      }
      HK_mls[CM(row, col, m)] = sum;
    }
  }

  //
  // ------------------------------------------------------------
  // 2. Resi = (I - HK_mls) * delta_z
  // ------------------------------------------------------------
  // This follows MATLAB more closely than:
  // state_correction = K*delta_z;
  // Resi = delta_z - H*state_correction;
  //

  for (size_t row = 0; row < m; ++row) {
    double sum = 0.0;
    for (size_t col = 0; col < m; ++col) {
      const double I_minus_HK = ((row == col) ? 1.0 : 0.0) - HK_mls[CM(row, col, m)];
      sum += I_minus_HK * delta_z[col];
    }
    Resi[row] = sum;
  }

  //
  // ------------------------------------------------------------
  // 3. Build S = R_matrix + H_matrix * P_matrix_propagated * H_matrix'
  // MATLAB operation order:
  // S = R_matrix + (H_matrix * P_matrix_propagated) * H_matrix'
  // ------------------------------------------------------------
  // HP_mls = H_matrix * P_matrix_propagated
  //

  for (size_t col = 0; col < 17; ++col) {
    for (size_t row = 0; row < m; ++row) {
      double sum = 0.0;
      for (size_t k = 0; k < 17; ++k) {
        sum += H_matrix[CM(row, k, H_ld)] * P_matrix_propagated[CM(k, col, 17)];
      }
      HP_mls[CM(row, col, m)] = sum;
    }
  }

  //
  // KG_S = R_matrix + HP_mls * H_matrix'
  // KG_S size = m x m
  //

  for (size_t col = 0; col < m; ++col) {
    for (size_t row = 0; row < m; ++row) {
      double sum = R_matrix[CM(row, col, R_ld)];
      for (size_t k = 0; k < 17; ++k) {
        sum += HP_mls[CM(row, k, m)] * H_matrix[CM(col, k, H_ld)];
      }
      KG_S[CM(row, col, MLS_HROWS_MAX)] = sum;
    }
  }

  //
  // Do NOT call Symmetrize_S(m) here if your target is MATLAB matching.
  // MATLAB uses the matrix result as produced by the expression.
  //

  //
  // ------------------------------------------------------------
  // 4. UnitWeightVar = delta_z' * inv(S) * delta_z / m
  // MATLAB:
  // UnitWeightVar = (delta_z' / S) * delta_z / m;
  // ------------------------------------------------------------
  // Solve:
  // S * x = delta_z
  //
  // Then:
  // UnitWeightVar = delta_z' * x / m
  //

  for (size_t i = 0; i < m; ++i) {
    KG_rhs[i] = delta_z[i];
  }

  bool solve_ok = false;

  //
  // Try Cholesky first, as S should normally be SPD.
  // If it fails, use LU fallback.
  //

  if (Cholesky_Lower_InPlace(KG_S, m, MLS_HROWS_MAX)) {
    Cholesky_Lower_Solve_InPlace(KG_S, m, MLS_HROWS_MAX, KG_rhs);
    solve_ok = true;
  } else {

    Serial.println("MLS residual error: Cholesky failed for UnitWeightVar. GNSS update rejected.");
    return false;
    //
    // Rebuild S because Cholesky modifies KG_S
    //

    for (size_t col = 0; col < m; ++col) {
      for (size_t row = 0; row < m; ++row) {
        double sum = R_matrix[CM(row, col, R_ld)];
        for (size_t k = 0; k < 17; ++k) {
          sum += HP_mls[CM(row, k, m)] * H_matrix[CM(col, k, H_ld)];
        }
        KG_S[CM(row, col, MLS_HROWS_MAX)] = sum;
      }
    }

    for (size_t i = 0; i < m; ++i) {
      KG_rhs[i] = delta_z[i];
    }

    if (LU_Factor_InPlace(KG_S, m, MLS_HROWS_MAX, KG_piv)) {
      solve_ok = LU_Solve_InPlace(KG_S, m, MLS_HROWS_MAX, KG_piv, KG_rhs);
    }
  }

  if (!solve_ok) {
    Serial.println("MLS residual update error: failed to solve S for UnitWeightVar");
    return false;
  }

  double UnitWeightVar = 0.0;
  for (size_t i = 0; i < m; ++i) {
    UnitWeightVar += delta_z[i] * KG_rhs[i];
  }
  UnitWeightVar /= static_cast<double>(m);

  //
  // For MATLAB matching, do not replace invalid values with a floor.
  // But keep this check to avoid embedded crash.
  //

  if (!isfinite(UnitWeightVar)) {
    Serial.println("MLS residual update error: UnitWeightVar is not finite");
    return false;
  }

  //
  // ------------------------------------------------------------
  // 5. CovarianceOfResi = R_matrix - H_matrix * P_matrix_new * H_matrix'
  // MATLAB operation order:
  // CovarianceOfResi = R_matrix - (H_matrix * P_matrix) * H_matrix'
  // ------------------------------------------------------------
  // Here P_matrix already contains P_matrix_new from Joseph update.
  //
  // First compute:
  // HP_mls = H_matrix * P_matrix
  //

  for (size_t col = 0; col < 17; ++col) {
    for (size_t row = 0; row < m; ++row) {
      double sum = 0.0;
      for (size_t k = 0; k < 17; ++k) {
        sum += H_matrix[CM(row, k, H_ld)] *  P_matrix[CM(k, col, 17)];
      }
      HP_mls[CM(row, col, m)] = sum;
    }
  }

  //
  // Only diagonal of CovarianceOfResi is needed:
  // CovarianceOfResiDiag(row) = R(row,row) - HP(row,:) * H(row,:)'
  //

  const double sqrt_UnitWeightVar = sqrt(UnitWeightVar);
  for (size_t row = 0; row < m; ++row) {
    double hPh = 0.0;
    for (size_t k = 0; k < 17; ++k) {
      hPh += HP_mls[CM(row, k, m)] * H_matrix[CM(row, k, H_ld)];
    }
    const double CovarianceOfResiDiag = R_matrix[CM(row, row, R_ld)] - hPh;
    //
    // MATLAB-like behavior:
    // If CovarianceOfResiDiag is negative, sqrt() gives NaN.
    // If it is zero, division may give Inf.
    // Do not force 1e-12 here if you want MATLAB matching.
    //
    StdResi[row] = Resi[row] / sqrt_UnitWeightVar / sqrt(CovarianceOfResiDiag);
  }

  //
  // ------------------------------------------------------------
  // 6. MGain = eye(length(delta_z))
  // ------------------------------------------------------------
  // We store only diagonal of MGain.
  //

  for (size_t i = 0; i < m; ++i) {
    MGainDiag[i] = 1.0;
  }

  size_t OutlierCount = 0;

  //
  // ------------------------------------------------------------
  // 7. index_outlier and index_buffer
  // ------------------------------------------------------------
  // MATLAB:
  //
  // index_outlier = find(abs(StdResi) > MLS.K1);
  //
  // index_buffer = find(abs(StdResi) > MLS.K0 & ...
  //                     abs(StdResi) <= MLS.K1);
  //

  for (size_t i = 0; i < m; ++i) {

    const double absStd = fabs(StdResi[i]);

    if (absStd > MLS_Config.K1) {

      MGainDiag[i] = 1.0e4;
      OutlierCount++;

    } else if ((absStd > MLS_Config.K0) && (absStd <= MLS_Config.K1)) {

      const double K0 = MLS_Config.K0;
      const double K1 = MLS_Config.K1;

      //
      // MATLAB exact formula:
      //
      // MGain(i,i) = (K1-K0)^2/K0 * ...
      //              absStd / ((K1-absStd)*(K1-absStd));
      //
      // No safety replacement here.
      //

      MGainDiag[i] = ((K1 - K0) * (K1 - K0) / K0) * (absStd / ((K1 - absStd) * (K1 - absStd)));

      OutlierCount++;
    }
  }

  //
  // Optional:
  // MLS_Config.OutlierCount = OutlierCount;
  //
  //
  // ------------------------------------------------------------
  // 8. R_matrix_NextLoop = MGain * R_matrix
  // ------------------------------------------------------------
  // Since MGain is diagonal:
  // R_next(row,col) = MGainDiag(row) * R_matrix(row,col)
  //

  // Clear R_matrix_NextLoop
  for (size_t col = 0; col < m; ++col) {
    for (size_t row = 0; row < m; ++row) {
      R_matrix_NextLoop[CM(row, col, R_ld)] = 0.0;
    }
  }

  for (size_t col = 0; col < m; ++col) {
    for (size_t row = 0; row < m; ++row) {
      R_matrix_NextLoop[CM(row, col, R_ld)] = MGainDiag[row] * R_matrix[CM(row, col, R_ld)];
    }
  }
  ///-----------------------------------------------Correct Attitude, Velocity & Position Using (14.7-9)--------------------------------------------------///////////////////
  // MATLAB = est_C_b_e_new = (eye(3) - Skew_symmetric(x_est_new(1:3,1))) * est_C_b_e_old;

  // x_est_new[0] = attitude error x
  // x_est_new[1] = attitude error y
  // x_est_new[2] = attitude error z

  est_C_b_e[0] = C_interpt[0] + x_est_new[2] * C_interpt[1] - x_est_new[1] * C_interpt[2];
  est_C_b_e[1] = -x_est_new[2] * C_interpt[0] + C_interpt[1] + x_est_new[0] * C_interpt[2];
  est_C_b_e[2] = x_est_new[1] * C_interpt[0] - x_est_new[0] * C_interpt[1] + C_interpt[2];

  est_C_b_e[3] = C_interpt[3] + x_est_new[2] * C_interpt[4] - x_est_new[1] * C_interpt[5];
  est_C_b_e[4] = -x_est_new[2] * C_interpt[3] + C_interpt[4] + x_est_new[0] * C_interpt[5];
  est_C_b_e[5] = x_est_new[1] * C_interpt[3] - x_est_new[0] * C_interpt[4] + C_interpt[5];

  est_C_b_e[6] = C_interpt[6] + x_est_new[2] * C_interpt[7] - x_est_new[1] * C_interpt[8];
  est_C_b_e[7] = -x_est_new[2] * C_interpt[6] + C_interpt[7] + x_est_new[0] * C_interpt[8];
  est_C_b_e[8] = x_est_new[1] * C_interpt[6] - x_est_new[0] * C_interpt[7] + C_interpt[8];

  // MATLAB = est_v_eb_e_new = est_v_eb_e_old - x_est_new(4:6,1);

  // Your C++:
  // est_v_eb_e_old = v_interpt
  // est_v_eb_e_new = est_v_eb_e

  est_v_eb_e[0] = v_interpt[0] - x_est_new[3];
  est_v_eb_e[1] = v_interpt[1] - x_est_new[4];
  est_v_eb_e[2] = v_interpt[2] - x_est_new[5];


  // MATLAB = est_r_eb_e_new = est_r_eb_e_old - x_est_new(7:9,1);

  // Your C++:
  // est_r_eb_e_old = r_interpt
  // est_r_eb_e_new = est_r_eb_e

  est_r_eb_e[0] = r_interpt[0] - x_est_new[6];
  est_r_eb_e[1] = r_interpt[1] - x_est_new[7];
  est_r_eb_e[2] = r_interpt[2] - x_est_new[8];

  ///--------------------------------------------------Update IMU Bias & GNSS Receiver Clock Estimates----------------------------------------------------///////////////////
  // MATLAB = est_IMU_bias_new = est_IMU_bias_old + x_est_new(10:15,1);
  est_IMU_Bias[0] = est_IMU_Bias[0]+ x_est_new[9];
  est_IMU_Bias[1] = est_IMU_Bias[1]+ x_est_new[10];
  est_IMU_Bias[2] = est_IMU_Bias[2]+ x_est_new[11];
  est_IMU_Bias[3] = est_IMU_Bias[3]+ x_est_new[12];
  est_IMU_Bias[4] = est_IMU_Bias[4]+ x_est_new[13];
  est_IMU_Bias[5] = est_IMU_Bias[5]+ x_est_new[14];

  // MATLAB = est_clock_new = est_clock_old + x_est_new(16:17,1)';
  gnssConfig.est_clock[0] = gnssConfig.est_clock[0] + x_est_new[15];
  gnssConfig.est_clock[1] = gnssConfig.est_clock[1] + x_est_new[16];
  return true;
}
