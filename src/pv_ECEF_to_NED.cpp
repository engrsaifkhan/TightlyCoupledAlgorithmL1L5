#include "pv_ECEF_to_NED.h"
#include <cmath>

namespace {
  inline double signum_double(const double x)
  {
    return (x > 0.0) ? 1.0 : ((x < 0.0) ? -1.0 : 0.0);
  }

  inline bool isFinite3(const double a, const double b, const double c)
  {
    return std::isfinite(a) && std::isfinite(b) && std::isfinite(c);
  }
}

bool pv_ECEF_to_NED(const double r_eb_e[3], const double v_eb_e[3], double &old_est_L_b, double &old_est_lambda_b, double &old_est_h_b, double old_est_v_eb_n[3])
{
  if (!r_eb_e || !v_eb_e || !old_est_v_eb_n) {
    return false;
  }

  const double x = r_eb_e[0];
  const double y = r_eb_e[1];
  const double z = r_eb_e[2];

  if (!isFinite3(x, y, z) || !isFinite3(v_eb_e[0], v_eb_e[1], v_eb_e[2])) {
    return false;
  }

  // Original MATLAB constants
  const double R_0 = 6378137.0;          // WGS84 equatorial radius (m)
  const double e   = 0.0818191908425;    // WGS84 eccentricity

  const double e2 = e * e;
  const double sqrt_1_minus_e2 = std::sqrt(1.0 - e2);

  // From (2.113)
  old_est_lambda_b = std::atan2(y, x);

  // From (C.29) and (C.30)
  const double abs_z = std::fabs(z);
  const double beta  = std::sqrt(x * x + y * y);

  // Pole-safe handling
  if (beta == 0.0) {
    return false;
  }

  const double k1 = sqrt_1_minus_e2 * abs_z;
  const double k2 = e2 * R_0;
  const double E  = (k1 - k2) / beta;
  const double F  = (k1 + k2) / beta;

  // From (C.31)
  const double P = (4.0 / 3.0) * (E * F + 1.0);

  // From (C.32)
  const double Q = 2.0 * (E * E - F * F);

  // From (C.33)
  double D = P * P * P + Q * Q;
  if (D < 0.0) {
    D = 0.0;
  }

  // From (C.34)
  const double sqrt_D = std::sqrt(D);
  const double V = std::cbrt(sqrt_D - Q) - std::cbrt(sqrt_D + Q);

  // From (C.35)
  double G_term = E * E + V;
  if (G_term < 0.0) {
    G_term = 0.0;
  }
  const double G = 0.5 * (std::sqrt(G_term) + E);

  // From (C.36)
  const double denom = 2.0 * G - E;
  if (denom == 0.0) {
    return false;
  }

  double T_term = G * G + (F - V * G) / denom;
  if (T_term < 0.0) {
    T_term = 0.0;
  }
  const double T = std::sqrt(T_term) - G;

  if (T == 0.0) {
    return false;
  }

  // From (C.37)
  const double sign_z = signum_double(z);
  old_est_L_b = sign_z * std::atan((1.0 - T * T) / (2.0 * T * sqrt_1_minus_e2));

  // From (C.38)
  const double cos_lat  = std::cos(old_est_L_b);
  const double sin_lat  = std::sin(old_est_L_b);
  const double cos_long = std::cos(old_est_lambda_b);
  const double sin_long = std::sin(old_est_lambda_b);

  old_est_h_b = (beta - R_0 * T) * cos_lat + (z - sign_z * R_0 * sqrt_1_minus_e2) * sin_lat;

  // Direct ECEF to NED velocity transform
  old_est_v_eb_n[0] = (-sin_lat * cos_long) * v_eb_e[0] + (-sin_lat * sin_long) * v_eb_e[1] + ( cos_lat) * v_eb_e[2];

  old_est_v_eb_n[1] = (-sin_long) * v_eb_e[0] + ( cos_long) * v_eb_e[1];

  old_est_v_eb_n[2] = (-cos_lat * cos_long) * v_eb_e[0] + (-cos_lat * sin_long) * v_eb_e[1] + (-sin_lat) * v_eb_e[2];

  return isFinite3(old_est_L_b, old_est_lambda_b, old_est_h_b) && isFinite3(old_est_v_eb_n[0], old_est_v_eb_n[1], old_est_v_eb_n[2]);
}