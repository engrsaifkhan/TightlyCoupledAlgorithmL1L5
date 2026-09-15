#include "NED_to_ECEF.h"

namespace {
  constexpr double WGS84_A  = 6378137.0;
  constexpr double WGS84_E2 = 0.006694379990121435;
  constexpr double ONE_MINUS_E2 = 1.0 - WGS84_E2;
}

bool NED_to_ECEF(double L_b, double lambda_b, double h_b, const double v_eb_n[3], const double C_b_n[9], double r_eb_e[3], double v_eb_e[3], double C_b_e[9])
{
  if (!std::isfinite(L_b) || !std::isfinite(lambda_b) || !std::isfinite(h_b)) {
    return false;
  }

  const double sLat = std::sin(L_b);
  const double cLat = std::cos(L_b);
  const double sLon = std::sin(lambda_b);
  const double cLon = std::cos(lambda_b);

  const double denom = std::sqrt(1.0 - WGS84_E2 * sLat * sLat);
  if (denom <= 0.0 || !std::isfinite(denom)) {
    return false;
  }

  const double R_E = WGS84_A / denom;

  const double Rc = (R_E + h_b) * cLat;
  r_eb_e[0] = Rc * cLon;
  r_eb_e[1] = Rc * sLon;
  r_eb_e[2] = (ONE_MINUS_E2 * R_E + h_b) * sLat;

  const double c11 = -sLat * cLon;
  const double c12 = -sLon;
  const double c13 = -cLat * cLon;

  const double c21 = -sLat * sLon;
  const double c22 =  cLon;
  const double c23 = -cLat * sLon;

  const double c31 =  cLat;
  const double c32 =  0.0;
  const double c33 = -sLat;

  const double vn = v_eb_n[0];
  const double ve = v_eb_n[1];
  const double vd = v_eb_n[2];

  v_eb_e[0] = c11 * vn + c12 * ve + c13 * vd;
  v_eb_e[1] = c21 * vn + c22 * ve + c23 * vd;
  v_eb_e[2] = c31 * vn + c32 * ve + c33 * vd;

  C_b_e[0] = c11 * C_b_n[0] + c12 * C_b_n[1] + c13 * C_b_n[2];
  C_b_e[1] = c21 * C_b_n[0] + c22 * C_b_n[1] + c23 * C_b_n[2];
  C_b_e[2] = c31 * C_b_n[0] + c32 * C_b_n[1] + c33 * C_b_n[2];

  C_b_e[3] = c11 * C_b_n[3] + c12 * C_b_n[4] + c13 * C_b_n[5];
  C_b_e[4] = c21 * C_b_n[3] + c22 * C_b_n[4] + c23 * C_b_n[5];
  C_b_e[5] = c31 * C_b_n[3] + c32 * C_b_n[4] + c33 * C_b_n[5];

  C_b_e[6] = c11 * C_b_n[6] + c12 * C_b_n[7] + c13 * C_b_n[8];
  C_b_e[7] = c21 * C_b_n[6] + c22 * C_b_n[7] + c23 * C_b_n[8];
  C_b_e[8] = c31 * C_b_n[6] + c32 * C_b_n[7] + c33 * C_b_n[8];

  return true;
}

bool NED_to_ECEF(double L_b, double lambda_b, double h_b, const double v_eb_n[3], const double C_b_n[9], double C_b_e[9])
{
  double temp_r_eb_e[3];
  double temp_v_eb_e[3];

  return NED_to_ECEF(L_b, lambda_b, h_b, v_eb_n, C_b_n, temp_r_eb_e, temp_v_eb_e, C_b_e);
}

bool NED_to_ECEF(double L_b, double lambda_b, double h_b, const double v_eb_n[3], const double C_b_n[9])
{
  double temp_r_eb_e[3];
  double temp_v_eb_e[3];
  double temp_C_b_e[9];

  return NED_to_ECEF(L_b, lambda_b, h_b, v_eb_n, C_b_n, temp_r_eb_e, temp_v_eb_e, temp_C_b_e);
}