#include "Euler_to_CTM.h"         

// Precomputed yaw constants from your generated code
// yaw = atan2(-0.55919290347074668, -0.82903757255504185)
// static constexpr double CY = -0.82903757255504185;  // cos(yaw)
// static constexpr double SY = -0.55919290347074668;  // sin(yaw)

bool Euler_to_CTM(const double eul[3], double C[9])
{
  const double phi   = eul[0];
  const double theta = eul[1];
  const double psi = eul[2];

  const double sphi   = sin(phi);
  const double cphi   = cos(phi);
  const double stheta = sin(theta);
  const double ctheta = cos(theta);
  const double spsi = sin(psi);
  const double cpsi = cos(psi);


  const double sphi_stheta = sphi * stheta;
  const double cphi_stheta = cphi * stheta;

  // Column-major layout kept same as MATLAB-generated output
  C[0] = ctheta * cpsi;
  C[3] = ctheta * spsi;
  C[6] = -stheta;

  C[1] = -(cphi * spsi) + (sphi_stheta * cpsi);
  C[4] =  (cphi * cpsi) + (sphi_stheta * spsi);
  C[7] =  sphi * ctheta;

  C[2] =  (sphi * spsi) + (cphi_stheta * cpsi);
  C[5] = -(sphi * cpsi) + (cphi_stheta * spsi);
  C[8] =  cphi * ctheta;
  

  return true;
}