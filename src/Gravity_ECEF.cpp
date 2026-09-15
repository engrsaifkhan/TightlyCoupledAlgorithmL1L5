#include "Gravity_ECEF.h"
#include <math.h>

void Gravity_ECEF(const double r_eb_e[3], double g[3])
{
    // WGS84 constants
    const double R_0 = 6378137.0;          // Equatorial radius (m)
    const double mu = 3.986004418e14;      // Earth gravitational constant (m^3/s^2)
    const double J_2 = 1.082627e-3;        // Second zonal harmonic
    const double omega_ie = 7.292115e-5;   // Earth rotation rate (rad/s)

    // Magnitude of position vector
    const double x = r_eb_e[0];
    const double y = r_eb_e[1];
    const double z = r_eb_e[2];

    // Calculate Distance from Center of the Earth
    const double mag_r = sqrt(x * x + y * y + z * z);

    // Dummy output if position is zero
    if (mag_r == 0.0)
    {
      g[0] = 0.0;
      g[1] = 0.0;
      g[2] = 0.0;
      return;
    }

    // Common terms
    const double mag_r2 = mag_r * mag_r;
    const double mag_r3 = mag_r2 * mag_r;
    const double z_ratio = z / mag_r;
    const double z_scale = 5.0 * z_ratio * z_ratio;
    const double r0_over_r_sq = (R_0 / mag_r) * (R_0 / mag_r);

    const double factor = -mu / mag_r3;
    const double J2_term = 1.5 * J_2 * r0_over_r_sq;

    // Gravitational part gamma
    const double gamma_x = factor * (x + J2_term * (1.0 - z_scale) * x);
    const double gamma_y = factor * (y + J2_term * (1.0 - z_scale) * y);
    const double gamma_z = factor * (z + J2_term * (3.0 - z_scale) * z);

    // Centrifugal acceleration
    // -cross(omega_vec, cross(omega_vec, r_eb_e))
    // for omega_vec = [0, 0, omega_ie]
    const double omega2 = omega_ie * omega_ie;
    const double centrifugal_x = omega2 * x;
    const double centrifugal_y = omega2 * y;
    const double centrifugal_z = 0.0;

    // Total gravity
    g[0] = gamma_x + centrifugal_x;
    g[1] = gamma_y + centrifugal_y;
    g[2] = gamma_z + centrifugal_z;
}