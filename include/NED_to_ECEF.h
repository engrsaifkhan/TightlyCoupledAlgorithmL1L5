#ifndef NED_TO_ECEF_H
#define NED_TO_ECEF_H

#include <cmath>
#include <cstdint>

// MATLAB-compatible column-major 3x3 matrices
// C_b_n and C_b_e are stored as:
// [0] r1c1, [1] r2c1, [2] r3c1,
// [3] r1c2, [4] r2c2, [5] r3c2,
// [6] r1c3, [7] r2c3, [8] r3c3

bool NED_to_ECEF(
    double L_b,                  // latitude (rad)
    double lambda_b,             // longitude (rad)
    double h_b,                  // height (m)
    const double v_eb_n[3],      // NED velocity [vn, ve, vd]
    const double C_b_n[9],       // body-to-NED DCM (column-major)
    double r_eb_e[3],            // output ECEF position
    double v_eb_e[3],            // output ECEF velocity
    double C_b_e[9]              // output body-to-ECEF DCM
);

bool NED_to_ECEF(
    double L_b,
    double lambda_b,
    double h_b,
    const double v_eb_n[3],
    const double C_b_n[9],
    double C_b_e[9]              // output body-to-ECEF DCM
);

#endif