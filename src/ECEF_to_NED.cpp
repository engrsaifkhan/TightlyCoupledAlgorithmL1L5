#include "ECEF_to_NED.h"
#include <cmath>

bool ECEF_to_NED(const double est_C_b_e[9], const double est_v_eb_e[3], const double est_r_eb_e[3], double &L_b, double &Lambda_b, double &h_b, double v_eb_n[3], double C_b_n[9]) {
    const double R_0  = 6378137.0;
    const double ecc2 = 0.006694380022900;

    const double x = est_r_eb_e[0];
    const double y = est_r_eb_e[1];
    const double z = est_r_eb_e[2];

    double sign_z = 0.0;
    if (z > 0.0) {
        sign_z = 1.0;
    } else if (z < 0.0) {
        sign_z = -1.0;
    }

    //----- Convert position using Borkowski closed-form exact solution From (2.113)
    Lambda_b = atan2(y, x);

    double beta = sqrt(x*x + y*y);
    if (beta == 0.0) return false;
    
    //-- From (C.29) and (C.30)-----------
    double k1 = sqrt(1.0 - ecc2) * fabs(z);
    double k2 = ecc2 * R_0;
    double E  = (k1 - k2) / beta;
    double F  = (k1 + k2) / beta;

    //--- From (C.31) ------
    double P = (4.0 / 3.0) * (E * F + 1.0);
    //--- From (C.32) ------
    double Q = 2.0 * (E * E - F * F);
    //--- From (C.33) ------
    double D = P * P * P + Q * Q;
    //--- From (C.34) ------
    double V = cbrt(sqrt(D) - Q) - cbrt(sqrt(D) + Q);
    //--- From (C.35) ------
    double G = 0.5 * (sqrt(E * E + V) + E);
    //--- From (C.36) ------
    double T = sqrt(G * G + (F - V * G) / (2.0 * G - E)) - G;
    //--- From (C.37) ------
    L_b = sign_z * atan((1.0 - T * T) / (2.0 * T * sqrt(1.0 - ecc2)));
    //--- From (C.38) ------
    h_b = (beta - R_0 * T) * cos(L_b)
        + (z - sign_z * R_0 * sqrt(1.0 - ecc2)) * sin(L_b);

    // Calculate ECEF to NED coordinate transformation matrix using (2.150)
    double cos_lat  = cos(L_b);
    double sin_lat  = sin(L_b);
    double cos_long = cos(Lambda_b);
    double sin_long = sin(Lambda_b);

    // C_e_n in column-major
    double C_e_n[9];
    C_e_n[0] = -sin_lat * cos_long;
    C_e_n[1] = -sin_long;
    C_e_n[2] = -cos_lat * cos_long;

    C_e_n[3] = -sin_lat * sin_long;
    C_e_n[4] =  cos_long;
    C_e_n[5] = -cos_lat * sin_long;

    C_e_n[6] =  cos_lat;
    C_e_n[7] =  0.0;
    C_e_n[8] = -sin_lat;

    // Transform velocity using (2.73)
    v_eb_n[0] = C_e_n[0]*est_v_eb_e[0] + C_e_n[3]*est_v_eb_e[1] + C_e_n[6]*est_v_eb_e[2];
    v_eb_n[1] = C_e_n[1]*est_v_eb_e[0] + C_e_n[4]*est_v_eb_e[1] + C_e_n[7]*est_v_eb_e[2];
    v_eb_n[2] = C_e_n[2]*est_v_eb_e[0] + C_e_n[5]*est_v_eb_e[1] + C_e_n[8]*est_v_eb_e[2];

    // Transform attitude using (2.15)
    for (int col = 0; col < 3; col++) {
        C_b_n[0 + 3*col] = C_e_n[0]*est_C_b_e[0 + 3*col] + C_e_n[3]*est_C_b_e[1 + 3*col] + C_e_n[6]*est_C_b_e[2 + 3*col];
        C_b_n[1 + 3*col] = C_e_n[1]*est_C_b_e[0 + 3*col] + C_e_n[4]*est_C_b_e[1 + 3*col] + C_e_n[7]*est_C_b_e[2 + 3*col];
        C_b_n[2 + 3*col] = C_e_n[2]*est_C_b_e[0 + 3*col] + C_e_n[5]*est_C_b_e[1 + 3*col] + C_e_n[8]*est_C_b_e[2 + 3*col];
    }

    return true;
}