#ifndef ECEF_TO_NED_H
#define ECEF_TO_NED_H

#include <stddef.h>
#include <stdint.h>

bool ECEF_to_NED(const double est_C_b_e[9], const double est_v_eb_e[3], const double est_r_eb_e[3], double &L_b, double &Lambda_b, double &h_b, double v_eb_n[3], double C_b_n[9]);

#endif