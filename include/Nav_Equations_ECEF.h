#ifndef NAV_EQUATIONS_ECEF_H
#define NAV_EQUATIONS_ECEF_H

#include <stdbool.h>

bool Nav_Equations_ECEF(
  double tor_i,
  const double old_est_r_eb_e[3],
  const double old_est_v_eb_e[3],
  const double old_C_b_e[9],
  const double f_ib_b[3],
  const double omega_ib_b[3],
  double est_r_eb_e[3],
  double est_v_eb_e[3],
  double est_C_b_e[9]
);

bool Nav_Equations_ECEF(
  double tor_i,
  double old_est_r_eb_e[3],
  double old_est_v_eb_e[3],
  double old_C_b_e[9],
  const double f_ib_b[3],
  const double omega_ib_b[3]
);

#endif