#ifndef MLS_KF_H
#define MLS_KF_H

#include <Arduino.h>
#include <stddef.h>
#include <stdint.h>
#include "All_Configs.h"
#include "Globals.h"

bool MLS_KF_predict(
  const double tor_s,
  double C_interpt[9],
  double v_interpt[3],
  double r_interpt[3],
  double meas_f_ib_b[3],
  TC_KFConfig &TC_KF_Config,
  double meas_omega_ib_b[3],
  MLSConfig &MLS_Config,
  double P_matrix[324]);

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
  double P_matrix[324],
  double *R_matrix,
  double *R_matrix_NextLoop,
  size_t R_dim,
  GNSSConfig &gnssConfig,
  double est_C_b_e[9],
  double est_v_eb_e[3],
  double est_r_eb_e[3],
  double est_IMU_Bias[6]);

#endif
