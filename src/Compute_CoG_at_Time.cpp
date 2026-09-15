#include "Compute_CoG_at_Time.h"
#include <Arduino.h>
#include <math.h>
#include <algorithm>  

bool Compute_CoG_at_Time(double t_current, const CoGConfig& CoG_config, const double L_imu_b[3], double L_imu_b_current[3]) {
  // ---------------- Time since burn start ----------------
  double t_since_burn = t_current - CoG_config.t_burn_start;

  // ---------------- Initial CoG of structure + propellant ----------------
  double L_f0 = CoG_config.L_f0; //Initial Length of fuel
  double V_f0 = PI * CoG_config.R * CoG_config.R * L_f0;  //Initial volume of fuel
  double m_f0 = CoG_config.rho * V_f0;
  double z_f0 = CoG_config.z_f0 + (L_f0 / 2.0); //initial bottom position of fuel

  // Initial propellant CoG vector
  double rf0[3] = {
    CoG_config.x_offset,
    CoG_config.y_offset,
    z_f0
  };

  // ---------------- Total burn time ----------------
  double T_burn = L_f0 / CoG_config.r_b;

  // ---------------- Current propellant CoG ----------------
  double rf_current[3] = {0.0, 0.0, 0.0};
  double m_f_current = 0.0;
  double L_f = 0.0;
  double V_f = 0.0;
  double z_f = 0.0;

  if (t_since_burn < 0.0) {
    // Burn has not started yet
    L_f = CoG_config.L_f0;
    V_f = PI * CoG_config.R * CoG_config.R * L_f;
    m_f_current = CoG_config.rho * V_f;
    z_f = CoG_config.z_f0 + (L_f / 2.0); //shift in fuel CoG in Z axis

    rf_current[0] = CoG_config.x_offset;
    rf_current[1] = CoG_config.y_offset;
    rf_current[2] = z_f;
  }
  else if (t_since_burn > T_burn) {
    // Burn complete
    L_f = 0.0;
    m_f_current = 0.0;

    // Not really used when fuel mass is zero, but kept safe
    rf_current[0] = CoG_config.rs[0];
    rf_current[1] = CoG_config.rs[1];
    rf_current[2] = CoG_config.rs[2];
  }
  else {
    // During burn
    double x_burn = CoG_config.r_b * t_since_burn;
    L_f = std::max(CoG_config.L_f0 - x_burn, 0.0);
    V_f = PI * CoG_config.R * CoG_config.R * L_f;
    m_f_current = CoG_config.rho * V_f;
    z_f = CoG_config.z_f0 + (L_f / 2.0); //shift in fuel CoG in Z axis

    rf_current[0] = CoG_config.x_offset;
    rf_current[1] = CoG_config.y_offset;
    rf_current[2] = z_f;
  }

  // ---------------- Current total CoG ----------------
  double r_cg_total[3] = {0.0, 0.0, 0.0};

  if (m_f_current > 0.0) {
    double denom = CoG_config.m_s + m_f_current;

    if (denom <= 0.0) {
      return false;
    }

    for (int i = 0; i < 3; i++) {
      r_cg_total[i] =
        (CoG_config.m_s * CoG_config.rs[i] + m_f_current * rf_current[i]) / denom;
    }
  }
  else {
    for (int i = 0; i < 3; i++) {
      r_cg_total[i] = CoG_config.rs[i];
    }
  }

  double initial_fuel_cog_z =
    CoG_config.z_f0 + CoG_config.L_f0 / 2.0;

  L_imu_b_current[0] =
      r_cg_total[0] - CoG_config.rs[0] + L_imu_b[0];

  L_imu_b_current[1] =
      r_cg_total[1] - CoG_config.rs[1] + L_imu_b[1];

  L_imu_b_current[2] =
      r_cg_total[2]
      - (CoG_config.rs[2] + initial_fuel_cog_z)
      + L_imu_b[2];

  return true;
}  
