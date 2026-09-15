#ifndef PV_ECEF_TO_NED_H
#define PV_ECEF_TO_NED_H

bool pv_ECEF_to_NED(const double r_eb_e[3], const double v_eb_e[3], double &old_est_L_b, double &old_est_lambda_b, double &old_est_h_b, double old_est_v_eb_n[3]);

#endif  // PV_ECEF_TO_NED_H