#include "Nav_Equations_ECEF.h"
#include "ReOrthonormalizeDCM_matlab.h"
#include "Gravity_ECEF.h"
#include "Skew_Symmetric.h"
#include <Arduino.h>
#include <math.h>

bool Nav_Equations_ECEF(double tor_i, const double old_r_eb_e[3], const double old_v_eb_e[3], const double old_C_b_e[9], 
                        const double f_ib_b[3], const double omega_ib_b[3],
                        double est_r_eb_e[3], double est_v_eb_e[3], double est_C_b_e[9]) {

  // Nav_equations_ECEF - Runs precision ECEF-frame inertial navigation equations

  // Inputs:
  //   tor_i         time interval between epochs (s)
  //   old_r_eb_e    previous Cartesian position of body frame w.r.t. ECEF frame, resolved along ECEF-frame axes (m)
  //   old_C_b_e     previous body-to-ECEF-frame coordinate transformation matrix
  //   old_v_eb_e    previous velocity of body frame w.r.t. ECEF frame, resolved along ECEF-frame axes (m/s)
  //   f_ib_b        specific force of body frame w.r.t. ECEF frame, resolved along body-frame axes, averaged over time interval (m/s^2)
  //   omega_ib_b    angular rate of body frame w.r.t. ECEF frame, resolved about body-frame axes, averaged over time interval (rad/s)
    
  // Outputs:
  //   est_r_eb_e        Cartesian position of body frame w.r.t. ECEF frame, resolved along ECEF-frame axes (m)
  //   est_v_eb_e        velocity of body frame w.r.t. ECEF frame, resolved along ECEF-frame axes (m/s)
  //   est_ C_b_e         body-to-ECEF-frame coordinate transformation matrix

  // Parameters
  double omega_ie = 7.292115E-5;                 // Earth rotation rate (rad/s)


  // ATTITUDE UPDATE
  // From (2.145) determine the Earth rotation over the update interval
  // C_Earth = C_e_i' * old_C_e_i
  double alpha_ie = omega_ie * tor_i;

  // //Debug Print
  // Serial.println(F("-----------------------------alpha_ie------------------------------"));
  // Serial.print("alpha_ie:");
  // Serial.println(alpha_ie, 7);
  
  double C_Earth[9] = {cos(alpha_ie), -sin(alpha_ie), 0, sin(alpha_ie), cos(alpha_ie), 0, 0, 0, 1};

  // //Debug Print
  // Serial.println(F("-----------------------------C_Earth------------------------------"));
  // Serial.print("C_Earth[0]:");
  // Serial.println(C_Earth[0], 7);
  // Serial.print("C_Earth[1]:");
  // Serial.println(C_Earth[1], 7);
  // Serial.print("C_Earth[2]:");
  // Serial.println(C_Earth[2], 7);
  // Serial.print("C_Earth[3]:");
  // Serial.println(C_Earth[3], 7);
  // Serial.print("C_Earth[4]:");
  // Serial.println(C_Earth[4], 7);
  // Serial.print("C_Earth[5]:");
  // Serial.println(C_Earth[5], 7);
  // Serial.print("C_Earth[6]:");
  // Serial.println(C_Earth[6], 7);
  // Serial.print("C_Earth[7]:");
  // Serial.println(C_Earth[7], 7);
  // Serial.print("C_Earth[8]:");
  // Serial.println(C_Earth[8], 7);

  double alpha_ib_b[3] = {omega_ib_b[0] * tor_i,
                          omega_ib_b[1] * tor_i,
                          omega_ib_b[2] * tor_i};

  // //Debug Print
  // Serial.println(F("-----------------------------alpha_ib_b------------------------------"));
  // Serial.print("alpha_ib_b[0]:");
  // Serial.println(alpha_ib_b[0], 7);
  // Serial.print("alpha_ib_b[1]:");
  // Serial.println(alpha_ib_b[1], 7);
  // Serial.print("alpha_ib_b[2]:");
  // Serial.println(alpha_ib_b[2], 7);

  double mag_alpha = sqrt(alpha_ib_b[0] * alpha_ib_b[0] + alpha_ib_b[1] * alpha_ib_b[1] + alpha_ib_b[2] * alpha_ib_b[2]);

  //Debug Print
  // Serial.println(F("-----------------------------mag_alpha------------------------------"));
  // Serial.print("mag_alpha:");
  // Serial.println(mag_alpha, 7);

  double Alpha_ib_b[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  Skew_Symmetric(alpha_ib_b, Alpha_ib_b);

  // //Debug Print
  // Serial.println(F("-----------------------------Alpha_ib_b------------------------------"));
  // Serial.print("Alpha_ib_b[0]:");
  // Serial.println(Alpha_ib_b[0], 7);
  // Serial.print("Alpha_ib_b[1]:");
  // Serial.println(Alpha_ib_b[1], 7);
  // Serial.print("Alpha_ib_b[2]:");
  // Serial.println(Alpha_ib_b[2], 7);
  // Serial.print("Alpha_ib_b[3]:");
  // Serial.println(Alpha_ib_b[3], 7);
  // Serial.print("Alpha_ib_b[4]:");
  // Serial.println(Alpha_ib_b[4], 7);
  // Serial.print("Alpha_ib_b[5]:");
  // Serial.println(Alpha_ib_b[5], 7);
  // Serial.print("Alpha_ib_b[6]:");
  // Serial.println(Alpha_ib_b[6], 7);
  // Serial.print("Alpha_ib_b[7]:");
  // Serial.println(Alpha_ib_b[7], 7);
  // Serial.print("Alpha_ib_b[8]:");
  // Serial.println(Alpha_ib_b[8], 7);


  // Obtain Coordinate Transformation Matrix From The New Attitude W.R.T. An Inertial Frame To The Old Using Rodrigues' Formula, (5.73) Or 4th Order Approximation (5.72)
  double Alpha_sq[9];
  double C_new_old[9];

  // Alpha_sq = Alpha_ib_b * Alpha_ib_b
  for (int col = 0; col < 3; col++) {
    for (int row = 0; row < 3; row++) {
      Alpha_sq[row + 3 * col] = Alpha_ib_b[row + 3 * 0] * Alpha_ib_b[0 + 3 * col] + Alpha_ib_b[row + 3 * 1] * Alpha_ib_b[1 + 3 * col] +
      Alpha_ib_b[row + 3 * 2] * Alpha_ib_b[2 + 3 * col];
    }
  }

  if (mag_alpha > 1.0E-8) {
    double c1 = sin(mag_alpha) / mag_alpha;
    double c2 = (1.0 - cos(mag_alpha)) / (mag_alpha * mag_alpha);

    for (int i = 0; i < 9; i++) {
      C_new_old[i] = c1 * Alpha_ib_b[i] + c2 * Alpha_sq[i];
    }
  } else {
    double mag2 = mag_alpha * mag_alpha;
    double c1 = 1.0 - (mag2 / 6.0);
    double c2 = 0.5 - (mag2 / 24.0);

    for (int i = 0; i < 9; i++) {
      C_new_old[i] = c1 * Alpha_ib_b[i] + c2 * Alpha_sq[i];
    }
  }

  // Add eye(3)
  C_new_old[0] += 1.0;   // (1,1)
  C_new_old[4] += 1.0;   // (2,2)
  C_new_old[8] += 1.0;   // (3,3)

  // //Debug Print
  // Serial.println(F("-----------------------------C_new_old------------------------------"));
  // Serial.print("C_new_old[0]:");
  // Serial.println(C_new_old[0], 7);
  // Serial.print("C_new_old[1]:");
  // Serial.println(C_new_old[1], 7);
  // Serial.print("C_new_old[2]:");
  // Serial.println(C_new_old[2], 7);
  // Serial.print("C_new_old[3]:");
  // Serial.println(C_new_old[3], 7);
  // Serial.print("C_new_old[4]:");
  // Serial.println(C_new_old[4], 7);
  // Serial.print("C_new_old[5]:");
  // Serial.println(C_new_old[5], 7);
  // Serial.print("C_new_old[6]:");
  // Serial.println(C_new_old[6], 7);
  // Serial.print("C_new_old[7]:");
  // Serial.println(C_new_old[7], 7);
  // Serial.print("C_new_old[8]:");
  // Serial.println(C_new_old[8], 7);


  double temp[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  // temp = C_Earth * old_C_b_e
  for (int col = 0; col < 3; col++) {
    for (int row = 0; row < 3; row++) {
      temp[row + 3 * col] = C_Earth[row + 3 * 0] * old_C_b_e[0 + 3 * col] + C_Earth[row + 3 * 1] * old_C_b_e[1 + 3 * col] + C_Earth[row + 3 * 2] * old_C_b_e[2 + 3 * col];
    }
  }

  // est_C_b_e = temp * C_new_old
  for (int col = 0; col < 3; col++) {
    for (int row = 0; row < 3; row++) {
      est_C_b_e[row + 3 * col] = temp[row + 3 * 0] * C_new_old[0 + 3 * col] + temp[row + 3 * 1] * C_new_old[1 + 3 * col] + temp[row + 3 * 2] * C_new_old[2 + 3 * col];
    }
  }

  // //Debug Print
  // Serial.println(F("-----------------------------est_C_b_e------------------------------"));
  // Serial.print("est_C_b_e[0]:");
  // Serial.println(est_C_b_e[0], 7);
  // Serial.print("est_C_b_e[1]:");
  // Serial.println(est_C_b_e[1], 7);
  // Serial.print("est_C_b_e[2]:");
  // Serial.println(est_C_b_e[2], 7);
  // Serial.print("est_C_b_e[3]:");
  // Serial.println(est_C_b_e[3], 7);
  // Serial.print("est_C_b_e[4]:");
  // Serial.println(est_C_b_e[4], 7);
  // Serial.print("est_C_b_e[5]:");
  // Serial.println(est_C_b_e[5], 7);
  // Serial.print("est_C_b_e[6]:");
  // Serial.println(est_C_b_e[6], 7);
  // Serial.print("est_C_b_e[7]:");
  // Serial.println(est_C_b_e[7], 7);
  // Serial.print("est_C_b_e[8]:");
  // Serial.println(est_C_b_e[8], 7);

  // Re-Orthonormalize DCM
  ReOrthonormalizeDCM_matlab(est_C_b_e);

  // //Debug Print
  // Serial.println(F("-----------------------------est_C_b_e------------------------------"));
  // Serial.print("est_C_b_e[0]:");
  // Serial.println(est_C_b_e[0], 7);
  // Serial.print("est_C_b_e[1]:");
  // Serial.println(est_C_b_e[1], 7);
  // Serial.print("est_C_b_e[2]:");
  // Serial.println(est_C_b_e[2], 7);
  // Serial.print("est_C_b_e[3]:");
  // Serial.println(est_C_b_e[3], 7);
  // Serial.print("est_C_b_e[4]:");
  // Serial.println(est_C_b_e[4], 7);
  // Serial.print("est_C_b_e[5]:");
  // Serial.println(est_C_b_e[5], 7);
  // Serial.print("est_C_b_e[6]:");
  // Serial.println(est_C_b_e[6], 7);
  // Serial.print("est_C_b_e[7]:");
  // Serial.println(est_C_b_e[7], 7);
  // Serial.print("est_C_b_e[8]:");
  // Serial.println(est_C_b_e[8], 7);

  // SPECIFIC FORCE FRAME TRANSFORMATION
  // Calculate the average body-to-ECEF-frame coordinate transformation
  // matrix over the update interval using (5.84) 
  double A2[9];
  double temp1[9];
  double temp2[9];
  double temp3[9];
  double earth_skew_vec[3] = {0.0, 0.0, omega_ie};
  double Earth_skew[9];
  double ave_C_b_e[9];

  // A2 = Alpha_ib_b * Alpha_ib_b
  for (int col = 0; col < 3; col++) {
    for (int row = 0; row < 3; row++) {
      A2[row + 3 * col] = Alpha_ib_b[row + 3 * 0] * Alpha_ib_b[0 + 3 * col] + Alpha_ib_b[row + 3 * 1] * Alpha_ib_b[1 + 3 * col] + Alpha_ib_b[row + 3 * 2] * Alpha_ib_b[2 + 3 * col];
    }
  }

  // Earth_skew = Skew_symmetric([0;0;omega_ie])
  Skew_Symmetric(earth_skew_vec, Earth_skew);

  if (mag_alpha > 1.0E-8) {
    // temp1 = I + 0.5*Alpha_ib_b + (1/6)*A2
    for (int i = 0; i < 9; i++) {
      temp1[i] = 0.5 * Alpha_ib_b[i] + (1.0 / 6.0) * A2[i];
    }
    temp1[0] += 1.0;
    temp1[4] += 1.0;
    temp1[8] += 1.0;

    // temp2 = old_C_b_e * temp1
    for (int col = 0; col < 3; col++) {
      for (int row = 0; row < 3; row++) {
        temp2[row + 3 * col] = old_C_b_e[row + 3 * 0] * temp1[0 + 3 * col] + old_C_b_e[row + 3 * 1] * temp1[1 + 3 * col] + old_C_b_e[row + 3 * 2] * temp1[2 + 3 * col];
      }
    }

    // temp3 = Earth_skew * old_C_b_e
    for (int col = 0; col < 3; col++) {
        for (int row = 0; row < 3; row++) {
          temp3[row + 3 * col] = Earth_skew[row + 3 * 0] * old_C_b_e[0 + 3 * col] + Earth_skew[row + 3 * 1] * old_C_b_e[1 + 3 * col] + Earth_skew[row + 3 * 2] * old_C_b_e[2 + 3 * col];
        }
    }

    // ave_C_b_e = temp2 - 0.5 * temp3 * tor_i
    for (int i = 0; i < 9; i++) {
      ave_C_b_e[i] = temp2[i] - 0.5 * temp3[i] * tor_i;
    }
  }
  else {
    double mag2 = mag_alpha * mag_alpha;

    // temp1 = I + (1-mag^2/6)*Alpha_ib_b + (0.5-mag^2/24)*A2
    double c1 = 1.0 - (mag2 / 6.0);
    double c2 = 0.5 - (mag2 / 24.0);

    for (int i = 0; i < 9; i++) {
      temp1[i] = c1 * Alpha_ib_b[i] + c2 * A2[i];
    }
    temp1[0] += 1.0;
    temp1[4] += 1.0;
    temp1[8] += 1.0;

    // temp2 = old_C_b_e * temp1
    for (int col = 0; col < 3; col++) {
      for (int row = 0; row < 3; row++) {
        temp2[row + 3 * col] = old_C_b_e[row + 3 * 0] * temp1[0 + 3 * col] + old_C_b_e[row + 3 * 1] * temp1[1 + 3 * col] + old_C_b_e[row + 3 * 2] * temp1[2 + 3 * col];
      }
    }

    // temp3 = Earth_skew * old_C_b_e
    for (int col = 0; col < 3; col++) {
      for (int row = 0; row < 3; row++) {
        temp3[row + 3 * col] = Earth_skew[row + 3 * 0] * old_C_b_e[0 + 3 * col] + Earth_skew[row + 3 * 1] * old_C_b_e[1 + 3 * col] + Earth_skew[row + 3 * 2] * old_C_b_e[2 + 3 * col];
      }
    }

    // ave_C_b_e = temp2 - 0.5 * temp3 * tor_i
    for (int i = 0; i < 9; i++) {
      ave_C_b_e[i] = temp2[i] - 0.5 * temp3[i] * tor_i;
    }
  }

  // //Debug Print
  // Serial.println(F("-----------------------------ave_C_b_e------------------------------"));
  // Serial.print("ave_C_b_e[0]:");
  // Serial.println(ave_C_b_e[0], 7);
  // Serial.print("ave_C_b_e[1]:");
  // Serial.println(ave_C_b_e[1], 7);
  // Serial.print("ave_C_b_e[2]:");
  // Serial.println(ave_C_b_e[2], 7);
  // Serial.print("ave_C_b_e[3]:");
  // Serial.println(ave_C_b_e[3], 7);
  // Serial.print("ave_C_b_e[4]:");
  // Serial.println(ave_C_b_e[4], 7);
  // Serial.print("ave_C_b_e[5]:");
  // Serial.println(ave_C_b_e[5], 7);
  // Serial.print("ave_C_b_e[6]:");
  // Serial.println(ave_C_b_e[6], 7);
  // Serial.print("ave_C_b_e[7]:");
  // Serial.println(ave_C_b_e[7], 7);
  // Serial.print("ave_C_b_e[8]:");
  // Serial.println(ave_C_b_e[8], 7);

  // // Transform specific force to ECEF-frame resolving axes using (5.85)
  double f_ib_e[3] = {0.0, 0.0, 0.0};

  f_ib_e[0] = ave_C_b_e[0] * f_ib_b[0] + ave_C_b_e[3] * f_ib_b[1] + ave_C_b_e[6] * f_ib_b[2];
  f_ib_e[1] = ave_C_b_e[1] * f_ib_b[0] + ave_C_b_e[4] * f_ib_b[1] + ave_C_b_e[7] * f_ib_b[2];
  f_ib_e[2] = ave_C_b_e[2] * f_ib_b[0] + ave_C_b_e[5] * f_ib_b[1] + ave_C_b_e[8] * f_ib_b[2];

  // //Debug Print
  // Serial.println(F("-----------------------------f_ib_b------------------------------"));
  // Serial.print("f_ib_b[0]:");
  // Serial.println(f_ib_b[0], 7);
  // Serial.print("f_ib_b[1]:");
  // Serial.println(f_ib_b[1], 7);
  // Serial.print("f_ib_b[2]:");
  // Serial.println(f_ib_b[2], 7);
  
  // //Debug Print
  // Serial.println(F("-----------------------------f_ib_e------------------------------"));
  // Serial.print("f_ib_e[0]: ");
  // Serial.print(f_ib_e[0], 7);
  // Serial.print("f_ib_e[1]: ");
  // Serial.print(f_ib_e[1], 7);
  // Serial.print("f_ib_e[2]: ");
  // Serial.println(f_ib_e[2], 7);
  
  // UPDATE VELOCITY From (5.36)
  double g[3] = {0.0, 0.0, 0.0};
  double v_eb_e_pred[3] = {0.0, 0.0, 0.0};

  double omega_ie_vec[3] = {0.0, 0.0, omega_ie};
  double Omega_ie[9];
  double Omega_v[3];
  double v_sum[3];
  double Omega_vsum[3];

  Gravity_ECEF(old_r_eb_e, g);
  Skew_Symmetric(omega_ie_vec, Omega_ie);

  // Omega_v = Omega_ie * old_v_eb_e
  Omega_v[0] = Omega_ie[0] * old_v_eb_e[0] + Omega_ie[3] * old_v_eb_e[1] + Omega_ie[6] * old_v_eb_e[2];
  Omega_v[1] = Omega_ie[1] * old_v_eb_e[0] + Omega_ie[4] * old_v_eb_e[1] + Omega_ie[7] * old_v_eb_e[2];
  Omega_v[2] = Omega_ie[2] * old_v_eb_e[0] + Omega_ie[5] * old_v_eb_e[1] + Omega_ie[8] * old_v_eb_e[2];

  // Prediction
  v_eb_e_pred[0] = old_v_eb_e[0] + tor_i * (f_ib_e[0] + g[0] - 2.0 * Omega_v[0]);
  v_eb_e_pred[1] = old_v_eb_e[1] + tor_i * (f_ib_e[1] + g[1] - 2.0 * Omega_v[1]);
  v_eb_e_pred[2] = old_v_eb_e[2] + tor_i * (f_ib_e[2] + g[2] - 2.0 * Omega_v[2]);

  // //Debug Print
  // Serial.println(F("-----------------------------v_eb_e_pred------------------------------"));
  // Serial.print("v_eb_e_pred[0]:");
  // Serial.println(v_eb_e_pred[0], 7);
  // Serial.print("v_eb_e_pred[1]:");
  // Serial.println(v_eb_e_pred[1], 7);
  // Serial.print("v_eb_e_pred[2]:");
  // Serial.println(v_eb_e_pred[2], 7);

  // v_sum = old_v_eb_e + v_eb_e_pred
  v_sum[0] = old_v_eb_e[0] + v_eb_e_pred[0];
  v_sum[1] = old_v_eb_e[1] + v_eb_e_pred[1];
  v_sum[2] = old_v_eb_e[2] + v_eb_e_pred[2];

  // Omega_vsum = Omega_ie * (old_v_eb_e + v_eb_e_pred)
  Omega_vsum[0] = Omega_ie[0] * v_sum[0] + Omega_ie[3] * v_sum[1] + Omega_ie[6] * v_sum[2];
  Omega_vsum[1] = Omega_ie[1] * v_sum[0] + Omega_ie[4] * v_sum[1] + Omega_ie[7] * v_sum[2];
  Omega_vsum[2] = Omega_ie[2] * v_sum[0] + Omega_ie[5] * v_sum[1] + Omega_ie[8] * v_sum[2];

  // Final corrected update
  est_v_eb_e[0] = old_v_eb_e[0] + tor_i * (f_ib_e[0] + g[0] - 0.5 * 2.0 * Omega_vsum[0]);
  est_v_eb_e[1] = old_v_eb_e[1] + tor_i * (f_ib_e[1] + g[1] - 0.5 * 2.0 * Omega_vsum[1]);
  est_v_eb_e[2] = old_v_eb_e[2] + tor_i * (f_ib_e[2] + g[2] - 0.5 * 2.0 * Omega_vsum[2]);

  // //Debug Print
  // Serial.println(F("-----------------------------est_v_eb_e------------------------------"));
  // Serial.print("est_v_eb_e[0]:");
  // Serial.println(est_v_eb_e[0], 7);
  // Serial.print("est_v_eb_e[1]:");
  // Serial.println(est_v_eb_e[1], 7);
  // Serial.print("est_v_eb_e[2]:");
  // Serial.println(est_v_eb_e[2], 7);

  // UPDATE CARTESIAN POSITION From (5.38)
  for (int i = 0; i < 3; i++) {
    est_r_eb_e[i] = old_r_eb_e[i] + 0.5 * tor_i * (est_v_eb_e[i] + old_v_eb_e[i]);
  }

  // //Debug Print
  // Serial.println(F("-----------------------------est_r_eb_e------------------------------"));
  // Serial.print("est_r_eb_e[0]:");
  // Serial.println(est_r_eb_e[0], 7);
  // Serial.print("est_r_eb_e[1]:");
  // Serial.println(est_r_eb_e[1], 7);
  // Serial.print("est_r_eb_e[2]:");
  // Serial.println(est_r_eb_e[2], 7);

  return true;
}


