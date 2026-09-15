#ifndef COMPUTE_COG_AT_TIME_H
#define COMPUTE_COG_AT_TIME_H

#include "All_Configs.h"

bool Compute_CoG_at_Time(double t_current, const CoGConfig &CoG_config, const double L_imu_b[3], double L_imu_b_current[3]);
  
#endif