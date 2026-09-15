#ifndef GLOBALS_H
#define GLOBALS_H

#include <Arduino.h>
// #include <SD.h>
#include "All_Configs.h"

// Files
// extern File file_REF;
// extern File file_GNSS;
// extern File file_RLG;
// extern File Output_Profile;

// // CSV row buffers
// extern char row_REF[LINE_BUF_SIZE];
// extern char row_GNSS[LINE_BUF_SIZE];
// extern char row_RLG[LINE_BUF_SIZE];

// extern uint32_t rowPtr_REF;
// extern uint32_t rowPtr_GNSS;
// extern uint32_t rowPtr_RLG;

// extern bool okREF;
// extern bool okGNSS;
// extern bool okRLG;
// extern bool pREF;
// extern bool pGNSS;
// extern bool pRLG;

// extern bool Start_GNSS_READ;
extern bool GNSS_EPOCH_AVAILABLE;
extern double lastAcceptedGnssTime;
extern double latestGNSSTime;
// Time/constants
extern double old_time;
extern double tor_rem;
extern const float TOLERANCE;

// Initial attitude / navigation states
extern double attitude_init[3];
extern double psi_nb_deg;

extern double old_est_C_b_n[9];
extern double old_est_C_b_e[9];
extern double old_est_r_eb_e[3];
extern double old_est_v_eb_e[3];
extern double range_rate_Bias;

extern double old_est_L_b;
extern double old_est_lambda_b;
extern double old_est_h_b;
extern double old_est_v_eb_n[3];

// Large matrices
extern double P_matrix[P_DIM];
extern double R_Matrix[R_DIM];
extern double R_Matrix_NextLoop[R_DIM];

// Lever arms / state variables
extern double STA[3];

extern double prev_w_ib_b_raw[3];
extern double est_IMU_Bias[6];

extern const double L_ia_b[3];
extern const double L_imu_b[3];
extern double L_ba_b[3];

extern double L_imu_b_current[3];
extern double alpha[3];

extern double est_r_eb_e[3];
extern double est_v_eb_e[3];
extern double est_C_b_e[9];
extern double est_r_eb_e_temp_old[3];

//newly added
static constexpr size_t STATE_DIM = 18;
static constexpr size_t BIAS_DIM  = 6;
static constexpr size_t CLOCK_DIM = 2;
// INS circular buffer
extern double INS_time_store[INS_BUFFER_SIZE];
extern double INS_r_store[INS_BUFFER_SIZE][3];
extern double INS_v_store[INS_BUFFER_SIZE][3];
extern double INS_C_store[INS_BUFFER_SIZE][9];
extern double INS_f_store[INS_BUFFER_SIZE][3];
extern double INS_w_store[INS_BUFFER_SIZE][3];

//Newly added INS buffer
extern double INS_f_raw_store[INS_BUFFER_SIZE][3];
extern double INS_w_raw_store[INS_BUFFER_SIZE][3];
extern float INS_P_store[INS_BUFFER_SIZE][P_DIM];
extern double INS_bias_store[INS_BUFFER_SIZE][BIAS_DIM];
extern double INS_clock_store[INS_BUFFER_SIZE][CLOCK_DIM];
extern double INS_dopp_store[INS_BUFFER_SIZE];
extern bool INS_valid_store[INS_BUFFER_SIZE];

extern unsigned long IMU_epoch_counter;
extern int latest_INS_store_idx;

extern double matchedGnssTime;
extern double time_last_GNSS;

extern double current_DGyro;
extern double stationaryperiod;

// Satellite/debug variables
extern int tot_sat[MAX_GNSS_ROWS_PER_EPOCH];
extern int sat_cnt_epoch;

extern const double micro_g_to_meters_per_second_squared;

// // Manual outage parameters
// extern double end_time;
// extern const int first_gap;
// extern int outage_duration;
// extern double total_time;
// extern double spacing;

// extern double Outage_start[NUM_OUTAGES];
// extern double Outage_end[NUM_OUTAGES];
// extern double Outage_satCount[NUM_OUTAGES];

// Packet and row objects
extern EpochPacket packet;
extern REFRow refData; //pocketSDR
extern REFRow neo7m; //neo7m
extern GNSSRow gnssData;
extern RLGRow rlgData;

//===================================================//
//                Global Variables                   //
//===================================================//
// extern String rawBuffer;
// extern double gps_obs;
// extern TcaData currentTca;

// Config objects
extern CoGConfig CoG_config;
extern GNSSConfig GNSS_config;
extern TC_KFConfig TC_KF_config;
extern MLSConfig MLS_config;
extern LLASIMconfig LLASIM_config;

#endif