#include "All_Configs.h"
#include "Globals.h"


// Files
// File file_REF;
// File file_GNSS;
// File file_RLG;
// File Output_Profile;

// CSV row buffers
// char row_REF[LINE_BUF_SIZE];
// char row_GNSS[LINE_BUF_SIZE];
// char row_RLG[LINE_BUF_SIZE];

// uint32_t rowPtr_REF  = 1;
// uint32_t rowPtr_GNSS = 1;
// uint32_t rowPtr_RLG  = 1;

// bool okREF  = false;
// bool okGNSS = false;
// bool okRLG  = false;

// bool pREF  = false;
// bool pGNSS = false;
// bool pRLG  = false;

// bool Start_GNSS_READ = false;
bool GNSS_EPOCH_AVAILABLE = false;
double lastAcceptedGnssTime = -1.0;
double latestGNSSTime =0.0;

// Time/constants
double old_time =0.0;
double tor_rem = 0.0;
const float TOLERANCE = 1.5f;

// Initial attitude / navigation states
double attitude_init[3] = {0.0, 0.0, 0.0};
double psi_nb_deg = 0.0;

double old_est_C_b_n[9] = {0.0};
double old_est_C_b_e[9] = {0.0};
double old_est_r_eb_e[3] = {0.0};
double old_est_v_eb_e[3] = {0.0};
double range_rate_Bias ={0.0};

double old_est_L_b = 0.0;
double old_est_lambda_b = 0.0;
double old_est_h_b = 0.0;
double old_est_v_eb_n[3] = {0.0};

// Large matrices in PSRAM
double P_matrix[P_DIM]={0.0};                  //  P_DIM = 17*17
double R_Matrix[R_DIM]={0.0};                  //  R-DIM = 15x15
double R_Matrix_NextLoop[R_DIM]={0.0};

// Lever arms / state variables
double STA[3] = {0.0, 0.0, 0.0};

double prev_w_ib_b_raw[3] = {0.0, 0.0, 0.0};
double est_IMU_Bias[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

extern const double L_ia_b[3] = {0, 0, 0}; // IMU antenna levaer arm
extern const double L_imu_b[3] = {0, 0, 0}; // IMU body COG lever arm
double L_ba_b[3] = {0.0, 0.0, 0.0}; // Body COG to Antenna lever arm

double L_imu_b_current[3] = {0.0, 0.0, 0.0};
double alpha[3] = {0.0, 0.0, 0.0};

double est_r_eb_e[3] = {0.0, 0.0, 0.0};
double est_v_eb_e[3] = {0.0, 0.0, 0.0};
double est_C_b_e[9] = {0.0};

double est_r_eb_e_temp_old[3] = {0.0, 0.0, 0.0};

// INS circular buffer
DMAMEM double INS_time_store[INS_BUFFER_SIZE]= {0.0};
DMAMEM double INS_r_store[INS_BUFFER_SIZE][3]= {0.0};
DMAMEM double INS_v_store[INS_BUFFER_SIZE][3]= {0.0};
DMAMEM double INS_C_store[INS_BUFFER_SIZE][9]= {0.0};
DMAMEM double INS_f_store[INS_BUFFER_SIZE][3]= {0.0};
DMAMEM double INS_w_store[INS_BUFFER_SIZE][3]= {0.0};

//Newly added INS buffer
DMAMEM double INS_f_raw_store[INS_BUFFER_SIZE][3]= {0.0};
DMAMEM double INS_w_raw_store[INS_BUFFER_SIZE][3]= {0.0};
DMAMEM float INS_P_store[INS_BUFFER_SIZE][P_DIM]= {0.0};
DMAMEM double INS_bias_store[INS_BUFFER_SIZE][BIAS_DIM];
DMAMEM double INS_clock_store[INS_BUFFER_SIZE][CLOCK_DIM]= {0.0};
DMAMEM double INS_dopp_store[INS_BUFFER_SIZE]= {0.0};
bool INS_valid_store[INS_BUFFER_SIZE]= {0.0};

unsigned long IMU_epoch_counter = 0;
int latest_INS_store_idx = -1;

double matchedGnssTime = 0.0;
double time_last_GNSS= 0.0;

double current_DGyro = 0.0;
double stationaryperiod = 0.0;

// Satellite/debug variables
int tot_sat[MAX_GNSS_ROWS_PER_EPOCH]= {0.0};
int sat_cnt_epoch = 1;

extern const double micro_g_to_meters_per_second_squared = 9.80665e-6;

// // Manual outage parameters
// double end_time = 0.0;

// extern const int first_gap = 100;
// int outage_duration = 20;

// double total_time = 0.0;
// double spacing = 0.0;


// double Outage_start[NUM_OUTAGES];
// double Outage_end[NUM_OUTAGES];
// double Outage_satCount[NUM_OUTAGES];

// Packet and row objects
EpochPacket packet;
REFRow refData;
REFRow neo7m;
GNSSRow gnssData;
RLGRow rlgData;

//===================================================//
//                Global Variables   (Dummy code)                //
//===================================================//
// String rawBuffer;
// double gps_obs = 0.0;
// TcaData currentTca;

// Config objects
CoGConfig CoG_config;
GNSSConfig GNSS_config;
TC_KFConfig TC_KF_config;
MLSConfig MLS_config;
LLASIMconfig LLASIM_config;