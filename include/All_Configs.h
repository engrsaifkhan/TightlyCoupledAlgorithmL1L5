#ifndef ALL_CONFIGS_H
#define ALL_CONFIGS_H

#include <Arduino.h>
#include <stdint.h>
#include <stddef.h>


enum INSIntervalStatus {
  INS_INTERVAL_OK = 0,
  INS_INTERVAL_INVALID_GNSS_TIME,
  INS_INTERVAL_NOT_ENOUGH_INS,
  INS_INTERVAL_INVALID_BUFFER_TIME,
  INS_INTERVAL_GNSS_AHEAD_OF_INS,
  INS_INTERVAL_GNSS_OLDER_THAN_BUFFER,
  INS_INTERVAL_NOT_FOUND
};

static INSIntervalStatus lastINSIntervalStatus = INS_INTERVAL_NOT_FOUND;
// =====================================================
// Global compile-time constants
// =====================================================
static constexpr size_t MAX_GNSS_ROWS_PER_EPOCH = 15;

static constexpr size_t P_DIM = 17 * 17;

static constexpr size_t R_ROWS_MAX = 2 * MAX_GNSS_ROWS_PER_EPOCH + 1;
static constexpr size_t R_DIM = R_ROWS_MAX * R_ROWS_MAX;

// static constexpr size_t LINE_BUF_SIZE = 512;
static constexpr size_t INS_BUFFER_SIZE = 200;
// static constexpr size_t NUM_OUTAGES = 7;

#define GNSS_SIGNAL_LEN    20
// #define MAX_SATS        40 // dummy code
#define MAX_GNSS_ROWS_PER_EPOCH 16
#define GNSS_LINE_MAX 300

enum GNSSSignalID : uint8_t {
  SIG_UNKNOWN = 0,

  SIG_GPS_L1_CA        = 1,   // GPS L1 C/A
  SIG_GALILEO_E1B      = 2,   // Galileo E1B
  SIG_GALILEO_E1C      = 3,   // Galileo E1C
  SIG_GALILEO_E5A_I    = 4,   // Galileo E5a I
  SIG_GALILEO_E5B_I    = 5,   // Galileo E5b I
  SIG_BEIDOU_B1I       = 6,   // BeiDou B1I
  SIG_BEIDOU_B1CD      = 7,   // BeiDou B1CD
  SIG_BEIDOU_B2A_D     = 8,   // BeiDou B2a D
  SIG_BEIDOU_B2B       = 9,   // BeiDou B2B
  SIG_BEIDOU_B2I       = 10,  // BeiDou B2I
  SIG_BEIDOU_B3I       = 11,  // BeiDou B3I
  SIG_GLONASS_G1_CA    = 12   // GLONASS G1 C/A
};
// ==============================
// REFRow struct
// ==============================
struct REFRow {
  double lat;
  double lon;
  double geoid_sep = 0.0;
  double hei_msl = 0.0;
  double hei_ellipsoid = 0.0;
  bool geoid_valid = false;
  float hei;

  float speed;
  float course;

  float vel_north;
  float vel_east;
  float vel_down;

  float hdop;

  int year;
  int month;
  int day;

  int hour;
  int minute;
  int second;
  int centisecond;

  int satellites;

  double timestamp;
};

// struct REFRow {
//   double utc;
//   double timestamp;
//   double lat;
//   double lon;
//   double hei;
//   double vel_north;
//   double vel_east;
//   double vel_down;
// };

// ==============================
// RLGRow struct (old IMU version)
// ==============================
// struct RLGRow {
//   float gx, gy, gz;
//   float ax, ay, az;
//   float roll, pitch, yaw;
//   double lat, lon;
//   float height;
//   float vel_east, vel_north, vel_up;
//   float temperature;        //Temperature
//   double timestamp;
// };

// ==============================
// RLGRow struct (new IMU version)
// ==============================
// struct RLGRow {
//   // IMU angular velocity, deg/s
//   double gx, gy, gz;

//   // IMU acceleration, m/s^2
//   double ax, ay, az;

//   // Attitude, deg
//   int roll, pitch, yaw;

//   // Integrated INS/GNSS position
//   double lat, lon;
//   float height;

//   // Navigation velocity, converted to ENU style
//   float vel_east;
//   float vel_north;
//   float vel_up;

//   // GPS time from IMU packet
//   double timestamp;
//   double timeOffset;

//   // Frame information
//   uint8_t frameCount;
//   uint8_t gpsStatus;

//   // Temperatures, deg C
//   int gyro_x_temp;
//   int gyro_y_temp;
//   int gyro_z_temp;

//   int acc_x_temp;
//   int acc_y_temp;
//   int acc_z_temp;

//   int if_temp;

//   // Optional GPS-only solution
//   double gps_lat;
//   double gps_lon;
//   float gps_alt;

//   int gps_heading;
//   float gps_vel_north;
//   float gps_vel_east;
//   float gps_vel_up;
// };

// ==============================
// RLGRow struct
// WT61P data mapped into the existing RLG interface
// ==============================
struct RLGRow {
  // Angular velocity, deg/s
  double gx;
  double gy;
  double gz;

  // Acceleration, m/s^2
  double ax;
  double ay;
  double az;

  // Attitude, deg
  double roll;
  double pitch;
  double yaw;

  // Not supplied by WT61P.
  // Retained to avoid changing downstream interfaces.
  double lat;
  double lon;
  float height;

  float vel_east;
  float vel_north;
  float vel_up;

  // Local controller timestamp when the WT61P sample begins
  double timestamp;

  // WT61P has no GPS time-offset field
  double timeOffset;

  // Software sample counter for WT61P
  uint8_t frameCount;

  // WT61P has no GPS status
  uint8_t gpsStatus;

  // WT61P reports one internal temperature.
  // The same temperature is mapped to the existing fields.
  double gyro_x_temp;
  double gyro_y_temp;
  double gyro_z_temp;

  double acc_x_temp;
  double acc_y_temp;
  double acc_z_temp;

  double if_temp;

  // Not supplied by WT61P
  double gps_lat;
  double gps_lon;
  float gps_alt;

  double gps_heading;
  float gps_vel_north;
  float gps_vel_east;
  float gps_vel_up;
};


// struct RLGRow {
//   double timestamp;
//   double ax;
//   double ay;
//   double az;
//   double gx;
//   double gy;
//   double gz;
//   double lat;
//   double lon;
//   double height;
//   double vel_north;
//   double vel_east;
// };

//===================================================//
//                GNSS parsing Data Structures  (Dummy code)                  //
//===================================================//
// struct SatelliteData {
//     double pr_corr;
//     double prrate_corr;
//     double x, y, z;
//     double vx, vy, vz;
//     double clk_rate;
//     double elev;
//     double cn0;
//     String signalName;
// };

// struct TcaData {
//     double obssec;
//     int num_sats;
//     SatelliteData sats[MAX_SATS];
//     bool isValid;
// };

// ==============================
// GNSSRow struct
// ==============================
struct GNSSRow {
  double pr_corr;       // corrected pseudorange, m
  double prrate_corr;   // corrected pseudorange rate, m/s

  double sat_x;         // satellite ECEF X, m
  double sat_y;         // satellite ECEF Y, m
  double sat_z;         // satellite ECEF Z, m

  double sat_vx;        // satellite ECEF VX, m/s
  double sat_vy;        // satellite ECEF VY, m/s
  double sat_vz;        // satellite ECEF VZ, m/s

  double sat_clk_rate;      // satellite clock rate correction
  float el_deg;           // elevation angle, deg
  float cn0;            // carrier-to-noise ratio, dB-Hz

  char signalName[GNSS_SIGNAL_LEN];
  uint8_t signal;
};


struct EpochPacket { 
  size_t gnssCount;                 //  // number of Obeservations in Epoch             
  double gnss_obssec;
  GNSSRow gnssRows[MAX_GNSS_ROWS_PER_EPOCH];    // satellite measurements
  bool valid;
};



// ==============================
// CoG Configuration
// ==============================
struct CoGConfig {
  double m_s;              // structural mass (kg)
  double rs[3];            // structure CoG [x y z]
  double R;                // solid propellant radius (m)
  double L_f0;             // initial propellant length (m)
  double z_f0;             // propellant start with nose reference (m)
  double x_offset;         // offset in X-axis
  double y_offset;         // offset in Y-axis
  double rho;              // propellant density (kg/m^3)
  double r_b;              // burn rate (m/s)
  double t_burn_start;     // burn start time
};

// ==============================
// GNSS Configuration
// ==============================
struct GNSSConfig {
  double rangeInnovationThreshold;     // range residual threshold
  double rateInnovationThreshold;      // rate residual threshold
  double mask_angle;              // elevation mask angle (deg)
  double mask_SignalStrength;     // C/N0 mask (dB-Hz)

  int intend_no_GNSS_meas;        // desired number of GNSS measurements

  uint8_t omit[12];                // signals to omit
  size_t omit_count;              // number of valid entries in omit[]

  double est_clock[2];            // [clock offset (m), clock drift (m/s)]
};

// ==============================
// Tight-Coupled KF Configuration
// ==============================
struct TC_KFConfig {
  // Initial state covariance
  double init_att_unc;
  double init_vel_unc;
  double init_pos_unc;
  double init_b_a_unc;
  double init_b_g_unc;
  double init_clock_offset_unc;
  double init_clock_drift_unc;
  double rangerateBias;
  // Propagation noise covariance
  double gyro_noise_PSD;
  double accel_noise_PSD;
  double accel_bias_PSD;
  double gyro_bias_PSD;
  double clock_phase_PSD;
  double clock_freq_PSD;
  double rangerateBias_PSD;

  // Other KF parameters
  double ZARU_DGrySD;
  double StationaryDetectionHorizontalSpeed;
  double GYRO_STATIC_THRESH;
  double ACC_STATIC_THRESH;
  double StationaryDetectionAccelerometerThreshold;
  double StationaryDetectionHSTimeWindow;

  bool StationaryFlag;

};

// ==============================
// MLS Configuration
// ==============================
struct MLSConfig {
  int MAXloop;
  double DeltaPosi;
  double K0;
  double K1;
};

struct LLASIMconfig {
  double Lat_init;
  double Lon_init;
  double Alt_init;
  double vn_init;
  double ve_init;
  double vd_init;
};

// ==============================
// Init function prototypes
// ==============================
void InitCoGConfig(CoGConfig &cfg, double burn_start_time);
void InitGNSSConfig(GNSSConfig &cfg);
void InitTCKFConfig(TC_KFConfig &cfg);
void InitMLSConfig(MLSConfig &cfg);
void InitLLASIMconfig(LLASIMconfig &cfg);
#endif