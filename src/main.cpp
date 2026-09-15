// This code under test with sensor integration

// #include <SPI.h>
// #include <cmath>
// #include <vector>
// #include <cstring>
// #include <cstdlib>
// #include <stdint.h>
#include <Arduino.h>
#include <TinyGPSPlus.h>

#include "All_Configs.h"
#include "Globals.h"

#include "ch_LLA2ECEF.h"
#include "NED_to_ECEF.h"
#include "MLS_KF.h"
#include "Euler_to_CTM.h"
#include "CTM_to_Euler.h"
#include "pv_ECEF_to_NED.h"
#include "ECEF_to_NED.h"
#include "Nav_Equations_ECEF.h"
#include "Compute_CoG_at_Time.h"


// Forward declarations (required for .cpp; Arduino .ino auto-generated these)
void InitGNSSConfig(GNSSConfig& GNSS_Config);
void InitTCKFConfig(TC_KFConfig& TC_KF_Config);
void InitMLSConfig(MLSConfig& MLS_Config);
void InitLLASIMConfig(LLASIMconfig &LLASIM_Config);
void InitCoGConfig(CoGConfig& CoG_Config, double burn_start_time);
void clearSerialInput(Stream &s, uint32_t quiet_ms, uint32_t max_ms);
void setup();
void loop();
static inline void Transpose3x3CM_InPlace(double C[9]);
bool Initial_Initializations();
bool ComputeInitialAligmentMean(double& time_last_GNSS, double& old_time, double& mean_gx, double& mean_gy, double& mean_gz);
bool ReadNEO7M_PVT_NearestTime(double target_time, uint32_t timeout_ms, double max_allowed_dt_sec);
bool Get_ComputeThetaPhi(double mean_ax, double mean_ay, double mean_az, double mean_gx, double mean_gy, double mean_gz);
bool PVTAlignmentAtOldTime(double mean_gx, double mean_gy, double mean_gz, double& lat0, double& lon0, double& h0);
void Initialize_P_Matrix(const TC_KFConfig& TC_KF_Config);
bool FindINSIntervalForGNSSTime(double gnssTime, int &idx1, int &idx2);
bool waitForInitialGNSS(int requiredPackets, uint32_t timeoutMs);
bool read_NEO7M_PVT();
void ppsISR();
double getPreciseTime();
uint32_t daysSince1970(int Y, int M, int D);
double utcToGpsSeconds(int y, int m, int d, int h, int min, double sec);
static bool parsePythonDouble(const char* text, double& value);
static bool processPythonIMULine(char* line, double packetStartTime);
bool waitForInitialIMU(int requiredPackets, uint32_t timeoutMs);
bool readIMU();
bool IsGnssObssecAtIntegerSecond(double gps_obs);
bool waitForIntegerSecondPocketRawForSync(uint32_t timeoutMs);
static void trimString(char* str);
static bool parseTCAHeader(char* line);
uint8_t convertSignalNameToID(const char* signalName);
static bool parseGNSSObsRow(char* line, GNSSRow& row);
bool readGNSSData();
bool isEndLine(const char *line);
void resetGNSSParser();
bool processGNSSLine(char *line);
bool SyncPreciseTimeFromPocketRaw(double gps_obs);
void printSignalNameSafe(const char* name);
void printGNSSPacketCompact(const EpochPacket& packet);
void ClearINSBuffer();
void drainSerialInput(Stream &s, uint32_t drain_ms = 300);
void resetAllInputParserStates();
void drainAllSerialInputs(uint32_t drain_ms);
bool WaitForFreshRawGNSSAtOrAfterTime(double target_time, double& matchedGnssTime, uint32_t timeout_ms, double max_delay_sec);
void ClearPendingGNSS();


//===================================================//
//                POCKET SDR  (Simulated GNSS)                       //
//===================================================//
#define GNSSRAW_SERIAL Serial1
static uint8_t gnssrawRxBuffer[16384];
#define GNSSRAW_BAUDRATE 500000
#define GNSS_LINE_MAX 300
enum GNSSParserState {
  WAIT_TCA_HEADER,
  READ_TCA_ROWS
};

static GNSSParserState gnssState = WAIT_TCA_HEADER;

static EpochPacket tempGNSS;
static uint16_t expectedSats = 0;

double gnssLatencyEstimate = 0.0;
bool gnssLatencyValid = false;
//===================================================//
//                NEO7M-GPS                         //
//===================================================//
#define GNSS_SERIAL Serial4
#define GPS_BAUDRATE 9600
TinyGPSPlus gnss;
TinyGPSCustom gga_geoid_GPGGA(gnss, "GPGGA", 11);
TinyGPSCustom gga_geoid_GNGGA(gnss, "GNGGA", 11);


//===================================================//
//                IMU_SIMULATED (RLG specs)                          //
//===================================================//
#define IMU_SERIAL Serial2

static constexpr uint32_t IMU_BAUDRATE = 921600;  // Must match Python --baud

static constexpr size_t IMU_FRAME_LEN = 512;  // ASCII line buffer; 

// Additional UART memory for Teensy Serial2
static uint8_t imuRxBuffer[16384] = {0};

// Current frame parser state
static uint8_t imuFrame[IMU_FRAME_LEN] = {0};
static size_t imuFrameIndex = 0;
static double imuFrameStartTime = 0.0;

// ===================================================
// IMUSIM -> Teensy IMU receive diagnostics
// ===================================================

// When an oversized line is detected, ignore all remaining characters
// until its newline so the tail cannot be mistaken for a new packet.
static bool IMUsimDiscardOversizedLine = false;

//===================================================//
//                PPS Time Management                //
//===================================================//
#define PPS_PIN 22
int pps_count = 0;
uint32_t ppsCountAtSync = 0;
bool preciseTimeValid = false;
bool pps_seen = false;
//===================================================//
//                Precise Time Value                 //
//===================================================//
double starting_time = 0;
static unsigned long last_pps_time = 0;
uint32_t start_micros = 0;

//===================================================//
//                  UTC → GPS seconds                //
//===================================================//
#define GPS_UTC_LEAP 18

//===================================================//
//             IMUSIM-PARSING            //
//===================================================//

double lastTimeIMU = 0.0;
bool IMUSIM_Available = false;

// Throttle TCA output. INS still runs at IMU rate; only the serial output is limited.
static uint32_t lastTcaOutputMs = 0;


int y = 0;
int count = 0;


// -------------------------------- GNSS CONFIG --------------------------------
void InitGNSSConfig(GNSSConfig& GNSS_Config) {
  GNSS_Config.mask_angle = 15;           // Elevation Angle (Signals Below Than That Declare as an Outlier) (deg)
  GNSS_Config.mask_SignalStrength = 35;  // CNo (Signals Below than that Declare as an Outlier) (db-Hz)
  GNSS_Config.est_clock[0] = 0.0;                 // receiver clock offset (m)
  GNSS_Config.est_clock[1] = 0.0;                 // receiver clock drift (m/s)
}

// -------------------------------- TC KF CONFIG --------------------------------
void InitTCKFConfig(TC_KFConfig& TC_KF_Config) {
  // -------- State covariance (P-matrix) parameters --------

  TC_KF_Config.init_att_unc = 2.1 * DEG_TO_RAD;
  TC_KF_Config.init_vel_unc = 0.01;
  TC_KF_Config.init_pos_unc = 0.5;
  TC_KF_Config.init_b_a_unc = 30.0 * micro_g_to_meters_per_second_squared;
  TC_KF_Config.init_b_g_unc = 3.1 * DEG_TO_RAD / 3600.0;
  TC_KF_Config.init_clock_offset_unc = 1.1;
  TC_KF_Config.init_clock_drift_unc = 0.05;
  
  // -------- Propagation noise covariance (Q-matrix) parameters --------
   TC_KF_config.gyro_noise_PSD = 1.50e-7;   // rad^2/s
   TC_KF_config.accel_noise_PSD = 9.62e-6;  // m^2/s^3
   TC_KF_config.accel_bias_PSD = 6.0e-6;    // m^2/s^5, provisional
   TC_KF_config.gyro_bias_PSD = 2.12e-11;   // rad^2/s^3
   TC_KF_Config.clock_phase_PSD = 100.1;
   TC_KF_Config.clock_freq_PSD = 10.5;

  // -------- Other KF parameters --------
  TC_KF_Config.ZARU_DGrySD = 4.8481e-7;
  TC_KF_Config.StationaryDetectionHorizontalSpeed = 0.05;
  TC_KF_Config.GYRO_STATIC_THRESH = 0.02 * DEG_TO_RAD;
  TC_KF_Config.ACC_STATIC_THRESH = 0.08;
  TC_KF_Config.StationaryDetectionAccelerometerThreshold = 0.01;
  TC_KF_Config.StationaryDetectionHSTimeWindow = 0.2;
  TC_KF_Config.StationaryFlag = false;
}

// -------------------------------- MLS CONFIG --------------------------------
void InitMLSConfig(MLSConfig& MLS_Config) {
  MLS_Config.MAXloop = 0;
  MLS_Config.DeltaPosi = 0.0001;
  MLS_Config.K0 = 3.0;
  MLS_Config.K1 = 5.0;
}

void InitLLASIMConfig(LLASIMconfig &LLASIM_Config){
  LLASIM_Config.Lat_init = -121.915773;
  LLASIM_Config.Lon_init = 37.352721;
  LLASIM_Config.Alt_init = 20;
  LLASIM_Config.vn_init = 0;
  LLASIM_Config.ve_init = 0;
  LLASIM_Config.vd_init = 0;
}
// -------------------------------- CoG CONFIG --------------------------------///
void InitCoGConfig(CoGConfig& CoG_Config, double burn_start_time) {
  CoG_Config.m_s = 1.0;  // Structural mass (kg)

  CoG_Config.rs[0] = 0.0;  // Structure CoG x
  CoG_Config.rs[1] = 0.0;  // Structure CoG y
  CoG_Config.rs[2] = 0.0;  // Structure CoG z

  CoG_Config.R = 0.0;         // Solid propellant radius (m)
  CoG_Config.L_f0 = 0.0;      // Initial propellant length (m)
  CoG_Config.z_f0 = 0.0;      // Propellant start with nose reference (m)
  CoG_Config.x_offset = 0.0;  // Offset in X-axis
  CoG_Config.y_offset = 0.0;  // Offset in Y-axis
  CoG_Config.rho = 0.0;       // Propellant density (kg/m^3)
  CoG_Config.r_b = 0.0;       // Burn rate (m/s)
  CoG_Config.t_burn_start = burn_start_time;
}

void clearSerialInput(Stream &s, uint32_t quiet_ms, uint32_t max_ms)
{
  uint32_t start = millis();
  uint32_t lastByteTime = millis();

  while ((millis() - start) < max_ms) {
    while (s.available() > 0) {
      s.read();
      lastByteTime = millis();
    }

    if ((millis() - lastByteTime) >= quiet_ms) {
      break;  // no more pending burst data
    }
  }
}

void setup() {
  Serial.begin(921600);  // Teensy USB debug serial
  // ---------------------------------------------------------------------------
  // Pocket SDR raw-GNSS UART
  // ---------------------------------------------------------------------------
  GNSSRAW_SERIAL.addMemoryForRead(
      gnssrawRxBuffer,
      sizeof(gnssrawRxBuffer)
  );

  GNSSRAW_SERIAL.begin(GNSSRAW_BAUDRATE);
  GNSSRAW_SERIAL.setTimeout(5);

  Serial.println(F("Pocket SDR UART started"));

  // ---------------------------------------------------------------------------
  // NEO-7M PVT UART
  // ---------------------------------------------------------------------------
  GNSS_SERIAL.begin(GPS_BAUDRATE);
  GNSS_SERIAL.setTimeout(5);

  Serial.println(F("Initial NEO7M GPS UART started."));

  // ---------------------------------------------------------------------------
  // Python IMU UART
  // ---------------------------------------------------------------------------
  IMU_SERIAL.addMemoryForRead(
      imuRxBuffer,
      sizeof(imuRxBuffer)
  );

  IMU_SERIAL.begin(IMU_BAUDRATE);
  IMU_SERIAL.setTimeout(5);

  imuFrameIndex = 0;
  IMUsimDiscardOversizedLine = false;

  IMUSIM_Available = false;

  Serial.print(F("Python IMU UART started at "));
  Serial.print(IMU_BAUDRATE);
  Serial.println(F(" baud."));

  // ---------------------------------------------------------------------------
  // Drain stale queued data from every input port
  //
  // This removes old bytes already present in UART queues before initialization.
  // After this point, waitForInitialIMU(), waitForInitialGNSS(), and
  // waitForInitialRawGNSS() will use fresh packets only.
  // ---------------------------------------------------------------------------
  Serial.println(F("Draining queued serial data from all ports..."));

  drainAllSerialInputs(500);

  Serial.println(F("Serial input queues drained."));
  Serial.println(F("Waiting for fresh IMU, NEO7M, and Pocket SDR packets..."));

  // ---------------------------------------------------------------------------
  // Config initialization
  // ---------------------------------------------------------------------------
  packet.valid = false;

  InitGNSSConfig(GNSS_config);
  InitTCKFConfig(TC_KF_config);
  InitMLSConfig(MLS_config);

  // ---------------------------------------------------------------------------
  // PPS setup
  // ---------------------------------------------------------------------------
  pinMode(PPS_PIN, INPUT);

  last_pps_time = 0;
  pps_count = 0;
  pps_seen = false;
  preciseTimeValid = false;

  attachInterrupt(
      digitalPinToInterrupt(PPS_PIN),
      ppsISR,
      RISING
  );

  // ---------------------------------------------------------------------------
  // Wait for fresh IMU packet
  // ---------------------------------------------------------------------------
  const bool imuReady = waitForInitialIMU(1, 5000);

  if (!imuReady) {
    Serial.println(F("ERROR: Python IMU handshake failed."));
    Serial.println(F(
        "Check converter COM port, crossed TX/RX, common GND, and baud."
    ));
  } else {
    Serial.println(F(
        "Python IMU handshake completed and first fresh sample received."
    ));
  }

  // ---------------------------------------------------------------------------
  // Wait for fresh NEO7M PVT packet
  // ---------------------------------------------------------------------------
  bool gnssReady = waitForInitialGNSS(1, 5000);

  if (!gnssReady) {
    Serial.println(F("ERROR: Initial GNSS PVT packet not detected."));
  } else {
    Serial.println(F("Initial fresh GNSS PVT packet detected."));
  }

  // ---------------------------------------------------------------------------
  // Drain Pocket SDR queue once more before raw GNSS sync
  //
  // This avoids synchronizing time from an old raw-GNSS epoch that arrived
  // during IMU/NEO waiting.
  // ---------------------------------------------------------------------------
  clearSerialInput(GNSSRAW_SERIAL, 50, 500);
  resetGNSSParser();

  packet.valid = false;
  packet.gnssCount = 0;
  packet.gnss_obssec = 0.0;

  // ---------------------------------------------------------------------------
  // Wait for fresh Pocket SDR raw GNSS packet
  // ---------------------------------------------------------------------------
  bool rawGnssReady = waitForIntegerSecondPocketRawForSync(5000);

  if (!rawGnssReady) {
    Serial.println(F("ERROR: No integer-second Pocket GNSS Raw packet detected."));
  } else {
    Serial.println(F("Integer-second Pocket GNSS Raw packet detected."));

    if (!SyncPreciseTimeFromPocketRaw(packet.gnss_obssec)) {
      Serial.println(F("ERROR: Could not synchronize precise time from PocketRaw."));
      rawGnssReady = false;
    }
  }
  
  // ---------------------------------------------------------------------------
  // If any input source is missing, do not initialize TCA
  // ---------------------------------------------------------------------------
  if (!imuReady || !gnssReady || !rawGnssReady) {
    Serial.println(F(
        "TCA-INS initialization paused because required input data is missing."
    ));
    Serial.println(F(
        "UART parsers will continue running in loop()."
    ));

    
    return;
  }

  // ---------------------------------------------------------------------------
  // Initial alignment and PVT alignment
  // ---------------------------------------------------------------------------
  if (!Initial_Initializations()) {
    Serial.println(F("Failed to initialize TCA-INS algorithm."));

    
    return;
  }

  Serial.println(F("TCA-INS Algorithm initialized successfully."));

  // ---------------------------------------------------------------------------
  // IMPORTANT:
  // Do NOT overwrite start_micros here.
  //
  // start_micros should already be set inside SyncPreciseTimeFromPocketRaw()
  // to the PPS edge corresponding to starting_time.
  //
  // Do NOT use:
  // start_micros = micros();
  // ---------------------------------------------------------------------------

  InitCoGConfig(CoG_config, old_time);

  // ---------------------------------------------------------------------------
  // Clear INS circular buffer after successful alignment
  // ---------------------------------------------------------------------------
  ClearINSBuffer();

  // ---------------------------------------------------------------------------
  // Reset pending GNSS after clearing INS buffer
  //
  // Any GNSS epoch received during alignment is stale relative to the now-empty
  // INS buffer and cannot be interpolated.
  // ---------------------------------------------------------------------------
  GNSS_EPOCH_AVAILABLE = false;
  packet.valid = false;
  packet.gnssCount = 0;
  packet.gnss_obssec = 0.0;
  matchedGnssTime = 0.0;
  lastAcceptedGnssTime = 0.0;

  resetGNSSParser();

  // ---------------------------------------------------------------------------
  // Drain all serial input queues again so loop() starts from fresh data
  // ---------------------------------------------------------------------------
  Serial.println(F("Final serial drain before starting TCA loop..."));

  drainAllSerialInputs(500);

  // Final GNSS parser reset after drain
  resetGNSSParser();

  GNSS_EPOCH_AVAILABLE = false;
  packet.valid = false;
  packet.gnssCount = 0;
  packet.gnss_obssec = 0.0;
  matchedGnssTime = 0.0;
  lastAcceptedGnssTime = 0.0;
  readIMU();
  old_time = lastTimeIMU;
  Serial.println(F("Fresh-input TCA startup complete."));
}

void loop() {
  
  // -----------------------------------------------------------------------------
  // Raw GNSS receive/accept logic
  // Do not overwrite a pending GNSS epoch.
  // Accept only usable fresh epochs.
  // -----------------------------------------------------------------------------
  if (!GNSS_EPOCH_AVAILABLE) {

    if (GNSSRAW_SERIAL.available()) {

      if (readGNSSData()) {

        if (packet.valid &&
            packet.gnssCount > 0 &&
            packet.gnss_obssec > 0.0 &&
            isfinite(packet.gnss_obssec)) {

          double gnssTime = packet.gnss_obssec;

          // ---------------------------------------------------------------------
          // If INS buffer is not ready yet, discard this GNSS epoch.
          // ---------------------------------------------------------------------
          if (IMU_epoch_counter < 3 || latest_INS_store_idx < 0) {

            packet.valid = false;
            packet.gnssCount = 0;
            packet.gnss_obssec = 0.0;
            GNSS_EPOCH_AVAILABLE = false;

          } else {
            
              double tNow = old_time;
              double measuredLatency = tNow - gnssTime;
              Serial.print("measured latencyb/w GNSS sim and GNSS real = ");
              Serial.println(measuredLatency);
              if (!gnssLatencyValid) {
                gnssLatencyEstimate = measuredLatency;
                gnssLatencyValid = true;
              } else {
                // low-pass filter to avoid jitter
                gnssLatencyEstimate =
                    0.95 * gnssLatencyEstimate +
                    0.05 * measuredLatency;
              }
              double gnssTimeCorrected = gnssTime + gnssLatencyEstimate;
              matchedGnssTime = gnssTimeCorrected;
              GNSS_EPOCH_AVAILABLE = true;
              lastAcceptedGnssTime = matchedGnssTime;
            
          }
        }
      }
    }
  }
  
  if (readIMU()) {
  
    ////-----------------------------------------------------------------------INERTIAL NAVIGATION SYSTEM---------------------------------------///////
    // Split IMU Time, Accelerometer & Gyroscope Values
    double IMU_time = lastTimeIMU;

    double meas_f_ib_b_raw[3] = { rlgData.ax, rlgData.ay, rlgData.az };
    double meas_omega_ib_b_raw[3] = { (rlgData.gx * DEG_TO_RAD), (rlgData.gy * DEG_TO_RAD), (rlgData.gz * DEG_TO_RAD) };
    // Serial.print("Ax= ");
    // Serial.print(meas_f_ib_b_raw[0],7);
    // Serial.print("Ay= ");
    // Serial.print(meas_f_ib_b_raw[1],7);
    // Serial.print("Az= ");
    // Serial.print(meas_f_ib_b_raw[2],7);
    // Avoid 200 Hz USB debug printing; it can block the loop and create
    // IMU/GNSS replay gaps. Re-enable only for short diagnostics.
    // Serial.print("Gx= ");
    // Serial.print(rlgData.gx,7);
    // Serial.print("Gy= ");
    // Serial.print(rlgData.gy,7);
    // Serial.print("Gz= ");
    // Serial.println(rlgData.gz,7);
    //Time Interval Between Two IMU update
    double tor_i = IMU_time - old_time;
    Serial.print(" tor_i= ");
      Serial.println(tor_i,4);

    // 200 Hz IMU should give about 0.005 s. If the loop was blocked or
    // an input gap occurred, skip this epoch instead of propagating with
    // an unrealistic IMU integration interval.
    if (tor_i <= 0.0 || tor_i > 0.02) {
      Serial.print(F("WARNING: bad tor_i = "));
      Serial.print(tor_i, 7);
      Serial.println(F(" -> skipping this IMU epoch."));
      old_time = IMU_time;
      return;
    }

    GNSS_config.est_clock[0] = GNSS_config.est_clock[0] + GNSS_config.est_clock[1] * tor_i;
    //---------------------IMU BIAS CORRECTION - GYROSCOPE-----------------------------///
    // Direct gyro bias correction only.
    // No simulator Earth-rate / transport-rate conversion is applied here.
    double prev_w_ib_b_corr[3] = {
      prev_w_ib_b_raw[0] - est_IMU_Bias[3],
      prev_w_ib_b_raw[1] - est_IMU_Bias[4],
      prev_w_ib_b_raw[2] - est_IMU_Bias[5]
    };

    double meas_omega_ib_b[3] = {
      meas_omega_ib_b_raw[0] - est_IMU_Bias[3],
      meas_omega_ib_b_raw[1] - est_IMU_Bias[4],
      meas_omega_ib_b_raw[2] - est_IMU_Bias[5]
    };
    //---------------------IMU BIAS CORRECTION - GYROSCOPE-----------------------------///


    //---------------------IMU BIAS CORRECTION - ACCELEROMETER-----------------------------///
    //-------------------- Current Lever Arm Between IMU & Test Structure CoG -------------------///
    if (!Compute_CoG_at_Time(IMU_time, CoG_config, L_imu_b, L_imu_b_current)) {
      Serial.println("Compute_COG_at_Time failed");
      while (1) {}
    }

    //-------------------- Lever Arm Produces Centripetal & Tangential Acceleration -------------------///
    double prev_alpha[3] = { alpha[0], alpha[1], alpha[2] };

    // Define & Set alpha_inst
    double alpha_inst[3]={0.0};
    for (int i = 0; i < 3; i++) {
      alpha_inst[i] = (meas_omega_ib_b[i] - prev_w_ib_b_corr[i]) / tor_i;
    }

    //Update alpha Values
    for (int i = 0; i < 3; i++) {
      alpha[i] = 0.5 * (alpha_inst[i] + prev_alpha[i]);
    }

    double omega_mid[3] = { 0.5 * (meas_omega_ib_b[0] + prev_w_ib_b_corr[0]),
                            0.5 * (meas_omega_ib_b[1] + prev_w_ib_b_corr[1]),
                            0.5 * (meas_omega_ib_b[2] + prev_w_ib_b_corr[2]) };

    double centrip[3] = { 0.0, 0.0, 0.0 };
    // outer cross: cross(omega_mid, omega_cross_L)
    centrip[0] = omega_mid[1] * (omega_mid[0] * L_imu_b_current[1] - omega_mid[1] * L_imu_b_current[0]) - omega_mid[2] * (omega_mid[2] * L_imu_b_current[0] - omega_mid[0] * L_imu_b_current[2]);
    centrip[1] = omega_mid[2] * (omega_mid[1] * L_imu_b_current[2] - omega_mid[2] * L_imu_b_current[1]) - omega_mid[0] * (omega_mid[0] * L_imu_b_current[1] - omega_mid[1] * L_imu_b_current[0]);
    centrip[2] = omega_mid[0] * (omega_mid[2] * L_imu_b_current[0] - omega_mid[0] * L_imu_b_current[2]) - omega_mid[1] * (omega_mid[1] * L_imu_b_current[2] - omega_mid[2] * L_imu_b_current[1]);

    double tang[3] = { 0.0, 0.0, 0.0 };
    // tang = cross(alpha, L_imu_b_current)
    tang[0] = alpha[1] * L_imu_b_current[2] - alpha[2] * L_imu_b_current[1];
    tang[1] = alpha[2] * L_imu_b_current[0] - alpha[0] * L_imu_b_current[2];
    tang[2] = alpha[0] * L_imu_b_current[1] - alpha[1] * L_imu_b_current[0];

    //--------------------- Applying Bias and lever arm correction --------------------///
    double meas_f_ib_b[3] = { meas_f_ib_b_raw[0] - est_IMU_Bias[0] - centrip[0] - tang[0],
                              meas_f_ib_b_raw[1] - est_IMU_Bias[1] - centrip[1] - tang[1],
                              meas_f_ib_b_raw[2] - est_IMU_Bias[2] - centrip[2] - tang[2] };

    // Serial.print("corrected Ax= ");
    // Serial.print(meas_f_ib_b[0],7);
    // Serial.print("corrected Ay= ");
    // Serial.print(meas_f_ib_b[1],7);
    // Serial.print("corrected Az= ");
    // Serial.print(meas_f_ib_b[2],7);
    // ///---------------------------------------------------UPADTE INS PVA (NAVIGATION) SOLUTION IN ECEF FRAME-------------------------------------////////
    if (!Nav_Equations_ECEF(tor_i, old_est_r_eb_e, old_est_v_eb_e, old_est_C_b_e, meas_f_ib_b, meas_omega_ib_b, est_r_eb_e, est_v_eb_e, est_C_b_e)) {
      Serial.println("Nav_Equations_ECEF Failed");
      while (1) {}
    }
    // Serial.print(est_r_eb_e[0]);
    // Serial.print(", ");
    // Serial.print(est_r_eb_e[1]);
    // Serial.print(", ");
    // Serial.print(est_r_eb_e[2]);
    // Serial.print(", ");
    // Serial.print(est_v_eb_e[0]);
    // Serial.print(", ");
    // Serial.print(est_v_eb_e[1]);
    // Serial.print(", ");
    // Serial.println(est_v_eb_e[2]);

       // KF prediction at IMU rate
      if (!MLS_KF_predict(tor_i,
                          old_est_C_b_e,
                          old_est_v_eb_e,
                          old_est_r_eb_e,
                          meas_f_ib_b,
                          TC_KF_config,
                          meas_omega_ib_b,
                          MLS_config,
                          P_matrix)) {
        Serial.println(F("Kalman prediction failed"));
        return;
      }

      IMU_epoch_counter++;
      int store_idx = (IMU_epoch_counter - 1) % INS_BUFFER_SIZE;

      INS_time_store[store_idx] = IMU_time;

      for (int i = 0; i < 3; i++) {
        INS_r_store[store_idx][i] = est_r_eb_e[i];
        INS_v_store[store_idx][i] = est_v_eb_e[i];
        INS_f_store[store_idx][i] = meas_f_ib_b[i];
        INS_w_store[store_idx][i] = meas_omega_ib_b[i];
        INS_f_raw_store[store_idx][i] = meas_f_ib_b_raw[i];
        INS_w_raw_store[store_idx][i] = meas_omega_ib_b_raw[i];
      }
        
      //3D Matrix
      for (int i = 0; i < 9; i++) {
        INS_C_store[store_idx][i] = est_C_b_e[i];
      }

      
      for (size_t i = 0; i < P_DIM; i++) {
        INS_P_store[store_idx][i] = P_matrix[i];
      }

      for (size_t i = 0; i < BIAS_DIM; i++) {
        INS_bias_store[store_idx][i] = est_IMU_Bias[i];
      }

      for (size_t i = 0; i < CLOCK_DIM; i++) {
        INS_clock_store[store_idx][i] = GNSS_config.est_clock[i];
      }

      latest_INS_store_idx = store_idx;
      INS_valid_store[store_idx] = true;
    //*=====> Robust ZARU / stationary detection.
    //  This detector requires small
    // acceleration-norm error, small gyro norm, and low horizontal speed.
    const double g0_static = 9.80665;
    
    const double acc_norm_static = sqrt(
        meas_f_ib_b_raw[0] * meas_f_ib_b_raw[0] +
        meas_f_ib_b_raw[1] * meas_f_ib_b_raw[1] +
        meas_f_ib_b_raw[2] * meas_f_ib_b_raw[2]
    );

    const double gyro_norm_static = sqrt(
        meas_omega_ib_b[0] * meas_omega_ib_b[0] +
        meas_omega_ib_b[1] * meas_omega_ib_b[1] +
        meas_omega_ib_b[2] * meas_omega_ib_b[2]
    );

    const double horizontal_speed_static = sqrt(
        old_est_v_eb_n[0] * old_est_v_eb_n[0] +
        old_est_v_eb_n[1] * old_est_v_eb_n[1]
    );

    const bool imuLooksStatic =
        (fabs(acc_norm_static - g0_static) < TC_KF_config.ACC_STATIC_THRESH) &&
        (gyro_norm_static < TC_KF_config.GYRO_STATIC_THRESH) &&
        (horizontal_speed_static < TC_KF_config.StationaryDetectionHorizontalSpeed);

    if (imuLooksStatic) {
      stationaryperiod += tor_i;
      TC_KF_config.StationaryFlag =
          (stationaryperiod > TC_KF_config.StationaryDetectionHSTimeWindow);
    } else {
      TC_KF_config.StationaryFlag = false;
      stationaryperiod = 0.0;
    }
    ////-----------------Robust ZARU / stationary detection---------------------------/////
    ////-----------------------------------------------------------------------INERTIAL NAVIGATION SYSTEM---------------------------------------///////


    ///-----------------------------------------------------------------------TIGHTLY COUPLED ALGORITHM---------------------------------------///////
    ///---------------------------------Update GNSS Observation & Run Kalman Filter-----------------------------------/////////////
    if (GNSS_EPOCH_AVAILABLE) {
      
      if (matchedGnssTime <= 0.0) {
      Serial.println(F("ERROR: GNSS_EPOCH_AVAILABLE true but matchedGnssTime is zero."));
      GNSS_EPOCH_AVAILABLE = false;
      old_time = IMU_time;
      return;
      }

          ///------------------------------Retrieve Matched INS Data For Coupling----------------------------------------/////////////
      int idx2 = -1;

      int idx1 = -1;
        
      if (!FindINSIntervalForGNSSTime(matchedGnssTime, idx1, idx2)) {

        if (lastINSIntervalStatus == INS_INTERVAL_GNSS_AHEAD_OF_INS) {
          // GNSS is newer than latest INS.
          // Keep it pending and wait for next IMU sample.
          // Do not clear GNSS.
          // Do not treat this as an error.
          goto FINISH_IMU_STEP;
        }

        // GNSS is stale, invalid, or no valid INS bracket exists.
        // Discard it.
        ClearPendingGNSS();

        goto FINISH_IMU_STEP;
      }
      
      double t1 = INS_time_store[idx1];
      double t2 = INS_time_store[idx2];
      double dt = t2 - t1; //two successive INS samples time interval

      double tor_g = matchedGnssTime - t1;  //Where matchedGnssTime = GNSS_time
      // Serial.print("tor_g= ");
      // Serial.println(tor_g,4);
      if (tor_g < 0.0 || tor_g > dt) {
        Serial.println("GNSS time is outside current INS interval.");
        
        GNSS_EPOCH_AVAILABLE = false;
        old_time = IMU_time;
        return;
      }       
        
      double r1[3]; //INS buffered position
      double v1[3]; //INS buffered velocity
      double f2[3]; //INS buffered lever arm compensated acceleration
      double w2[3]; //INS buffered lever arm compensated angular velocity
      double C1[9]; //INS bufferedbody to ecef transformation matrix
      
      double P1[P_DIM]; //Kalman predicted state  covariance buffered 
      double bias1[BIAS_DIM]; //IMU bias buffered
      double clock1[CLOCK_DIM];

      double f2_raw[3]; //INS buffered raw accelerration
      double w2_raw[3]; //INS buffered raw angular velcoity
    
      //newly added INS buffer time synchronization
      for (int i = 0; i < 3; i++) {
        r1[i] = INS_r_store[idx1][i];
        v1[i] = INS_v_store[idx1][i];

        f2_raw[i] = INS_f_raw_store[idx2][i];
        w2_raw[i] = INS_w_raw_store[idx2][i];
      }

      for (int i = 0; i < 9; i++) {
        C1[i] = INS_C_store[idx1][i];
      }

      for (size_t i = 0; i < P_DIM; i++) {
        P1[i] = INS_P_store[idx1][i];
      }

      for (size_t i = 0; i < BIAS_DIM; i++) {
        bias1[i] = INS_bias_store[idx1][i];
      }

      for (size_t i = 0; i < CLOCK_DIM; i++) {
        clock1[i] = INS_clock_store[idx1][i];
      }

       for (int i = 0; i < 3; i++) {
          f2[i] = f2_raw[i] - bias1[i];
          w2[i] = w2_raw[i] - bias1[i + 3];
        }

        // Restore covariance, bias, clock, and Doppler-bias state at t1
        for (size_t i = 0; i < P_DIM; i++) {
          P_matrix[i] = P1[i];
        }

        for (size_t i = 0; i < BIAS_DIM; i++) {
          est_IMU_Bias[i] = bias1[i];
        }

        for (size_t i = 0; i < CLOCK_DIM; i++) {
          GNSS_config.est_clock[i] = clock1[i];
        }
        GNSS_config.est_clock[0] = GNSS_config.est_clock[0] + GNSS_config.est_clock[1] * tor_g;
        double r_interpt[3] = {0.0, 0.0, 0.0};
        double v_interpt[3] = {0.0, 0.0, 0.0};
        double C_interpt[9] = {0.0};
        if (tor_g <= 1e-9) {

          for (int i = 0; i < 3; i++) {
            r_interpt[i] = r1[i];
            v_interpt[i] = v1[i];
          }

          for (int i = 0; i < 9; i++) {
            C_interpt[i] = C1[i];
          }

        } else {

          if (!Nav_Equations_ECEF(tor_g,
                                  r1,
                                  v1,
                                  C1,
                                  f2,
                                  w2,
                                  r_interpt,
                                  v_interpt,
                                  C_interpt)) {
            Serial.println(F("Nav_Equations_ECEF Failed during GNSS rollback."));
            GNSS_EPOCH_AVAILABLE = false;
            return;
          }

          // Propagate covariance from t1 to exact GNSS time.
          // Use the interval-start state C1/v1/r1 for the covariance linearization.
          if (!MLS_KF_predict(tor_g,
                              C1,
                              v1,
                              r1,
                              f2,
                              TC_KF_config,
                              w2,
                              MLS_config,
                              P_matrix)) {
            Serial.println(F("MLS_KF_predict failed during GNSS rollback."));
            GNSS_EPOCH_AVAILABLE = false;
            return;
          }
        }

      
      ///---------------------------------------------Update GNSS Observation---------------------//////
  
          ////-----------------------------Determine Whether Enough GNSS Signal Available -------------------------------///////////
          if (packet.gnssCount >= 3) {
            ///---------------------------Update GNSS Observation Based on Available Satellites For TCA--------------------------///////
              time_last_GNSS = matchedGnssTime;  //matchedGnssTime = GNSS_time
            ///-----AT THIS STAGE WE UPDATE PRRATE_CORR VALUE BASED ON THIS-----------------//
            ///-----MATLAB SYNTAX "rangerate + rateofsatelliteclock * c"--------------------//
            ///-----TO AVOID CREATING NEW VECTOR (I.E. GNSS_MEASUREMENTS)-------------------//

            // Serial.println(F("========== GNSS RAW PACKET AFter picksubset=========="));
            // printGNSSPacketCompact(packet);
            // Serial.println(F("====================================="));
            // GNSS_EPOCH_AVAILABLE = false;
            for (size_t i = 0; i < packet.gnssCount; ++i) {
              double temp_value = packet.gnssRows[i].prrate_corr + (packet.gnssRows[i].sat_clk_rate * 299792458.0);  // Where Constant value "299792458.0" is value c that is Speed of Light
              packet.gnssRows[i].prrate_corr = 0;
              packet.gnssRows[i].prrate_corr = temp_value;
            }
            
            ///---------------------------------Measuerment Noise Covariance 'R' Matrix Calculation----------------------------------////

              // ------------------------------------------------------------
              // Declaration of R-Matrix
              // MATLAB:
              // R_matrix = eye(no_GNSS_meas * 2);
              // if TC_KF_config.StationarityFlag == 1
              //     R_matrix = eye(no_GNSS_meas * 2 + 1);
              // end
              // ------------------------------------------------------------

              size_t R_dim = 2 * packet.gnssCount;

              if (TC_KF_config.StationaryFlag) {
                R_dim += 1;
              }

              if (R_dim > R_ROWS_MAX) {
                Serial.println(F("ERROR: R_dim exceeds R_ROWS_MAX. Discarding epoch."));

                packet.valid = false;
                packet.gnssCount = 0;
                GNSS_EPOCH_AVAILABLE = false;
                return;
              }

              // ------------------------------------------------------------
              // Clear R matrix
              // Column-major indexing: R_Matrix[row + R_dim * col]
              // ------------------------------------------------------------

              for (size_t col = 0; col < R_dim; col++) {
                for (size_t row = 0; row < R_dim; row++) {
                  R_Matrix[row + R_dim * col] = 0.0;
                }
              }

              // ------------------------------------------------------------
              // Initialize as identity
              // ------------------------------------------------------------

              for (size_t i = 0; i < R_dim; i++) {
                R_Matrix[i + R_dim * i] = 1.0;
              }

              // ------------------------------------------------------------
              // Urban pseudorange weighting parameters
              // MATLAB:
              // a = 1e9;
              // b = 20;
              // ------------------------------------------------------------

              const double a = 0.08;
              const double b = 30.0;

              // ------------------------------------------------------------
              // Preallocate variance vectors
              // MATLAB:
              // sigma2_vec_pr  = zeros(1,no_GNSS_meas);
              // sigma2_vec_prr = zeros(1,no_GNSS_meas);
              // ------------------------------------------------------------

              double sigma2_vec_pr[MAX_GNSS_ROWS_PER_EPOCH]  = {0.0};
              double sigma2_vec_prr[MAX_GNSS_ROWS_PER_EPOCH] = {0.0};

              // ------------------------------------------------------------
              // Calculate pseudorange and pseudorange-rate variances
              // ------------------------------------------------------------

              for (size_t k = 0; k < packet.gnssCount; ++k) {

                const double elev_deg = packet.gnssRows[k].el_deg;
                const double CNo      = packet.gnssRows[k].cn0;

                // ============================================================
                // Pseudorange variance
                // MATLAB:
                // if elev_deg(k) < GNSS_config.mask_angle ||
                //    CNo(k) < GNSS_config.mask_SignalStrenth
                //
                //     sigma2_vec_pr(k)  = 1e8;
                //     sigma2_vec_prr(k) = 1e8;
                //
                // else
                //     sin_elev = sind(elev_deg(k));
                //     if sin_elev < 1e-3
                //         sin_elev = 1e-3;
                //     end
                //
                //     sigma2_vec_pr(k) =
                //       (a + b * 10^(-CNo(k)/10.0)) / sin_elev;
                // end
                // ============================================================

                if ((elev_deg < GNSS_config.mask_angle) ||
                    (CNo < GNSS_config.mask_SignalStrength)) {

                  sigma2_vec_pr[k]  = 1e8;
                  sigma2_vec_prr[k] = 1e8;

                } else {

                  double sin_elev = sin(elev_deg * DEG_TO_RAD);

                  if (sin_elev < 1e-3) {
                    sin_elev = 1e-3;
                  }

                  sigma2_vec_pr[k] =
                      (a + b * pow(10.0, -CNo / 10.0)) / sin_elev;
                }

                // ============================================================
                // Pseudorange-rate variance
                //
                // MATLAB:
                // sigma_prr = 0.02;
                // sin_elev_prr = sind(max(elev_deg(k), GNSS_config.mask_angle));
                // sigma2_vec_prr(k) = sigma_prr^2 / sin_elev_prr;
                //
                // NOTE:
                // This follows MATLAB exactly and overwrites sigma2_vec_prr(k),
                // even if the satellite failed elevation/CN0 mask above.
                // ============================================================

                const double sigma_prr = 0.02;   // m/s

                double elev_for_prr = elev_deg;

                if (elev_for_prr < GNSS_config.mask_angle) {
                  elev_for_prr = GNSS_config.mask_angle;
                }

                double sin_elev_prr = sin(elev_for_prr * DEG_TO_RAD);

                if (sin_elev_prr < 1e-3) {
                  sin_elev_prr = 1e-3;
                }

                sigma2_vec_prr[k] =
                    (sigma_prr * sigma_prr) / sin_elev_prr;
              }

              // ------------------------------------------------------------
              // Fill R matrix
              // MATLAB:
              // R_matrix(1:no_GNSS_meas,1:no_GNSS_meas)
              //     = diag(sigma2_vec_pr);
              //
              // R_matrix(no_GNSS_meas+1:2*no_GNSS_meas,
              //          no_GNSS_meas+1:2*no_GNSS_meas)
              //     = diag(sigma2_vec_prr);
              // ------------------------------------------------------------

              // Pseudorange block
              for (size_t k = 0; k < packet.gnssCount; ++k) {
                R_Matrix[k + R_dim * k] = sigma2_vec_pr[k];
              }

              // Pseudorange-rate block
              for (size_t k = 0; k < packet.gnssCount; ++k) {
                size_t idx = packet.gnssCount + k;
                R_Matrix[idx + R_dim * idx] = sigma2_vec_prr[k];
              }

              // ------------------------------------------------------------
              // ZARU measurement if stationary
              // MATLAB:
              // R_matrix(2*no_GNSS_meas+1,2*no_GNSS_meas+1)
              //     = TC_KF_config.ZARU_DGrySD^2;
              // ------------------------------------------------------------

              if (TC_KF_config.StationaryFlag) {
                size_t idx = 2 * packet.gnssCount;
                R_Matrix[idx + R_dim * idx] =
                    TC_KF_config.ZARU_DGrySD * TC_KF_config.ZARU_DGrySD;
              }

              ///-------------------------------Measuerment Noise Covariance R-Matrix Calculation---------------------------------------///////////////

            // ///---------------------------------------------------Kalman Filter Running-------------------------------------------------/////////////
            bool MLS_LoopFlag = true;
            int MLS_LoopCount = 0;

            while (MLS_LoopFlag) {
             if (!MLS_KF_update(packet.gnssRows,
                                packet.gnssCount,
                                C_interpt,
                                v_interpt,
                                r_interpt,
                                meas_f_ib_b,
                                TC_KF_config,
                                L_ba_b,
                                meas_omega_ib_b,
                                MLS_config,
                                P_matrix,
                                R_Matrix,
                                R_Matrix_NextLoop,
                                R_dim,
                                GNSS_config,
                                est_C_b_e,
                                est_v_eb_e,
                                est_r_eb_e,
                                est_IMU_Bias)) {
                Serial.println(F("Kalman update failed"));
                GNSS_EPOCH_AVAILABLE = false;
                return;
              }
              
              MLS_LoopCount = MLS_LoopCount + 1;

              const double dx0 = est_r_eb_e_temp_old[0] - est_r_eb_e[0];
              const double dx1 = est_r_eb_e_temp_old[1] - est_r_eb_e[1];
              const double dx2 = est_r_eb_e_temp_old[2] - est_r_eb_e[2];

              const double pos_diff_norm = sqrt(dx0 * dx0 + dx1 * dx1 + dx2 * dx2);

              if ((pos_diff_norm < MLS_config.DeltaPosi) || (MLS_LoopCount >= MLS_config.MAXloop)) {
                MLS_LoopFlag = false;
                
              } else {
                est_r_eb_e_temp_old[0] = r_interpt[0];
                est_r_eb_e_temp_old[1] = r_interpt[1];
                est_r_eb_e_temp_old[2] = r_interpt[2];

                memcpy(R_Matrix, R_Matrix_NextLoop, R_dim * R_dim * sizeof(double));
              }
              // MATLAB: est_r_eb_e_temp_old = est_r_eb_e_temp;
                for (size_t i = 0; i < 3; ++i) {
                  est_r_eb_e_temp_old[i] = est_r_eb_e[i];
                }

              
            }

            

            ///--------------------------------------------REPROPAGATE FROM GNSS TIME TO CURRENT IMU TIME----------------------------////

            double tor_rem = IMU_time - matchedGnssTime;
              // Serial.print("tor_rem= ");
              // Serial.println(tor_rem);
              
            if (tor_rem > 0.0) {

              // After MLS_KF_update(), solution is at GNSS time.
              double r_cur[3];
              double v_cur[3];
              double C_cur[9];

              static double P_cur[P_DIM];

              double bias_cur[BIAS_DIM];
              double clock_cur[CLOCK_DIM];

              for (int i = 0; i < 3; i++) {
                r_cur[i] = est_r_eb_e[i];
                v_cur[i] = est_v_eb_e[i];
              }

              for (int i = 0; i < 9; i++) {
                C_cur[i] = est_C_b_e[i];
              }

              for (size_t i = 0; i < P_DIM; i++) {
                P_cur[i] = P_matrix[i];
              }

              for (size_t i = 0; i < BIAS_DIM; i++) {
                bias_cur[i] = est_IMU_Bias[i];
              }

              for (size_t i = 0; i < CLOCK_DIM; i++) {
                clock_cur[i] = GNSS_config.est_clock[i];
              }

              double last_time = matchedGnssTime;

              double prev_w_raw_replay[3];

              for (int i = 0; i < 3; i++) {
                prev_w_raw_replay[i] = INS_w_raw_store[idx2][i];
              }

              // Replay buffered IMU samples after GNSS time up to current IMU time.
              // Start from idx2 because idx2 is the first INS sample at/after GNSS time.
              size_t replay_idx = idx2;

              int nStored = (IMU_epoch_counter < INS_BUFFER_SIZE)
                              ? IMU_epoch_counter
                              : INS_BUFFER_SIZE;
              
              for (int replay_count = 0; replay_count < nStored; replay_count++) {
                
                if (replay_idx < 0 || replay_idx >= INS_BUFFER_SIZE) {
                  Serial.println(F("ERROR: invalid replay_idx during INS repropagation."));
                  break;
                }
                
                if (!INS_valid_store[replay_idx]) {
                  replay_idx = (replay_idx + 1) % INS_BUFFER_SIZE;
                  continue;
                }
                double this_time = INS_time_store[replay_idx];

                // Replay only samples after GNSS time and up to current IMU time.
                if (this_time > matchedGnssTime + 1e-9 && this_time <= IMU_time + 1e-9) {

                  double dt_replay = this_time - last_time;
                  // Serial.print(F("dt_replay = "));
                  //   Serial.println(dt_replay, 7);
                  if (dt_replay <= 0.0) {
                    replay_idx = (replay_idx + 1) % INS_BUFFER_SIZE;
                    continue;
                  }

                  // Safety check for 200 Hz IMU.
                  // Normal dt_replay should be around 0.005 s.
                  if (dt_replay > 0.05) {
                    Serial.print(F("WARNING: large dt_replay = "));
                    Serial.print(dt_replay, 7);
                    // Serial.println(F(" -> discarding this GNSS epoch."));
                    // ClearPendingGNSS();
                    goto FINISH_IMU_STEP;
                  }

                  double f_raw[3];
                  double w_raw[3];

                  for (int i = 0; i < 3; i++) {
                    f_raw[i] = INS_f_raw_store[replay_idx][i];
                    w_raw[i] = INS_w_raw_store[replay_idx][i];
                  }

                  // MATLAB:
                  // prev_w_corr = prev_w_raw_replay - bias_cur(4:6);
                  // w_corr      = w_raw - bias_cur(4:6);
                  // Direct gyro bias correction only. No simulator Earth-rate /
                  // transport-rate conversion is applied in replay.
                  double prev_w_corr[3];
                  double w_corr[3];

                  for (int i = 0; i < 3; i++) {
                    prev_w_corr[i] = prev_w_raw_replay[i] - bias_cur[i + 3];
                    w_corr[i]      = w_raw[i]              - bias_cur[i + 3];
                  }

                  // MATLAB:
                  // L_bi_b_current_replay = compute_CoG_at_time(...)
                  //
                  // If you have time-varying CoG, replace this assignment with your
                  // compute_CoG_at_time() equivalent.
                  double L_imu_b_current_replay[3];

                  for (int i = 0; i < 3; i++) {
                    L_imu_b_current_replay[i] = L_imu_b_current[i];
                  }

                  // MATLAB:
                  // alpha_replay = (w_corr - prev_w_corr) / dt_replay;
                  // omega_mid_replay = 0.5 * (w_corr + prev_w_corr);
                  double alpha_replay[3];
                  double omega_mid_replay[3];

                  for (int i = 0; i < 3; i++) {
                    alpha_replay[i] = (w_corr[i] - prev_w_corr[i]) / dt_replay;
                    omega_mid_replay[i] = 0.5 * (w_corr[i] + prev_w_corr[i]);
                  }

                  // centrip_replay = cross(omega_mid_replay, cross(omega_mid_replay, L_imu_b_current_replay))
                  double temp_cross[3];
                  double centrip_replay[3];
                  double tang_replay[3];
                  double f_corr[3];

                  temp_cross[0] = omega_mid_replay[1] * L_imu_b_current_replay[2]
                                - omega_mid_replay[2] * L_imu_b_current_replay[1];

                  temp_cross[1] = omega_mid_replay[2] * L_imu_b_current_replay[0]
                                - omega_mid_replay[0] * L_imu_b_current_replay[2];

                  temp_cross[2] = omega_mid_replay[0] * L_imu_b_current_replay[1]
                                - omega_mid_replay[1] * L_imu_b_current_replay[0];

                  centrip_replay[0] = omega_mid_replay[1] * temp_cross[2]
                                    - omega_mid_replay[2] * temp_cross[1];

                  centrip_replay[1] = omega_mid_replay[2] * temp_cross[0]
                                    - omega_mid_replay[0] * temp_cross[2];

                  centrip_replay[2] = omega_mid_replay[0] * temp_cross[1]
                                    - omega_mid_replay[1] * temp_cross[0];

                  // tang_replay = cross(alpha_replay, L_imu_b_current_replay)
                  tang_replay[0] = alpha_replay[1] * L_imu_b_current_replay[2]
                                - alpha_replay[2] * L_imu_b_current_replay[1];

                  tang_replay[1] = alpha_replay[2] * L_imu_b_current_replay[0]
                                - alpha_replay[0] * L_imu_b_current_replay[2];

                  tang_replay[2] = alpha_replay[0] * L_imu_b_current_replay[1]
                                - alpha_replay[1] * L_imu_b_current_replay[0];

                  // MATLAB:
                  // f_corr = f_raw - bias_cur(1:3) - centrip_replay - tang_replay;
                  for (int i = 0; i < 3; i++) {
                    f_corr[i] = f_raw[i] - bias_cur[i] - centrip_replay[i] - tang_replay[i];
                  }

                  // Propagate navigation state.
                  // Save interval-start state for MLS_KF_predict().
                  double r_replay_old[3];
                  double v_replay_old[3];
                  double C_replay_old[9];

                  for (int i = 0; i < 3; i++) {
                    r_replay_old[i] = r_cur[i];
                    v_replay_old[i] = v_cur[i];
                  }

                  for (int i = 0; i < 9; i++) {
                    C_replay_old[i] = C_cur[i];
                  }

                  double r_next[3];
                  double v_next[3];
                  double C_next[9];

                  if (!Nav_Equations_ECEF(dt_replay,
                                          r_replay_old,
                                          v_replay_old,
                                          C_replay_old,
                                          f_corr,
                                          w_corr,
                                          r_next,
                                          v_next,
                                          C_next)) {
                    Serial.println(F("Nav_Equations_ECEF Failed during INS replay."));
                    GNSS_EPOCH_AVAILABLE = false;
                    return;
                  }

                  for (int i = 0; i < 3; i++) {
                    r_cur[i] = r_next[i];
                    v_cur[i] = v_next[i];
                  }

                  for (int i = 0; i < 9; i++) {
                    C_cur[i] = C_next[i];
                  }

                  // Propagate covariance using interval-start state.
                  if (!MLS_KF_predict(dt_replay,
                                      C_replay_old,
                                      v_replay_old,
                                      r_replay_old,
                                      f_corr,
                                      TC_KF_config,
                                      w_corr,
                                      MLS_config,
                                      P_cur)) {
                    Serial.println(F("MLS_KF_predict Failed during INS replay."));
                    GNSS_EPOCH_AVAILABLE = false;
                    return;
                  }

                  // Update corrected buffer entries
                  for (int i = 0; i < 3; i++) {
                    INS_r_store[replay_idx][i] = r_cur[i];
                    INS_v_store[replay_idx][i] = v_cur[i];

                    INS_f_store[replay_idx][i] = f_corr[i];
                    INS_w_store[replay_idx][i] = w_corr[i];
                  }
                      
                  for (int i = 0; i < 9; i++) {
                    INS_C_store[replay_idx][i] = C_cur[i];
                  }

                  for (size_t i = 0; i < P_DIM; i++) {
                    INS_P_store[replay_idx][i] = P_cur[i];
                  }

                  for (size_t i = 0; i < BIAS_DIM; i++) {
                    INS_bias_store[replay_idx][i] = bias_cur[i];
                  }

                  for (size_t i = 0; i < CLOCK_DIM; i++) {
                    INS_clock_store[replay_idx][i] = clock_cur[i];
                  }

                  // Prepare next replay step
                  for (int i = 0; i < 3; i++) {
                    prev_w_raw_replay[i] = w_raw[i];
                  }

                  last_time = this_time;
                }

                if (replay_idx == latest_INS_store_idx) {
                  break;
                }

                replay_idx = (replay_idx + 1) % INS_BUFFER_SIZE;
              }
                
              // Assign corrected current solution
              for (int i = 0; i < 3; i++) {
                est_r_eb_e[i] = r_cur[i];
                est_v_eb_e[i] = v_cur[i];

                est_IMU_Bias[i]     = bias_cur[i];
                est_IMU_Bias[i + 3] = bias_cur[i + 3];
              }

              for (int i = 0; i < 9; i++) {
                est_C_b_e[i] = C_cur[i];
              }

              for (size_t i = 0; i < P_DIM; i++) {
                P_matrix[i] = P_cur[i];
              }

              for (size_t i = 0; i < CLOCK_DIM; i++) {
                GNSS_config.est_clock[i] = clock_cur[i];
              }

            }
            ClearPendingGNSS();
            ///------------------------------------------------------------END REPROPAGATION---------------------------------------------////
            
          }
          ////-----------------------------Determine Whether Enough GNSS Signal Available -------------------------------///////////  
    }

    ///---------------------------------Update GNSS Observation & Run Kalman Filter-----------------------------------/////////////

    ////-----------------------------------------------------------------------TIGHTLY COUPLED ALGORITHM---------------------------------------///////

    ///--------------------------------Update INS-Cirular Buffer at Current-Interval----------------------------------/////////////
        for (int i = 0; i < 3; i++) {
          INS_r_store[store_idx][i] = est_r_eb_e[i];
          INS_v_store[store_idx][i] = est_v_eb_e[i];
          INS_f_store[store_idx][i] = meas_f_ib_b[i];
          INS_w_store[store_idx][i] = meas_omega_ib_b[i];
          
          INS_f_raw_store[store_idx][i] = meas_f_ib_b_raw[i];
          INS_w_raw_store[store_idx][i] = meas_omega_ib_b_raw[i];
        }

        //3D Matrix
        for (int i = 0; i < 9; i++) {
          INS_C_store[store_idx][i] = est_C_b_e[i];
        }

        
        for (size_t i = 0; i < P_DIM; i++) {
          INS_P_store[store_idx][i] = P_matrix[i];
        }

        for (size_t i = 0; i < BIAS_DIM; i++) {
          INS_bias_store[store_idx][i] = est_IMU_Bias[i];
        }

        for (size_t i = 0; i < CLOCK_DIM; i++) {
          INS_clock_store[store_idx][i] = GNSS_config.est_clock[i];
        }

    
    ///-------------------------------Reset ECEF IMU PVA to Directly Use in INS Mechanization of Next Epoch -------------------------/////////////


    ///--------------------------------Output Profile----------------------------------/////////////
    FINISH_IMU_STEP:
    // Store raw current gyro for the next IMU epoch.
    prev_w_ib_b_raw[0] = meas_omega_ib_b_raw[0];
    prev_w_ib_b_raw[1] = meas_omega_ib_b_raw[1];
    prev_w_ib_b_raw[2] = meas_omega_ib_b_raw[2];

    old_time = IMU_time;

    for (size_t i = 0; i < 3; ++i) {
      old_est_r_eb_e[i] = est_r_eb_e[i];
      old_est_v_eb_e[i] = est_v_eb_e[i];
    }

    for (size_t i = 0; i < 9; ++i) {
      old_est_C_b_e[i] = est_C_b_e[i];
    }

    double est_r_ea_e[3];

    // User-requested legacy output lever-arm multiplication retained.
    est_r_ea_e[0] = est_r_eb_e[0]
                    + est_C_b_e[0] * L_ba_b[0]
                    + est_C_b_e[1] * L_ba_b[1]
                    + est_C_b_e[2] * L_ba_b[2];

    est_r_ea_e[1] = est_r_eb_e[1]
                    + est_C_b_e[3] * L_ba_b[0]
                    + est_C_b_e[4] * L_ba_b[1]
                    + est_C_b_e[5] * L_ba_b[2];

    est_r_ea_e[2] = est_r_eb_e[2]
                    + est_C_b_e[6] * L_ba_b[0]
                    + est_C_b_e[7] * L_ba_b[1]
                    + est_C_b_e[8] * L_ba_b[2];

    // Skew_symmetric(meas_omega_ib_b) * L_ba_b
    // This is equivalent to: meas_omega_ib_b cross L_ba_b

    double omega_cross_L_b[3];

    omega_cross_L_b[0] =
      meas_omega_ib_b[1] * L_ba_b[2]
      - meas_omega_ib_b[2] * L_ba_b[1];

    omega_cross_L_b[1] =
      meas_omega_ib_b[2] * L_ba_b[0]
      - meas_omega_ib_b[0] * L_ba_b[2];

    omega_cross_L_b[2] =
      meas_omega_ib_b[0] * L_ba_b[1]
      - meas_omega_ib_b[1] * L_ba_b[0];

    double est_v_ea_e[3];
    // User-requested legacy output lever-arm multiplication retained.
    est_v_ea_e[0] = est_v_eb_e[0]
                    + est_C_b_e[0] * omega_cross_L_b[0]
                    + est_C_b_e[1] * omega_cross_L_b[1]
                    + est_C_b_e[2] * omega_cross_L_b[2];

    est_v_ea_e[1] = est_v_eb_e[1]
                    + est_C_b_e[3] * omega_cross_L_b[0]
                    + est_C_b_e[4] * omega_cross_L_b[1]
                    + est_C_b_e[5] * omega_cross_L_b[2];

    est_v_ea_e[2] = est_v_eb_e[2]
                    + est_C_b_e[6] * omega_cross_L_b[0]
                    + est_C_b_e[7] * omega_cross_L_b[1]
                    + est_C_b_e[8] * omega_cross_L_b[2];
    //------------Testing mode
    if(!ECEF_to_NED(est_C_b_e, est_v_eb_e, est_r_eb_e, old_est_L_b, old_est_lambda_b, old_est_h_b, old_est_v_eb_n, old_est_C_b_n)){
      Serial.println(F("Failed to convert DCM ECEF to NED"));
    }

    if(!CTM_to_Euler(old_est_C_b_n, attitude_init)){
      Serial.println(F("Failed to compute Attitude (RPY)"));
    }
    // Serial.print(attitude_init[2]*RAD_TO_DEG);
    //--------------Testing mode
    if (!pv_ECEF_to_NED(est_r_ea_e, est_v_ea_e, old_est_L_b, old_est_lambda_b, old_est_h_b, old_est_v_eb_n)) {
      Serial.println(F("Failed to convert PV ECEF to NED."));
      while (1) {}
    }
    // double tca_h_msl = old_est_h_b - neo7m.geoid_sep;
    // Serial.println(neo7m.geoid_sep);
    // old_est_h_b = tca_h_msl;
    // Serial.println(F("-----------------------------PVT------------------------------"));

    // Serial.println(GNSS_config.est_clock[1]); //print rx clk drift error
    // Serial.printf("TCA,%0.3f,%0.7f,%0.7f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f\n",
    //             IMU_time, old_est_L_b * RAD_TO_DEG, old_est_lambda_b * RAD_TO_DEG,
    //              old_est_h_b,
    //              old_est_v_eb_n[0], old_est_v_eb_n[1], old_est_v_eb_n[2],
    //              attitude_init[0]*RAD_TO_DEG, attitude_init[1]*RAD_TO_DEG, attitude_init[2]*RAD_TO_DEG,
    //              GNSS_config.est_clock[0],GNSS_config.est_clock[1]);

    // Throttle TCA serial output to 10 Hz. The INS propagation remains 200 Hz.
    // if ((millis() - lastTcaOutputMs) >= 100U) {
    //   lastTcaOutputMs = millis();

      Serial1.printf("TCA,%0.3f,%0.7f,%0.7f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,%0.3f,0.0000,0.0000\n",
                  IMU_time, old_est_L_b * RAD_TO_DEG, old_est_lambda_b * RAD_TO_DEG,
                   old_est_h_b,
                   old_est_v_eb_n[0], old_est_v_eb_n[1], old_est_v_eb_n[2],
                   attitude_init[0]*RAD_TO_DEG, attitude_init[1]*RAD_TO_DEG, attitude_init[2]*RAD_TO_DEG);
    // }             
    
  }

}


static inline void Transpose3x3CM_InPlace(double C[9]){
  double t;

  t = C[1]; C[1] = C[3]; C[3] = t;
  t = C[2]; C[2] = C[6]; C[6] = t;
  t = C[5]; C[5] = C[7]; C[7] = t;
}

//----------------------------------------------------------------INITIAL INITIALIZATIONS-----------------------------------------------------------/////
bool Initial_Initializations() {
  Serial.println(F("-------------------Attitude/Coarse Alignment Started--------------------"));
  double mean_gx, mean_gy, mean_gz;
  double lat0, lon0, h0;

  if (!ComputeInitialAligmentMean(time_last_GNSS, old_time, mean_gx, mean_gy, mean_gz)) {
    Serial.println(F("Failed to Compute Initial Alignment Mean."));
    return false;
  }
  Serial.println(F("------------------Attitude/Coarse Alignment Completed-------------------"));

  Serial.println(F("-------------------------PVT Alignment Started--------------------------"));
  if (!PVTAlignmentAtOldTime(mean_gx, mean_gy, mean_gz, lat0, lon0, h0)) {
    Serial.println(F("Failed to Align PVT."));
    return false;
  }
  Serial.println(F("-------------------------PVT Alignment Completed--------------------------"));

  //Serial.println(F("-------------------------Origin of Receiver Position LLA to ECEF Started--------------------------"));
  ch_LLA2ECEF(lat0, lon0, h0, STA);

  // Serial.println(F("-------------------------Origin of Receiver Position Completed--------------------------"));

  //Serial.println(F("-------------------------Kalman Filter State Covariance P-MATRIX [17 x 17] Initialization Started--------------------------"));
  Initialize_P_Matrix(TC_KF_config);
  //Serial.println(F("-------------------------Kalman Filter State Covariance P-MATRIX [17 x 17] Initialization Completed--------------------------"));

  return true;
}
//----------------------------------------------------------------INITIAL INITIALIZATIONS-----------------------------------------------------------/////

///---------------------------COMPUTE INITIAL ALIGNMENT MEAN USING RLG DATA UPTILL OLD_TIME (DON NOT INCLUDE OLD_TIME)---------------------------------//

//---------------------------------------MEAN BASED ON RLG DATA WITH RESPECT TO NO OF ROWS--------------------------------------------//
bool ComputeInitialAligmentMean(double& time_last_GNSS,
                                double& old_time,
                                double& mean_gx,
                                double& mean_gy,
                                double& mean_gz)
{
  const uint32_t IMU_RATE_HZ = 200;
  const uint32_t ALIGNMENT_TIME_SEC = 60;
  const uint32_t REQUIRED_IMU_EPOCHS = IMU_RATE_HZ * ALIGNMENT_TIME_SEC;
  const uint32_t TIMEOUT_MS = (ALIGNMENT_TIME_SEC + 10) * 1000UL;

  const double MAX_GNSS_ALIGN_DELAY_SEC = 0.105;  // 10 Hz GNSS + tolerance

  double sum_ax = 0.0, sum_ay = 0.0, sum_az = 0.0;
  double sum_gx = 0.0, sum_gy = 0.0, sum_gz = 0.0;

  uint32_t count = 0;
  uint32_t startMs = millis();

  Serial.println(F("Computing initial RLG alignment mean for 60 seconds..."));

  // Clear stale GNSS state before alignment
  matchedGnssTime = 0.0;
  packet.valid = false;
  packet.gnssCount = 0;
  packet.gnss_obssec = 0.0;
  resetGNSSParser();

  while (count < REQUIRED_IMU_EPOCHS) {

    if ((millis() - startMs) > TIMEOUT_MS) {
      Serial.println(F("ERROR: IMU alignment mean timeout."));
      return false;
    }

    // Service raw GNSS parser, but do not use packet.gnss_obssec
    // as matchedGnssTime inside every IMU sample.
    if (readGNSSData()) {
      // Optional diagnostic only
      if (packet.valid && packet.gnss_obssec > 0.0) {
        // Do not assign matchedGnssTime here for alignment ending.
        // This packet may become stale before the 60 s IMU mean completes.
      }
    }

    read_NEO7M_PVT();

    if (readIMU()) {

      if (lastTimeIMU <= 0.0 || !isfinite(lastTimeIMU)) {
        Serial.println(F("ERROR: Invalid IMU timestamp during alignment."));
        return false;
      }

      sum_ax += rlgData.ax;
      sum_ay += rlgData.ay;
      sum_az += rlgData.az;

      sum_gx += rlgData.gx;
      sum_gy += rlgData.gy;
      sum_gz += rlgData.gz;

      // This is the final IMU time of the alignment window.
      old_time = lastTimeIMU;

      count++;
    }
  }

  if (count == 0) {
    Serial.println(F("ERROR: No IMU samples collected."));
    return false;
  }

  double mean_ax = sum_ax / (double)count;
  double mean_ay = sum_ay / (double)count;
  double mean_az = sum_az / (double)count;

  mean_gx = sum_gx / (double)count;
  mean_gy = sum_gy / (double)count;
  mean_gz = sum_gz / (double)count;

  Serial.println(F("IMU 60-second mean completed."));

  // ------------------------------------------------------------------
  // Now align the GNSS raw epoch to old_time.
  //
  // We require the first raw GNSS epoch at or after old_time.
  // For 10 Hz GNSS:
  // matchedGnssTime - old_time should be approximately 0.000 to 0.100 s.
  // ------------------------------------------------------------------
  if (!WaitForFreshRawGNSSAtOrAfterTime(old_time,
                                        matchedGnssTime,
                                        3000,
                                        MAX_GNSS_ALIGN_DELAY_SEC)) {
    Serial.println(F("ERROR: Could not find fresh GNSS raw epoch near old_time."));
    return false;
  }

  double align_dt = matchedGnssTime - old_time;

  Serial.printf("PRT,%0.6f,%0.6f,%0.6f\n",
                old_time,
                matchedGnssTime,
                align_dt);

  if (align_dt < -1.0e-6 || align_dt > MAX_GNSS_ALIGN_DELAY_SEC) {
    Serial.println(F("ERROR: matchedGnssTime is not close to old_time."));
    return false;
  }

  time_last_GNSS = matchedGnssTime;

  /*
   * Now read NEO7M PVT aligned to old_time.
   * But in your simulator architecture, neo7m.timestamp is currently assigned
   * from packet.gnss_obssec. Therefore this should use matchedGnssTime
   * or forced/manual PVT timestamp, not a stale packet.
   */
  if (!ReadNEO7M_PVT_NearestTime(matchedGnssTime, 3000, 1)) {
    Serial.println(F("ERROR: No NEO7M PVT close to matchedGnssTime."));
    return false;
  }

  return Get_ComputeThetaPhi(mean_ax,
                             mean_ay,
                             mean_az,
                             mean_gx,
                             mean_gy,
                             mean_gz);
}
//---------------------------------------MEAN BASED ON RLG DATA WITH RESPECT TO NO OF ROWS--------------------------------------------//

bool ReadNEO7M_PVT_NearestTime(double target_time,
                               uint32_t timeout_ms,
                               double max_allowed_dt_sec) {
  uint32_t startMs = millis();

  bool haveBest = false;
  double bestAbsDt = 1e9;

  // Store best PVT values manually
  double best_timestamp = 0.0;
  double best_lat = 0.0;
  double best_lon = 0.0;
  double best_hei = 0.0;
  double best_speed = 0.0;
  double best_course = 0.0;

  while ((millis() - startMs) < timeout_ms) {

    if (read_NEO7M_PVT()) {

      double absDt = fabs(neo7m.timestamp - target_time);
      Serial.println(neo7m.timestamp, 7);
      Serial.println(target_time, 7);
      Serial.println(absDt);
      if (absDt < bestAbsDt) {
        bestAbsDt = absDt;
        haveBest = true;

        best_timestamp = neo7m.timestamp;
        best_lat = neo7m.lat;
        best_lon = neo7m.lon;
        best_hei = neo7m.hei;
        best_speed = neo7m.speed;
        best_course = neo7m.course;
      }
      readGNSSData();
      readIMU();

      // Stop after passing target time
      if (neo7m.timestamp >= target_time && haveBest) {
        break;
      }
    }
  }

  if (!haveBest || bestAbsDt > max_allowed_dt_sec) {
    return false;
  }

  // Copy best PVT back to global neo7m
  neo7m.timestamp = best_timestamp;
  neo7m.lat = best_lat;
  neo7m.lon = best_lon;
  neo7m.hei = best_hei;
  neo7m.speed = best_speed;
  neo7m.course = best_course;

  return true;
}
bool Get_ComputeThetaPhi(double mean_ax,
                         double mean_ay,
                         double mean_az,
                         double mean_gx,
                         double mean_gy,
                         double mean_gz) {
  /*
    Assumption:
    - accelerometer units: m/s^2
    - gyro units from RLG: deg/s
    - body/static level convention gives mean_az ≈ -9.81 m/s^2
  */

  const double DEG_TO_RAD_LOCAL = PI / 180.0;

  // Convert gyro from deg/s to rad/s
  double wx = mean_gx * DEG_TO_RAD_LOCAL;
  double wy = mean_gy * DEG_TO_RAD_LOCAL;
  double wz = mean_gz * DEG_TO_RAD_LOCAL;

  // Roll and pitch leveling
  double phi_nb = atan2(-mean_ay, -mean_az);

  double horizontal_g = sqrt(mean_ay * mean_ay + mean_az * mean_az);
  double theta_nb = atan2(mean_ax, horizontal_g);

  // Direct gyrocompassing
  double sin_phi = sin(phi_nb);
  double cos_phi = cos(phi_nb);
  double sin_theta = sin(theta_nb);
  double cos_theta = cos(theta_nb);

  double sin_psi = -wy * cos_phi + wz * sin_phi;

  double cos_psi =  wx * cos_theta
                  + wy * sin_phi * sin_theta
                  + wz * cos_phi * sin_theta;

  double psi_nb = atan2(sin_psi, cos_psi);

  attitude_init[0] = phi_nb;
  attitude_init[1] = theta_nb;
  attitude_init[2] = psi_nb;

  Serial.println(F("-----------------------Roll-Pitch-Yaw Estimation------------------------"));
  Serial.print(F("phi_nb [deg]: "));
  Serial.println(attitude_init[0] * RAD_TO_DEG, 7);

  Serial.print(F("theta_nb [deg]: "));
  Serial.println(attitude_init[1] * RAD_TO_DEG, 7);

  Serial.print(F("psi_nb [deg]: "));
  Serial.println(attitude_init[2] * RAD_TO_DEG, 7);

  if (!Euler_to_CTM(attitude_init, old_est_C_b_n)) {
    Serial.println(F("Failed to Compute Body to NED Transformation Matrix."));
    return false;
  }

  Transpose3x3CM_InPlace(old_est_C_b_n);
  return true;
}
///---------------------------COMPUTE INITIAL ALIGNMENT MEAN USING RLG DATA UPTILL OLD_TIME (DON NOT INCLUDE OLD_TIME)---------------------------------//


///-------------------------------------------------------COMPUTE PVT ALIGNMENT-----------------------------------------------------------//
bool PVTAlignmentAtOldTime(double mean_gx, double mean_gy, double mean_gz, double& lat0, double& lon0, double& h0) {

  // Step 4: use REF lat/lon/h for conversion
  lat0 = neo7m.lat * DEG_TO_RAD;
  lon0 = neo7m.lon * DEG_TO_RAD;
  h0 = neo7m.hei;

  // double inti_v_eb_n[3]   = {refData.vel_north, refData.vel_east, refData.vel_down}; //ned vel from pocketsdr
  double inti_v_eb_n[3] = { neo7m.speed * cos(neo7m.course * DEG_TO_RAD), neo7m.speed * sin(neo7m.course * DEG_TO_RAD), 0 };
  double old_est_r_ea_e[3] = { 0.0, 0.0, 0.0 };  // [1x3] ECEF position of receiver XYZ (m)
  double old_est_v_ea_e[3] = { 0.0, 0.0, 0.0 };  // [1x3] ECEF velocity of receiver XYZ (m/s)

  Serial.println(F("---------------------Initial PVT at old_time------------------------------ "));
  Serial.print(F("lat0:"));
  Serial.println(lat0 * RAD_TO_DEG, 7);
  Serial.print(F("lon0:"));
  Serial.println(lon0 * RAD_TO_DEG, 7);
  Serial.print(F("h0:"));
  Serial.println(h0, 7);
  Serial.print(F("inti_v_eb_n[0]:"));
  Serial.println(inti_v_eb_n[0], 7);
  Serial.print(F("inti_v_eb_n[1]:"));
  Serial.println(inti_v_eb_n[1], 7);
  Serial.print(F("inti_v_eb_n[2]:"));
  Serial.println(inti_v_eb_n[2], 7);
  if (!NED_to_ECEF(lat0, lon0, h0, inti_v_eb_n, old_est_C_b_n, old_est_r_ea_e, old_est_v_ea_e, old_est_C_b_e)) {
    Serial.println(F("DCM NED_to_ECEF Failed"));
    return false;
  }
  for(int i =0; i<3; i++){
    L_ba_b[i] = L_ia_b[i] + L_imu_b[i];
  }
  
  ///------------------- PVA of IMU(IMU to Receiver L_ba_b Compensated) (a) w.r.t Earth Center(e) in ECEF(e) -------------------------/////
  for (int i = 0; i < 3; i++) {
    old_est_r_eb_e[i] = old_est_r_ea_e[i] - (old_est_C_b_e[i] * L_ba_b[0] + old_est_C_b_e[i + 3] * L_ba_b[1] + old_est_C_b_e[i + 6] * L_ba_b[2]);
  }

  const double omega_ie = 7.292115e-5;
  
  double wx = mean_gx * DEG_TO_RAD;
  double wy = mean_gy * DEG_TO_RAD;
  double wz = mean_gz * DEG_TO_RAD;

  double omega_ie_b[3];

  omega_ie_b[0] = old_est_C_b_e[2] * omega_ie;  
  omega_ie_b[1] = old_est_C_b_e[5] * omega_ie;  
  omega_ie_b[2] = old_est_C_b_e[8] * omega_ie;  


  // omega_eb_b = omega_ib_b - C_e_b * omega_ie_e
  double omega_eb_b[3];

  omega_eb_b[0] = wx - omega_ie_b[0];
  omega_eb_b[1] = wy - omega_ie_b[1];
  omega_eb_b[2] = wz - omega_ie_b[2];

  // Compute Skew(w) * L without making the skew matrix
  // This is Skew_symmetric([wx;wy;wz]) * L_ba_b
  const double t0 = omega_eb_b[1] * L_ba_b[2] - omega_eb_b[2] * L_ba_b[1];
  const double t1 = omega_eb_b[2] * L_ba_b[0] - omega_eb_b[0] * L_ba_b[2];
  const double t2 = omega_eb_b[0] * L_ba_b[1] - omega_eb_b[1] * L_ba_b[0];

  // ECEF velocity of IMU from antenna velocity
  for (int i = 0; i < 3; i++) {
    old_est_v_eb_e[i] = old_est_v_ea_e[i] - (old_est_C_b_e[i] * t0 + old_est_C_b_e[i + 3] * t1 + old_est_C_b_e[i + 6] * t2);
  }

  ///------------------- PVA of IMU(IMU to Receiver L_ba_b Compensated) (a) w.r.t Earth Center(e) in ECEF(e) -------------------------/////

  if (!pv_ECEF_to_NED(old_est_r_eb_e, old_est_v_eb_e, old_est_L_b, old_est_lambda_b, old_est_h_b, old_est_v_eb_n)) {
    Serial.println(F("Failed to convert pv_ECEF to NED."));
    while (1) {}
  }

  if (!NED_to_ECEF(old_est_L_b, old_est_lambda_b, old_est_h_b, old_est_v_eb_n, old_est_C_b_n, old_est_C_b_e)) {
    Serial.println(F("Failed to convert NED to ECEF."));
    while (1) {}
  }
  
 
    //------------Testing mode
    if(!ECEF_to_NED(old_est_C_b_e, old_est_v_eb_e, old_est_r_eb_e, old_est_L_b, old_est_lambda_b, old_est_h_b, old_est_v_eb_n, old_est_C_b_n)){
      Serial.println(F("Failed to convert DCM ECEF to NED"));
    }

    Transpose3x3CM_InPlace(old_est_C_b_n);

    if(!CTM_to_Euler(old_est_C_b_n, attitude_init)){
      Serial.println(F("Failed to compute Attitude (RPY)"));
    }
    Serial.print(attitude_init[2]*RAD_TO_DEG);
    
    Transpose3x3CM_InPlace(old_est_C_b_n);

  Serial.println(F("---------------------Lever arm compensated Initial PVAT------------------------------ "));
  Serial.println(F("-----------------------------Old_est_L_b------------------------------"));
  Serial.print(F("old_est_L_b:"));
  Serial.println(old_est_L_b * RAD_TO_DEG, 7);
  Serial.println(F("-----------------------------Old_est_lambda_b------------------------------"));
  Serial.print(F("old_est_lambda_b:"));
  Serial.println(old_est_lambda_b * RAD_TO_DEG, 7);
  Serial.println(F("-----------------------------Old_est_h_b------------------------------"));
  Serial.print(F("old_est_h_b:"));
  Serial.println(old_est_h_b, 7);
  Serial.println(F("-----------------------------Old_est_v_eb_n------------------------------"));
  Serial.print(F("old_est_v_eb_n[0]:"));
  Serial.println(old_est_v_eb_n[0], 7);
  Serial.print(F("old_est_v_eb_n[1]:"));
  Serial.println(old_est_v_eb_n[1], 7);
  Serial.print(F("old_est_v_eb_n[2]:"));
  Serial.println(old_est_v_eb_n[2], 7);
  Serial.println(F("-----------------------------Old_est_C_b_e------------------------------"));
  Serial.print(F("old_est_C_b_e[0]:"));
  Serial.println(old_est_C_b_e[0], 7);
  Serial.print(F("old_est_C_b_e[1]:"));
  Serial.println(old_est_C_b_e[1], 7);
  Serial.print(F("old_est_C_b_e[2]:"));
  Serial.println(old_est_C_b_e[2], 7);
  Serial.print(F("old_est_C_b_e[3]:"));
  Serial.println(old_est_C_b_e[3], 7);
  Serial.print(F("old_est_C_b_e[4]:"));
  Serial.println(old_est_C_b_e[4], 7);
  Serial.print(F("old_est_C_b_e[5]:"));
  Serial.println(old_est_C_b_e[5], 7);
  Serial.print(F("old_est_C_b_e[6]:"));
  Serial.println(old_est_C_b_e[6], 7);
  Serial.print(F("old_est_C_b_e[7]:"));
  Serial.println(old_est_C_b_e[7], 7);
  Serial.print(F("old_est_C_b_e[8]:"));
  Serial.println(old_est_C_b_e[8], 7);

  return true;
}
///-------------------------------------------------------PVT ALIGNMENT-----------------------------------------------------------//
///-------------------------------------------------------INITIALIZE KALMAN FILTER STATE COVARIANCE P-MATRIX [17 x 17]-----------------------------------------------------------//
void Initialize_P_Matrix(const TC_KFConfig& TC_KF_Config) {
  memset(P_matrix, 0, sizeof(P_matrix));

  P_matrix[0] = TC_KF_Config.init_att_unc * TC_KF_Config.init_att_unc;
  P_matrix[18] = TC_KF_Config.init_att_unc * TC_KF_Config.init_att_unc;
  P_matrix[36] = TC_KF_Config.init_att_unc * TC_KF_Config.init_att_unc;

  P_matrix[54] = TC_KF_Config.init_vel_unc * TC_KF_Config.init_vel_unc;
  P_matrix[72] = TC_KF_Config.init_vel_unc * TC_KF_Config.init_vel_unc;
  P_matrix[90] = TC_KF_Config.init_vel_unc * TC_KF_Config.init_vel_unc;

  P_matrix[108] = TC_KF_Config.init_pos_unc * TC_KF_Config.init_pos_unc;
  P_matrix[126] = TC_KF_Config.init_pos_unc * TC_KF_Config.init_pos_unc;
  P_matrix[144] = TC_KF_Config.init_pos_unc * TC_KF_Config.init_pos_unc;

  P_matrix[162] = TC_KF_Config.init_b_a_unc * TC_KF_Config.init_b_a_unc;
  P_matrix[180] = TC_KF_Config.init_b_a_unc * TC_KF_Config.init_b_a_unc;
  P_matrix[198] = TC_KF_Config.init_b_a_unc * TC_KF_Config.init_b_a_unc;

  P_matrix[216] = TC_KF_Config.init_b_g_unc * TC_KF_Config.init_b_g_unc;
  P_matrix[234] = TC_KF_Config.init_b_g_unc * TC_KF_Config.init_b_g_unc;
  P_matrix[252] = TC_KF_Config.init_b_g_unc * TC_KF_Config.init_b_g_unc;

  P_matrix[270] = TC_KF_Config.init_clock_offset_unc * TC_KF_Config.init_clock_offset_unc;
  P_matrix[288] = TC_KF_Config.init_clock_drift_unc * TC_KF_Config.init_clock_drift_unc;

}
///-------------------------------------------------------INITIALIZE KALMAN FILTER STATE COVARIANCE P-MATRIX [17 x 17]-----------------------------------------------------------//

///-------------------INS current time retrieval-----------------------------
bool FindINSIntervalForGNSSTime(double gnssTime, int &idx1, int &idx2) {
  idx1 = -1;
  idx2 = -1;

  lastINSIntervalStatus = INS_INTERVAL_NOT_FOUND;

  // Small tolerance for floating-point timestamp comparison
  const double TIME_TOL = 1.0e-6;

  if (gnssTime <= 0.0 || isnan(gnssTime) || !isfinite(gnssTime)) {
    Serial.println(F("ERROR: Invalid GNSS time."));
    lastINSIntervalStatus = INS_INTERVAL_INVALID_GNSS_TIME;
    return false;
  }

  if (IMU_epoch_counter < 2 || latest_INS_store_idx < 0) {
    Serial.println(F("ERROR: Not enough INS samples stored."));
    lastINSIntervalStatus = INS_INTERVAL_NOT_ENOUGH_INS;
    return false;
  }

  int nStored = (IMU_epoch_counter < INS_BUFFER_SIZE)
                  ? IMU_epoch_counter
                  : INS_BUFFER_SIZE;

  int latestIdx = latest_INS_store_idx;

  if (latestIdx < 0 || latestIdx >= INS_BUFFER_SIZE) {
    Serial.println(F("ERROR: Invalid latest_INS_store_idx."));
    lastINSIntervalStatus = INS_INTERVAL_INVALID_BUFFER_TIME;
    return false;
  }

  double latestTime = INS_time_store[latestIdx];

  int oldestIdx = (IMU_epoch_counter >= INS_BUFFER_SIZE)
                    ? ((latest_INS_store_idx + 1) % INS_BUFFER_SIZE)
                    : 0;

  if (oldestIdx < 0 || oldestIdx >= INS_BUFFER_SIZE) {
    Serial.println(F("ERROR: Invalid oldest INS index."));
    lastINSIntervalStatus = INS_INTERVAL_INVALID_BUFFER_TIME;
    return false;
  }

  double oldestTime = INS_time_store[oldestIdx];

  if (latestTime <= 0.0 || oldestTime <= 0.0 ||
      !isfinite(latestTime) || !isfinite(oldestTime)) {
    Serial.println(F("ERROR: INS buffer has invalid oldest/latest time."));
    lastINSIntervalStatus = INS_INTERVAL_INVALID_BUFFER_TIME;
    return false;
  }


  // ------------------------------------------------------------
  // GNSS epoch is newer than latest INS.
  // This is NOT an error.
  // Caller should HOLD this GNSS epoch pending.
  // ------------------------------------------------------------
  if (gnssTime > latestTime + TIME_TOL) {
    lastINSIntervalStatus = INS_INTERVAL_GNSS_AHEAD_OF_INS;
    return false;
  }

  // ------------------------------------------------------------
  // GNSS epoch is older than available INS history.
  // Caller should DISCARD this GNSS epoch.
  // ------------------------------------------------------------
  if (gnssTime < oldestTime - TIME_TOL) {
    Serial.println(F("GNSS time is older than oldest INS buffer time."));
    lastINSIntervalStatus = INS_INTERVAL_GNSS_OLDER_THAN_BUFFER;
    return false;
  }

  // ------------------------------------------------------------
  // If GNSS is very close to latest INS, use the interval ending
  // at latestIdx.
  // ------------------------------------------------------------
  if (fabs(gnssTime - latestTime) <= TIME_TOL) {
    int older = (latestIdx - 1 + INS_BUFFER_SIZE) % INS_BUFFER_SIZE;

    if (INS_valid_store[older] &&
        INS_time_store[older] > 0.0 &&
        INS_time_store[older] < latestTime) {
      idx1 = older;
      idx2 = latestIdx;
      lastINSIntervalStatus = INS_INTERVAL_OK;
      return true;
    }
  }

  // ------------------------------------------------------------
  // Search backward from latest sample until the bracketing
  // interval is found.
  // ------------------------------------------------------------
  for (int k = 0; k < nStored - 1; k++) {
    int newer = (latestIdx - k + INS_BUFFER_SIZE) % INS_BUFFER_SIZE;
    int older = (newer - 1 + INS_BUFFER_SIZE) % INS_BUFFER_SIZE;

    if (!INS_valid_store[older] || !INS_valid_store[newer]) {
      continue;
    }

    double t_old = INS_time_store[older];
    double t_new = INS_time_store[newer];

    if (t_old <= 0.0 || t_new <= 0.0 ||
        !isfinite(t_old) || !isfinite(t_new)) {
      continue;
    }

    // Timestamp must increase from older to newer.
    // If this fails, there is a time jump or corrupted buffer entry.
    if (t_new <= t_old) {
      Serial.print(F("WARNING: non-monotonic INS buffer interval: t_old="));
      Serial.print(t_old, 6);
      Serial.print(F(" t_new="));
      Serial.println(t_new, 6);
      continue;
    }

    if ((t_old - TIME_TOL) <= gnssTime &&
        gnssTime <= (t_new + TIME_TOL)) {
      idx1 = older;
      idx2 = newer;
      lastINSIntervalStatus = INS_INTERVAL_OK;
      return true;
    }
  }

  Serial.println(F("ERROR: GNSS time within oldest/latest range but no bracketing INS interval found."));
  lastINSIntervalStatus = INS_INTERVAL_NOT_FOUND;
  return false;
} 



bool waitForInitialGNSS(int requiredPackets, uint32_t timeoutMs) {
  int validCount = 0;
  uint32_t startTime = millis();

  while (validCount < requiredPackets) {

    readIMU();
    readGNSSData();
    // if (readGNSSData()) {  
    //     gps_obs = packet.gnss_obssec;   // keep latest raw GNSS time     
    // }
    if (read_NEO7M_PVT()) {
      validCount++;

      Serial.print(F("Valid Initial GNSS packet received: "));
      Serial.println(validCount);
    }

    if (millis() - startTime > timeoutMs) {
      return false;
    }
  }

  return true;
}

//===================================================//
//     NEO7m and Pocket SDR PVT Parsing function      //
//===================================================//

bool read_NEO7M_PVT() {
  while (GNSS_SERIAL.available() > 0) {
    char c = GNSS_SERIAL.read();
    gnss.encode(c);
    // Serial.write(c);
  }

  if (gnss.location.isUpdated() && gnss.altitude.isUpdated() && gnss.speed.isUpdated() && gnss.course.isUpdated() && gnss.date.isUpdated() && gnss.time.isUpdated() && gnss.satellites.isUpdated()) {
    // neo7m.lat = gnss.location.lat();
    // neo7m.lon = gnss.location.lng();
    // TinyGPS++ altitude from NMEA GGA is usually MSL / orthometric height
    neo7m.hei_msl = gnss.altitude.meters();
    neo7m.lat = 37.352721;
    neo7m.lon = -121.915773;
    neo7m.hei_msl = 20;
    // Parse GGA geoid separation field
    neo7m.geoid_valid = false;

    if (gga_geoid_GPGGA.isValid()) {
      neo7m.geoid_sep = atof(gga_geoid_GPGGA.value());
      neo7m.geoid_valid = true;
    } else if (gga_geoid_GNGGA.isValid()) {
      neo7m.geoid_sep = atof(gga_geoid_GNGGA.value());
      neo7m.geoid_valid = true;
    }

    // Convert MSL height to ellipsoid height for ECEF initialization
    if (neo7m.geoid_valid) {
      neo7m.hei_ellipsoid = neo7m.hei_msl + neo7m.geoid_sep;
    } else {
      // fallback only; not ideal for ECEF initialization
      neo7m.hei_ellipsoid = neo7m.hei_msl;
    }

    // Keep old variable name as ellipsoid height
    // neo7m.hei = neo7m.hei_ellipsoid;
    neo7m.hei = 20;
    // neo7m.speed = 0.0;
    // neo7m.course = 318.91;
    
    neo7m.speed = gnss.speed.mps();
    neo7m.course = gnss.course.deg();

    neo7m.hdop = gnss.hdop.hdop();

    neo7m.year = gnss.date.year();
    neo7m.month = gnss.date.month();
    neo7m.day = gnss.date.day();

    neo7m.hour = gnss.time.hour();
    neo7m.minute = gnss.time.minute();
    neo7m.second = gnss.time.second();
    neo7m.centisecond = gnss.time.centisecond();

    // neo7m.timestamp = utcToGpsSeconds(neo7m.year, neo7m.month, neo7m.day, neo7m.hour, neo7m.minute, float(neo7m.second) + float(neo7m.centisecond) / 100.0);
    neo7m.timestamp = packet.gnss_obssec;
    neo7m.satellites = gnss.satellites.value();

    return true;
  }

  return false;
}

//===================================================//
//                PPS Time Management                //
//===================================================//

volatile uint32_t prev_pps_time = 0;
volatile uint32_t pps_interval_us = 1000000;
volatile bool pps_interval_valid = false;

void ppsISR()
{
  const uint32_t now_us = micros();

  if (pps_seen) {
    prev_pps_time = last_pps_time;

    // unsigned subtraction is safe across micros() rollover
    pps_interval_us = now_us - last_pps_time;

    pps_interval_valid = true;
  }

  last_pps_time = now_us;
  pps_count++;
  pps_seen = true;
}

//===================================================//
//                Precise Time Value                 //
//===================================================//

double getPreciseTime() {
  
  noInterrupts();
  uint32_t lastPpsUsCopy = last_pps_time;
  uint32_t ppsCountCopy  = pps_count;
  uint32_t ppsCountAtSyncCopy = ppsCountAtSync;
  double startingTimeCopy = starting_time;
  interrupts();

  uint32_t nowUs = micros();
  double secSinceLastPps = (double)(nowUs - lastPpsUsCopy) * 1e-6;

  uint32_t ppsDelta = ppsCountCopy - ppsCountAtSyncCopy;
  if (ppsCountCopy == 0){
    return startingTimeCopy + (double(nowUs - start_micros) * 1e-6);
  }
  
  return startingTimeCopy + (double)ppsDelta + secSinceLastPps;
}

//===================================================//
//                  UTC → GPS seconds                //
//===================================================//

uint32_t daysSince1970(int Y, int M, int D) {
  if (M <= 2) {
    Y--;
    M += 12;
  }
  return 365UL * Y + Y / 4 - Y / 100 + Y / 400 + (153 * (M - 3) + 2) / 5 + D - 719469UL;
}

double utcToGpsSeconds(int y, int m, int d, int h, int min, double sec) {
  uint32_t days = daysSince1970(y, m, d);
  double unixSec = days * 86400.0 + h * 3600.0 + min * 60.0 + sec;

  double gpsSec = unixSec + GPS_UTC_LEAP - 315964800.0;
  return fmod(gpsSec, 604800.0);
}




//===================================================//
//     IMU  parser and data type conversion functions                         //
//===================================================//

static bool parsePythonDouble(
    const char* text,
    double& value) {

  if (text == nullptr || *text == '\0') {
    return false;
  }

  char* endPointer = nullptr;
  value = strtod(text, &endPointer);

  if (endPointer == text ||
      *endPointer != '\0' ||
      !isfinite(value)) {
    return false;
  }

  return true;
}

// Returns true only when one complete IMU sample has
// been assigned to the existing rlgData variables.
static bool processPythonIMULine(
    char* line,
    double packetStartTime)
{
    char* savePointer = nullptr;

    // ===========================================================
    // First token = message type
    // ===========================================================
    char* messageType =
        strtok_r(
            line,
            ",",
            &savePointer);

    if (messageType == nullptr)
    {
        // pythonParseErrors++;
        return false;
    }

    // ===========================================================
    // Python may still transmit:
    //
    // IMU_STREAM_BEGIN,count,sample_rate
    //
    // We simply ignore it.
    //
    // NO READY response is transmitted.
    // ===========================================================
    if (strcmp(
            messageType,
            "IMU_STREAM_BEGIN") == 0)
    {
        // Optional fresh counters at beginning of simulation.
        // pythonReceivedSamples = 0;
        // pythonParseErrors = 0;
        // pythonOversizedLines = 0;
        // pythonMaxRxBacklog = 0;

        // IMPORTANT:
        // No READY command.
        return false;
    }

    // ===========================================================
    // Python may still transmit:
    //
    // IMU_STREAM_END,count
    //
    // Ignore it.
    //
    // NO DONE response is transmitted.
    // ===========================================================
    if (strcmp(
            messageType,
            "IMU_STREAM_END") == 0)
    {
        return false;
    }

    // ===========================================================
    // We are interested only in IMU packets.
    // ===========================================================
    if (strcmp(messageType, "IMU") != 0)
    {
        return false;
    }

    // ===========================================================
    // Current Python packet:
    //
    // IMU,index,time,ax,ay,az,gx,gy,gz
    //
    // Ignore index and time.
    // ===========================================================

    char* ignoredIndex =
        strtok_r(
            nullptr,
            ",",
            &savePointer);

    char* ignoredTime =
        strtok_r(
            nullptr,
            ",",
            &savePointer);

    // They must exist because the Python format contains them,
    // but their numerical values are deliberately not parsed.
    if ((ignoredIndex == nullptr) ||
        (ignoredTime == nullptr))
    {
        // pythonParseErrors++;
        return false;
    }

    // Avoid compiler unused-variable warnings.
    (void)ignoredIndex;
    (void)ignoredTime;

    // ===========================================================
    // Extract ONLY:
    //
    // ax
    // ay
    // az
    // gx
    // gy
    // gz
    // ===========================================================

    char* axToken =
        strtok_r(
            nullptr,
            ",",
            &savePointer);

    char* ayToken =
        strtok_r(
            nullptr,
            ",",
            &savePointer);

    char* azToken =
        strtok_r(
            nullptr,
            ",",
            &savePointer);

    char* gxToken =
        strtok_r(
            nullptr,
            ",",
            &savePointer);

    char* gyToken =
        strtok_r(
            nullptr,
            ",",
            &savePointer);

    char* gzToken =
        strtok_r(
            nullptr,
            ",",
            &savePointer);

    // ===========================================================
    // All six sensor fields must exist.
    // ===========================================================
    if ((axToken == nullptr) ||
        (ayToken == nullptr) ||
        (azToken == nullptr) ||
        (gxToken == nullptr) ||
        (gyToken == nullptr) ||
        (gzToken == nullptr))
    {
        // pythonParseErrors++;
        return false;
    }

    // ===========================================================
    // Reject unexpected additional fields.
    // ===========================================================
    if (strtok_r(
            nullptr,
            ",",
            &savePointer) != nullptr)
    {
        // pythonParseErrors++;
        return false;
    }

    // ===========================================================
    // Convert six ASCII values -> double.
    // ===========================================================

    double ax = 0.0;
    double ay = 0.0;
    double az = 0.0;

    double gx = 0.0;
    double gy = 0.0;
    double gz = 0.0;

    if (!parsePythonDouble(axToken, ax) ||
        !parsePythonDouble(ayToken, ay) ||
        !parsePythonDouble(azToken, az) ||
        !parsePythonDouble(gxToken, gx) ||
        !parsePythonDouble(gyToken, gy) ||
        !parsePythonDouble(gzToken, gz))
    {
        // pythonParseErrors++;
        return false;
    }

    // ===========================================================
    // ACCELEROMETER
    //
    // gnss-ins-sim:
    // m/s^2
    // ===========================================================

    rlgData.ax = ax;
    rlgData.ay = ay;
    rlgData.az = az;

    // ===========================================================
    // GYROSCOPE
    //
    // gnss-ins-sim sends rad/s.
    //
    // Your existing INS later multiplies gyro by DEG_TO_RAD,
    // therefore convert rad/s -> deg/s here.
    // ===========================================================

    constexpr double PYTHON_GYRO_RAD_TO_DEG =
        180.0 / PI;

    rlgData.gx =
        gx * PYTHON_GYRO_RAD_TO_DEG;

    rlgData.gy =
        gy * PYTHON_GYRO_RAD_TO_DEG;

    rlgData.gz =
        gz * PYTHON_GYRO_RAD_TO_DEG;

    // ===========================================================
    // Timestamp
    //
    // Python simulator timestamp is NOT parsed.
    // Use Teensy/PPS timestamp instead.
    // ===========================================================

    rlgData.timestamp =
        packetStartTime;

    lastTimeIMU =
        rlgData.timestamp;

    // ===========================================================
    // Fields unavailable from the simulator
    // ===========================================================

    rlgData.roll = 0.0;
    rlgData.pitch = 0.0;
    rlgData.yaw = 0.0;

    rlgData.acc_x_temp = 0.0;
    rlgData.acc_y_temp = 0.0;
    rlgData.acc_z_temp = 0.0;

    rlgData.gyro_x_temp = 0.0;
    rlgData.gyro_y_temp = 0.0;
    rlgData.gyro_z_temp = 0.0;

    rlgData.if_temp = 0.0;
    rlgData.timeOffset = 0.0;
    rlgData.gpsStatus = 0;

    // Existing legacy frame counter.
    rlgData.frameCount =
        static_cast<uint8_t>(
            rlgData.frameCount + 1U);

    // ===========================================================
    // Inform existing INS code that a new sample is available.
    // ===========================================================

    IMUSIM_Available = true;

    // pythonReceivedSamples++;

    return true;
}

bool waitForInitialIMU(
    int requiredPackets,
    uint32_t timeoutMs)
{
    if (requiredPackets <= 0)
    {
        return true;
    }

    int validCount = 0;

    const uint32_t startTime =
        millis();

    while (validCount < requiredPackets)
    {
        // --------------------------------------------------------
        // Continue servicing other serial interfaces.
        // --------------------------------------------------------
        readGNSSData();
        read_NEO7M_PVT();

        // --------------------------------------------------------
        // Drain all currently available valid IMU packets.
        // --------------------------------------------------------
        if (readIMU()) {
      validCount++;

      Serial.print(F("Valid RLG packet received: "));
      Serial.println(validCount);
    }

    if (millis() - startTime > timeoutMs) {
      return false;
    }
  }

  return true;
}

bool readIMU()
{
    // ===========================================================
    // Used only to discard a partial packet left by a disconnect.
    //
    // At 200 Hz, normal spacing is 5 ms.
    // A 100 ms gap is therefore safely considered a long idle
    // period for purposes of clearing a PARTIAL ASCII line.
    // ===========================================================
    static uint32_t lastRxByteMs = 0;
    static bool haveReceivedByte = false;

    uint32_t currentMs = millis();

    int availableNow =
        IMU_SERIAL.available();

    // ===========================================================
    // If UART has been idle and a partial packet remains,
    // discard only that partial packet.
    //
    // Do NOT reset valid rlgData.
    // ===========================================================
    if (availableNow == 0)
    {
        if (haveReceivedByte &&
            ((uint32_t)(currentMs - lastRxByteMs) > 100U))
        {
            if ((imuFrameIndex != 0) ||
                IMUsimDiscardOversizedLine)
            {
                imuFrameIndex = 0;
                IMUsimDiscardOversizedLine = false;
            }
        }

        return false;
    }

    // ===========================================================
    // Track maximum UART backlog.
    //
    // Do not print from here.
    // ===========================================================
    // if ((availableNow > 0) &&
    //     (static_cast<uint32_t>(availableNow) >
    //      pythonMaxRxBacklog))
    // {
    //     pythonMaxRxBacklog =
    //         static_cast<uint32_t>(availableNow);
    // }

    // ===========================================================
    // Consume Serial2 bytes.
    //
    // Return after exactly ONE valid IMU packet so rlgData isn't
    // overwritten before the INS consumes it.
    // ===========================================================
    while (IMU_SERIAL.available() > 0)
    {
        availableNow =
            IMU_SERIAL.available();

        // if ((availableNow > 0) &&
        //     (static_cast<uint32_t>(availableNow) >
        //      pythonMaxRxBacklog))
        // {
        //     pythonMaxRxBacklog =
        //         static_cast<uint32_t>(availableNow);
        // }

        const int receivedByte =
            IMU_SERIAL.read();

        if (receivedByte < 0)
        {
            break;
        }

        currentMs = millis();

        // =======================================================
        // Detect a long physical gap before this new byte.
        //
        // If the old connection disappeared in the middle of a
        // packet, discard that old partial packet before accepting
        // bytes from the new connection.
        // =======================================================
        if (haveReceivedByte &&
            ((uint32_t)(currentMs - lastRxByteMs) > 100U))
        {
            if ((imuFrameIndex != 0) ||
                IMUsimDiscardOversizedLine)
            {
                imuFrameIndex = 0;
                IMUsimDiscardOversizedLine = false;
            }
        }

        lastRxByteMs = currentMs;
        haveReceivedByte = true;

        const char incomingCharacter =
            static_cast<char>(receivedByte);

        // =======================================================
        // Previously detected oversized/corrupted line.
        //
        // Discard until newline and then resynchronize.
        // =======================================================
        if (IMUsimDiscardOversizedLine)
        {
            if (incomingCharacter == '\n')
            {
                IMUsimDiscardOversizedLine = false;
                imuFrameIndex = 0;
            }

            continue;
        }

        // Ignore carriage return from CRLF.
        if (incomingCharacter == '\r')
        {
            continue;
        }

        // =======================================================
        // Newline = complete Python packet.
        // =======================================================
        if (incomingCharacter == '\n')
        {
            // Ignore blank line.
            if (imuFrameIndex == 0)
            {
                continue;
            }

            // Null terminate ASCII line.
            imuFrame[imuFrameIndex] = '\0';

            const double completedPacketStartTime =
                imuFrameStartTime;

            // Immediately prepare for next packet.
            imuFrameIndex = 0;

            // ===================================================
            // true means:
            //
            // ax ay az gx gy gz
            //
            // have successfully been copied into rlgData.
            // ===================================================
            if (processPythonIMULine(
                    reinterpret_cast<char*>(imuFrame),
                    completedPacketStartTime))
            {
                return true;
            }

            // BEGIN, END, bad line, etc.
            // Continue looking for a valid IMU line.
            continue;
        }

        // =======================================================
        // First character of packet.
        //
        // This creates the local Teensy timestamp. The Python
        // timestamp is deliberately ignored.
        // =======================================================
        if (imuFrameIndex == 0)
        {
            imuFrameStartTime =
                getPreciseTime();
        }

        // =======================================================
        // Store received ASCII character.
        // =======================================================
        if (imuFrameIndex <
            (IMU_FRAME_LEN - 1))
        {
            imuFrame[imuFrameIndex++] =
                static_cast<uint8_t>(
                    incomingCharacter);
        }
        else
        {
            // ===================================================
            // Line is too large or corrupted.
            //
            // Discard everything until next newline.
            // ===================================================
            // pythonOversizedLines++;

            IMUsimDiscardOversizedLine = true;

            imuFrameIndex = 0;
        }
    }

    return false;
}




//===================================================//
//     Pocket SDR RAW DATA Parsing functions    //
//===================================================//
bool IsGnssObssecAtIntegerSecond(double gps_obs)
{
  if (gps_obs <= 0.0 || !isfinite(gps_obs)) {
    return false;
  }

  double nearestInteger = round(gps_obs);
  double error = fabs(gps_obs - nearestInteger);

  // For 10 Hz GNSS, integer epochs should be .000.
  // This tolerance handles floating-point/parser small errors.
  return (error < 0.02);
}
bool waitForIntegerSecondPocketRawForSync(uint32_t timeoutMs)
{
  uint32_t startMs = millis();

  while ((millis() - startMs) < timeoutMs) {

    read_NEO7M_PVT();

    if (readGNSSData()) {

      if (packet.valid &&
          packet.gnssCount > 0 &&
          packet.gnss_obssec > 0.0 &&
          isfinite(packet.gnss_obssec)) {

        Serial.print(F("SYNC_CANDIDATE obssec="));
        Serial.println(packet.gnss_obssec, 6);

        if (IsGnssObssecAtIntegerSecond(packet.gnss_obssec)) {
          Serial.print(F("SYNC_ACCEPT obssec="));
          Serial.println(packet.gnss_obssec, 6);
          return true;
        } else {
          Serial.print(F("SYNC_SKIP obssec="));
          Serial.println(packet.gnss_obssec, 6);
        }
      }
    }
  }

  return false;
}

// --------------------------------------------------
// Trim leading and trailing spaces
// --------------------------------------------------
static void trimString(char* str) {
  if (str == nullptr) return;

  // Trim leading spaces
  while (*str == ' ' || *str == '\t') {
    memmove(str, str + 1, strlen(str));
  }

  // Trim trailing spaces
  int len = strlen(str);
  while (len > 0 && (str[len - 1] == ' ' || str[len - 1] == '\t' || str[len - 1] == '\r')) {
    str[len - 1] = '\0';
    len--;
  }
}

// --------------------------------------------------
// Parse TCA header:
// TCA,<num_sats>,<obssec>
// --------------------------------------------------
static bool parseTCAHeader(char* line) {
  char* savePtr = nullptr;

  char* token = strtok_r(line, ",", &savePtr);
  if (token == nullptr) return false;

  if (strcmp(token, "TCA") != 0) return false;

  token = strtok_r(nullptr, ",", &savePtr);
  if (token == nullptr) return false;

  int sats = atoi(token);

  if (sats <= 0 || sats > static_cast<int>(MAX_GNSS_ROWS_PER_EPOCH)) {
    return false;
  }

  token = strtok_r(nullptr, ",", &savePtr);
  if (token == nullptr) return false;

  double obssec = atof(token);

  expectedSats = static_cast<uint16_t>(sats);

  memset(&tempGNSS, 0, sizeof(tempGNSS));
  tempGNSS.gnssCount = 0;
  tempGNSS.gnss_obssec = obssec;
  tempGNSS.valid = false;

  return true;
}
//--------------------------------------
//      Signal name to number conversion
//--------------------------------------
uint8_t convertSignalNameToID(const char* signalName) {
  if (strcmp(signalName, "GPS L1 C/A") == 0) return SIG_GPS_L1_CA;
  if (strcmp(signalName, "Galileo E1B") == 0) return SIG_GALILEO_E1B;
  if (strcmp(signalName, "Galileo E1C") == 0) return SIG_GALILEO_E1C;
  if (strcmp(signalName, "Galileo E5a I") == 0) return SIG_GALILEO_E5A_I;
  if (strcmp(signalName, "Galileo E5b I") == 0) return SIG_GALILEO_E5B_I;
  if (strcmp(signalName, "BeiDou B1I") == 0) return SIG_BEIDOU_B1I;
  if (strcmp(signalName, "BeiDou B1CD") == 0) return SIG_BEIDOU_B1CD;
  if (strcmp(signalName, "BeiDou B2a D") == 0) return SIG_BEIDOU_B2A_D;
  if (strcmp(signalName, "BeiDou B2B") == 0) return SIG_BEIDOU_B2B;
  if (strcmp(signalName, "BeiDou B2I") == 0) return SIG_BEIDOU_B2I;
  if (strcmp(signalName, "BeiDou B3I") == 0) return SIG_BEIDOU_B3I;
  if (strcmp(signalName, "GLONASS G1 C/A") == 0) return SIG_GLONASS_G1_CA;

  return SIG_UNKNOWN;
}


// --------------------------------------------------
// Parse one satellite observation row:
//
// <pr_corr>,<prrate_corr>,<x>,<y>,<z>,
// <vx>,<vy>,<vz>,<clk_rate>,<elev>,<CN0>,<Signal Name>
// --------------------------------------------------
static bool parseGNSSObsRow(char* line, GNSSRow& row) {
  char* savePtr = nullptr;
  char* token = nullptr;

  token = strtok_r(line, ",", &savePtr);
  if (token == nullptr) return false;
  row.pr_corr = atof(token);

  token = strtok_r(nullptr, ",", &savePtr);
  if (token == nullptr) return false;
  row.prrate_corr = atof(token);

  token = strtok_r(nullptr, ",", &savePtr);
  if (token == nullptr) return false;
  row.sat_x = atof(token);

  token = strtok_r(nullptr, ",", &savePtr);
  if (token == nullptr) return false;
  row.sat_y = atof(token);

  token = strtok_r(nullptr, ",", &savePtr);
  if (token == nullptr) return false;
  row.sat_z = atof(token);

  token = strtok_r(nullptr, ",", &savePtr);
  if (token == nullptr) return false;
  row.sat_vx = atof(token);

  token = strtok_r(nullptr, ",", &savePtr);
  if (token == nullptr) return false;
  row.sat_vy = atof(token);

  token = strtok_r(nullptr, ",", &savePtr);
  if (token == nullptr) return false;
  row.sat_vz = atof(token);

  token = strtok_r(nullptr, ",", &savePtr);
  if (token == nullptr) return false;
  row.sat_clk_rate = atof(token);

  token = strtok_r(nullptr, ",", &savePtr);
  if (token == nullptr) return false;
  row.el_deg = atof(token);

  token = strtok_r(nullptr, ",", &savePtr);
  if (token == nullptr) return false;
  row.cn0 = atof(token);

  token = strtok_r(nullptr, ",", &savePtr);
  if (token == nullptr) return false;

  trimString(token);

  strncpy(row.signalName, token, GNSS_SIGNAL_LEN - 1);
  row.signalName[GNSS_SIGNAL_LEN - 1] = '\0';
  row.signal = convertSignalNameToID(row.signalName);

  return true;
}

static char gnssLine[GNSS_LINE_MAX];
static uint16_t gnssLineLen = 0;

bool readGNSSData() {
  while (GNSSRAW_SERIAL.available() > 0) {
    
    char ch = (char)GNSSRAW_SERIAL.read();
    // Serial.println(ch);
    if (ch == '\r') {
      continue;
    }

    if (ch == '\n') {
      gnssLine[gnssLineLen] = '\0';

      bool epochComplete = processGNSSLine(gnssLine);

      gnssLineLen = 0;

      if (epochComplete) {
        return true;
      }
    } 
    else {
      if (gnssLineLen < GNSS_LINE_MAX - 1) {
        gnssLine[gnssLineLen++] = ch;
      } else {
        Serial.println(F("ERROR: GNSS line overflow."));
        gnssLineLen = 0;
        resetGNSSParser();
      }
    }
  }

  return false;
}

bool isEndLine(const char *line) {
  return (strcmp(line, "END") == 0);
}

void resetGNSSParser() {
  gnssState = WAIT_TCA_HEADER;
  expectedSats = 0;

  tempGNSS.gnssCount = 0;
  tempGNSS.gnss_obssec = 0.0;
  tempGNSS.valid = false;
}

bool processGNSSLine(char *line) {
  trimString(line);

  if (strlen(line) == 0) {
    return false;
  }

  switch (gnssState) {

    case WAIT_TCA_HEADER:
    {
      if (strncmp(line, "TCA,", 4) == 0) {

        resetGNSSParser();

        if (parseTCAHeader(line)) {
          tempGNSS.gnssCount = 0;
          tempGNSS.valid = false;
          gnssState = READ_TCA_ROWS;
        } else {
          resetGNSSParser();
        }
      }

      return false;
    }

    case READ_TCA_ROWS:
    {
      if (isEndLine(line)) {

        if (tempGNSS.gnssCount == expectedSats && expectedSats > 0 &&
      tempGNSS.gnss_obssec > 0.0) {
          tempGNSS.valid = true;

          packet = tempGNSS;

          resetGNSSParser();

          return true;
        }

        Serial.print(F("ERROR: GNSS row mismatch. expected="));
        Serial.print(expectedSats);
        Serial.print(F(", parsed="));
        Serial.println(tempGNSS.gnssCount);

        resetGNSSParser();
        return false;
      }

      if (tempGNSS.gnssCount < expectedSats &&
          tempGNSS.gnssCount < MAX_GNSS_ROWS_PER_EPOCH) {

        GNSSRow row;

        if (parseGNSSObsRow(line, row)) {
          tempGNSS.gnssRows[tempGNSS.gnssCount] = row;
          tempGNSS.gnssCount++;
        } else {
          Serial.println(F("ERROR: Bad GNSS observation row."));
          // Serial.println(line);
          resetGNSSParser();
        }

        return false;
      }

      Serial.println(F("ERROR: Too many GNSS rows."));
      resetGNSSParser();
      return false;
    }
  }

  return false;
}

bool SyncPreciseTimeFromPocketRaw(double gps_obs) {
  if (gps_obs <= 0.0 || !isfinite(gps_obs)) {
    Serial.println(F("ERROR: Invalid gps_obs."));
    return false;
  }
  

  // if (!pps_seen) {
  //   Serial.println(F("ERROR: PPS not seen yet."));
  //   return false;
  // }

  noInterrupts();
  uint32_t lastPpsUsCopy = last_pps_time;
  uint32_t ppsCountCopy  = pps_count;
  interrupts();

  uint32_t nowUs = micros();
  double secSinceLastPps = (double)(nowUs - lastPpsUsCopy) * 1e-6;

  if (secSinceLastPps < 0.0 || secSinceLastPps > 1.2) {
    Serial.println(F("ERROR: PPS age invalid."));
    return false;
  }

  double gpsFloor = floor(gps_obs);
  double gpsFrac  = gps_obs - gpsFloor;

  /*
    starting_time is the GPS TOW at the last PPS edge.

    Example 1:
      gps_obs = 554773.201
      secSinceLastPps = 0.350
      last PPS corresponds to 554773.000

    Example 2:
      gps_obs = 554773.201
      secSinceLastPps = 0.050
      last PPS likely corresponds to 554774.000
  */
  if (secSinceLastPps < gpsFrac) {
    starting_time = gpsFloor + 1.0;
  } else {
    starting_time = gpsFloor;
  }

  ppsCountAtSync = ppsCountCopy;
  preciseTimeValid = true;

  
  return true;
}


// ==============================
// Debug print helper functions
// Pocet SDR RAW Data
// ==============================

void printSignalNameSafe(const char* name) {
  for (int i = 0; i < GNSS_SIGNAL_LEN; ++i) {
    if (name[i] == '\0') {
      break;
    }
    Serial.print(name[i]);
  }
}

void printGNSSPacketCompact(const EpochPacket& packet) {

  Serial.print(F(", count: "));
  Serial.print((unsigned long)packet.gnssCount);

  Serial.print(F(",GNSS obssec: "));
  Serial.print(packet.gnss_obssec, 6);

  Serial.print(F(", valid: "));
  Serial.println(packet.valid ? F("true") : F("false"));

  for (size_t i = 0; i < packet.gnssCount; ++i) {
    const GNSSRow& row = packet.gnssRows[i];

    // Serial.print(F("i="));
    // Serial.print((unsigned long)i);

    Serial.print(F(", sigName="));
    printSignalNameSafe(row.signalName);

    Serial.print(F(", sig="));
    Serial.print(row.signal);

    Serial.print(F(", PR="));
    Serial.print(row.pr_corr, 3);

    Serial.print(F(", PRR="));
    Serial.print(row.prrate_corr, 6);

    Serial.print(F(", el="));
    Serial.print(row.el_deg, 2);

    Serial.print(F(", cn0="));
    Serial.print(row.cn0, 2);

    Serial.print(F(", sat_clk_rate="));
    Serial.print(row.sat_clk_rate, 2);
    
    Serial.print(F(", sat_xyz=["));
    Serial.print(row.sat_x, 3);
    Serial.print(F(", "));
    Serial.print(row.sat_y, 3);
    Serial.print(F(", "));
    Serial.print(row.sat_z, 3);
    Serial.print(F("]"));

    Serial.print(F(", sat_vxyz=["));
    Serial.print(row.sat_vx, 6);
    Serial.print(F(", "));
    Serial.print(row.sat_vy, 6);
    Serial.print(F(", "));
    Serial.print(row.sat_vz, 6);
    Serial.println(F("]"));
  }
}

void ClearINSBuffer() {
  for (size_t i = 0; i < INS_BUFFER_SIZE; i++) {
    INS_time_store[i] = 0.0;
    INS_valid_store[i] = false;
  }

  IMU_epoch_counter = 0;
  latest_INS_store_idx = -1;
}

// -----------------------------------------------------------------------------
// Drain queued serial input
// -----------------------------------------------------------------------------
void drainSerialInput(Stream &s,
                      uint32_t drain_ms)
{
  uint32_t t0 = millis();

  while ((millis() - t0) < drain_ms) {
    while (s.available() > 0) {
      s.read();
    }
    yield();
  }
}


// -----------------------------------------------------------------------------
// Reset parser states after draining stale data
// -----------------------------------------------------------------------------
void resetAllInputParserStates()
{
  // Pocket SDR raw GNSS parser
  resetGNSSParser();

  // Clear current GNSS packet
  packet.valid = false;
  packet.gnssCount = 0;
  packet.gnss_obssec = 0.0;

  GNSS_EPOCH_AVAILABLE = false;
  matchedGnssTime = 0.0;
  lastAcceptedGnssTime = 0.0;

  // IMU ASCII parser state
  imuFrameIndex = 0;
  IMUsimDiscardOversizedLine = false;

  IMUSIM_Available = false;
}


// -----------------------------------------------------------------------------
// Drain all serial inputs used in this project
// -----------------------------------------------------------------------------
void drainAllSerialInputs(uint32_t drain_ms)
{
  // USB debug serial input, if any data is being sent from PC to Teensy
  drainSerialInput(Serial, drain_ms);

  // Pocket SDR raw GNSS UART
  drainSerialInput(GNSSRAW_SERIAL, drain_ms);

  // NEO7M PVT UART
  drainSerialInput(GNSS_SERIAL, drain_ms);

  // Python IMU UART
  drainSerialInput(IMU_SERIAL, drain_ms);

  resetAllInputParserStates();
}

bool WaitForFreshRawGNSSAtOrAfterTime(double target_time,
                                      double& matchedGnssTime,
                                      uint32_t timeout_ms,
                                      double max_delay_sec)
{
  uint32_t startMs = millis();

  matchedGnssTime = 0.0;

  while ((millis() - startMs) < timeout_ms) {

    // Keep NEO parser serviced if required
    read_NEO7M_PVT();

    if (readGNSSData()) {

      if (!packet.valid || packet.gnssCount == 0 || packet.gnss_obssec <= 0.0) {
        continue;
      }

      double gnssTime = packet.gnss_obssec;
      double dt = gnssTime - target_time;

      Serial.print(F("RAW_ALIGN_CHECK gnssTime="));
      Serial.print(gnssTime, 6);
      Serial.print(F(" old_time="));
      Serial.print(target_time, 6);
      Serial.print(F(" gnss-old="));
      Serial.println(dt, 6);

      // GNSS is still older than alignment end time.
      // Discard it.
      if (dt < -1.0e-6) {
        continue;
      }

      // First GNSS epoch at or after old_time.
      if (dt <= max_delay_sec) {
        matchedGnssTime = gnssTime;
        return true;
      }

      // If the first future GNSS epoch is too far away,
      // then GNSS was missing or parser was delayed.
      Serial.println(F("ERROR: Fresh GNSS epoch too far after old_time."));
      return false;
    }
  }

  return false;
}

void ClearPendingGNSS()
{
  GNSS_EPOCH_AVAILABLE = false;
  packet.valid = false;
  packet.gnssCount = 0;
  packet.gnss_obssec = 0.0;
  matchedGnssTime = 0.0;
}