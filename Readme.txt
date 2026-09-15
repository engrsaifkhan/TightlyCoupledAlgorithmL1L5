This code is designed for realtime application with sensor integration
Update:
Fused PVT is transmitting but accuracy improvement work is on going...
Last static position accuracy achieved less than 10m RMS in both horizontal and vertical and max is 21m in horizontal
Dynamic test is hold on due to working of GNSS data corection 
Last test date: 3-july-2026

--- VS Code / PlatformIO (Teensy 4.1) ---
Project is converted from Arduino .ino to C++ PlatformIO layout:
  src/main.cpp   - main sketch (setup/loop)
  include/       - headers
  src/*.cpp      - supporting sources
  platformio.ini - Teensy 4.1 target
Original sketch kept as TCA_IMUGNSS_SIM_15_09_2026.ino.bak

Build/upload in VS Code or Cursor:
  1) Install PlatformIO IDE extension
  2) Open this folder
  3) PlatformIO: Build  (checkmark) then Upload (arrow)
Serial monitor baud: 921600