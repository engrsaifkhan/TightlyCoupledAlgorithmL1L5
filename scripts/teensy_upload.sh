#!/usr/bin/env bash
# Reliable Teensy 4.1 upload for PlatformIO on Linux.
# HalfKay is flaky for ~1-2s after soft reboot; settle + retry fixes "error writing".
set -uo pipefail

HEX="${1:-}"
if [[ -z "$HEX" || ! -f "$HEX" ]]; then
  echo "Usage: $0 <firmware.hex>"
  exit 1
fi

TOOL_DIR="${HOME}/.platformio/packages/tool-teensy"
export LD_LIBRARY_PATH="${TOOL_DIR}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
REBOOT="${TOOL_DIR}/teensy_reboot"
LOADER="${TOOL_DIR}/teensy_loader_cli"

halfkay_present() {
  lsusb -d 16c0:0478 >/dev/null 2>&1
}

serial_present() {
  lsusb -d 16c0:0483 >/dev/null 2>&1
}

if command -v fuser >/dev/null 2>&1; then
  if fuser /dev/ttyACM0 >/dev/null 2>&1; then
    echo "WARNING: /dev/ttyACM0 is open — close Serial Monitor if upload fails."
  fi
fi

# Enter bootloader if currently running serial firmware
if serial_present && ! halfkay_present; then
  echo "Soft-rebooting Teensy into HalfKay..."
  if [[ -x "$REBOOT" ]]; then
    "$REBOOT" -s 2>/dev/null || true
  fi
fi

# Wait for HalfKay (press physical reset button if soft reboot fails)
echo "Waiting for HalfKay bootloader (press Teensy PROGRAM/reset button if needed)..."
for _ in $(seq 1 60); do
  if halfkay_present; then
    break
  fi
  sleep 0.5
done

if ! halfkay_present; then
  echo "ERROR: HalfKay not found. Press the Teensy reset button and retry."
  exit 1
fi

# Critical: let USB/HID finish enumerating before the first write attempt
echo "HalfKay found — settling USB..."
sleep 2

# Retry writes; first attempt after reboot often fails with "error writing to Teensy"
for attempt in 1 2 3 4 5; do
  echo "Programming attempt ${attempt}/5..."
  if "$LOADER" --mcu=TEENSY41 -v "$HEX"; then
    echo "Upload OK."
    exit 0
  fi
  echo "Write failed; re-checking HalfKay..."
  sleep 1
  if ! halfkay_present; then
    echo "HalfKay disappeared — trying soft reboot again..."
    [[ -x "$REBOOT" ]] && "$REBOOT" -s 2>/dev/null || true
    sleep 2
  else
    sleep 1
  fi
done

echo "ERROR: upload failed after retries. Unplug/replug Teensy, press reset, try again."
exit 1
