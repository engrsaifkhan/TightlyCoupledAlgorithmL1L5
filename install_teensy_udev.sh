#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
RULES_SRC="$ROOT/00-teensy.rules"
RULES_DST="/etc/udev/rules.d/00-teensy.rules"

if [[ ! -f "$RULES_SRC" ]]; then
  echo "Missing $RULES_SRC"
  exit 1
fi

echo "Installing Teensy udev rules (needs sudo)..."
sudo cp "$RULES_SRC" "$RULES_DST"
sudo chmod 644 "$RULES_DST"
sudo udevadm control --reload-rules
sudo udevadm trigger

# Immediate fix for currently attached Teensy (bootloader or serial)
echo "Relaxing permissions on currently attached Teensy USB nodes..."
for node in /dev/bus/usb/*/*; do
  if udevadm info -q property -n "$node" 2>/dev/null | grep -q 'ID_VENDOR_ID=16c0'; then
    sudo chmod 666 "$node" || true
    ls -l "$node"
  fi
done
for h in /dev/hidraw* /dev/ttyACM*; do
  [[ -e "$h" ]] || continue
  if udevadm info -q property -n "$h" 2>/dev/null | grep -q 'ID_VENDOR_ID=16c0'; then
    sudo chmod 666 "$h" || true
    ls -l "$h"
  fi
done

echo
echo "Done. Unplug/replug the Teensy (or press the reset button), then upload again:"
echo "  platformio run --target upload --environment teensy41"
