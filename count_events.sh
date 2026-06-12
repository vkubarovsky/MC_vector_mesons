#!/usr/bin/env bash
# Count RC events processed so far across all job directories.
# Each event has 2 proton lines (" 2212 ") in the lund file.
#
# Usage: ./count_events.sh {rho|phi|jpsi}

LUND_BASE=~/Downloads/DIFFRAD_lund

usage() { echo "Usage: $0 {rho|phi|jpsi}"; exit 1; }

[[ $# -eq 1 ]] || usage
MESON=$1
MESON_DIR="$LUND_BASE/$MESON"
[[ -d "$MESON_DIR" ]] || { echo "Directory not found: $MESON_DIR"; exit 1; }

jobs=( $(ls -d "$MESON_DIR"/job_* 2>/dev/null | sort -V) )
[[ ${#jobs[@]} -gt 0 ]] || { echo "No job_* directories found in $MESON_DIR"; exit 1; }

total_lines=0
printf "%-10s  %10s\n" "Job" "Events"
printf "%-10s  %10s\n" "----------" "----------"

for job in "${jobs[@]}"; do
    name=$(basename "$job")
    lund_file=$(ls "$job"/*rc*.lund 2>/dev/null | head -1)
    if [[ -z "$lund_file" ]]; then
        printf "%-10s  %10s\n" "$name" "no file"
        continue
    fi
    lines=$(grep -c " 2212 " "$lund_file" 2>/dev/null || echo 0)
    events=$(( lines / 2 ))
    total_lines=$(( total_lines + lines ))
    printf "%-10s  %10d\n" "$name" "$events"
done

total=$(( total_lines / 2 ))
printf "%-10s  %10s\n" "----------" "----------"
printf "%-10s  %10d\n" "TOTAL" "$total"
