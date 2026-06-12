#!/bin/bash
# run_diffrad.sh -- Compile and run Born-only + full-RC samples
#
# Usage:
#   ./run_diffrad.sh [options]
#
# Options:
#   -src         FILE   Fortran source file    (default: ../diffrad_gen.f90)
#   -born-input  FILE   Born input file        (default: gen_input_born.dat)
#   -rc-input    FILE   RC input file          (default: gen_input_rc.dat)
#   -born-lund   FILE   Born LUND output       (default: born_events.lund)
#   -rc-lund     FILE   RC LUND output         (default: rc_events.lund)
#   -no-compile         Skip compilation

set -e

SRC="../diffrad_akushevich.f90"
BORN_INPUT="gen_input_born.dat"
RC_INPUT="gen_input_rc.dat"
BORN_LUND="born_events.lund"
RC_LUND="rc_events.lund"
COMPILE=true

while [[ $# -gt 0 ]]; do
    case $1 in
        -src)         SRC="$2";         shift 2 ;;
        -born-input)  BORN_INPUT="$2";  shift 2 ;;
        -rc-input)    RC_INPUT="$2";    shift 2 ;;
        -born-lund)   BORN_LUND="$2";   shift 2 ;;
        -rc-lund)     RC_LUND="$2";     shift 2 ;;
        -no-compile)  COMPILE=false;    shift   ;;
        *) echo "Unknown argument: $1"; exit 1  ;;
    esac
done

if $COMPILE; then
    echo "Compiling $SRC ..."
    gfortran -ffree-line-length-none -std=legacy -fwrapv -w -O2 \
             "$SRC" -o diffrad_gen.exe
    echo "Compilation successful."
    echo ""
fi

echo "=== Run 1: Born only ==="
./diffrad_gen.exe -input "$BORN_INPUT" -lund "$BORN_LUND"
echo ""

echo "=== Run 2: Full RC ==="
./diffrad_gen.exe -input "$RC_INPUT" -lund "$RC_LUND"
echo ""

echo "Both runs complete."
echo "Now run: python plot_eta.py -born $BORN_LUND -rc $RC_LUND"
