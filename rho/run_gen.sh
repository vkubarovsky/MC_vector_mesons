#!/bin/bash
# ============================================================
# run_gen.sh  --  Compile and run diffrad_gen
# ============================================================

set -e

INPUT="gen_input.dat"
LUND="events.lund"
COMPILE=true

while [[ $# -gt 0 ]]; do
    case $1 in
        -input)      INPUT="$2";      shift 2 ;;
        -lund)       LUND="$2";       shift 2 ;;
        -no-compile) COMPILE=false;   shift   ;;
        -help|--help)
            echo "Usage: ./run_gen.sh [options]"
            echo ""
            echo "Options:"
            echo "  -input  FILE  Generator input file  (default: $INPUT)"
            echo "  -lund   FILE  Output LUND file      (default: $LUND)"
            echo "  -no-compile   Skip compilation"
            echo "  -help         Show this help"
            exit 0 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if $COMPILE; then
    echo "Compiling diffrad_gen.f90 ..."
    gfortran -ffree-line-length-none -std=legacy -fwrapv -w -O2 \
             diffrad_gen.f90 -o diffrad_gen.exe
    echo "Compilation successful."
    echo ""
fi

echo "Running generator ..."
echo "Input: $INPUT"
echo "-------------------------------"
./diffrad_gen.exe -input "$INPUT" -lund "$LUND"
echo "-------------------------------"

echo ""
echo "Output files:"
echo "  $LUND : $(wc -l < "$LUND") lines"
echo ""
echo "To plot results:"
echo "  python plot_gen.py -lund $LUND"
