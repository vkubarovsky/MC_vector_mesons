#!/bin/bash
# ============================================================
# run_eta.sh  --  Generate Born and RC samples, compute eta
# ============================================================
set -e

echo "Compiling diffrad_gen.f90 ..."
gfortran -ffree-line-length-none -std=legacy -fwrapv -w -O2 \
         diffrad_gen.f90 -o diffrad_gen.exe
echo "Compilation successful."
echo ""

# Run 1: Born only
echo "=== Run 1: Born only ==="
cp gen_input_born.dat gen_input.dat
./diffrad_gen.exe
mv born_events.lund born_events.lund
echo "Born events: $(wc -l < born_events.lund) lines"
echo ""

# Run 2: RC with vmax=1.2
echo "=== Run 2: RC vmax=1.2 ==="
cp gen_input_rc12.dat gen_input.dat
./diffrad_gen.exe
mv rc_events.lund rc_events.lund
echo "RC events: $(wc -l < rc_events.lund) lines"
echo ""

echo "Now run: python plot_eta.py"
