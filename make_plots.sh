#!/bin/bash
# make_plots.sh -- Regenerate all plots for rho, phi, jpsi and copy to report/
#
# Usage:  bash make_plots.sh [-nev N] [-no-open]
#   -nev N      Max events to read per file (default: 0 = all)
#   -no-open    Skip opening PNGs in Preview
#
# Run from MC_vector_mesons/ directory.

set -e

BASE="$(cd "$(dirname "$0")" && pwd)"
NEV=0
OPEN=true

while [[ $# -gt 0 ]]; do
    case $1 in
        -nev)     NEV="$2";   shift 2 ;;
        -no-open) OPEN=false; shift   ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

NEV_ARG=""
[ "$NEV" -gt 0 ] && NEV_ARG="-nev $NEV"

ok()   { echo "  [OK]  $1"; }
warn() { echo "  [--]  $1 (skipped — no lund files)"; }

run_plots() {
    local meson="$1"
    local DIR="$BASE/$meson"
    local LUND="$DIR/lund"

    echo ""
    echo "========================================"
    echo "  $meson"
    echo "========================================"

    if [ ! -f "$LUND/born_events.lund" ] || [ ! -f "$LUND/rc_events.lund" ]; then
        warn "$meson: lund/born_events.lund or rc_events.lund not found"
        return
    fi

    cd "$DIR"

    echo "  plot_compare.py ..."
    python3 plot_compare.py \
        -born lund/born_events.lund \
        -rc   lund/rc_events.lund   \
        $NEV_ARG 2>&1 | grep -E "Saved|ERROR|error" | sed 's/^/    /'
    ok "plot_compare done"

    echo "  plot_eta.py ..."
    python3 plot_eta.py \
        -born lund/born_events.lund \
        -rc   lund/rc_events.lund   \
        $NEV_ARG 2>&1 | grep -E "Saved|sigma|eta|ERROR" | sed 's/^/    /'
    ok "plot_eta done"

    echo "  plot_rc_vcut.py ..."
    python3 "$BASE/plot_rc_vcut.py" \
        -born lund/born_events.lund \
        -rc   lund/rc_events.lund   \
        -label "$meson" \
        -out  "rc_vcut.png" \
        $NEV_ARG 2>&1 | grep -E "Saved|eta|ERROR" | sed 's/^/    /'
    ok "plot_rc_vcut done"
}

run_plots rho
run_plots phi
run_plots jpsi

echo ""
echo "========================================"
echo "  Copying figures to report/"
echo "========================================"
bash "$BASE/report/copy_figures.sh" | sed 's/^/  /'

if $OPEN; then
    echo ""
    echo "Opening figures in Preview ..."
    PNGS=$(ls \
        "$BASE/report/rho_born_validation.png" \
        "$BASE/report/rho_rc_validation.png" \
        "$BASE/report/rho_kinematics.png" \
        "$BASE/report/rho_eta.png" \
        "$BASE/report/phi_born_validation.png" \
        "$BASE/report/phi_rc_validation.png" \
        "$BASE/report/phi_kinematics.png" \
        "$BASE/report/phi_eta.png" \
        "$BASE/report/jpsi_born_validation.png" \
        "$BASE/report/jpsi_rc_validation.png" \
        "$BASE/report/jpsi_kinematics.png" \
        "$BASE/report/jpsi_eta.png" \
        2>/dev/null)
    [ -n "$PNGS" ] && open $PNGS
fi

echo ""
echo "All done.  Now compile:  cd report && pdflatex report.tex"
