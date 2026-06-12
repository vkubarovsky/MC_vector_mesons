#!/bin/bash
# run_all.sh -- Generate Born+RC for rho, jpsi, phi simultaneously
#
# Usage:
#   ./run_all.sh -version NAME [-nev N] [-njobs N] [-born-only] [-no-compile]
#
# All three species start at the same time.
# Output logs: logs/rho.log, logs/jpsi.log, logs/phi.log

ROOT="$(cd "$(dirname "$0")" && pwd)"
LOGDIR="$ROOT/logs"

NJOBS=10
NEV=10000
VERSION="vm"
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case $1 in
        -nev)        NEV="$2";            shift 2 ;;
        -njobs)      NJOBS="$2";          shift 2 ;;
        -version)    VERSION="$2";        shift 2 ;;
        -born-only)  EXTRA_ARGS+=(-born-only); shift ;;
        -no-compile) EXTRA_ARGS+=(-no-compile); shift ;;
        -help|--help)
            echo "Usage: ./run_all.sh -version NAME [-nev N] [-njobs N] [-born-only] [-no-compile]"
            echo ""
            echo "  -version    NAME  Model version: akushevich, harut, ...  (required)"
            echo "                    Compiles diffrad_NAME.f90, writes to DIFFRAD_lund/NAME/"
            echo "  -nev        N     Events per job per species  (default: $NEV)"
            echo "  -njobs      N     Parallel jobs per species   (default: $NJOBS)"
            echo "  -born-only        Skip RC generation"
            echo "  -no-compile       Skip recompilation"
            exit 0 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

if [[ -z "$VERSION" ]]; then
    echo "ERROR: -version NAME is required (e.g. -version akushevich)"
    echo "       Run ./run_all.sh -help for usage."
    exit 1
fi

SRC_FILE="$ROOT/diffrad_${VERSION}.f90"
if [[ ! -f "$SRC_FILE" ]]; then
    echo "ERROR: source file not found: $SRC_FILE"
    exit 1
fi

LUND_BASE="$HOME/Downloads/DIFFRAD_lund/${VERSION}"
mkdir -p "$LOGDIR"

# Write lund_VERSION pointer files (plain text with the target path).
# Plain text works on OneDrive; Python scripts read the path via resolve_lund_dir().
for SP in rho phi jpsi; do
    LINK="$ROOT/$SP/lund_${VERSION}"
    TARGET="$LUND_BASE/$SP"
    mkdir -p "$TARGET"
    echo "$TARGET" > "$LINK"
    echo "  Pointer: $SP/lund_${VERSION} -> $TARGET"
done

echo "======================================================"
echo " DIFFRAD generator — all species"
echo "======================================================"
echo " Version        : $VERSION  (diffrad_${VERSION}.f90)"
echo " Events per job : $NEV"
echo " Jobs per species: $NJOBS"
echo " Total events   : $((NEV * NJOBS)) per species"
echo " LUND output    : $LUND_BASE/"
echo " Species        : rho, jpsi, phi"
echo " Logs           : $LOGDIR/"
echo "======================================================"
echo ""
echo "Launching all three species now ..."
echo ""

START=$(date +%s)

for SP in rho jpsi phi; do
    LOG="$LOGDIR/${SP}.log"
    (
        cd "$ROOT/$SP"
        ./run_parallel.sh \
            -nev    "$NEV" \
            -njobs  "$NJOBS" \
            -src    "$SRC_FILE" \
            -lund   "$LUND_BASE/$SP/born_events.lund" \
            "${EXTRA_ARGS[@]}" \
            > "$LOG" 2>&1
    ) &
    eval "PID_${SP}=$!"
    eval "echo \"  $SP  started  (PID \$PID_${SP}, log: logs/${SP}.log)\""
done

echo ""
echo "Waiting for all jobs to finish ..."
echo ""

ALL_OK=true
for SP in rho jpsi phi; do
    eval "pid=\$PID_${SP}"
    if wait "$pid"; then
        echo "  $SP  DONE"
    else
        echo "  $SP  FAILED  (check logs/${SP}.log)"
        ALL_OK=false
    fi
done

END=$(date +%s)
ELAPSED=$((END - START))
MINS=$((ELAPSED / 60))
SECS=$((ELAPSED % 60))

echo ""
echo "======================================================"
echo " Elapsed: ${MINS}m ${SECS}s"
echo "======================================================"

if $ALL_OK; then
    echo ""
    echo " Output files:"
    for SP in rho jpsi phi; do
        DIR="$LUND_BASE/$SP"
        for F in born_events.lund rc_events.lund; do
            if [[ -f "$DIR/$F" ]]; then
                LINES=$(wc -l < "$DIR/$F")
                printf "   %-8s %-20s  %d lines\n" "$SP" "$F" "$LINES"
            fi
        done
    done
    echo ""
    echo " Sigma summary:"
    for SP in rho jpsi phi; do
        DIR="$LUND_BASE/$SP"
        for STAT in born_events_stat.dat rc_events_stat.dat; do
            if [[ -f "$DIR/$STAT" ]]; then
                SIG=$(grep sigma_nb "$DIR/$STAT" | awk '{print $3}')
                printf "   %-8s %-28s  sigma = %s nb\n" "$SP" "$STAT" "$SIG"
            fi
        done
    done
    echo ""
    echo " To plot, run in each species dir:"
    echo "   python plot_compare.py -version $VERSION"
    echo "   python plot_eta.py     -version $VERSION"
fi

echo ""
echo "======================================================"
