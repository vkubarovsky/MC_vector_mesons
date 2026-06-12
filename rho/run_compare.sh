#!/bin/bash
# ============================================================
# run_compare.sh  --  Plot comparison of Born and RC samples
# ============================================================

BORN="lund/born_events.lund"
RC="lund/rc_events.lund"
INPUT=""
NEV=0

while [[ $# -gt 0 ]]; do
    case $1 in
        -born)  BORN="$2";  shift 2 ;;
        -rc)    RC="$2";    shift 2 ;;
        -input) INPUT="$2"; shift 2 ;;
        -nev)   NEV="$2";   shift 2 ;;
        -help|--help)
            echo "Usage: ./run_compare.sh [options]"
            echo ""
            echo "Options:"
            echo "  -born  FILE  Born LUND file                  (default: $BORN)"
            echo "  -rc    FILE  RC LUND file                    (default: $RC)"
            echo "  -input FILE  Input file for kinematic cuts   (cuts off if absent)"
            echo "  -nev   N     Max events to read (0 = all)    (default: $NEV)"
            echo "  -help        Show this help"
            exit 0 ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

ARGS="-born $BORN -rc $RC"
[[ -n "$INPUT" ]] && ARGS="$ARGS -input $INPUT"
[[ "$NEV" -gt 0 ]] && ARGS="$ARGS -nev $NEV"

echo "Running: python plot_compare.py $ARGS"
python3 plot_compare.py $ARGS
