#!/bin/bash
# run_parallel.sh -- Run parallel Born+RC generations, then combine

set -e

NJOBS=10
NEV=10000
BORN_INPUT="gen_input_born.txt"
RC_INPUT="gen_input_rc.txt"
COMPILE=true
STATS_ONLY=false
BORN_ONLY=false
VERSION="akushevich"
SRC_OVERRIDE=""
LUND_OVERRIDE=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -njobs)       NJOBS="$2";         shift 2 ;;
        -nev)         NEV="$2";           shift 2 ;;
        -born-input)  BORN_INPUT="$2";    shift 2 ;;
        -rc-input)    RC_INPUT="$2";      shift 2 ;;
        -lund)        LUND_OVERRIDE="$2"; shift 2 ;;
        -src)         SRC_OVERRIDE="$2";  shift 2 ;;
        -version)     VERSION="$2";       shift 2 ;;
        -no-compile)  COMPILE=false;      shift   ;;
        -stats-only)  STATS_ONLY=true;    shift   ;;
        -born-only)   BORN_ONLY=true;     shift   ;;
        -help|--help)
            echo "Usage: ./run_parallel.sh [options]"
            echo ""
            echo "Options:"
            echo "  -version    NAME  Generator version: akushevich, harut, vpk, ...  (default: akushevich)"
            echo "                    Sets source file ../diffrad_NAME.f90 and"
            echo "                    output dir ~/Downloads/DIFFRAD_lund/NAME/PARTICLE/"
            echo "  -njobs      N     Number of parallel jobs        (default: $NJOBS)"
            echo "  -nev        N     Events per job                 (default: $NEV)"
            echo "  -born-input FILE  Born input template            (default: $BORN_INPUT)"
            echo "  -rc-input   FILE  RC input template              (default: $RC_INPUT)"
            echo "  -lund       FILE  Override output born LUND file"
            echo "  -src        FILE  Override source Fortran file"
            echo "  -no-compile       Skip compilation"
            echo "  -born-only        Run Born only (skip RC)"
            echo "  -stats-only       Collect stats from existing job_* dirs only"
            echo "  -help             Show this help"
            echo ""
            echo "Examples:"
            echo "  ./run_parallel.sh -version vpk -nev 1000"
            echo "  ./run_parallel.sh -version akushevich -njobs 10 -nev 10000"
            exit 0 ;;
        *) echo "Unknown argument: $1"; exit 1  ;;
    esac
done

# Derive SRC and LUND_BORN from VERSION (overridden by explicit -src / -lund)
PARTICLE="$(basename "$(pwd)")"
SRC="${SRC_OVERRIDE:-../diffrad_${VERSION}.f90}"
LUND_BORN="${LUND_OVERRIDE:-$HOME/Downloads/DIFFRAD_lund/${VERSION}/${PARTICLE}/born_events.lund}"

LUND_DIR="$(dirname "$LUND_BORN")"
LUND_BASE="$(basename "$LUND_BORN" .lund)"
LUND_RC="$LUND_DIR/rc_events.lund"

WORKDIR="$(pwd)"
mkdir -p "$LUND_DIR"
[[ -d "$LUND_BORN" ]] && rm -rf "$LUND_BORN"

if $STATS_ONLY; then
    echo "=== Stats-only mode: reading existing job_* directories ==="
    echo "LUND directory: $LUND_DIR"
    echo ""
else
    echo "Born output: $LUND_BORN"
    echo "RC   output: $LUND_RC"
    echo ""

    if $COMPILE; then
        echo "Compiling $SRC ..."
        gfortran -ffree-line-length-none -std=legacy -fwrapv -w -O2 \
                 "$SRC" -o diffrad_gen.exe
        echo "Compilation successful."
        echo ""
    fi

    echo "Launching $NJOBS parallel jobs ($NEV events each) ..."
    pids=()
    for i in $(seq 1 $NJOBS); do
        DIR="$LUND_DIR/job_${i}"
        rm -rf "$DIR"
        mkdir -p "$DIR"
        cp diffrad_gen.exe "$DIR/"
        SEED=$((333522 + i * 100003))

        (
            cd "$DIR"

            awk -v nev=$NEV -v seed=$SEED \
                '/^nev[[:space:]]/{printf "nev    %d\n", nev;  next}
                 /^iy[[:space:]]/ {printf "iy     %d\n", seed; next}
                 {print}' "$WORKDIR/$BORN_INPUT" > born_input.dat

            ./diffrad_gen.exe -input born_input.dat -lund born_events.lund \
                > born_out.txt 2>&1

            if [ "$BORN_ONLY" = false ]; then
                awk -v nev=$NEV -v seed=$SEED \
                    '/^nev[[:space:]]/{printf "nev    %d\n", nev;  next}
                     /^iy[[:space:]]/ {printf "iy     %d\n", seed; next}
                     {print}' "$WORKDIR/$RC_INPUT" > rc_input.dat

                ./diffrad_gen.exe -input rc_input.dat -lund rc_events.lund \
                    > rc_out.txt 2>&1
            fi

        ) &
        pids+=($!)
        echo "  Job $i started (PID $!), seed=$SEED"
    done

    echo ""
    echo "Waiting for all $NJOBS jobs ..."
    for pid in "${pids[@]}"; do
        wait "$pid"
    done
    echo "All jobs complete."
    echo ""

    echo "Combining LUND files ..."
    cat "$LUND_DIR"/job_*/born_events.lund > "$LUND_BORN"
    echo "  $LUND_BORN : $(wc -l < "$LUND_BORN") lines"
    if [ "$BORN_ONLY" = false ]; then
        cat "$LUND_DIR"/job_*/rc_events.lund > "$LUND_RC"
        echo "  $LUND_RC : $(wc -l < "$LUND_RC") lines"
    fi
    echo ""
fi

python3 - << PYEOF
import glob, os, sys

lund_dir  = "$LUND_DIR"
lund_born = "$LUND_BORN"
lund_rc   = "$LUND_RC"
born_only = "$BORN_ONLY" == "true"

def read_stat(fname):
    d = {}
    with open(fname) as f:
        for line in f:
            if '=' in line:
                k, v = line.split('=', 1)
                d[k.strip()] = v.split()[0]
    return d

def write_stat(fname, ngen, nsoft, nhard, sigma, sigma_err):
    hard_frac = nhard / max(1, ngen)
    with open(fname, 'w') as f:
        f.write(f'ngen      = {ngen:8d}\n')
        f.write(f'nsoft     = {nsoft:8d}\n')
        f.write(f'nhard     = {nhard:8d}\n')
        f.write(f'hard_frac =  {hard_frac:.4f}\n')
        f.write(f'sigma_nb  =  {sigma:.4e}\n')
        f.write(f'sigma_err =  {sigma_err:.4e}\n')

all_dirs = sorted(glob.glob(f'{lund_dir}/job_*/'))
born_ok, rc_ok = [], []
for d in all_dirs:
    bf = os.path.join(d, 'born_events_stat.dat')
    rf = os.path.join(d, 'rc_events_stat.dat')
    if os.path.exists(bf) and os.path.getsize(bf) > 0:
        born_ok.append(bf)
    if os.path.exists(rf) and os.path.getsize(rf) > 0:
        rc_ok.append(rf)

print(f"  Jobs with Born stat : {len(born_ok)} / {len(all_dirs)}")
print(f"  Jobs with RC stat   : {len(rc_ok)} / {len(all_dirs)}")

if not born_ok:
    print("  No completed jobs found — nothing to summarise.")
    sys.exit(0)

born_stats  = [read_stat(f) for f in born_ok]
ngen_born   = sum(int(s['ngen'])  for s in born_stats)
nsoft_born  = sum(int(s.get('nsoft','0')) for s in born_stats)
nhard_born  = sum(int(s.get('nhard','0')) for s in born_stats)
sigmas_born = [float(s['sigma_nb']) for s in born_stats]
sigma_born  = sum(sigmas_born) / len(sigmas_born)
err_born    = sigma_born / ngen_born ** 0.5
stat_born   = lund_born.replace('.lund', '_stat.dat')
write_stat(stat_born, ngen_born, nsoft_born, nhard_born, sigma_born, err_born)
print()
print(f"  Born : ngen={ngen_born:8d}  sigma={sigma_born:.4e} +/- {err_born:.4e} nb")

if not born_only and rc_ok:
    rc_stats  = [read_stat(f) for f in rc_ok]
    ngen_rc   = sum(int(s['ngen'])  for s in rc_stats)
    nsoft_rc  = sum(int(s.get('nsoft','0')) for s in rc_stats)
    nhard_rc  = sum(int(s.get('nhard','0')) for s in rc_stats)
    sigmas_rc = [float(s['sigma_nb']) for s in rc_stats]
    sigma_rc  = sum(sigmas_rc) / len(sigmas_rc)
    err_rc    = sigma_rc / ngen_rc ** 0.5
    stat_rc   = lund_rc.replace('.lund', '_stat.dat')
    write_stat(stat_rc, ngen_rc, nsoft_rc, nhard_rc, sigma_rc, err_rc)

    vfiles = sorted(glob.glob(f'{lund_dir}/job_*/rc_events_vdist.dat'))
    nv = 0
    vdist_out = lund_rc.replace('.lund', '_vdist.dat')
    with open(vdist_out, 'w') as out:
        for vf in vfiles:
            with open(vf) as inp:
                content = inp.read()
                out.write(content)
                nv += content.count('\n')

    print(f"  RC   : ngen={ngen_rc:8d}  sigma={sigma_rc:.4e} +/- {err_rc:.4e} nb")
    print(f"         nhard={nhard_rc} ({100*nhard_rc/max(1,ngen_rc):.2f}%)")
    print(f"  eta  = sigma_RC / sigma_Born = {sigma_rc/sigma_born:.4f}")
    print(f"  vdist: {nv} hard events from {len(vfiles)} jobs")
PYEOF

echo ""
echo "Done. Now run:  python plot_compare.py -born \"$LUND_BORN\" -rc \"$LUND_RC\""
echo "                python plot_eta.py     -born \"$LUND_BORN\" -rc \"$LUND_RC\""
