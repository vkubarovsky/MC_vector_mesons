# Diffractive Vector Meson MC Generator with QED Radiative Corrections

Monte Carlo generator for exclusive electroproduction of vector mesons
(ρ, ω, φ, J/ψ) with QED radiative corrections based on the Akushevich
formalism (Eur. Phys. J. C8, 457, 1999).

## Source versions

| File | Status |
|---|---|
| `diffrad_vm.f90` | **The supported production code.** Combined generator for all mesons (`ivec` in input file); tuned June 2026 cross-section defaults (phi: exponential t, bt=1.284; J/psi: dipole t, mg2=3.112), parameters overridable via input-file keys `alf1 alf2 alf3 nuT cR bt mg2`; `sig_hard_fix` RC scheme; SDME code; 64-bit counters. |
| `diffrad_akushevich.f90` | Alternative physics model (Akushevich cross sections), kept for comparison via `run_all.sh -version akushevich`. |
| `diffrad_harut.f90` | Alternative physics model (Harut's cross sections), for comparison. |
| `diffrad_gen_akushevich.f90` | **Original Akushevich code — reference, do not edit or delete.** |

The tuning that produced the defaults lives in
[vector_mesons_generator_tuning](https://github.com/vkubarovsky/vector_mesons_generator_tuning);
its `combined/diffrad_vm.f90` is a frozen record of the tuning result —
development continues only here.

## Validation of the sig_hard_fix RC scheme (2026-06-12)

The hard-photon fix (signed Bardin–Shumeiko remainder `sig_F` without the
`max(0,·)` clamp; hard-event probability `(sig_F + ana_soft)/sig_total` with
the analytic soft logarithm over `[vcut_ir, vmax]`, `vcut_ir = 1e-2 GeV²`;
photon energy sampled by `sample_v_soft` from the `1/v` kernel weighted by
the Born suppression at the shifted `W'² = W² − v`) was validated against
independent runs (14 jobs × 1000 events, paired seeds):

**J/ψ** (Ebeam 10.6, W 4.05–4.55, alf2=4.105, mg2=3.106): this code gives
σ_Born = 269.2 ± 2.3 nb, σ_RC = 245.8 ± 2.1 nb, hard fraction 17.8% —
reproducing the standalone reference generator
(`jpsi/rc_test/sighard_fix`: 267.0 ± 1.9, 247.4 ± 1.7 nb, 17.9%) within
statistics. η(W) of the fixed generator independently matches the exact
`idiffrad` integration. Before the fix the J/ψ hard fraction was 0%
(negative qqt remainder clamped to zero ⇒ no radiated photons).

**φ** (Ebeam 10.6, RGA kinematics Q² 0.5–7): paired same-seed comparison
pre-fix vs post-fix gives Δσ_RC = −0.34% ± 1.2% — the total rate is
unchanged (the small shift is the removed clamp bias). The hard-event
fraction rises 11% → 40% by construction: events with photon energy
ω ≳ 5 MeV now carry an explicit radiated photon instead of being folded
into the non-radiated sample; the photon spectrum follows the analytic
`1/v` shape with threshold suppression.

Raw validation outputs (scripts, per-job statistics) live in the tuning
repo's local scratch `combined/validation/` (gitignored; conclusions here).

---

## Directory Structure

```
MC_vector_mesons/
├── diffrad_vm.f90          # Main Fortran source (shared by all mesons)
├── README.md                # This file
├── rho/                     # ρ meson working directory
│   ├── gen_input_born.dat   # Born-only input (ivec=1, cutv=0)
│   ├── gen_input_rc12.dat   # Full RC input (ivec=1, cutv=1.2 GeV²)
│   ├── run_diffrad.sh       # Single compile+run (quick test)
│   ├── run_parallel.sh      # Parallel production run
│   ├── run_eta.sh           # Legacy single-run script
│   ├── run_gen.sh           # Legacy compile+run script
│   ├── plot_compare.py      # Born vs RC validation plots
│   ├── plot_eta.py          # RC correction factor η vs Q²
│   ├── plot_gen.py          # Single LUND file diagnostic plots
│   └── lund -> ~/Downloads/DIFFRAD_lund/rho/   # symlink (not in OneDrive)
│
└── phi/                     # φ meson working directory
    ├── gen_input_born.dat   # Born-only input (ivec=3, cutv=0)
    ├── gen_input_rc.dat     # Full RC input (ivec=3, cutv=1.2 GeV²)
    ├── run_diffrad.sh       # Single compile+run (quick test)
    ├── run_parallel.sh      # Parallel production run
    ├── plot_compare.py      # Born vs RC validation plots
    ├── plot_eta.py          # RC correction factor η vs Q²
    ├── plot_gen.py          # Single LUND file diagnostic plots
    └── lund -> ~/Downloads/DIFFRAD_lund/phi/   # symlink (not in OneDrive)
```

LUND files are large and stored **outside OneDrive** at:
- `~/Downloads/DIFFRAD_lund/rho/`
- `~/Downloads/DIFFRAD_lund/phi/`

The `lund/` symlink in each meson directory points there transparently.

---

## Generator Source: `diffrad_vm.f90`

Single Fortran 90 source file implementing the full RC generator.
Compiled with:
```bash
gfortran -ffree-line-length-none -std=legacy -fwrapv -w -O2 diffrad_vm.f90 -o diffrad_gen.exe
```

### Key subroutines

| Subroutine | Purpose |
|---|---|
| `diffrad_gen` (main) | Accept/reject loop; Born + RC event generation |
| `build_4vectors` | Constructs all 4-momenta from (Ebeam, x, y, t, φ) |
| `sample_born` | Samples (Q², y, t, φ) from Born cross section |
| `sample_vrad` | Samples radiated v from podinl-weighted distribution |
| `podinl` | Akushevich radiator integrand; IR-divergent at small v |
| `sample_bw` | Breit-Wigner mass sampling for the vector meson |
| `decay_rho` | Isotropic decay V → π⁺π⁻ (or K⁺K⁻ for φ) in CM frame |
| `write_lund` | Writes one event in LUND format |
| `conkin` | Computes kinematic invariants (S, W², t_min, t_max, …) |
| `build_atm` | Builds the Akushevich tensor T_μν for a given (v, t, φ) |

### Input file format (14 lines)

```
27.5          ! bmom   — beam momentum (GeV)
0.0           ! tmom   — target momentum (0 = fixed target)
1             ! lepton — 1=electron
1             ! ivec   — 1=ρ, 2=ω, 3=φ, 4=J/ψ
1.2           ! cutv   — v_cut (GeV²); 0 = Born only
1000          ! nev    — number of events to generate
333522        ! seed   — random seed
0.40          ! q2min  (GeV²)
4.50          ! q2max  (GeV²)
0.54          ! ymin
0.56          ! ymax
0.109         ! tmin   (GeV²)
0.111         ! tmax   (GeV²)
0             ! iborn  — 0=full RC, 1=Born only
```

### Output files (written to same directory as exe)

| File | Content |
|---|---|
| `<lund>.lund` | Events in LUND format (6 particles Born, 7 with photon) |
| `<lund>_stat.dat` | Summary: ngen, nsoft, nhard, hard_frac, σ_Born, σ_err |
| `<lund>_vdist.dat` | Sampled v values for each hard event (GeV²) |

### LUND format

Each event: one header line + one line per particle.

Header: `npart  1  1  0  0  0  0  0  0  iev`

Particle columns: `idx  PID  px  py  pz  E  mass`

Particle order: e⁻(beam), p(target), e⁻(scattered), π⁺, π⁻, p(recoil), [γ]

PID codes: 11=e⁻, 2212=p, 211=π⁺, -211=π⁻, 321=K⁺, -321=K⁻, 22=γ

---

## Physics: Radiative Corrections

### Inelasticity variable v

In the Akushevich formalism the inelasticity variable is:
```
v = M_X² − M_p²  =  2·M_p·E_γ  (lab frame, proton at rest)
```

This gives the three equivalent representations:
```
E_γ  = v / (2·M_p)
M_X  = sqrt(v + M_p²)
v    = M_X² − M_p²
```

For cutv = 1.2 GeV²:
- E_γ,max = 1.2 / (2 × 0.938) ≈ 640 MeV
- M_X,max = sqrt(0.938² + 1.2) ≈ 1.44 GeV

### ISR vs FSR

Hard events are classified as ISR (photon from beam electron) or FSR
(photon from scattered electron) with probability:

```
P(ISR) = E₂² / (E₁² + E₂²)
```

**ISR**: photon collinear with beam, k₁ restored to full beam energy in LUND output.
**FSR**: photon collinear with scattered electron, k₂ scaled by (E₂−ω)/E₂.

4-momentum is conserved exactly in both cases.

### Bugs fixed (April 2026)

**Bug #1 — Wrong photon energy formula**
- Wrong:  `omega = vrad / (2·Ebeam)`  (beam energy in denominator)
- Fixed:  `omega = vrad / (2·amp)`    (proton mass in denominator)
- Impact: photon energy was wrong by factor Ebeam/Mp ≈ 29

**Bug #2 — Wrong vrad sampling distribution**
- Wrong: log-flat sampling (∝ 1/v), then accept/reject on full weight
- Fixed: dedicated `sample_vrad` subroutine using podinl-weighted
  accept/reject for efficient sampling of the physical distribution

**Bug #3 — Fortran continuation style**
- Several lines in `sample_vrad` used old fixed-form `&` at column 6
  instead of free-form `&` at end of line; caused compilation errors

---

## Shell Scripts

### `run_diffrad.sh`

Quick single-job compile and run. Use for testing and validation.

```bash
./run_diffrad.sh [options]

Options:
  -born-input FILE   Born input file   (default: gen_input_born.dat)
  -rc-input   FILE   RC input file     (default: gen_input_rc12.dat / gen_input_rc.dat)
  -born-lund  FILE   Born LUND output  (default: born_events.lund)
  -rc-lund    FILE   RC LUND output    (default: rc_events.lund)
  -no-compile        Skip compilation
```

Compiles `../diffrad_vm.f90`, runs Born-only then full-RC, prints statistics.

**Quick test:**
```bash
./run_diffrad.sh       # uses defaults from gen_input_*.dat (1000 events)
```

---

### `run_parallel.sh`

Production parallel runner. Launches N independent jobs, then combines output.

```bash
./run_parallel.sh [options]

Options:
  -njobs       N    Number of parallel jobs   (default: 14)
  -nev         N    Events per job            (default: 100000)
  -born-input FILE  Born input template       (default: gen_input_born.dat)
  -rc-input   FILE  RC input template         (rho default: gen_input_rc12.dat)
                                              (phi default: gen_input_rc.dat)
  -no-compile       Skip compilation
  -stats-only       Only collect statistics from existing job_* dirs —
                    no compilation, no jobs launched, no LUND combining.
                    Works with however many jobs have already finished.
```

**What it does (normal run):**
1. Creates `~/Downloads/DIFFRAD_lund/{rho,phi}/` and `lund/` symlink
2. Compiles `../diffrad_vm.f90`
3. Spawns N background jobs, each in `~/Downloads/DIFFRAD_lund/{rho,phi}/job_N/`
4. Each job runs Born first, then RC sequentially (rc files appear after Born finishes)
5. Each job gets a unique seed: `seed = 333522 + N × 100003`
6. Waits for all jobs, then combines LUND files with `cat`
7. Scans finished `job_*/` dirs, sums ngen/nsoft/nhard, averages σ, writes
   combined `born_events_stat.dat`, `rc_events_stat.dat`, `rc_events_vdist.dat`

**Combined stat file fields:**
```
ngen      — total events generated (sum across jobs)
nsoft     — non-radiated events
nhard     — hard-radiated events
hard_frac — nhard / ngen
sigma_nb  — cross section (average across jobs)
sigma_err — statistical error
```

**Quick test (1000 events, 1 job):**
```bash
./run_parallel.sh -njobs 1 -nev 1000
```

**Production run (1.4M events, 14 jobs):**
```bash
./run_parallel.sh
```

**Check progress while jobs are running:**
```bash
./run_parallel.sh -stats-only
```
Prints how many jobs have finished and the current summed statistics.
Safe to run at any time — does not interfere with running jobs.

**Recalculate statistics after all jobs finish (e.g. after fixing combining bug):**
```bash
./run_parallel.sh -stats-only
```

**Running on remote machine without losing jobs if SSH drops:**
```bash
ssh user@machine
screen -S rho
cd /path/to/rho && ./run_parallel.sh
# Ctrl+A D to detach; screen -r rho to reattach
```

---

### `run_eta.sh` / `run_gen.sh` (rho only, legacy)

Older single-job scripts kept for reference. Superseded by `run_diffrad.sh`
and `run_parallel.sh`. Use the newer scripts for all new work.

---

## Python Scripts

All scripts accept command-line arguments. Defaults point to `lund/` (the
symlink to `~/Downloads/DIFFRAD_lund/{rho,phi}/`).

---

### `plot_compare.py`

Validation plots comparing Born and RC samples side by side.
Each sample produces a 3×3 panel figure saved as `*_validation.png`.

```bash
python plot_compare.py [-born lund/born_events.lund] [-rc lund/rc_events.lund] [-nev N]
```

**Panels:**
- Q², y, −t distributions
- W, ρ/φ invariant mass, missing mass M_X = M(p_tar + k_γ)
- E_γ (hard events only), missing energy E_miss (should be 0), v distribution

**Key checks:**
- E_miss = 0 for all events → 4-momentum conservation
- M_X, E_γ, v mutually consistent: v = 2·M_p·E_γ = M_X² − M_p²

---

### `plot_eta.py`

Computes and plots the RC correction factor η = σ_RC / σ_Born vs Q².

```bash
python plot_eta.py [-born lund/born_events.lund] [-rc lund/rc_events.lund] [-nev N]
```

Reads LUND files and the corresponding `*_stat.dat` for absolute normalization.
Output: `eta_mc_vs_exact.png`

---

### `plot_gen.py`

Diagnostic plots for a single LUND file (typically the RC sample).

```bash
python plot_gen.py [-lund lund/rc_events.lund] [-nev N]
```

Output: `gen_validation.png`

---

## Typical Workflow

```bash
# 1. Quick test (rho)
cd rho
./run_parallel.sh -njobs 1 -nev 1000
python plot_compare.py         # check lund/born_events_validation.png, lund/rc_events_validation.png

# 2. Production run
./run_parallel.sh              # 14 × 100000 = 1.4M events
python plot_eta.py             # compute eta = sigma_RC / sigma_Born
python plot_compare.py         # full validation plots

# 3. Same for phi
cd ../phi
./run_parallel.sh -njobs 1 -nev 1000   # test
./run_parallel.sh                       # production
python plot_compare.py
```

---

## Physics Parameters (current runs)

| Parameter | Value |
|---|---|
| Beam energy | 27.5 GeV |
| Target | Fixed proton |
| ivec | 1 (ρ), 3 (φ) |
| Q² range | 0.4 – 4.5 GeV² |
| y range | 0.54 – 0.56 |
| \|t\| range | 0.109 – 0.111 GeV² |
| cutv | 1.2 GeV² |
| E_γ,max | ~640 MeV |
| M_X,max | ~1.44 GeV |
