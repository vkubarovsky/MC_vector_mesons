#!/usr/bin/env python3
"""
check_integral.py  --  Verify LUND ar_weight integral vs stat file sigma_nb,
                       and compute sigma for different v_max cuts.

The generator uses accept/reject (A/R) sampling:
    sigma = wmax * ngen / ntry   (stored in stat file as sigma_nb)

Each accepted event writes ar_weight = w(x), the full weight before the A/R test.
For accepted events w(x) is biased high (high-w events accepted more often), so:
    mean(ar_weight) != sigma_nb   in general.

However the RATIO of weight sums IS meaningful:
    sigma(v_max=V') / sigma(v_max=V) = wsum(v<V') / wsum(v<V)

because the ntry cancels in the ratio when it's the same run.

Usage:
    python check_integral.py [-born lund/born_events.lund] [-rc lund/rc_events.lund]
    python check_integral.py -rc lund/rc_events.lund
"""

import numpy as np
import os, argparse

parser = argparse.ArgumentParser()
parser.add_argument('-born', default='lund/born_events.lund')
parser.add_argument('-rc',   default='lund/rc_events.lund')
parser.add_argument('-nev',  default=0, type=int)
args = parser.parse_args()

def stat_name(lund):
    return lund[:-5]+'_stat.dat' if lund.endswith('.lund') else lund+'_stat.dat'

def vdist_name(lund):
    return lund[:-5]+'_vdist.dat' if lund.endswith('.lund') else lund+'_vdist.dat'

def read_stat(fname):
    d = {}
    if not os.path.exists(fname):
        return d
    with open(fname) as f:
        for line in f:
            if '=' in line:
                k, v = line.split('=', 1)
                d[k.strip()] = float(v.split()[0])
    return d

def read_lund_weights(fname, max_ev=0):
    """Read ar_weight (field 10) and npart from all header lines.
    Returns arrays: weights, has_photon (bool)."""
    weights = []
    hard    = []
    with open(fname) as f:
        for line in f:
            if max_ev > 0 and len(weights) >= max_ev:
                break
            tok = line.split()
            if not tok: continue
            try:
                npart = int(tok[0])
                if npart not in (6, 7): raise ValueError
                if len(tok) < 15: raise ValueError
                w = float(tok[9])
            except:
                continue
            weights.append(w)
            hard.append(npart == 7)
            for _ in range(npart):
                next(f, '')
    return np.array(weights), np.array(hard)

def read_vdist(fname):
    if os.path.exists(fname) and os.path.getsize(fname) > 0:
        return np.loadtxt(fname)
    return np.array([])

def analyse(lund_file, label):
    stat_file  = stat_name(lund_file)
    vdist_file = vdist_name(lund_file)

    if not os.path.exists(lund_file):
        print(f"\n[{label}] {lund_file} not found — skipping.\n")
        return

    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"  LUND : {lund_file}")
    print(f"  Stat : {stat_file}")
    print(f"{'='*60}")

    # ── stat file ──────────────────────────────────────────────
    stat = read_stat(stat_file)
    if stat:
        sigma_stat = stat.get('sigma_nb', float('nan'))
        ngen_stat  = int(stat.get('ngen',  0))
        nhard_stat = int(stat.get('nhard', 0))
        print(f"\n  From stat file:")
        print(f"    ngen       = {ngen_stat:8d}")
        print(f"    nhard      = {nhard_stat:8d}  ({100*nhard_stat/max(1,ngen_stat):.2f}%)")
        print(f"    sigma_nb   = {sigma_stat:.6e} nb")
    else:
        sigma_stat = float('nan')
        ngen_stat  = 0
        print(f"\n  Stat file not found or empty.")

    # ── LUND weights ───────────────────────────────────────────
    print(f"\n  Reading LUND weights ...")
    weights, hard = read_lund_weights(lund_file, args.nev)
    N      = len(weights)
    N_hard = hard.sum()
    N_soft = N - N_hard

    wsum       = weights.sum()
    wsum_soft  = weights[~hard].sum()
    wsum_hard  = weights[hard].sum()
    wmean      = wsum / max(1, N)
    wmax_obs   = weights.max()
    wmin_obs   = weights[weights > 0].min() if (weights > 0).any() else 0

    print(f"\n  From LUND file ({N} events):")
    print(f"    N_soft           = {N_soft:8d}")
    print(f"    N_hard           = {N_hard:8d}  ({100*N_hard/max(1,N):.2f}%)")
    print(f"    sum(ar_weight)   = {wsum:.6e} nb")
    print(f"    mean(ar_weight)  = {wmean:.6e} nb")
    print(f"    min(ar_weight)   = {wmin_obs:.6e} nb")
    print(f"    max(ar_weight)   = {wmax_obs:.6e} nb")

    print(f"\n  Comparison with stat file:")
    print(f"    sigma_nb  (stat) = {sigma_stat:.6e} nb")
    print(f"    mean(w)  (LUND)  = {wmean:.6e} nb  (ratio = {wmean/sigma_stat:.4f})")
    print(f"    sum(w)/N  = mean = same as above")

    # The A/R estimator: sigma = wmax_gen * ngen/ntry
    # We can estimate ntry from stat: ntry = wmax_gen * ngen / sigma_stat
    # But wmax_gen is unknown. We can infer: if sigma = mean(w)*ngen/ntry,
    # then ntry = mean(w)*ngen/sigma_stat, efficiency = sigma_stat/mean(w)
    eff_implied = sigma_stat / wmean if wmean > 0 else float('nan')
    ntry_implied = int(round(N / eff_implied)) if np.isfinite(eff_implied) else 0
    print(f"\n  Implied A/R efficiency = sigma_stat / mean(w) = {eff_implied:.4f}")
    print(f"  Implied ntry           = ngen / efficiency   = {ntry_implied:d}")
    print(f"  Check: wsum/ntry       = {wsum/max(1,ntry_implied):.6e} nb  "
          f"(should = sigma_stat = {sigma_stat:.6e} nb)")

    # ── v_max dependence from hard events ──────────────────────
    vv = read_vdist(vdist_file)
    if len(vv) > 0 and N_hard > 0:
        print(f"\n  v distribution: {len(vv)} hard events from vdist file")
        vmax_gen = vv.max()
        print(f"    v range: {vv.min():.4f} -- {vmax_gen:.4f} GeV²")

        # Cross section for different v_max cuts
        # sigma(v_max) = sigma_soft + sigma_hard(v < v_max)
        # Ratio = (wsum_soft + wsum_hard(v<cut)) / (wsum_soft + wsum_hard)
        # Scale by sigma_stat to get absolute value
        cuts = [0.2, 0.3, 0.5, 0.8, 1.0, 1.2, vmax_gen]
        print(f"\n  sigma(v_max) from LUND weights (normalised to stat sigma):")
        print(f"  {'v_max':>8}  {'N_hard':>7}  {'wsum_hard':>14}  "
              f"{'sigma_RC':>14}  {'ratio_to_full':>14}")

        w_hard = weights[hard]
        # Sort vv and w_hard together  (vv has same length as N_hard)
        if len(vv) == N_hard:
            idx = np.argsort(vv)
            vv_s = vv[idx]
            w_hard_s = w_hard[idx]
        else:
            # length mismatch — vdist may be from parallel jobs, use histogram approach
            vv_s, w_hard_s = None, None

        sigma_full = wsum_soft + wsum_hard  # unnormalised total

        for vcut in cuts:
            if vv_s is not None:
                mask = vv_s <= vcut
                nh_cut   = mask.sum()
                wh_cut   = w_hard_s[mask].sum()
            else:
                # Estimate from fraction of vdist histogram
                nh_cut  = (vv <= vcut).sum()
                frac_h  = nh_cut / max(1, len(vv))
                wh_cut  = wsum_hard * frac_h
            w_cut    = wsum_soft + wh_cut
            sigma_cut = sigma_stat * w_cut / sigma_full
            ratio    = w_cut / sigma_full
            print(f"  {vcut:8.2f}  {nh_cut:7d}  {wh_cut:14.4e}  "
                  f"{sigma_cut:14.6e}  {ratio:14.4f}")
    else:
        print(f"\n  No vdist file or no hard events — skipping v_max scan.")

analyse(args.born, "BORN sample")
analyse(args.rc,   "RC sample")
