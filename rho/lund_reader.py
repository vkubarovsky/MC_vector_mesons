"""
lund_reader.py  --  Universal LUND file reader.

Particle identification by PID among GEANT-tracked particles only (N3==1).

Logic:
  Only particles with column p[2]==1 (GEANT tracking flag) are used.
  This automatically excludes the beam e- and the vector meson (both N3=0).

  Scattered e-   : pid=11, lowest energy among N3=1 electrons
  Recoil proton  : pid=2212, highest |p| among N3=1 protons
  pi+ / K+       : pid=211 / 321
  pi- / K-       : pid=-211 / -321
  J/psi decay e+ : pid=-11
  J/psi decay e- : pid=11, the N3=1 electron that is not the scattered one
  Photon         : pid=22

Beam energy from header tok[6].  Particle columns:
  p[0]=idx  p[1]=type  p[2]=N3(geant)  p[3]=PID  p[4]=parent  p[5]=parent2
  p[6]=px   p[7]=py    p[8]=pz         p[9]=E    p[10]=mass
"""

import numpy as np
import os

Mp = 0.938272


def resolve_lund_dir(name):
    """Return the actual directory for a lund_VERSION pointer.
    Works for a real directory, a symlink, or an OneDrive text file
    containing the target path."""
    if os.path.isdir(name):
        return name
    if os.path.isfile(name):
        return open(name).read().strip()
    return name


def _v4(d):
    return np.array([d['E'], d['px'], d['py'], d['pz']])


def read_lund_events(fname, max_ev=0):
    """
    Read all events from a LUND file.  Returns a list of event dicts:
        ebeam, parts
        Q2, xB, nu, y, t, W, Eprime
        k1, k2, pp, pip, pim, kgam   (numpy [E,px,py,pz])
        has_photon
    """
    events = []
    with open(fname) as f:
        for line in f:
            if max_ev > 0 and len(events) >= max_ev:
                break
            tok = line.split()
            if not tok:
                continue
            try:
                npart = int(tok[0])
                if npart < 4:
                    raise ValueError
                ebeam = float(tok[6])
            except Exception:
                continue

            parts = []
            for _ in range(npart):
                p = next(f, '').split()
                if len(p) >= 10:
                    try:
                        parts.append(dict(
                            pid   =int(p[3]),
                            geant =int(p[2]),
                            parent=int(p[4]),
                            E  =float(p[9]),
                            px =float(p[6]),
                            py =float(p[7]),
                            pz =float(p[8])))
                    except Exception:
                        pass
            if len(parts) < 4:
                continue

            # ── only use GEANT-tracked particles (N3==1) ──────────────────────
            # This excludes beam e- and vector meson (both have N3=0)
            active    = [p for p in parts if p['geant'] == 1]
            electrons = [p for p in active if p['pid'] ==   11]
            protons   = [p for p in active if p['pid'] == 2212]
            photons   = [p for p in active if p['pid'] ==   22]
            pips      = [p for p in active if p['pid'] ==  211]
            pims      = [p for p in active if p['pid'] == -211]
            kps       = [p for p in active if p['pid'] ==  321]
            kms       = [p for p in active if p['pid'] == -321]
            positrons = [p for p in active if p['pid'] ==  -11]

            # beam: always synthesized from header
            k1 = np.array([ebeam, 0., 0., ebeam])

            # scattered e-: N3=1, PID=11, parent==0 (no parent)
            scat_cands = [p for p in electrons if p['parent'] == 0]
            if not scat_cands:
                scat_cands = electrons  # fallback
            if not scat_cands:
                continue
            scat_e = min(scat_cands, key=lambda p: p['E'])
            k2 = _v4(scat_e)

            # recoil proton: highest |p| (target proton is always at rest)
            if not protons:
                continue
            recoil_p = max(protons, key=lambda p: p['px']**2 + p['py']**2 + p['pz']**2)
            pp = _v4(recoil_p)

            # meson daughters: pi+/K+ and pi-/K- (or e+/e- for J/psi)
            pos_d = pips or kps or positrons
            neg_d = pims or kms
            # J/psi -> e+e-: decay e- has parent != 0 (parent = J/psi index)
            if not neg_d:
                neg_d = [p for p in electrons if p['parent'] != 0]
            if not pos_d or not neg_d:
                continue

            pip  = _v4(pos_d[0])
            pim  = _v4(neg_d[0])
            kgam = _v4(photons[0]) if photons else np.zeros(4)

            # ── kinematics from 4-vectors ──────────────────────────────────────
            q  = k1 - k2
            Q2 = -(q[0]**2 - q[1]**2 - q[2]**2 - q[3]**2)
            nu = q[0]
            y  = nu / ebeam
            xB = Q2 / (2*Mp*nu) if nu > 0 else np.nan
            # t from recoil proton: t = (pp - ptar)^2, ptar=[Mp,0,0,0]
            t  = 2*Mp*(Mp - pp[0])
            ptar = np.array([Mp, 0., 0., 0.])
            Wv   = ptar + q
            W2   = Wv[0]**2 - Wv[1]**2 - Wv[2]**2 - Wv[3]**2
            W    = np.sqrt(max(0., W2))

            events.append(dict(
                ebeam=ebeam, parts=parts,
                Q2=Q2, xB=xB, nu=nu, y=y, t=t, W=W,
                Eprime=scat_e['E'],
                k1=k1, k2=k2, pp=pp, pip=pip, pim=pim, kgam=kgam,
                has_photon=bool(photons),
            ))
    return events


def read_lund(fname, max_ev=0):
    """Lightweight reader — returns (Q2, xB, t, W, Eprime) arrays."""
    evs = read_lund_events(fname, max_ev=max_ev)
    if not evs:
        return tuple(np.array([]) for _ in range(5))
    return (np.array([e['Q2']     for e in evs]),
            np.array([e['xB']     for e in evs]),
            np.array([e['t']      for e in evs]),
            np.array([e['W']      for e in evs]),
            np.array([e['Eprime'] for e in evs]))
