"""
lund_reader.py  --  Universal LUND file reader for DIFFRAD / gagrho formats.

Particle identification by PID (column 3) and kinematics — independent of
particle ordering in the file.

PIDs used:
    11    electron (beam or scattered)
    2212  proton   (target at rest, or recoil)
    211   pi+
   -211   pi-
    321   K+
   -321   K-
     22   photon
    113   rho0 (DIFFRAD, type=0 — may be absent from GEANT4 tracking)
    111   rho0 (gagrho)
    333   phi
    443   J/psi

Particle line columns (same in all known formats):
    p[0]=idx  p[1]=status  p[2]=status2  p[3]=PID
    p[4]=parent  p[5]=parent2
    p[6]=px  p[7]=py  p[8]=pz  p[9]=E  p[10]=mass

Beam energy is at tok[6] in the event header.
"""

import numpy as np

Mp = 0.938272


def _v4(d):
    """Return [E, px, py, pz] numpy array from particle dict."""
    return np.array([d['E'], d['px'], d['py'], d['pz']])


def read_lund_events(fname, max_ev=0):
    """
    Read all events from a LUND file.  Returns a list of event dicts, each with:
        ebeam   : float  beam energy (GeV)
        parts   : list of dicts {pid, E, px, py, pz}
        -- computed kinematics --
        Q2, xB, nu, y, t, W, Eprime
        -- identified 4-vectors (np arrays [E,px,py,pz]) --
        k1      : beam electron
        k2      : scattered electron
        pp      : recoil proton
        pip     : pi+  (or positron for J/psi)
        pim     : pi-  (or electron for J/psi)
        kgam    : radiated photon (zeros if absent)
        has_photon : bool
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
            except:
                continue

            parts = []
            for _ in range(npart):
                p = next(f, '').split()
                if len(p) >= 10:
                    try:
                        parts.append(dict(pid   =int(p[3]),
                                          parent=int(p[4]),
                                          E  =float(p[9]),
                                          px =float(p[6]),
                                          py =float(p[7]),
                                          pz =float(p[8])))
                    except:
                        pass
            if len(parts) < 4:
                continue

            # ── identify particles ────────────────────────────────────────────
            electrons = [p for p in parts if p['pid'] ==   11]
            positrons = [p for p in parts if p['pid'] ==  -11]
            protons   = [p for p in parts if p['pid'] == 2212]
            photons   = [p for p in parts if p['pid'] ==   22]
            pips      = [p for p in parts if p['pid'] ==  211]
            pims      = [p for p in parts if p['pid'] == -211]
            kps       = [p for p in parts if p['pid'] ==  321]
            kms       = [p for p in parts if p['pid'] == -321]

            if len(electrons) < 1 or len(protons) < 1:
                continue

            # ── beam / scattered electron ─────────────────────────────────────
            # Use parent==0 to separate primary (beam/scattered) from decay e-
            primary_em = [p for p in electrons if p['parent'] == 0]
            decay_em   = [p for p in electrons if p['parent'] != 0]

            if len(primary_em) >= 2:
                # Both beam and scattered e- are in the file (DIFFRAD/gagrho)
                beam_e = min(primary_em, key=lambda p: abs(p['E'] - ebeam))
                scat_es = [p for p in primary_em if p is not beam_e]
                scat_e  = max(scat_es, key=lambda p: p['E'])
            elif len(primary_em) == 1:
                # Only scattered e- in file; synthesize beam along z-axis
                scat_e = primary_em[0]
                beam_e = dict(pid=11, parent=-1, E=ebeam, px=0., py=0., pz=ebeam)
            else:
                # No parent info — fall back to energy proximity
                if len(electrons) >= 2:
                    beam_e  = min(electrons, key=lambda p: abs(p['E'] - ebeam))
                    scat_es = [p for p in electrons if p is not beam_e]
                    scat_e  = max(scat_es, key=lambda p: p['E'])
                else:
                    scat_e = electrons[0]
                    beam_e = dict(pid=11, parent=-1, E=ebeam, px=0., py=0., pz=ebeam)

            # recoil proton: highest |p| (target is at rest)
            recoil_p = max(protons, key=lambda p: p['px']**2+p['py']**2+p['pz']**2)

            # ── meson daughters ───────────────────────────────────────────────
            # pi+/K+ for rho/phi;  positron (e+) for J/psi -> e+e-
            pos_d = pips or kps or positrons
            # pi-/K- for rho/phi;  decay e- for J/psi -> e+e-
            neg_d = pims or kms or decay_em
            if not pos_d or not neg_d:
                continue

            k1   = _v4(beam_e)
            k2   = _v4(scat_e)
            pp   = _v4(recoil_p)
            pip  = _v4(pos_d[0])
            pim  = _v4(neg_d[0])
            kgam = _v4(photons[0]) if photons else np.zeros(4)

            # ── kinematics ────────────────────────────────────────────────────
            q  = k1 - k2
            Q2 = -(q[0]**2 - q[1]**2 - q[2]**2 - q[3]**2)
            nu = q[0]
            y  = nu / ebeam
            xB = Q2 / (2*Mp*nu) if nu > 0 else np.nan
            t  = Mp**2 + pp[0]**2 - pp[1]**2 - pp[2]**2 - pp[3]**2 - 2*Mp*pp[0]
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
    """
    Lightweight reader — returns (Q2, xB, t, W, Eprime) arrays only.
    """
    evs = read_lund_events(fname, max_ev=max_ev)
    if not evs:
        return tuple(np.array([]) for _ in range(5))
    Q2s    = np.array([e['Q2']    for e in evs])
    xBs    = np.array([e['xB']    for e in evs])
    ts     = np.array([e['t']     for e in evs])
    Ws     = np.array([e['W']     for e in evs])
    Eprimes= np.array([e['Eprime']for e in evs])
    return Q2s, xBs, ts, Ws, Eprimes
