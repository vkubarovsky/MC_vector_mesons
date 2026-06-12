#!/usr/bin/env python3
"""
check_wuu_norm.py
Numerically integrate W_UU over all decay angles (phi, kap, theta) for a
specific kinematic point and verify that
   int W_UU dOmega / (4*pi^2) == sigma_T + eps * sigma_L

using the loadPhysicsInputs parameterization from gagrho generator.h.
"""
import numpy as np

# ── kinematic point ───────────────────────────────────────────────────────────
Q2 = 2.0   # GeV^2
xB = 0.20
t  = -0.30 # GeV^2

# epsilon (virtual photon polarisation)
Mp  = 0.938272
E   = 10.6   # beam energy  GeV
nu  = Q2 / (2*Mp*xB)
y   = nu / E
gam2 = 4*Mp**2*xB**2/Q2
eps_num = 1 - y - 0.25*y**2*gam2
eps_den = 1 - y + 0.5*y**2 + 0.25*y**2*gam2
eps = eps_num / eps_den

# ── physics inputs (from loadPhysicsInputs) ───────────────────────────────────
bt, bl = 2.75, 4.25
Q = np.sqrt(Q2)
sqt = np.sqrt(max(-t, 0.0))

tdept = np.exp(bt*t)
tdepl = np.exp(bl*t)
xdept = xB / (1.0 - xB)
xdepl = xB

sigT = 30.0 * tdept * xdept / Q2    # nb/GeV^2
sigL = 30.0 * tdepl * xdepl          # nb/GeV^2

# helicity matrix (only non-zero elements set in generator.h)
# convention: u[nu][nup][mu][mup]  index: + -> 0, 0 -> 1, - -> 2
u = np.zeros((3,3,3,3), dtype=complex)

def h(c):
    return 0 if c=='+' else (1 if c=='0' else 2)

# LL part
u[h('0')][h('0')][h('0')][h('0')] =  30.0*tdepl*xdepl + 0j              # sigma_L
u[h('0')][h('0')][h('0')][h('+')] =  0 + 80.0*tdept*xdepl*sqt/Q * 1j   # sin u^00_0+
# LT part
u[h('0')][h('+')][h('+')][h('+')] = (25.0*tdepl*xdepl*sqt/Q +
                                      25.0*tdepl*xdepl*sqt/Q * 1j)
# TT part
u[h('+')][h('+')][h('+')][h('+')] =  30.0*tdept*xdept/Q2 + 0j           # sigma_T
u[h('+')][h('+')][h('0')][h('+')] =  0 - 10.0*tdept*xdept*sqt/Q * 1j   # sin u^++_0+

def U(nu, nup, mu, mup):
    return u[h(nu)][h(nup)][h(mu)][h(mup)]

def W_UU(phi, kap, theta, eps):
    """W_UU kernel from w_kernels.hpp Eq.(15), returning the LL/LT/TT components."""
    ce = np.sqrt(eps*(1.0+eps))
    cphi  = np.cos(phi);    c2phi  = np.cos(2*phi)
    cpk   = np.cos(kap)
    cphi_k  = np.cos(phi+kap);  cphi_mk = np.cos(phi-kap)
    c2phi_k = np.cos(2*phi+kap); c2phi_mk= np.cos(2*phi-kap)
    c2phi2k = np.cos(2*phi+2*kap)
    cm2k    = np.cos(2*kap)

    LL  = ( np.real(U('0','0','+','+'))
          + eps * np.real(U('0','0','0','0'))
          - 2.0*cphi*ce * np.real(U('0','0','0','+'))
          - eps*c2phi   * np.real(U('0','0','-','+')))

    LT  = ( cphi_k*ce * np.real(U('0','+','0','+') - U('-','0','0','+'))
          - cpk * np.real(U('0','+','+','+') - U('-','0','+','+')
                          + 2.0*eps*U('0','+','0','0'))
          + c2phi_k*eps * np.real(U('0','+','-','+'))
          - cphi_mk*ce * np.real(U('0','-','0','+') - U('+','0','0','+'))
          + c2phi_mk*eps* np.real(U('+','0','-','+')))

    cphi2k_p = np.cos(phi+2*kap)
    cphi2k_m = np.cos(phi-2*kap)
    TT  = ( 0.5*(np.real(U('+','+','+','+')) + np.real(U('-','-','+','+')))
          + 0.5*c2phi2k*eps * np.real(U('-','+','-','+'))
          - cphi*ce * np.real(U('+','+','0','+') + U('-','-','0','+'))
          + cphi2k_p*ce * np.real(U('-','+','0','+'))
          - cm2k * np.real(U('-','+','+','+') + eps*U('-','+','0','0'))
          - c2phi*eps * np.real(U('+','+','-','+'))
          + cphi2k_m*ce * np.real(U('+','-','0','+'))
          + 0.5*c2phi_mk*eps * np.real(U('+','-','-','+')))

    ct = np.cos(theta); st = np.sin(theta)
    return ct**2 * LL + np.sqrt(2)*ct*st * LT + st**2 * TT

# ── numerical integration over (phi, kap, theta) ─────────────────────────────
N = 200   # grid points per angle (N^3 total ~ 8M)
phi_arr = np.linspace(0, 2*np.pi, N, endpoint=False)
kap_arr = np.linspace(0, 2*np.pi, N, endpoint=False)
th_arr  = np.linspace(0, np.pi,   N, endpoint=False)

# vectorize over phi and kap for a fixed theta slice
integral = 0.0
dphi = 2*np.pi/N
dkap = 2*np.pi/N
dth  = np.pi/N

print(f"Integrating W_UU over {N}^3 = {N**3:,} points ...")
PHI, KAP = np.meshgrid(phi_arr, kap_arr, indexing='ij')
for theta in th_arr:
    w = W_UU(PHI, KAP, theta, eps)
    integral += np.sum(w) * np.sin(theta) * dphi * dkap * dth

# ── compare ───────────────────────────────────────────────────────────────────
target = sigT + eps * sigL
integral_norm = integral / (4*np.pi**2)

print()
print(f"Kinematic point: Q2={Q2}, xB={xB}, t={t}")
print(f"  nu={nu:.3f} GeV, y={y:.4f}, eps={eps:.4f}")
print(f"  sigT = {sigT:.6f}  nb/GeV^2")
print(f"  sigL = {sigL:.6f}  nb/GeV^2")
print(f"  sigT + eps*sigL = {target:.6f}  nb/GeV^2")
print()
print(f"  int W_UU dOmega              = {integral:.6f}")
print(f"  int W_UU dOmega / (4*pi^2)  = {integral_norm:.6f}")
print()
ratio = integral_norm / target if target > 0 else float('nan')
print(f"  [int W_UU dOmega/(4pi^2)] / (sigT + eps*sigL) = {ratio:.6f}")
if abs(ratio - 1.0) < 0.01:
    print("  --> NORMALIZED CORRECTLY (ratio = 1.00)")
else:
    print(f"  --> OFF by factor {ratio:.4f}  (expected 1.00)")
