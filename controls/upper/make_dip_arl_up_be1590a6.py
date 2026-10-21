# AREAL version (check_arl_up / tpmahler-areal-upper-v1). "Narrow dip" negative control for an upper certificate: add a tiny asymmetric tent of mass eps to G (nodes x0 - w, x0,
# x0 + 3w; w = 1e-6: the tent's slope jumps eps/(2w^2) ~ 1e5 stay of the size of the certificate's own, so ball arithmetic does not
# blow up through cancellation (w = 1e-20 made every window fail for the wrong reason), and the tent (width 4e-6) still sits strictly
# between two check-grid points 2e-5 apart) at a point x0
# of the LP's 2e-5 check grid; x0 is chosen among the low-h grid points so that every (U2) integral stays positive (if x0 > 1 the
# claim v is raised by 1e-6, as in the POSITIVE control, so that only (U3) can fail). Its potential adds
# eps*(log w + phi((x - x0)/w)) near x0 (phi = potential of the unit tent on [-1, 3], peak at 0) and about eps*log|x - x0| far away.
# eps is chosen so that h stays > 0 at the tent's nodes AND at every point of the check grid, but h < 0 strictly inside (x0, x0 + 3w)
# (the minimum of phi is right of the peak). A grid pre-screen therefore passes; only an interval-wide argument can reject it.
# Usage: make_dip_up.py CERT.json OUT.json
import os, json, sys, math
from fractions import Fraction as Fr
import numpy as np, mpmath as mp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'checker', 'upper'))
import check_arl_up_a17af685 as CU
from flint import fmpz_poly
mp.mp.dps = 60
src, out = sys.argv[1], sys.argv[2]
base = json.load(open(src)); cert = CU.Cert(src)
# unit-mass asymmetric tent on [-1, 3], peak 1/2 at 0: slope jumps +1/2 at -1, -2/3 at 0, +1/6 at 3
Ct = [(mp.mpf(-1), mp.mpf(1)/2), (mp.mpf(0), -mp.mpf(2)/3), (mp.mpf(3), mp.mpf(1)/6)]
H = lambda t: (t*t/2*mp.log(abs(t)) - 3*t*t/4) if t != 0 else mp.mpf(0)
phi = lambda s: sum(c*H(e - s) for e, c in Ct)
ss = [mp.mpf(k)/1000 for k in range(1, 3000)]
vals = [phi(s) for s in ss]; kmin = min(range(len(vals)), key=lambda k: vals[k])
phi0, phimin, smin = phi(mp.mpf(0)), vals[kmin], ss[kmin]
print(f"phi(0) = {mp.nstr(phi0, 10)}, min phi = {mp.nstr(phimin, 10)} at s = {mp.nstr(smin, 6)}, phi(-1) = {mp.nstr(phi(-1), 10)}, phi(3) = {mp.nstr(phi(3), 10)}")
# float h on a 2e-5 grid of [e_0 + 0.01, 0.95]: choose the grid point (away from roots of weighted polys) with the smallest h
e = np.array([float(x) for x in cert.e]); G = np.array([float(g) for g in cert.G])
sl = np.diff(G)/np.diff(e); C = np.zeros(len(e)); C[0] = sl[0]; C[-1] = -sl[-1]; C[1:-1] = sl[1:] - sl[:-1]
def Uf(x, chunk=4000):                   # chunked: a dense (len(x), len(e)) matrix at once does not fit in memory
    out = np.empty(len(x))
    for i in range(0, len(x), chunk):
        t = e[None, :] - x[i:i + chunk, None]; at = np.abs(t)
        with np.errstate(divide='ignore', invalid='ignore'):
            Hm = np.where(at > 0, t*t/2*np.log(np.where(at > 0, at, 1.0)) - 0.75*t*t, 0.0)
        out[i:i + chunk] = Hm @ C
    return out
roots = np.concatenate([np.array([complex(z) for z, m in fmpz_poly(q[::-1]).complex_roots()]) for q, _, _, _ in cert.active])
rr = roots[np.abs(roots.imag) < 1e-12].real
xg = np.arange(float(cert.lo) + 0.01, 0.95, 2e-5)
xg = xg[np.min(np.abs(xg[:, None] - rr[None, :]), axis=1) > 2e-3]
def hf(x):
    v = Uf(x)
    for q, l, _, _ in cert.active:
        v -= float(l.mid())*np.log(np.abs(np.polyval(np.array(q, float), x)))   # low degree region check only; exact h below
    return v
# root products for log|Q| (Horner is unsafe for the high-degree Q): use FLINT roots
def logQ(q, x, chunk=50000):
    z = np.array([complex(zz) for zz, m in fmpz_poly(q[::-1]).complex_roots()])
    out = np.empty(len(x))
    for i in range(0, len(x), chunk):
        out[i:i + chunk] = np.log(abs(q[0])) + np.sum(np.log(np.abs(x[i:i + chunk, None] - z[None, :])), axis=1)
    return out
hg = Uf(xg) - sum(float(l.mid())*logQ(q, xg) for q, l, _, _ in cert.active)
# (U2) integrals of the base certificate (float, root products): s log|a| + sum_rho U(rho)
def Uz(z):
    t = e.astype(complex) - z
    with np.errstate(divide='ignore', invalid='ignore'):
        Hm = np.where(np.abs(t) > 0, t*t/2*np.log(np.where(np.abs(t) > 0, t, 1.0)) - 0.75*t*t, 0.0)
    return float(np.real(Hm @ C))
S = float(cert.s)
U2 = [S*np.log(abs(q[0])) + sum(Uz(complex(zz)) for zz, m in fmpz_poly(q[::-1]).complex_roots()) for q, _, _, _ in cert.active]
print('smallest (U2) integrals:', sorted(U2)[:4])
xg2 = np.arange(float(cert.lo) + 0.01, float(cert.hi) - 0.01, 2e-5)
xg2 = xg2[np.min(np.abs(xg2[:, None] - rr[None, :]), axis=1) > 2e-3]
hg2 = Uf(xg2) - sum(float(l.mid())*logQ(q, xg2) for q, l, _, _ in cert.active)
best = None
for kk in np.argsort(hg2)[:400]:
    xx = xg2[kk]; epsf = hg2[kk]/(13.8155 + 0.71)
    marg = min(U2[i] + epsf*logQ(q, np.array([xx]))[0] for i, (q, _, _, _) in enumerate(cert.active))
    if best is None or marg > best[0]: best = (marg, kk)
k = best[1]; xg, hg = xg2, hg2
print(f'chosen grid point {xg[k]:.5f} (h {hg[k]:.3e}); smallest (U2) after the tent ~ {best[0]:.3e}')
x0 = Fr(round(xg[k]*10**6) + 10, 10**6)      # 1e-5 to the right of a 1e-6-rational grid point
h0 = CU.h_exact(cert, x0); h0f = mp.mpf(str(float(h0.mid())))
print(f"x0 = {x0} = {float(x0)}: exact h(x0) in {h0}; grid min h = {hg[k]:.3e}")
w = Fr(1, 10**6)
lw = mp.log(mp.mpf(10)**-6)
eps = h0f/(-lw - (phi0 + phimin)/2)
epsq = Fr(str(mp.nstr(eps, 20)))
print(f"w = 1e-6, eps = {mp.nstr(eps, 6)}: expected h(x0) = {mp.nstr(h0f + eps*(lw + phi0), 4)}, min h = {mp.nstr(h0f + eps*(lw + phimin), 4)}, "
      f"grid drop ~ eps*log(1e-5) = {mp.nstr(eps*mp.log(mp.mpf('1e-5')), 4)}")
# insert nodes x0 - w, x0, x0 + 3w into the piecewise-linear G
nodes = [Fr(x) for x in base['nodes']]; Gs = [Fr(g) for g in base['G']]
j = max(i for i in range(len(nodes)) if nodes[i] <= x0 - w)
assert nodes[j] < x0 - w and x0 + 3*w < nodes[j + 1], "tent must sit inside one node interval"
def Gat(x): return Gs[j] + (Gs[j + 1] - Gs[j])*(x - nodes[j])/(nodes[j + 1] - nodes[j])
peak = epsq/(2*w)                       # mass of the tent = peak * (4w) / 2 = eps
newn = [x0 - w, x0, x0 + 3*w]; newg = [Gat(x0 - w), Gat(x0) + peak, Gat(x0 + 3*w)]
c = dict(base); c['nodes'] = [str(x) for x in nodes[:j + 1] + newn + nodes[j + 1:]]
c['G'] = [str(g) for g in Gs[:j + 1] + newg + Gs[j + 1:]]
c['definition'] = base['definition'] + ' | NARROW-DIP NEGATIVE CONTROL: tent of mass eps at x0 (must NOT verify)'
if x0 > 1:
    c['v'] = str(Fr(base['v']) + Fr(1, 10**6)); print('x0 > 1: claim v raised by 1e-6 (as in the POSITIVE control)')
json.dump(c, open(out, 'w'))
# float pre-screen of the perturbed h on the 2e-5 check grid near x0 (far-field tent potential eps*log|x - x0|)
xs = float(x0) - 1e-5 + 2e-5*np.arange(-500, 501)
hs = Uf(xs) - sum(float(l.mid())*logQ(q, xs) for q, l, _, _ in cert.active) + float(eps)*np.log(np.abs(xs - float(x0) - float(w)/3))
print(f"pre-screen on the check grid within +-0.01 of x0: min h = {hs.min():.3e} (> 0 means a grid test would pass); wrote {out}")
