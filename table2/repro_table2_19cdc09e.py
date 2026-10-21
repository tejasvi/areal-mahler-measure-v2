# Reproduce Flammang's areal constant 1.380047 (J. Number Theory 2015, Thm 2) from a transcription of her Table 2.
# Usage: repro_table2.py TRANSCRIPTION.json   (a list of {"j", "c": printed weight, "desc": integer coefficients, leading first})
# f(x) = g(x) - sum_j c_j log|Q_j(x)|, g(x) = log+ x + (min(1, x)^2 - 1)/2, on a float grid of (0, 70]; log|Q_j| by root products
# (flint roots). Prints (a) min f with her printed weights, (b) the LP max m over weights c_j >= 0 with sum c_j deg Q_j <= 1 and
# g - sum c_j log|Q_j| >= m on the grid (cutting planes, HiGHS). Float only: a consistency check of the transcription, not a proof.
# This is the code run inline at 00:35 PDT (text transcription) and 02:08 PDT (image transcription) on 2026-10-21, saved as a file;
# the last two output lines (whole-grid feasibility of the final weights) were added at freeze time.
import json, sys, hashlib, numpy as np
from flint import fmpz_poly
from scipy.optimize import linprog
src = sys.argv[1]
T = json.load(open(src))
print(f"input {src} sha256 {hashlib.sha256(open(src, 'rb').read()).hexdigest()}, {len(T)} rows")
polys = [r['desc'] for r in T]; c = np.array([float(r['c']) for r in T])
def roots(q): return np.array([complex(z) for z, m in fmpz_poly(q[::-1]).complex_roots() for _ in range(int(m))])
R = [roots(q) for q in polys]
rr = np.concatenate([z[np.abs(z.imag) < 1e-12].real for z in R]); rr = rr[rr > 0]
x = np.unique(np.concatenate([np.geomspace(1e-6, 0.05, 4000), np.arange(0.05, 70, 2e-4)] + [r + np.geomspace(1e-8, 1e-3, 30)*s for r in rr for s in (1, -1)]))
x = x[(x > 0)]; x = x[np.min(np.abs(x[:, None] - rr[None, :]), axis=1) > 1e-10]
L = np.array([np.log(abs(q[0])) + np.sum(np.log(np.abs(x[:, None] - z[None, :])), axis=1) for q, z in zip(polys, R)]).T
g = np.where(x >= 1, np.log(np.maximum(x, 1)), (np.minimum(x, 1)**2 - 1)/2)
f = g - L @ c
print(f"grid {len(x)} points on [{x[0]:.1e}, {x[-1]:.4f}]")
print(f"printed weights: min f = {f.min():.7f} at x = {x[np.argmin(f)]:.5f}  (N0 >= {np.exp(f.min()):.7f}; her constant 1.380047 -> log {np.log(1.380047):.7f}); sum c deg = {np.array([len(q) - 1 for q in polys], float) @ c:.6f}")
degs = np.array([len(q) - 1 for q in polys], float); J = len(polys)
sel = np.arange(0, len(x), 10)
for it in range(40):
    A = np.hstack([np.ones((len(sel), 1)), L[sel]]); A = np.vstack([A, np.concatenate([[0], degs])[None, :]])
    res = linprog(np.concatenate([[-1.0], np.zeros(J)]), A_ub=A, b_ub=np.concatenate([g[sel], [1.0]]), bounds=[(None, None)] + [(0, None)]*J, method='highs')
    m, cc = res.x[0], res.x[1:]; F = g - L @ cc - m
    if F.min() > -1e-11: break
    bad = np.nonzero(F < -1e-11)[0]; sel = np.unique(np.concatenate([sel, bad[::max(1, len(bad)//2000)]]))
print(f"re-optimized LP over the transcribed set: m = {m:.7f} (N0 >= {np.exp(m):.7f}), {it+1} rounds, sum c deg = {degs @ cc:.4f}, weight on the last row: {cc[-1]:.3e}")
mf = m + min(F.min(), 0.0)                   # the last weights are feasible on the whole grid for m + min(F, 0)
print(f"  last round: min F on the whole grid = {F.min():.3e} ({'converged' if F.min() > -1e-11 else 'round cap reached'}); "
      f"value valid on the whole grid with the last weights: {mf:.7f} (N0 >= {np.exp(mf):.7f})")
