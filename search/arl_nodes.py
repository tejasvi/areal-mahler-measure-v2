# Adaptive nodes for the areal lower-bound primal: base grid on [0, X] plus fine grids around every positive real root of the
# polynomials weighted by a previous dual (and a node at each root).  Usage: arl_nodes.py SET DUAL.npy X BASE FINE WIN OUT.npy
import os, numpy as np, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tpsets import load
from flint import fmpz_poly
setname, dualf, X, base, fine, win, out = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), sys.argv[7]
P = load(setname); c = np.load(dualf)[3:]
roots = []
for j, q in enumerate(P):
    if c[j] > 1e-6 and q != [1, 0]:
        for z, m in fmpz_poly(q[::-1]).complex_roots():
            if z.imag.is_zero():
                r = float(z.real.mid())
                if 0.02 < r < X: roots.append(r)
roots = sorted(set(np.round(roots, 10)))
pts = set(np.round(np.arange(0.0, X + 1e-12, base), 10))
for r in roots:
    pts |= set(np.round(np.arange(max(r - win, fine), min(r + win, X) + 1e-12, fine), 10)); pts.add(round(r, 10))
e = np.array(sorted(x for x in pts if 0 <= x <= X)); keep = [0]
for i in range(1, len(e)):
    if e[i] - e[keep[-1]] >= 0.4*fine: keep.append(i)
if keep[-1] != len(e) - 1: keep[-1] = len(e) - 1
e = e[keep]; np.save(out, e)
print(f"{len(roots)} weighted roots, {len(e)} nodes on [0, {X}]; min gap {np.diff(e).min():.2e}")
