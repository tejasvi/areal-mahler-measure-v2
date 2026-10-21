# "Narrow notch" negative control for an areal LOWER certificate (tpmahler-areal-lower-v1): remove an asymmetric tent of mass eps
# (nodes x0 - w, x0, x0 + 3w; w = 3e-6 by default; peak at x0) from G where G is large enough. Since F = g - m - sum c log|Q| - 2U + I/s, the
# removed mass raises U only locally and makes F dip inside (x0, x0 + 3w) (the potential of the tent is lowest right of its peak),
# while F stays > 0 at the tent's nodes and on a 2e-5 grid. eps is tuned from two exact trial evaluations of F (F is affine in eps up
# to O(eps^2)). x0 is taken in (0, 1) when possible, so the new term g'' = 1 of the areal checker is exercised.
# Usage: make_dip_arl_lo.py CERT.json OUT.json [X0 = a contact point of the dual, e.g. 23275/100000]
import os, json, sys, copy
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'checker', 'lower'))
import check_arl_lo_59a1c54c as CL
src, out = sys.argv[1], sys.argv[2]
base = json.load(open(src))
nodes = [Fr(x) for x in base['nodes']]; Gs = [Fr(g) for g in base['G']]
w = Fr(sys.argv[4]) if len(sys.argv) > 4 else Fr(3, 10**6)     # tent half-width scale (total width 4w < 2e-5)
def build(x0, eps):
    j = max(i for i in range(len(nodes)) if nodes[i] <= x0 - w)
    assert nodes[j] < x0 - w and x0 + 3*w < nodes[j + 1]
    Gat = lambda x: Gs[j] + (Gs[j + 1] - Gs[j])*(x - nodes[j])/(nodes[j + 1] - nodes[j])
    peak = eps/(2*w)
    newn = [x0 - w, x0, x0 + 3*w]; newg = [Gat(x0 - w), Gat(x0) - peak, Gat(x0 + 3*w)]
    assert min(newg) >= 0
    c = copy.deepcopy(base); c['nodes'] = [str(x) for x in nodes[:j + 1] + newn + nodes[j + 1:]]
    c['G'] = [str(g) for g in Gs[:j + 1] + newg + Gs[j + 1:]]
    return c
def Fvals(c, pts):
    tmp = out + '.tmp.json'; json.dump(c, open(tmp, 'w'))
    cert = CL.setup(CL.Cert(tmp))
    return [float(CL.F_exact(cert, p).mid()) for p in pts]
# candidate x0: a node interval midpoint in (0.2, 0.95) with large G (float scan of G only; F is evaluated exactly below)
# x0: a contact point of the dual (where F is near its minimum, ~ the margin) given on the command line; F must be small there,
# because the tent's peak eps/(2w) must stay below G (eps ~ F(x0)/27).
x0 = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(23275, 100000)
sint = Fr(1, 3)                                       # the tent potential is lowest near s = 1/3 (right of the peak)
pts = [x0, x0 + sint*w, x0 - Fr(1, 10**5), x0 + Fr(1, 10**5), x0 - Fr(2, 10**5), x0 + Fr(2, 10**5)]
F0 = Fvals(base, pts)
e1 = Fr(1, 10**8); F1 = Fvals(build(x0, e1), pts)
D = [(a - b)/float(e1) for a, b in zip(F1, F0)]       # dF/deps at the points (affine model)
jj = max(i for i in range(len(nodes)) if nodes[i] <= x0)
print('G at x0 ~', float(Gs[jj] + (Gs[jj + 1] - Gs[jj])*(x0 - nodes[jj])/(nodes[jj + 1] - nodes[jj])))
print('x0 =', x0, '; F_old at node/interior/grid:', ['%.3e' % v for v in F0], '; slopes:', ['%.3e' % v for v in D])
lo = F0[1]/(-D[1]); hi = F0[0]/(-D[0])               # need F0[1] + eps D[1] < 0 < F0[0] + eps D[0]
assert D[1] < D[0] < 0 and lo < hi, (lo, hi)
eps = Fr(str('%.6e' % ((lo + hi)/2)))
c = build(x0, eps)
F2 = Fvals(c, pts)
print('eps =', float(eps), '; F_new at node, interior, grid points:', ['%.3e' % v for v in F2])
c['definition'] = base['definition'] + ' | NARROW-NOTCH NEGATIVE CONTROL (must NOT verify)'
json.dump(c, open(out, 'w'))
import os; os.remove(out + '.tmp.json')
print('wrote', out)
