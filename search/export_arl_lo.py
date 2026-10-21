# AREAL version (format tpmahler-areal-lower-v1). Export primalN (nodes, hat weights) + fixed dual (m, lam, E, c on a polynomial set) to a tpmahler-cert-v1 JSON.
import os, numpy as np, sys, json, math
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tpsets import load
prim, dualf, dualset, margin, out = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5]
a = np.load(prim); ne = int(a[1]); e = a[2:2 + ne]; p = a[2 + ne:]
d = np.load(dualf); m, lam, E = d[:3]; c = d[3:]
if dualset == 'same':
    P = load('F27'); idx = json.load(open(prim.replace('.npy', '_idx.json'))); polys = [P[j] for j in idx]
else:
    polys = load(dualset)
assert len(polys) == len(c)
keep = [j for j in range(len(polys)) if c[j] > 0 or polys[j] == [1, 0]]
nodes = [Fr(str(round(float(x), 10))) for x in e]
assert all(nodes[i] < nodes[i + 1] for i in range(len(nodes) - 1)) and nodes[0] >= 0
G = [Fr(0)]*len(nodes)
for k in range(1, len(nodes) - 1):
    if p[k - 1] > 0:
        G[k] = Fr(float(lam*p[k - 1]*2.0/(e[k + 1] - e[k - 1])))
first = next(nodes[k - 1] for k in range(1, len(nodes) - 1) if G[k] > 0)
x_lo = min(Fr(1, 1000), Fr(math.floor(float(first)*100), 100))
s = sum((nodes[k + 1] - nodes[k])*(G[k] + G[k + 1])/2 for k in range(len(nodes) - 1))
degsum = sum(float(c[j])*(len(polys[j]) - 1) for j in keep)
kappa = 1 - degsum - 2*float(s)
beta = float(m) + sum(float(c[j])*math.log(sum(abs(v) for v in polys[j])) for j in keep) - float(E)/float(lam)
x_hi = Fr(max(math.ceil(float(nodes[-1])) + 1, math.ceil(math.exp(beta/kappa + 2.0))))
cert = {'format': 'tpmahler-areal-lower-v1',
        'definition': 'F(x) = g(x) - m - sum_j c_j log|Q_j(x)| - 2 U(x) + I/s >= 0 on (0, inf); rho = G(x) dx piecewise linear on nodes',
        'polys': [polys[j] for j in keep], 'c': [str(Fr(float(c[j]))) for j in keep],
        'nodes': [str(x) for x in nodes], 'G': [str(g) for g in G],
        'm': str(Fr(float(m)) - Fr(margin).limit_denominator(10**12)),
        'x_lo': str(x_lo), 'x_hi': str(x_hi), 'm_float_dual': float(m), 'margin': margin}
json.dump(cert, open(out, 'w'))
print(f"polys {len(keep)}, nodes {len(nodes)}, s = {float(s):.6f} (lam {lam:.6f}), budget sum c deg + 2s = {degsum + 2*float(s):.4f}, "
      f"kappa {kappa:.4f}, beta {beta:.4f}, x_lo {x_lo}, x_hi {x_hi}, m_cert {float(Fr(cert['m'])):.10f} (Omega {math.exp(float(Fr(cert['m']))):.7f})")
