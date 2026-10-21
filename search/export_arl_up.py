# AREAL version (objective g; format tpmahler-areal-upper-v1). Export an upper-bound LP solution (upper2_<tag>.npy = [val, len(e), e..., p..., b...]) to a tpmahler-upper-v1 JSON certificate.
# G_k = 2 p_k/(e_{k+1} - e_{k-1}) (unit-mass hats); nodes trimmed to the support (one zero at each end); only b_j > 0 kept;
# v = (1/s) int log+ dmu (float) + margin, rounded UP to a rational with denominator 10^12 (the checker recomputes it in arb).
import os, numpy as np, sys, json, math
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tpsets import load
src, setname, vmargin, out = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4]
a = np.load(src); ne = int(a[1]); e = a[2:2 + ne]; n = ne - 2; p = a[2 + ne:2 + ne + n]; b = a[2 + ne + n:]
polys = load(setname); assert len(polys) == len(b)
nodes = [Fr(str(round(float(x), 10))) for x in e]
assert all(nodes[i] < nodes[i + 1] for i in range(len(nodes) - 1)) and nodes[0] >= 0
G = [Fr(0)]*len(nodes)
for k in range(1, len(nodes) - 1):
    if p[k - 1] > 0:
        G[k] = Fr(float(p[k - 1]*2.0/(e[k + 1] - e[k - 1])))
pos = [k for k in range(len(G)) if G[k] > 0]
k0, k1 = max(pos[0] - 1, 0), min(pos[-1] + 1, len(G) - 1)
nodes, G = nodes[k0:k1 + 1], G[k0:k1 + 1]
keep = [j for j in range(len(polys)) if b[j] > 0]
s = sum((nodes[k + 1] - nodes[k])*(G[k] + G[k + 1])/2 for k in range(len(nodes) - 1))
lp = 0.0
for k in range(len(nodes) - 1):
    a_, b_ = float(nodes[k]), float(nodes[k + 1])
    beta = float((G[k + 1] - G[k])/(nodes[k + 1] - nodes[k])); alpha = float(G[k]) - beta*a_
    if a_ < 1:                                   # g = (x^2 - 1)/2 on [a, min(b, 1)]
        hi = min(b_, 1.0); P = lambda x: alpha*(x**3/6 - x/2) + beta*(x**4/8 - x**2/4)
        lp += P(hi) - P(a_)
    if b_ > 1:
        lo = max(a_, 1.0)
        F1 = lambda x: x*math.log(x) - x; F2 = lambda x: x*x/2*math.log(x) - x*x/4
        lp += alpha*(F1(b_) - F1(lo)) + beta*(F2(b_) - F2(lo))
v = lp/float(s) + vmargin
vq = Fr(math.ceil(v*10**12), 10**12)
X = max(Fr(5), Fr(math.ceil(float(nodes[-1]))))
cert = {'format': 'tpmahler-areal-upper-v1',
        'definition': 'mu = G(x) dx piecewise linear on nodes; h(x) = U_mu(x) - sum_j b_j log|Q_j(x)| >= 0 on [e_0, e_N]; Smith criterion (Lemma U1)',
        'polys': [polys[j] for j in keep], 'b': [str(Fr(float(b[j]))) for j in keep],
        'nodes': [str(x) for x in nodes], 'G': [str(g) for g in G], 'X': str(X), 'v': str(vq),
        'lp_value': float(a[0]), 'mass_float': float(s)}
json.dump(cert, open(out, 'w'))
print(f"{len(keep)} active polys, {len(nodes)} nodes on [{float(nodes[0]):.4f}, {float(nodes[-1]):.4f}], s = {float(s):.12f}, "
      f"sum b deg = {sum(float(b[j])*(len(polys[j]) - 1) for j in keep):.6f}, (1/s) int g = {lp/float(s):.10f}, v = {float(vq):.10f} "
      f"(lambda+ <= {math.ceil(math.exp(float(vq))*10**7)/10**7:.7f}, float display rounded up; check_up.py gives the rigorous value)")
