# Per-certificate summary: polynomials as Table-2 row numbers (image transcription), weights, mass s, degree budget,
# X / tails and node ranges.  Usage: cert_summary_arl.py CERT.json TRANSCRIPTION.json
import json, sys
from fractions import Fraction as Fr
c = json.load(open(sys.argv[1])); T = json.load(open(sys.argv[2]))
idx = {tuple(t['desc']): t['j'] for t in T}
nodes = [Fr(x) for x in c['nodes']]; G = [Fr(g) for g in c['G']]
s = sum((nodes[k + 1] - nodes[k])*(G[k] + G[k + 1])/2 for k in range(len(nodes) - 1))
W = c.get('c') if 'c' in c else c.get('b'); kind = 'lower' if 'c' in c else 'upper'
rows = sorted(((idx.get(tuple(q), 'NOT IN TABLE 2'), Fr(w), len(q) - 1) for q, w in zip(c['polys'], W)), key=lambda t: (isinstance(t[0], str), t[0]))
degsum = sum(w*d for _, w, d in rows)
first = next(nodes[k] for k in range(len(G)) if G[k] > 0); last = max(nodes[k] for k in range(len(G)) if G[k] > 0)
print(f"certificate {sys.argv[1]} ({c['format']}, {kind}): {len(rows)} polynomials, {len(nodes)} nodes on [{nodes[0]}, {nodes[-1]}] "
      f"(G > 0 on [{first}, {last}]), mass s = {float(s):.10f}")
if kind == 'lower':
    print(f"  claim log lambda_A+ >= m = {c['m']}; budget sum c_j deg Q_j + 2 s = {float(degsum + 2*s):.6f} < 1; x_lo = {c['x_lo']}, x_hi = {c['x_hi']}")
else:
    print(f"  claim log lambda_A+ <= v = {c['v']}; (U1) sum b_j deg Q_j = {float(degsum):.6f} < s = {float(s):.6f}; X = {c['X']}")
print("  Table-2 row j : weight (float) : degree")
for j, w, d in rows: print(f"    {j} : {float(w):.10e} : {d}")
