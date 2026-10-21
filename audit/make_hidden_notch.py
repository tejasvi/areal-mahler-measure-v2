#!/usr/bin/env python3
"""make_hidden_notch.py -- hidden negative control for the areal LOWER format: remove a flat-bottomed notch (trapezoid:
ramps of width wt, flat part of half-width a) from G at x0, so that F < 0 at the notch centre, while the certificate is
otherwise untouched (m, c, polynomials, x_lo, x_hi unchanged). The notch depth is tuned by bisection so that F_new(x0) < 0
but F_new > 0 at every node within 0.01 of x0, at the 4 inserted nodes and (option grid) at the points of the 2e-5 grid
within 1e-4. F_new = F + 2 U_T + (I_new/s_new - I/s), I_new = I - 2 sum C^G C^T J + sum C^T C^T J, s_new = s - mass
(exact closed forms).

    python3 -I make_hidden_notch.py CERT X0_NUM X0_DEN A_NUM A_DEN WT_NUM WT_DEN OUT.json [grid]"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_arl_second as M
from flint import arb, fmpq, ctx

SRC = sys.argv[1]
x0 = fmpq(int(sys.argv[2]), int(sys.argv[3]))
a = fmpq(int(sys.argv[4]), int(sys.argv[5]))
wt = fmpq(int(sys.argv[6]), int(sys.argv[7]))
OUT = sys.argv[8]
GRID = len(sys.argv) > 9 and sys.argv[9] == "grid"
M.setup(SRC)
e0, G0, N = list(M.D["e"]), list(M.D["G"]), M.D["N"]
s0 = M.D["s"]
ctx.prec = 192
I0 = M.energy_lower()
ctx.prec = 160


def interp(x):
    lo, hi = 0, N
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if e0[mid] <= x:
            lo = mid
        else:
            hi = mid
    k = min(lo, N - 1)
    return G0[k] + (G0[k + 1] - G0[k]) * (x - e0[k]) / (e0[k + 1] - e0[k])


pts = [x0 - a - wt, x0 - a, x0 + a, x0 + a + wt]
for p in pts:
    assert all(p != x for x in e0), "notch node coincides with a node"


def tent_jumps(ht):
    # T = trapezoid with value 0, ht, ht, 0 at pts: slope jumps ht/wt, -ht/wt, -ht/wt, ht/wt
    sl = ht / wt
    return [sl, -sl, -sl, sl]


def F_old(xq):
    x = arb(xq)
    val = -2 * M.U_real(x)
    for j in M.D["W"]:
        val -= M.D["wa"][j] * abs(M.D["Q"][j](x)).log()
    val += (x * x - 1) / 2 if xq <= 1 else x.log()
    return val - arb(M.D["m"])


def U_T(xq, CT):
    x = arb(xq)
    return sum((arb(c) * M.H_ball(arb(p) - x) for c, p in zip(CT, pts)), arb(0))


def model(ht):
    CT = tent_jumps(ht)
    mass = ht * (2 * a + wt)
    cross = arb(0)
    for Ci, Ei in zip(M.D["C"], e0):
        if Ci == 0:
            continue
        for c, p in zip(CT, pts):
            cross += arb(Ci) * arb(c) * M.J_ball(arb(Ei - p))
    tt = sum((arb(ci) * arb(cj) * M.J_ball(arb(pi - pj)) for ci, pi in zip(CT, pts) for cj, pj in zip(CT, pts)), arb(0))
    Inew = I0 - 2 * cross + tt
    snew = s0 - mass
    return CT, mass, Inew, snew


checks = [x for x in e0 if abs(x - x0) <= fmpq(1, 100)] + pts
if GRID:
    k0 = int((x0 * 50000).p // (x0 * 50000).q)
    checks += [fmpq(kk, 50000) for kk in range(k0 - 5, k0 + 6)]
Fold = {x: F_old(x) for x in checks + [x0]}


def F_new(x, CT, Inew, snew):
    return Fold[x] + 2 * U_T(x, CT) + Inew / arb(snew)


Gx0 = interp(x0)
lo, hi = fmpq(0), Gx0 * fmpq(99, 100)      # notch height bracket (G must stay >= 0)
best = None
for it in range(40):
    ht = (lo + hi) / 2
    CT, mass, Inew, snew = model(ht)
    fc = F_new(x0, CT, Inew, snew)
    mn = min(((x, F_new(x, CT, Inew, snew)) for x in checks), key=lambda t: float(t[1].mid()))
    if fc < 0 and mn[1] > 0:
        best = (ht, mass, fc, mn, Inew, snew)
        break
    if fc < 0:
        hi = ht
    else:
        lo = ht
if not best:
    print("no hidden notch found; last ht", float(ht), fc, mn)
    sys.exit(1)
ht, mass, fc, mn, Inew, snew = best
d = json.load(open(SRC))
newG = {x: interp(x) for x in pts}
newG[pts[1]] -= ht
newG[pts[2]] -= ht
pairs = sorted(list(zip(e0, G0)) + [(x, newG[x]) for x in pts], key=lambda t: t[0])
# nodes inside the flat part or the ramps (none: the notch is inside one node interval) keep their values
assert all(g >= 0 for _, g in pairs)
d["nodes"] = [str(p[0]) for p in pairs]
d["G"] = [str(p[1]) for p in pairs]
json.dump(d, open(OUT, "w"))
print(f"notch at x0 = {float(x0):.9f} (G = {float(Gx0):.4g}), flat half-width {float(a):.1e}, ramps {float(wt):.1e}, height {float(ht):.4g}, "
      f"mass {float(mass):.4e}: F_new(x0) = {fc.str(5)} < 0; min over {len(checks)} check points = {mn[1].str(5)} at {float(mn[0]):.7f} > 0; "
      f"I/s change {float((Inew / arb(snew) - I0 / arb(s0)).mid()):+.2e}; file {OUT}")
