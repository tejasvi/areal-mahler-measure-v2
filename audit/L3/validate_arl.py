#!/usr/bin/env python3
"""validate_arl.py -- tests the minorant of check_arl_second.py (both formats). The checker's interval test is run on
NSAMP random node intervals; on accepted sub-intervals (at most 1200, chosen at random) the function (F for the lower
format, h for the upper format) is evaluated by the FULL closed form (all nodes, no Taylor model) at 4 random points each
and must be >= the minorant (quadratic + chord) there.

    python3 -I validate_arl.py CERT.json [NSAMP [SEED]]      (defaults 300, 1)"""
import sys, os, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_arl_second as M
from flint import arb, fmpq, ctx

M.setup(sys.argv[1])
NSAMP = int(sys.argv[2]) if len(sys.argv) > 2 else 300
rng = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
e, N = M.D["e"], M.D["N"]
if M.D["mode"] == "lower":
    s = M.D["s"]; ctx.prec = 192; I = M.energy_lower(); ctx.prec = 128
    M.D["const"] = -arb(M.D["m"]) + I.lower() / arb(s)
M.D["brackets"] = []
for j in M.D["W"]:
    for (kind, r, lo, hi, mult) in M.D["roots"][j]:
        if kind == "r":
            M.D["brackets"].append((lo, hi, j))
M.D["brackets"].sort()

def full(xq):
    x = arb(xq)
    val = M.D["sigma"] * M.U_real(x)
    for j in M.D["W"]:
        val -= M.D["wa"][j] * abs(M.D["Q"][j](x)).log()
    if M.D["mode"] == "lower":
        val += (x * x - 1) / 2 if xq <= 1 else x.log()
        val += M.D["const"]
    return val

acc = []
orig = M.lb_sub
def rec(p, q, model):
    r = orig(p, q, model)
    if r[0]:
        acc.append((p, q, model))
    return r
M.lb_sub = rec
ks = rng.sample(range(N), min(NSAMP, N))
for k in ks:
    M.run_interval(e[k], e[k + 1], M.far_model(k), "v")
print(f"{len(acc)} accepted sub-intervals")
cap = {}
def quad(a0, a1, a2, w):
    cap["a"] = (a0, a1, a2); return True, 0.0
bad, worst, nt = 0, None, 0
for (p, q, model) in rng.sample(acc, min(1200, len(acc))):
    sv = M.quad_lb_positive; M.quad_lb_positive = quad
    try:
        orig(p, q, model)
    finally:
        M.quad_lb_positive = sv
    a0, a1, a2 = cap["a"]; mid = (p + q) / 2
    for _ in range(4):
        x = p + (q - p) * fmpq(rng.randrange(0, 10 ** 6 + 1), 10 ** 6)
        if any(lo <= x <= hi for (lo, hi, j) in M.D["brackets"]):
            continue
        u = arb(x - mid); mino = a0 + a1 * u + a2 * u * u
        gap = full(x) - mino; nt += 1
        if gap < 0:
            bad += 1; print("MINORANT ABOVE FUNCTION at", float(x), gap)
        g = float(gap.mid())
        if worst is None or g < worst[0]:
            worst = (g, float(x))
print(f"pointwise tests {nt}: minorant above the function at {bad} points; smallest gap {worst[0]:.3e} at x = {worst[1]:.6f}")
