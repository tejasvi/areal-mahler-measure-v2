#!/usr/bin/env python3
"""scan_contacts.py -- evaluates the areal function (F for the lower format, h for the upper format) at node-interval
midpoints and finer points and lists the smallest values (contact points), to choose the sites of the hidden controls.

    python3 -I scan_contacts.py CERT.json"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_arl_second as M
from flint import arb, fmpq, ctx
M.setup(sys.argv[1])
e, N = M.D["e"], M.D["N"]
if M.D["mode"] == "lower":
    ctx.prec = 192; I = M.energy_lower(); ctx.prec = 128
    M.D["const"] = -arb(M.D["m"]) + I.mid() / arb(M.D["s"])
def full(xq):
    x = arb(xq)
    val = M.D["sigma"] * M.U_real(x)
    for j in M.D["W"]:
        val -= M.D["wa"][j] * abs(M.D["Q"][j](x)).log()
    if M.D["mode"] == "lower":
        val += (x * x - 1) / 2 if xq <= 1 else x.log()
        val += M.D["const"]
    return float(val.mid())
vals = []
k0 = next(k for k in range(N + 1) if M.D["G"][k] > 0) - 1
for k in range(max(k0, 0), N):
    if M.D["G"][k] == 0 and M.D["G"][k + 1] == 0:
        continue
    for t in (1, 2, 3):
        x = e[k] + (e[k + 1] - e[k]) * fmpq(t, 4)
        vals.append((full(x), k, float(x)))
vals.sort()
for v, k, x in vals[:25]:
    print(f"{v:.4e} at x = {x:.6f} (node interval {k}, G_k = {float(M.D['G'][k]):.4g}, C_k = {float(M.D['C'][k]):.4g})")
