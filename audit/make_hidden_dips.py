#!/usr/bin/env python3
"""make_hidden_dips.py -- two hidden negative controls for the areal UPPER format: a flat-top bump of density (feet of width
wt, top of width 2a) is added (dip A) with centre XNUM/XDEN + 1/(3*10^7), a = 2e-6, wt = 1e-6, and (dip B) with its left
foot 1e-7 to the right of the node with index NODE, a = wt = 5e-7. Each bump's mass is tuned by bisection so that
h_new < 0 at the bump centre but h_new > 0 at every node (old and inserted) within 0.01 and at the 11 points of the 2e-5
grid around the bump. v is raised by 1e-6 so that only (U3) can fail. Writes hidden_up_A_contact.json and
hidden_up_B_nearnode.json to OUTDIR.

    python3 -I make_hidden_dips.py CERT XNUM XDEN NODE OUTDIR"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_arl_second as M
from flint import arb, fmpq

SRC = sys.argv[1]
M.setup(SRC)
e0, G0, N0 = list(M.D["e"]), list(M.D["G"]), M.D["N"]

def h_full(xq):
    x = arb(xq)
    val = M.U_real(x)
    for j in M.D["W"]:
        val -= M.D["wa"][j] * abs(M.D["Q"][j](x)).log()
    return val

def interp(x):
    lo, hi = 0, N0
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if e0[mid] <= x: lo = mid
        else: hi = mid
    k = min(lo, N0 - 1)
    return G0[k] + (G0[k + 1] - G0[k]) * (x - e0[k]) / (e0[k + 1] - e0[k])

def make(x0, a, wt, mt, tag):
    ht = mt / (2 * a + wt)
    pts = [x0 - a - wt, x0 - a, x0 + a, x0 + a + wt]
    add = [fmpq(0), ht, ht, fmpq(0)]
    pairs = sorted(list(zip(e0, G0)) + [(x, interp(x) + dv) for x, dv in zip(pts, add)], key=lambda t: t[0])
    d = json.load(open(SRC))
    d["nodes"] = [str(p[0]) for p in pairs]
    d["G"] = [str(p[1]) for p in pairs]
    v = fmpq(*[int(t) for t in d["v"].split("/")]) + fmpq(1, 10 ** 6)
    d["v"] = f"{v.p}/{v.q}"
    path = os.path.join(OUTDIR, f"hidden_up_{tag}.json")
    json.dump(d, open(path, "w"))
    return path, pts

def probe(path, x0, pts):
    M.D.clear(); M.setup(path)
    hc = h_full(x0)
    checks = list(pts) + [x for x in M.D["e"] if abs(x - x0) <= fmpq(1, 100)]
    g0 = (x0 * 50000).floor() if hasattr(x0, "floor") else None
    grid = []
    k0 = int((x0 * 50000).p // (x0 * 50000).q)
    for kk in range(k0 - 5, k0 + 6):
        grid.append(fmpq(kk, 50000))
    vals = [(float(x), h_full(x)) for x in checks + grid]
    M.D.clear(); M.setup(SRC)
    mn = min(vals, key=lambda t: float(t[1].mid()))
    return hc, mn, all(v > 0 for _, v in vals)

XA = fmpq(int(sys.argv[2]), int(sys.argv[3])); IB = int(sys.argv[4]); OUTDIR = sys.argv[5]
for tag, x0, a, wt in (("A_contact", XA + fmpq(1, 3 * 10 ** 7), fmpq(2, 10 ** 6), fmpq(1, 10 ** 6)),
                       ("B_nearnode", None, fmpq(1, 2 * 10 ** 6), fmpq(1, 2 * 10 ** 6))):
    if x0 is None:
        iB = IB
        x0 = e0[iB] + a + wt + fmpq(1, 10 ** 7)
    lo, hi = fmpq(1, 10 ** 8), fmpq(3, 10 ** 7)       # mass bracket: lo -> centre positive, hi -> nodes negative
    best = None
    for it in range(30):
        mt = (lo + hi) / 2
        path, pts = make(x0, a, wt, mt, tag)
        hc, mn, allpos = probe(path, x0, pts)
        if hc < 0 and allpos:
            best = (mt, hc, mn, path)
            break
        if hc < 0:          # centre negative but some node/grid point negative too: less mass
            hi = mt
        else:
            lo = mt
    if best:
        mt, hc, mn, path = best
        print(f"{tag}: x0 = {float(x0):.9f}, mass {float(mt):.4e}: h_new(centre) = {hc.str(5)} < 0; min over nodes within 0.01, "
              f"the 4 inserted nodes and 11 grid points of spacing 2e-5: {mn[1].str(5)} at {mn[0]:.7f} (> 0); file {path}")
    else:
        print(tag, "no hidden mass found; last", float(mt), hc, mn)
