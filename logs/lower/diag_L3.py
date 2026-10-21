# Diagnostic (not part of any proof): find the chord intervals of a lower certificate that need deep bisection in check_arl_lo.
# Same setup as check_arl_lo.main(); per-task bisection capped at CAP checks; chunksize 1.  Usage: diag_L3.py CERT.json NPROC CAP
import os, sys, time, json, multiprocessing as mp
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'checker', 'lower'))
import check_arl_lo_59a1c54c as C
from fractions import Fraction as Fr
path, nproc, CAP = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
t0 = time.time()
cert = C.setup(C.Cert(path))
wins = C.windows_for(cert)
pts = set([cert.x_lo, cert.x_hi]) | set(e for e in cert.e if cert.x_lo < e < cert.x_hi)
if cert.x_lo < 1 < cert.x_hi: pts.add(Fr(1))
for (L, R) in wins: pts |= {L, R}
pts = sorted(pts)
tasks = [(p, q) for p, q in zip(pts, pts[1:]) if not any(L <= p and q <= R for (L, R) in wins)]
print(f"windows {len(wins)}, tasks {len(tasks)} ({time.time() - t0:.0f}s)", flush=True)
def probe(i):
    a_, b_ = tasks[i]; s = time.time(); stack = [(a_, b_, 0)]; n = 0; mn = None; md = 0; worst = None
    while stack and n < CAP:
        p, q, dep = stack.pop(); ok, lo = C.check_interval(cert, p, q); n += 1; md = max(md, dep)
        if ok:
            mn = lo if mn is None else min(mn, lo); continue
        if dep >= 50 or not (lo > 0): return (i, 'FAIL', n, md, time.time() - s, float(p), float(q), lo)
        mid = (p + q)/2; m2 = Fr(round(mid*2**48), 2**48); mid = m2 if p < m2 < q else mid
        stack += [(p, mid, dep + 1), (mid, q, dep + 1)]
        worst = (float(p), float(q))
    return (i, 'CAP' if stack else 'OK', n, md, time.time() - s, float(a_), float(b_), len(stack), worst)
with mp.get_context('fork').Pool(nproc) as pool:
    res = pool.map(probe, range(len(tasks)), chunksize=1)
heavy = sorted([r for r in res if r[1] != 'OK' or r[2] > 200], key=lambda r: -r[2])
chunk = max(1, len(tasks)//(8*48))
print(f"probe done ({time.time() - t0:.0f}s); total checks {sum(r[2] for r in res)}; heavy/capped/failed: {len(heavy)}; checker chunksize at 48 procs = {chunk}")
for r in heavy[:40]: print('  ', r, 'chunk', r[0] // chunk)
