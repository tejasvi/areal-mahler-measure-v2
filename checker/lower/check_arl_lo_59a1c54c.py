#!/usr/bin/env python3
"""check_arl_lo.py -- rigorous checker for AREAL (Pritsker-measure) lower-bound certificates: a lower bound for the smallest limit
point lambda_A+ (and the essential minimum E_A+) of the absolute Pritsker measure N0 = ||P_alpha||_0^{1/deg} of totally positive
algebraic integers (also of all totally positive algebraic numbers). Derived from check_tp_61f4ef8f.py (the audited tpmahler lower-bound
checker) by replacing log+ x with g(x) = log+ x + (min(1, x)^2 - 1)/2 (Lemma TA1 = Lemma TP1 with g; tpmahler/PROGRESS.md).

Certificate (JSON, "tpmahler-areal-lower-v1"; same fields as "tpmahler-cert-v1"):
  polys : integer polynomials Q_j (coefficients highest degree first); must contain x = [1, 0]
  c     : rationals c_j >= 0 ("p/q"), one per polynomial
  nodes : rationals e_0 < e_1 < ... < e_N, e_0 >= 0
  G     : rationals G_0..G_N >= 0 with G_0 = G_N = 0 (values of a continuous piecewise-linear density at the nodes)
  m     : rational (claim: log lambda_A+ >= m)
  x_lo, x_hi : rationals (tail cut points)
Definitions (rebuilt here from scratch): rho = G(x) dx, s = rho(R) = sum_k (e_{k+1}-e_k)(G_k+G_{k+1})/2,
  U(x) = int log|x-y| G(y) dy, I = int int log|x-y| G(x) G(y) dx dy,
  F(x) = g(x) - m - sum_j c_j log|Q_j(x)| - 2 U(x) + I/s,  g(x) = log+ x + (min(1, x)^2 - 1)/2.
VERIFIED means: (i) data constraints; (ii) budget sum_j c_j deg Q_j + 2 s < 1 (strict); (iii) F(x) >= 0 for all x in (0, inf)
(F = +inf at positive roots of the Q_j with c_j > 0). With nu = rho/s and lambda = s this is (H0)-(H3) and (H2') of Lemma TP1
with g in place of log+ (Lemma TA1), so log lambda_A+ >= m, E_A+ >= e^m, also for all totally positive numbers.
Closed forms: G'' = sum_i C_i delta_{e_i} (exact slope jumps), U(x) = sum_i C_i H(e_i - x), H(t) = t^2/2 log|t| - 3t^2/4,
I = sum_{i,j} C_i C_j J(e_i - e_j), J(t) = t^4/24 log|t| - 25 t^4/288 (H'' = log|t|, J'' = H).
Positivity: tails (0, x_lo] (g >= -1/2 there) and [x_hi, inf) by monotone bounds; on [x_lo, x_hi] root windows (squarefree parts,
Sturm counts, exact sign-change brackets, ball evaluation with log of the upper bound of |Q|; g increasing, so g >= g(lower end)), a
breakpoint at x = 1 (g'' jumps there), and the chord bound F >= min(F(p),F(q)) - (h/4) int_p^q (F'')_+ with rigorous ball bounds;
g'' = 1 on (0, 1) is added to the regular part, g'' = -1/x^2 <= 0 on (1, inf) is dropped.
Usage: check_arl_lo.py CERT.json [NPROC]   (refuses python -O)
"""
import json, sys, math, time
from fractions import Fraction as Fr
from flint import arb, acb, fmpz_poly, ctx

ctx.prec = 200


class CheckError(Exception):
    pass


def require(cond, msg):
    if not cond:
        raise CheckError(msg)


def A(fr):
    fr = Fr(fr)
    return arb(fr.numerator)/arb(fr.denominator)


def ball(lo, hi):            # arb ball containing [lo, hi] for Fractions lo <= hi
    mid = (lo + hi)/2
    return A(mid) + arb(0, A((hi - lo)/2 + Fr(1, 10**40)).upper())


def logabs(t):
    return abs(t).log()


def Hf(t):                   # H(t) = t^2/2 log|t| - 3 t^2/4 (t arb, possibly a ball containing 0)
    if t.is_zero():
        return arb(0)
    if t.contains(0):
        r = abs(t).upper()
        require(r < arb('0.6'), "H: ball too wide around a node")
        B = r*r*((1/r).log()/2 + arb(3)/4)
        return arb(0, B.upper())
    return t*t/2*logabs(t) - arb(3)/4*t*t


def Jf(t):
    if t.is_zero():
        return arb(0)
    t2 = t*t
    return t2*t2/24*logabs(t) - arb(25)/288*t2*t2


def polyval(c, x):
    v = arb(0)
    for a in c:
        v = v*x + a
    return v


def polyval_exact(c, xf):
    v = Fr(0)
    for a in c:
        v = v*xf + a
    return v


def polyder(c):
    d = len(c) - 1
    return [c[i]*(d - i) for i in range(d)]


class Cert:
    def __init__(self, path):
        d = json.load(open(path))
        require(d.get('format') == 'tpmahler-areal-lower-v1', "unknown format")
        self.polys = [[int(a) for a in q] for q in d['polys']]
        require(all(isinstance(a, int) for q in d['polys'] for a in q), "non-integer coefficient")
        self.c = [Fr(x) for x in d['c']]
        self.e = [Fr(x) for x in d['nodes']]
        self.G = [Fr(x) for x in d['G']]
        self.m = Fr(d['m'])
        self.x_lo = Fr(d['x_lo'])
        self.x_hi = Fr(d['x_hi'])
        require(len(self.c) == len(self.polys) and all(x >= 0 for x in self.c), "c_j must be >= 0, one per polynomial")
        require(all(len(q) >= 1 and q[0] != 0 for q in self.polys), "polynomials must have a nonzero leading coefficient")
        require(len(set(tuple(q) for q in self.polys)) == len(self.polys), "duplicate polynomial")
        require(len(self.G) == len(self.e) >= 3, "G and nodes must have the same length >= 3")
        require(all(self.e[i] < self.e[i + 1] for i in range(len(self.e) - 1)), "nodes not strictly increasing")
        require(self.e[0] >= 0, "support must lie in [0, inf)")
        require(all(g >= 0 for g in self.G) and self.G[0] == 0 and self.G[-1] == 0, "G >= 0 with G_0 = G_N = 0 required")
        self.active = [(q, A(cj), cj) for q, cj in zip(self.polys, self.c) if cj != 0]
        # certified complex roots (with multiplicity) of each active polynomial: -(log|Q|)'' = sum_rho Re 1/(x - rho)^2 for real x
        self.roots = []
        for q, l, cj in self.active:
            f = fmpz_poly(list(reversed(q)))
            rs = []
            if f.degree() >= 1:
                for z, mult in f.complex_roots():
                    rs += [z]*int(mult)
            require(len(rs) == f.degree(), 'root count mismatch')
            self.roots.append(rs)
        lx = [cj for q, cj in zip(self.polys, self.c) if q == [1, 0]]
        require(len(lx) == 1 and lx[0] > 0, "the polynomial x must be present with a positive weight")
        self.cx = lx[0]
        N = len(self.e) - 1
        self.s = sum((self.e[k + 1] - self.e[k])*(self.G[k] + self.G[k + 1])/2 for k in range(N))
        require(self.s > 0, "zero mass")
        slope = [(self.G[k + 1] - self.G[k])/(self.e[k + 1] - self.e[k]) for k in range(N)]
        C = [Fr(0)]*(N + 1)
        C[0] = slope[0]
        C[N] = -slope[N - 1]
        for i in range(1, N):
            C[i] = slope[i] - slope[i - 1]
        require(sum(C) == 0 and sum(C[i]*self.e[i] for i in range(N + 1)) == 0, "G'' moments (compact support) violated")
        self.C = C
        self.nz = [i for i in range(N + 1) if C[i] != 0]
        self.Ca = {i: A(C[i]) for i in self.nz}
        self.ea = {i: A(self.e[i]) for i in self.nz}
        self.degsum = sum(cj*(len(q) - 1) for q, cj in zip(self.polys, self.c))
        self.first = min(self.e[i] for i in self.nz)
        self.last = max(self.e[i] for i in self.nz)

    def energy(self):
        tot = arb(0)
        for a_ in self.nz:
            row = arb(0)
            ea = self.ea[a_]
            for b_ in self.nz:
                if b_ != a_:
                    row += self.Ca[b_]*Jf(ea - self.ea[b_])
            tot += self.Ca[a_]*row
        return tot

    def U(self, x):                    # x arb (ball allowed)
        v = arb(0)
        for i in self.nz:
            v += self.Ca[i]*Hf(self.ea[i] - x)
        return v

    def U_exact(self, xf):             # xf Fraction: exact differences
        v = arb(0)
        for i in self.nz:
            t = self.e[i] - xf
            if t != 0:
                v += self.Ca[i]*Hf(A(t))
        return v


def g_lower(x):                        # a lower bound of g over the ball x (g is increasing on [0, inf))
    lo = x.lower()
    if lo >= 1:
        return lo.log()
    lo0 = lo if lo > 0 else arb(0)
    return (lo0*lo0 - 1)/2


def g_exact(xf):                       # g at a rational point xf > 0
    return A(xf).log() if xf >= 1 else A((xf*xf - 1)/2)


def setup(cert):
    cert.I = cert.energy()
    cert.const = cert.I/A(cert.s) - A(cert.m)       # F = g(x) + const - sum c log|Q| - 2U
    return cert


def F_exact(cert, xf):
    v = cert.const + g_exact(xf)
    for q, l, _ in cert.active:
        Qx = polyval_exact(q, xf)
        require(Qx != 0, "evaluation exactly at a root")
        v -= l*logabs(A(Qx))
    v -= 2*cert.U_exact(xf)
    return v


def Fpp_upper(cert, p, q):
    """regular and singular upper bounds of (F'')_+ data on [p, q] (Fractions, p < q, no node strictly inside, 1 not inside)."""
    X = ball(p, q)
    h = q - p
    require(not (p < 1 < q), "chord interval crosses x = 1")
    reg = arb(1) if q <= 1 else arb(0)          # g'' = 1 on (0, 1); g'' = -1/x^2 <= 0 on (1, inf) is dropped
    Xc = acb(X)
    for (qq, l, _), rs in zip(cert.active, cert.roots):
        tot = arb(0)
        for z in rs:                             # -(log|Q|)'' = sum over roots of Re 1/(x - rho)^2 (x real, Q != 0 on [p, q])
            d = Xc - z
            if d.real.contains(0) and d.imag.contains(0):
                return None, None
            tot += (1/(d*d)).real
        reg += l*tot
    sing = arb(0)
    for i in cert.nz:
        ei = cert.e[i]
        if ei == p or ei == q:
            c = cert.C[i]
            if c > 0:                            # -2c log|e-x| >= 0 near e: int over (0, h] of -2c log t = 2c h (1 + log(1/h))
                if h > 1:
                    return None, None
                hh = A(h)
                sing += 2*cert.Ca[i]*hh*(1 + (1/hh).log())
            else:                                # -2c log|e-x| = 2|c| log|e-x| <= 2|c| log h  (h <= 1)
                if h > 1:
                    return None, None
                reg += -2*cert.Ca[i]*A(h).log()
        else:
            d = cert.ea[i] - X
            if d.contains(0):
                return None, None
            reg += -2*cert.Ca[i]*logabs(d)
    return reg, sing


def check_interval(cert, p, q):
    Fp = F_exact(cert, p)
    Fq = F_exact(cert, q)
    lo = Fp.lower() if Fp.lower() < Fq.lower() else Fq.lower()
    if not (lo > 0):
        return False, float(lo)
    reg, sing = Fpp_upper(cert, p, q)
    if reg is None or not (reg.is_finite() and sing.is_finite()):
        return False, float(lo)
    regu = reg.upper()
    h = A(q - p)
    pen = h/4*(h*(regu if regu > 0 else arb(0)) + sing)
    return bool((arb(lo) - pen).lower() > 0), float(lo)


def sturm_count(P0, lo, hi):
    def pdiv(a, b):
        a = a[:]
        while len(a) >= len(b) and any(a):
            f = a[0]/b[0]
            for i in range(len(b)):
                a[i] -= f*b[i]
            a.pop(0)
        while a and a[0] == 0:
            a.pop(0)
        return a
    P1 = [c*(len(P0) - 1 - i) for i, c in enumerate(P0[:-1])]
    seq = [P0, P1]
    while True:
        r = pdiv(seq[-2], seq[-1])
        if not r:
            break
        seq.append([-c for c in r])

    def V(x):
        vals = [polyval_exact(p, x) for p in seq]
        vals = [v for v in vals if v != 0]
        return sum(1 for a, b in zip(vals, vals[1:]) if (a > 0) != (b > 0))
    return V(lo) - V(hi)


def arb_to_fr(x):                      # exact rational value of the midpoint of an arb
    mm, ee = x.mid().man_exp()
    mm, ee = int(mm), int(ee)
    return Fr(mm)*(Fr(2)**ee if ee >= 0 else Fr(1, 2**(-ee)))


def root_brackets(q, lo, hi):
    """disjoint exact brackets (a, b) with a sign change of the squarefree part of q, one per distinct root in (lo, hi).
    Root locations come from flint's certified complex root isolation; every bracket is re-verified by exact signs and the
    total is matched against a Sturm count over Q."""
    f = fmpz_poly(list(reversed(q)))
    g = f.gcd(f.derivative()) if f.degree() > 0 else fmpz_poly(1)
    sq = f // g if g.degree() > 0 else f
    P0 = [Fr(int(c)) for c in reversed(sq.coeffs())]          # highest first
    require(polyval_exact(P0, lo) != 0 and polyval_exact(P0, hi) != 0, "root at a tail cut point")
    count = sturm_count(P0, lo, hi) if len(P0) > 1 else 0
    br = []
    if sq.degree() >= 1:
        for z, mult in sq.complex_roots():
            if not z.imag.contains(0):
                continue
            a0 = arb_to_fr(z.real.lower()); b0 = arb_to_fr(z.real.upper())
            if b0 <= lo or a0 >= hi:
                continue
            ok = False
            for dl in [Fr(1, 10**30), Fr(1, 10**20), Fr(1, 10**14), Fr(1, 10**10)]:
                a_, b_ = a0 - dl, b0 + dl
                if polyval_exact(P0, a_)*polyval_exact(P0, b_) < 0:
                    ok = True
                    break
            if ok:
                require(lo < a_ and b_ < hi, "root bracket crosses a tail cut point")
                br.append((a_, b_))
    br.sort()
    require(all(br[i][1] < br[i + 1][0] for i in range(len(br) - 1)), "overlapping root brackets")
    require(len(br) == count, f"Sturm count mismatch for {q}: {count} vs {len(br)}")
    return br


def Uprime_ball(cert, X):                 # U'(x) = -sum_i C_i H'(e_i - x), H'(t) = t log|t| - t (continuous, 0 at t = 0)
    v = arb(0)
    for i in cert.nz:
        t = cert.ea[i] - X
        if t.contains(0):
            r = abs(t).upper()
            require(r < arb('0.3'), "U': ball too wide around a node")
            v += cert.Ca[i]*arb(0, (r*((1/r).log() + 1)).upper())     # |t log|t| - t| <= r (log(1/r) + 1) for |t| <= r < 1/e
        else:
            v -= cert.Ca[i]*(t*logabs(t) - t)
    return v


def window_ok(cert, L, R):
    """lower bound of F on [L, R] in mean-value form: exact values at the rational centre plus derivative balls times the radius."""
    X = ball(L, R)
    cm = (L + R)/2
    rad = A((R - L)/2)
    v = cert.const + g_lower(X)
    for qq, l, _ in cert.active:
        Qc = abs(A(polyval_exact(qq, cm)))
        Q1 = polyval(polyder(qq), X) if len(qq) > 1 else arb(0)
        ub = Qc + abs(Q1).upper()*rad                         # max |Q| on [L, R]
        v -= l*arb(ub.upper()).log()
    v -= 2*(cert.U_exact(cm) + abs(Uprime_ball(cert, X)).upper()*rad)
    return v.lower() > 0, v


def windows_for(cert):
    excl = []
    for q, l, _ in cert.active:
        for (a_, b_) in root_brackets(q, cert.x_lo, cert.x_hi):
            excl.append((a_, b_))
    excl.sort()
    wins = []
    for (a_, b_) in excl:
        d = Fr(1, 10**3)
        ok = False
        while d >= Fr(1, 10**14):
            L, R = a_ - d, b_ + d
            if L > cert.x_lo and R < cert.x_hi:
                good, v = window_ok(cert, L, R)
                if good:
                    ok = True
                    break
            d /= 4
        require(ok, f"root window failed near {float(a_)}")
        wins.append((a_ - d, b_ + d))
    # merge overlapping windows
    wins.sort()
    merged = []
    for w in wins:
        if merged and w[0] <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], w[1]))
        else:
            merged.append(w)
    for (L, R) in merged:                       # re-verify merged windows as a whole, splitting if needed
        pieces = [(L, R)]
        out_ok = True
        while pieces:
            a_, b_ = pieces.pop()
            good, _ = window_ok(cert, a_, b_)
            if good:
                continue
            if b_ - a_ < Fr(1, 10**15):
                out_ok = False
                break
            mid = (a_ + b_)/2
            pieces += [(a_, mid), (mid, b_)]
        require(out_ok, f"merged window [{float(L)}, {float(R)}] failed")
    return merged


def tail_left(cert):
    x_lo = cert.x_lo
    require(0 < x_lo <= cert.first and x_lo < 1, "x_lo must be in (0, 1) and at or left of the support")
    v = cert.const - arb(1)/2 - A(cert.cx)*A(x_lo).log()   # g >= -1/2 on (0, x_lo]; -c_x log x >= -c_x log x_lo
    K = 64                                            # max |Q| on [0, x_lo] from K sub-balls (upper bounds)
    subs = [ball(x_lo*k/K, x_lo*(k + 1)/K) for k in range(K)]
    for q, l, _ in cert.active:
        if q == [1, 0]:
            continue
        P0 = [Fr(a) for a in q]
        require(polyval_exact(P0, Fr(0)) != 0 and polyval_exact(P0, x_lo) != 0, "an active polynomial vanishes at 0 or x_lo")
        require(sturm_count(P0, Fr(0), x_lo) == 0 if len(P0) > 1 else True, "an active polynomial has a root in (0, x_lo]")
        mx = abs(polyval(q, subs[0])).upper()
        for X in subs[1:]:
            mx = mx.max(abs(polyval(q, X)).upper())     # arb max: a ball containing the larger value
        v -= l*arb(mx.upper()).log()
    v -= 2*cert.U_exact(Fr(0))                       # U(x) <= U(0) for 0 <= x <= first support point
    return v.lower() > 0, v


def tail_right(cert):
    x_hi = cert.x_hi
    require(x_hi >= 1 and x_hi >= cert.last, "x_hi must be >= 1 and >= the support")
    kappa = 1 - A(cert.degsum) - 2*A(cert.s)
    beta = A(cert.m) - cert.I/A(cert.s)
    for q, l, _ in cert.active:
        beta += l*arb(sum(abs(a) for a in q)).log()
    g = kappa*A(x_hi).log() - beta                    # F(x) >= kappa log x - beta for x >= x_hi (increasing since kappa > 0)
    return kappa.lower() > 0 and g.lower() > 0, g, kappa


_CERT = None


def _work(task):
    a_, b_ = task
    cert = _CERT
    stack = [(a_, b_, 0)]
    n = 0
    mn = None
    while stack:
        p, q, dep = stack.pop()
        ok, lo = check_interval(cert, p, q)
        n += 1
        if ok:
            mn = lo if mn is None else min(mn, lo)
            continue
        if dep >= 50 or not (lo > 0):
            return False, n, lo, (str(p), str(q))
        mid = (p + q)/2
        md = Fr(round(mid*2**48), 2**48)
        mid = md if p < md < q else mid
        stack += [(p, mid, dep + 1), (mid, q, dep + 1)]
    return True, n, mn, None


def main():
    import multiprocessing as mp
    global _CERT
    if sys.flags.optimize:
        print("refusing to run under python -O / -OO")
        sys.exit(2)
    t0 = time.time()
    path = sys.argv[1]
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    cert = setup(Cert(path))
    print(f"polys {len(cert.polys)} (active {len(cert.active)}), nodes {len(cert.e)}, s = {float(cert.s):.8f}, m = {float(cert.m):.10f}")
    budget = cert.degsum + 2*cert.s
    require(budget < 1, "budget sum c_j deg Q_j + 2 s < 1 fails")       # exact rational comparison
    print(f"budget sum c deg + 2 s = {float(budget):.6f} < 1  OK;  I = {cert.I}")
    okL, vL = tail_left(cert)
    print("left tail (0, x_lo]:", okL, float(vL.lower()))
    okR, gR, kap = tail_right(cert)
    print("right tail [x_hi, inf):", okR, float(gR.lower()), "kappa", float(kap.lower()))
    require(okL and okR, "tail check failed")
    wins = windows_for(cert)
    print(f"root windows: {len(wins)}  ({time.time() - t0:.0f}s)")
    pts = set([cert.x_lo, cert.x_hi]) | set(e for e in cert.e if cert.x_lo < e < cert.x_hi)
    if cert.x_lo < 1 < cert.x_hi:
        pts.add(Fr(1))
    for (L, R) in wins:
        pts |= {L, R}
    pts = sorted(pts)
    tasks = [(p, q) for p, q in zip(pts, pts[1:]) if not any(L <= p and q <= R for (L, R) in wins)]
    # coverage
    cover = sorted(tasks + wins)
    pos = cert.x_lo
    for (L, R) in cover:
        require(L <= pos, f"coverage gap at {float(pos)}")
        pos = max(pos, R)
    require(pos >= cert.x_hi, "coverage does not reach x_hi")
    print(f"chord intervals: {len(tasks)}", flush=True)
    _CERT = cert
    with mp.get_context('fork').Pool(nproc) as pool:
        res = pool.map(_work, tasks, chunksize=max(1, len(tasks)//(8*nproc)))
    n = sum(r[1] for r in res)
    bad = [(t, r) for t, r in zip(tasks, res) if not r[0]]
    mins = [r[2] for r in res if r[0] and r[2] is not None]
    print(f"interval checks: {n}; smallest endpoint lower bound {min(mins):.4e}" if mins else "no intervals")
    if bad:
        for t, r in bad[:10]:
            print(f"  FAILED near x = {float(t[0]):.6f}: lower {r[2]:.3e} at {r[3]}")
        print(f"NOT VERIFIED ({len(bad)} failed; {time.time() - t0:.0f}s)")
        sys.exit(1)
    E = A(cert.m).exp()
    lo7 = math.floor(arb_to_fr(E.lower())*10**7)
    print(f"VERIFIED: log lambda_A+ >= m = {cert.m}; N0 >= e^m >= {lo7 // 10**7}.{str(lo7 % 10**7).zfill(7)} (also E_A+)  ({time.time() - t0:.0f}s)")
    sys.exit(0)


if __name__ == '__main__':
    try:
        main()
    except CheckError as e:
        print("CheckError:", e)
        print("NOT VERIFIED")
        sys.exit(1)
