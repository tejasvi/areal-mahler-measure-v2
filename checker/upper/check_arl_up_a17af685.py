#!/usr/bin/env python3
"""check_arl_up.py -- rigorous checker for AREAL (Pritsker-measure) UPPER-bound certificates: lambda_A+ (smallest limit point of
N0(alpha) = ||P_alpha||_0^{1/deg} over totally positive algebraic integers) and E_A+, via Smith's criterion (arXiv:2111.12660, Thm 1.5)
in the checkable form of Lemma UA1 (= Lemma U1 with log+ replaced by g(x) = log+ x + (min(1, x)^2 - 1)/2; tpmahler/PROGRESS.md).
Derived from check_up_29783014.py: only the format name and (V) change ((1/s) int g dmu instead of (1/s) int log+ dmu).

Certificate (JSON, "tpmahler-areal-upper-v1"; same fields as "tpmahler-upper-v1"):
  polys : integer polynomials Q_j (coefficients highest degree first)
  b     : rationals b_j >= 0 ("p/q"), one per polynomial
  nodes : rationals 0 <= e_0 < e_1 < ... < e_N
  G     : rationals G_0..G_N >= 0, G_0 = G_N = 0 (values of a continuous piecewise-linear density at the nodes)
  X     : rational > 4 with e_N <= X  (Sigma = [0, X] has capacity X/4 > 1)
  v     : rational (claim: log lambda_A+ <= v and log E_A+ <= v)
Definitions (rebuilt from scratch): mu = G(x) dx, s = mu(R), U(x) = int log|x - y| G(y) dy = sum_i C_i H(e_i - x) with exact slope
jumps C_i, H(t) = t^2/2 log|t| - 3t^2/4, h(x) = U(x) - sum_j b_j log|Q_j(x)|.
VERIFIED means: (D) data: integer polynomials, each primitive and irreducible over Q, pairwise distinct; b_j >= 0; nodes increasing, G >= 0, G_0 = G_N = 0, s > 0; X > 4, e_N <= X;
(U1) sum_j b_j deg Q_j < s (exact);  (U2) int log|Q_j| dmu = s log|a_j| + sum_rho U(rho) > 0 for every j with b_j > 0 (arb; for a
complex root z, U(z) = Re sum_i C_i Hc(e_i - z), Hc(t) = t^2/2 Log t - 3t^2/4: the branch terms cancel since sum C_i e_i^2 is real);
(U3) h(x) >= 0 for all x in [e_0, e_N] (root windows around roots of the Q_j with b_j > 0, chord bound
h >= min(h(p), h(q)) - (len/4) int_p^q (h'')_+ elsewhere, h'' = sum_i C_i log|e_i - x| + sum_j b_j sum_rho Re 1/(x - rho)^2);
(V) (1/s) int g(x) dmu(x) <= v (closed form, arb).  Then Lemma UA1 gives lambda_A+ <= e^v and E_A+ <= e^v.

How (U3) is proved on the WHOLE interval [e_0, e_N] (no sampling; the LP's grids play no role in the verdict):
 * Only [e_0, e_N] is needed: h is superharmonic off it and tends to +inf at infinity by (U1), so h >= 0 there gives h >= 0 on C
   (minimum principle; Lemma U1 step (1)). No tail estimates are required.
 * [e_0, e_N] is covered by root windows and chord intervals with exact rational ends; the coverage is checked exactly.
 * Root windows: for each weighted Q_j and each real root in [e_0, e_N], an exact rational bracket with a sign change of Q_j (from
   FLINT's certified root isolation, re-verified by exact rational evaluation); the number of brackets must equal the Sturm count of
   Q_j on [e_0, e_N]. On a window [L, R] with centre c and radius r: h >= U(c) - r max|U'| - sum_j b_j log(|Q_j(c)| + r max|Q_j'|)
   (mean value theorem for U and for every Q_j; max|U'| and max|Q_j'| by ball arithmetic over [L, R], nodes inside the window via
   |t log|t| - t| <= r (log(1/r) + 1)). Windows shrink by factors of 4 (from 1e-3) until the bound is positive; overlapping windows
   are merged and re-verified as a whole (bisected if needed).
 * Chord intervals [p, q] (no node strictly inside, no root of a weighted Q_j in [p, q] -- enforced by ball arithmetic): h is C^1 on
   [p, q] with integrable h'', so h(x) = l(x) - int_p^q g(x, t) h''(t) dt with l the chord through (p, h(p)), (q, h(q)) and the Green's
   function 0 <= g <= (q - p)/4; hence h >= min(h(p), h(q)) - ((q - p)/4) int_p^q (h'')_+ on the whole of [p, q]. h(p), h(q): exact
   rational arguments, arb. h''(x) = sum_i C_i log|e_i - x| + sum_j b_j sum_rho Re 1/(x - rho)^2: interior nodes and all roots by
   ball bounds over [p, q] (times q - p); an endpoint node with C_i < 0 contributes exactly int_0^L |C_i| log(1/t) dt =
   |C_i| L (1 + log(1/L)) (L = q - p <= 1); with C_i > 0 the bound C_i log L. If the bound is not positive the interval is bisected
   (dyadic midpoints, depth <= 50); any non-positive endpoint value means NOT VERIFIED.
 * Negative control "narrow dip" (make_dip_up.py): a tent of mass ~1.4e-7 and width 4e-6 placed strictly between two points of the
   LP's 2e-5 check grid keeps h > 0 at every check-grid point and at the tent's nodes but makes h < 0 inside (x0, x0 + 3w); this
   checker rejects it by the bisection above (FAILED near x0, endpoint value -4.8e-9).
Usage: check_arl_up.py CERT.json [NPROC]
"""
import json, sys, math, time
from fractions import Fraction as Fr
from flint import arb, acb, fmpz_poly, fmpz, ctx

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


def real_root_bound(q):                # Cauchy bound: all roots satisfy |x| < 1 + max|a_i/a_0|
    return 1 + max(Fr(abs(a), abs(q[0])) for a in q[1:]) if len(q) > 1 else Fr(1)


class Cert:
    def __init__(self, path):
        d = json.load(open(path))
        require(d.get('format') == 'tpmahler-areal-upper-v1', "unknown format")
        require(all(isinstance(a, int) for q in d['polys'] for a in q), "non-integer coefficient")
        self.polys = [[int(a) for a in q] for q in d['polys']]
        self.b = [Fr(x) for x in d['b']]
        self.e = [Fr(x) for x in d['nodes']]
        self.G = [Fr(x) for x in d['G']]
        self.X = Fr(d['X'])
        self.v = Fr(d['v'])
        require(len(self.b) == len(self.polys) and all(x >= 0 for x in self.b), "b_j must be >= 0, one per polynomial")
        require(all(len(q) >= 2 and q[0] != 0 for q in self.polys), "polynomials must be nonconstant with a nonzero leading coefficient")
        require(len(set(tuple(q) for q in self.polys)) == len(self.polys), "duplicate polynomial")
        require(len(self.G) == len(self.e) >= 3, "G and nodes must have the same length >= 3")
        require(all(self.e[i] < self.e[i + 1] for i in range(len(self.e) - 1)), "nodes not strictly increasing")
        require(self.e[0] >= 0, "support must lie in [0, inf)")
        require(all(g >= 0 for g in self.G) and self.G[0] == 0 and self.G[-1] == 0, "G >= 0 with G_0 = G_N = 0 required")
        require(self.X > 4 and self.e[-1] <= self.X, "X > 4 and e_N <= X required (capacity of [0, X] > 1)")
        self.active = []                                   # (q, b arb, b Fraction, roots as arb balls)
        for q, bj in zip(self.polys, self.b):
            f = fmpz_poly(list(reversed(q)))
            cont = 0
            for a in q:
                cont = math.gcd(cont, abs(a))
            require(cont == 1, f"polynomial not primitive: {q}")
            fac = f.factor()[1]
            require(len(fac) == 1 and fac[0][1] == 1 and fac[0][0].degree() == f.degree(), f"polynomial not irreducible over Q: {q}")
            if bj == 0:
                continue
            rs = []
            for z, mult in f.complex_roots():
                require(int(mult) == 1, "repeated root in an irreducible polynomial?")
                rs.append(z)                              # acb; real roots come with an exactly zero imaginary part
            require(len(rs) == f.degree(), 'root count mismatch')
            B = real_root_bound(q)
            nreal = sturm_count([Fr(a) for a in q], -B, B)
            require(sum(1 for z in rs if z.imag.is_zero()) == nreal, f"real-root count mismatch (flint vs Sturm) for {q}")
            self.active.append((q, A(bj), bj, rs))
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
        self.degsum = sum(bj*(len(q) - 1) for q, bj in zip(self.polys, self.b))
        self.lo, self.hi = self.e[0], self.e[-1]

    def U(self, x):                    # x arb (ball allowed)
        v = arb(0)
        for i in self.nz:
            v += self.Ca[i]*Hf(self.ea[i] - x)
        return v

    def U_complex(self, z):            # z acb with nonzero imaginary part: U(z) = Re sum_i C_i Hc(e_i - z)
        require(not z.imag.contains(0), "U_complex: root ball meets the real axis")
        v = acb(0)
        for i in self.nz:
            t = acb(self.ea[i]) - z
            v += acb(self.Ca[i])*(t*t/2*t.log() - acb(3)/4*t*t)
        return v.real

    def U_exact(self, xf):             # xf Fraction: exact differences
        v = arb(0)
        for i in self.nz:
            t = self.e[i] - xf
            if t != 0:
                v += self.Ca[i]*Hf(A(t))
        return v


def int_logQ(cert, q, rs):
    """int log|Q| dmu = s log|a| + sum over the (real, simple) roots rho of U(rho)  (exact identity; arb enclosure)."""
    v = A(cert.s)*arb(abs(q[0])).log()
    for z in rs:
        if not z.imag.is_zero():
            v += cert.U_complex(z)
            continue
        r = z.real
        if r.rad() == 0:
            v += cert.U_exact(arb_to_fr(r))
        else:
            v += cert.U(r)
    return v


def int_logplus(cert):
    """int g(x) G(x) dx, g = log+ x + (min(1, x)^2 - 1)/2, closed form on each piece (G = alpha + beta x):
    on [a, min(b, 1)]: g = (x^2 - 1)/2, exact rational antiderivative alpha (x^3/6 - x/2) + beta (x^4/8 - x^2/4);
    on [max(a, 1), b]: g = log x, F1 = x log x - x, F2 = x^2/2 log x - x^2/4 (as for log+)."""
    tot = arb(0)
    exact = Fr(0)
    for k in range(len(cert.e) - 1):
        a_, b_ = cert.e[k], cert.e[k + 1]
        beta = (cert.G[k + 1] - cert.G[k])/(b_ - a_)
        alpha = cert.G[k] - beta*a_
        if a_ < 1:
            hi = min(b_, Fr(1))
            P = lambda x: alpha*(x**3/6 - x/2) + beta*(x**4/8 - x**2/4)
            exact += P(hi) - P(a_)
        if b_ > 1:
            lo = max(a_, Fr(1))
            def F1(x): xa = A(x); return xa*xa.log() - xa
            def F2(x): xa = A(x); return xa*xa/2*xa.log() - xa*xa/4
            tot += A(alpha)*(F1(b_) - F1(lo)) + A(beta)*(F2(b_) - F2(lo))
    return tot + A(exact)


def h_exact(cert, xf):
    v = cert.U_exact(xf)
    for q, l, _, _ in cert.active:
        Qx = polyval_exact(q, xf)
        require(Qx != 0, "evaluation exactly at a root")
        v -= l*logabs(A(Qx))
    return v


def hpp_upper(cert, p, q):
    """regular (pointwise) and singular (integrated) upper bounds for (h'')_+ on [p, q] (no node strictly inside, no root of an
    active polynomial in [p, q])."""
    X = ball(p, q)
    h = q - p
    reg = arb(0)
    Xc = acb(X)
    for (qq, l, _, rs) in cert.active:
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
            if h > 1:
                return None, None
            if c < 0:                            # c log|e - x| = |c| log(1/|e - x|) >= 0: int over (0, h] = |c| h (1 + log(1/h))
                hh = A(h)
                sing += -cert.Ca[i]*hh*(1 + (1/hh).log())
            else:                                # c log|e - x| <= c log h (h <= 1)
                reg += cert.Ca[i]*A(h).log()
        else:
            d = cert.ea[i] - X
            if d.contains(0):
                return None, None
            reg += cert.Ca[i]*logabs(d)
    return reg, sing


def check_interval(cert, p, q):
    hp = h_exact(cert, p)
    hq = h_exact(cert, q)
    lo = hp.lower() if hp.lower() < hq.lower() else hq.lower()
    if not (lo > 0):
        return False, float(lo)
    reg, sing = hpp_upper(cert, p, q)
    if reg is None or not (reg.is_finite() and sing.is_finite()):
        return False, float(lo)
    regu = reg.upper()
    L = A(q - p)
    pen = L/4*(L*(regu if regu > 0 else arb(0)) + sing)
    return bool((arb(lo) - pen).lower() > 0), float(lo)


def Uprime_ball(cert, X):                 # U'(x) = -sum_i C_i H'(e_i - x), H'(t) = t log|t| - t (continuous, 0 at t = 0)
    v = arb(0)
    for i in cert.nz:
        t = cert.ea[i] - X
        if t.contains(0):
            r = abs(t).upper()
            require(r < arb('0.3'), "U': ball too wide around a node")
            v += cert.Ca[i]*arb(0, (r*((1/r).log() + 1)).upper())
        else:
            v -= cert.Ca[i]*(t*logabs(t) - t)
    return v


def window_ok(cert, L, R):
    """lower bound of h on [L, R]: U(centre) - max|U'| rad - sum_j b_j log max|Q_j| (mean-value bounds)."""
    X = ball(L, R)
    cm = (L + R)/2
    rad = A((R - L)/2)
    v = cert.U_exact(cm) - abs(Uprime_ball(cert, X)).upper()*rad
    for qq, l, _, _ in cert.active:
        Qc = abs(A(polyval_exact(qq, cm)))
        Q1 = polyval(polyder(qq), X) if len(qq) > 1 else arb(0)
        ub = Qc + abs(Q1).upper()*rad
        v -= l*arb(ub.upper()).log()
    return v.lower() > 0, v


def root_brackets(q, lo, hi):
    """exact disjoint brackets (a, b), a sign change of q, one per root in [lo, hi] (q irreducible, so squarefree); roots from
    flint's certified isolation, re-verified by exact signs; the count is matched against a Sturm count over Q."""
    f = fmpz_poly(list(reversed(q)))
    P0 = [Fr(a) for a in q]
    if polyval_exact(P0, lo) == 0 or polyval_exact(P0, hi) == 0:
        raise CheckError(f"a root of an active polynomial is an end of the support: {q}")
    count = sturm_count(P0, lo, hi)
    br = []
    for z, mult in f.complex_roots():
        if not z.imag.is_zero():
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
            br.append((max(a_, lo), min(b_, hi)))
    br.sort()
    require(all(br[i][1] < br[i + 1][0] for i in range(len(br) - 1)), "overlapping root brackets")
    require(len(br) == count, f"Sturm count mismatch for {q}: {count} vs {len(br)}")
    return br


def windows_for(cert):
    excl = []
    for q, l, _, _ in cert.active:
        excl += root_brackets(q, cert.lo, cert.hi)
    excl.sort()
    wins = []
    for (a_, b_) in excl:
        d = Fr(1, 10**3)
        ok = False
        while d >= Fr(1, 10**14):
            L, R = max(a_ - d, cert.lo), min(b_ + d, cert.hi)
            good, v = window_ok(cert, L, R)
            if good:
                ok = True
                break
            d /= 4
        require(ok, f"root window failed near {float(a_)}")
        wins.append((L, R))
    wins.sort()
    merged = []
    for w in wins:
        if merged and w[0] <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], w[1]))
        else:
            merged.append(w)
    for (L, R) in merged:                       # re-verify merged windows as a whole, splitting if needed
        pieces = [(L, R)]
        while pieces:
            a_, b_ = pieces.pop()
            good, _ = window_ok(cert, a_, b_)
            if good:
                continue
            require(b_ - a_ >= Fr(1, 10**15), f"merged window [{float(L)}, {float(R)}] failed")
            mid = (a_ + b_)/2
            pieces += [(a_, mid), (mid, b_)]
    return merged


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


def up_decimal(x, digits):              # decimal string >= the arb x (directed rounding up)
    q = arb_to_fr(x.upper()) if not x.upper().is_zero() else Fr(0)
    k = math.ceil(q*10**digits)
    return f"{k // 10**digits}.{str(k % 10**digits).zfill(digits)}"


def main():
    import multiprocessing as mp
    global _CERT
    if sys.flags.optimize:                 # -O guard: this checker never relies on assert, but refuse optimized runs anyway
        print("refusing to run under python -O / -OO")
        sys.exit(2)
    t0 = time.time()
    path = sys.argv[1]
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    cert = Cert(path)
    print(f"polys {len(cert.polys)} (active {len(cert.active)}), nodes {len(cert.e)}, support [{float(cert.lo)}, {float(cert.hi)}], "
          f"X = {cert.X}, s = {float(cert.s):.8f}, claim v = {float(cert.v):.10f}")
    require(cert.degsum < cert.s, "(U1) sum_j b_j deg Q_j < s fails")      # exact rational comparison
    print(f"(U1) sum b deg = {float(cert.degsum):.6f} < s = {float(cert.s):.6f}  OK")
    worst = None
    for q, l, bj, rs in cert.active:
        val = int_logQ(cert, q, rs)
        require(val.lower() > 0, f"(U2) int log|Q| dmu > 0 fails for {q}: {val}")
        worst = val if worst is None or val.lower() < worst.lower() else worst
    print(f"(U2) int log|Q_j| dmu > 0 for all {len(cert.active)} active Q_j; smallest {float(worst.lower()):.3e}")
    VL = int_logplus(cert)/A(cert.s)
    require(arb_to_fr(VL.upper()) <= cert.v, "(V) claim v below (1/s) int g dmu")        # exact dyadic upper end <= v
    print(f"(V) (1/s) int g dmu = {VL} <= v  OK")
    wins = windows_for(cert)
    print(f"root windows: {len(wins)}  ({time.time() - t0:.0f}s)")
    pts = set([cert.lo, cert.hi]) | set(cert.e)
    for (L, R) in wins:
        pts |= {L, R}
    pts = sorted(p for p in pts if cert.lo <= p <= cert.hi)
    tasks = [(p, q) for p, q in zip(pts, pts[1:]) if not any(L <= p and q <= R for (L, R) in wins)]
    cover = sorted(tasks + wins)
    pos = cert.lo
    for (L, R) in cover:
        require(L <= pos, f"coverage gap at {float(pos)}")
        pos = max(pos, R)
    require(pos >= cert.hi, "coverage does not reach e_N")
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
    E = A(cert.v).exp()
    print(f"VERIFIED: log lambda_A+ <= v = {cert.v} <= {up_decimal(A(cert.v), 10)}; lambda_A+ <= e^v <= {up_decimal(E, 7)} (also E_A+)  "
          f"({time.time() - t0:.0f}s)")
    sys.exit(0)


if __name__ == '__main__':
    try:
        main()
    except CheckError as e:
        print("CheckError:", e)
        print("NOT VERIFIED")
        sys.exit(1)
