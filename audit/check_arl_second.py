#!/usr/bin/env python3
"""check_arl_second.py -- second, independent checker for the areal certificates (formats tpmahler-areal-lower-v1 and
tpmahler-areal-upper-v1). Written from the two certificate formats (with log+ replaced by g) and the two criteria of the
paper (the lower-bound lemma and the upper-bound lemma), without reading the first checkers. Its method for the upper
format is that of the second checker of the companion paper on the totally positive case.

    python3 -I check_arl_second.py CERT.json NPROC

g(x) = log+ x + (min(1, x)^2 - 1)/2 = (x^2 - 1)/2 on [0, 1], log x on [1, inf)   (Pritsker's areal potential on R+).

UPPER (tpmahler-areal-upper-v1): (U1) sum b_j deg Q_j < s; (U2) int log|Q_j| dmu > 0 for weighted Q_j; (U3) h = U - sum b_j
log|Q_j| >= 0 on [e_0, e_N]; (V) (1/s) int g dmu <= v (closed form: exact rational on [e_k, min(e_{k+1}, 1)], arb above 1).
LOWER (tpmahler-areal-lower-v1): budget sum c_j deg Q_j + 2 s < 1; I = sum_{i,j} C_i C_j J(e_i - e_j), J(t) = t^4/24 log|t|
- 25 t^4/288 (J'' = H, J'''' = log|t|), a rigorous lower bound is used; F(x) = g(x) - m - sum c_j log|Q_j(x)| - 2 U(x) + I/s
>= 0 on (0, inf): left tail (0, x_lo] (g >= -1/2, U <= U(0) left of the support, c_x log x, max |Q_j| by ball arithmetic),
right tail [x_hi, inf) (F >= (1 - sum c deg - 2 s) log x - m - sum c_j log sum_k |q_jk| + I/s), main range by rigorous
quadratic minorants + bisection: F(x) >= F(m) + F'(m)u + inf_J R''u^2/2 + chords of the concave singular node terms
(-2 C_i H(e_i - x) with C_i < 0 near e_i), g handled exactly on each side of 1 (g'' = 1 on (0,1), -1/x^2 on (1,inf),
intervals never cross 1), weighted real roots cut out by 1e-40 brackets with a log bound.
"""
import sys, json, time, math
from fractions import Fraction
from multiprocessing import Pool

if sys.flags.optimize:
    sys.exit("refusing to run under python -O")

from flint import arb, acb, fmpq, fmpz, fmpz_poly, ctx

PREC = 128
ctx.prec = PREC
DNEAR = fmpq(1, 20)
BRACKET = fmpq(1, 10 ** 40)
MAXTESTS = 300000
ONE = fmpq(1)


def q_of(s):
    if isinstance(s, bool) or not isinstance(s, (str, int)):
        raise ValueError(f"bad rational {s!r}")
    fr = Fraction(s)
    return fmpq(fr.numerator, fr.denominator)


def fq_of_exact_arb(x):
    m, e = x.mid().man_exp()
    m, e = int(m), int(e)
    return fmpq(m * 2 ** e) if e >= 0 else fmpq(m, 2 ** (-e))


def qabs(x):
    return x if x >= 0 else -x


def qmax(*xs):
    r = xs[0]
    for x in xs[1:]:
        if x > r:
            r = x
    return r


def qmin(*xs):
    r = xs[0]
    for x in xs[1:]:
        if x < r:
            r = x
    return r


def H_ball(t):
    """H(t) = t^2/2 log|t| - 3t^2/4 (also for balls containing 0)."""
    if not t.contains(0):
        L = abs(t).log()
        return t * t * (L / 2 - arb(3) / 4)
    r = t.abs_upper()
    if r == 0:
        return arb(0)
    if not (r < arb(1) / 3):
        return arb("nan")
    lo = -(r * r / 2) * abs(r.log()) - 3 * r * r / 4
    return (lo / 2) + arb(0, abs(lo / 2).abs_upper())


def Hc(t):
    return t * t * (t.log() / 2 - acb(arb(3) / 4))


def J_ball(t):
    """J(t) = t^4/24 log|t| - 25 t^4/288, J(0) = 0."""
    if t == 0:
        return arb(0)
    t2 = t * t
    return t2 * t2 * (abs(t).log() / 24 - arb(25) / 288)


D = {}


def setup(path):
    d = json.load(open(path))
    fmt = d.get("format")
    errs = []
    if fmt == "tpmahler-areal-upper-v1":
        mode, wkey = "upper", "b"
    elif fmt == "tpmahler-areal-lower-v1":
        mode, wkey = "lower", "c"
    else:
        raise SystemExit(f"unknown format {fmt!r}")
    polys = []
    for c in d["polys"]:
        if not isinstance(c, list) or not c or any(isinstance(a, bool) or not isinstance(a, int) for a in c):
            errs.append(f"poly coefficients not integers: {c}")
            continue
        polys.append(c)
    w = [q_of(x) for x in d[wkey]]
    e = [q_of(x) for x in d["nodes"]]
    G = [q_of(x) for x in d["G"]]
    if len(w) != len(polys):
        errs.append("weights/polys length mismatch")
    Q, seen = [], set()
    for j, c in enumerate(polys):
        if c[0] == 0:
            errs.append(f"Q_{j}: leading coefficient 0")
        P = fmpz_poly(list(reversed(c)))
        if P.degree() < 1:
            errs.append(f"Q_{j}: degree < 1")
        cont = 0
        for a in c:
            cont = math.gcd(cont, a)
        if mode == "upper" and cont != 1:
            errs.append(f"Q_{j}: not primitive")
        if mode == "upper":
            fc, fl = P.factor()
            if not (len(fl) == 1 and fl[0][1] == 1 and fl[0][0].degree() == P.degree()):
                errs.append(f"Q_{j}: reducible over Q")
        key = tuple(c) if c[0] > 0 else tuple(-a for a in c)
        if key in seen:
            errs.append(f"Q_{j}: duplicate (up to sign)")
        seen.add(key)
        Q.append(P)
    if any(x < 0 for x in w):
        errs.append("negative weight")
    if len(e) != len(G) or len(e) < 2 or e[0] < 0 or any(not e[k] < e[k + 1] for k in range(len(e) - 1)):
        errs.append("nodes")
    if any(gv < 0 for gv in G) or G[0] != 0 or G[-1] != 0:
        errs.append("G")
    s = fmpq(0)
    for k in range(len(e) - 1):
        s += (e[k + 1] - e[k]) * (G[k] + G[k + 1]) / 2
    if not s > 0:
        errs.append("mass s <= 0")
    N = len(e) - 1
    sig = [(G[k + 1] - G[k]) / (e[k + 1] - e[k]) for k in range(N)]
    C = [(sig[i] if i < N else fmpq(0)) - (sig[i - 1] if i > 0 else fmpq(0)) for i in range(N + 1)]
    sigma = fmpq(1) if mode == "upper" else fmpq(-2)
    D.update(d=d, mode=mode, polys=polys, w=w, e=e, G=G, Q=Q, s=s, sig=sig, C=C, N=N, errs=errs, sigma=sigma)
    D["Ea"] = [arb(x) for x in e]
    D["Ca"] = [arb(x) for x in C]
    D["Cs"] = [sigma * x for x in C]
    D["Csa"] = [arb(sigma * x) for x in C]
    D["W"] = [j for j in range(len(Q)) if w[j] > 0]
    D["wa"] = [arb(x) for x in w]
    D["dQ"] = [P.derivative() for P in Q]
    if mode == "upper":
        D["X"] = q_of(d["X"])
        D["v"] = q_of(d["v"])
        if not D["X"] > 4 or not e[-1] <= D["X"]:
            errs.append("X")
    else:
        D["m"] = q_of(d["m"])
        D["x_lo"] = q_of(d["x_lo"])
        D["x_hi"] = q_of(d["x_hi"])
        xi = [j for j, c in enumerate(polys) if c in ([1, 0], [-1, 0])]
        if not xi or not w[xi[0]] > 0:
            errs.append("the polynomial x must be present with positive weight")
        D["jx"] = xi[0] if xi else None
    ctx.prec = 256
    roots = {}
    for j in D["W"]:
        rr = []
        for (r, mult) in Q[j].complex_roots():
            if r.imag == 0:
                mid = fq_of_exact_arb(r.real)
                lo, hi = mid - BRACKET, mid + BRACKET
                if not (r.real - arb(lo) > 0 and arb(hi) - r.real > 0):
                    errs.append(f"Q_{j}: real root ball not inside bracket")
                rr.append(("r", r, lo, hi, mult))
            else:
                rr.append(("c", r, None, None, mult))
        roots[j] = rr
    ctx.prec = PREC
    D["roots"] = roots


# ------------------------------------------------------------------------------------------------ potential pieces
def U_real(x):
    tot = arb(0)
    for Ei, Ci in zip(D["Ea"], D["Ca"]):
        tot += Ci * H_ball(Ei - x)
    return tot


def U_complex(z):
    tot = acb(0)
    for Ei, Ci in zip(D["Ea"], D["Ca"]):
        tot += acb(Ci) * Hc(acb(Ei) - z)
    return tot.real


def int_g_dmu():
    """int g(x) G(x) dx: exact rational part on [0, 1], arb part on [1, inf)."""
    e, G, sig = D["e"], D["G"], D["sig"]
    exact = fmpq(0)
    logpart = arb(0)
    for k in range(len(e) - 1):
        a_, b_ = e[k], e[k + 1]
        beta = sig[k]
        alpha = G[k] - beta * a_
        lo, hi = a_, qmin(b_, ONE)
        if hi > lo:     # int (alpha + beta x)(x^2 - 1)/2 dx
            P = lambda x: (alpha * (x ** 3 / 3 - x) + beta * (x ** 4 / 4 - x ** 2 / 2)) / 2
            exact += P(hi) - P(lo)
        lo, hi = qmax(a_, ONE), b_
        if hi > lo:
            A, B = arb(alpha), arb(beta)

            def prim(x):
                xa = arb(x)
                L = xa.log()
                return A * (xa * L - xa) + B * (xa * xa / 2 * L - xa * xa / 4)
            logpart += prim(hi) - prim(lo)
    return exact, logpart


def energy_lower():
    """Rigorous ball for I = sum_{i,j} C_i C_j J(e_i - e_j) (symmetric: 2 sum_{i<j})."""
    e, Ca = D["e"], D["Ca"]
    N = len(e)
    tot = arb(0)
    for i in range(N):
        if D["C"][i] == 0:
            continue
        row = arb(0)
        for j in range(i + 1, N):
            if D["C"][j] == 0:
                continue
            row += Ca[j] * J_ball(arb(e[j] - e[i]))
        tot += Ca[i] * row
    return 2 * tot


# ------------------------------------------------------------------------------------------------ interval minorant
FACT = [1, 1, 2, 6, 24, 120, 720, 5040]


def far_model(k):
    """Taylor model (order 6) of sum_{far} Cs_i H(e_i - x) on node interval k."""
    e, N = D["e"], D["N"]
    c = (e[k] + e[k + 1]) / 2
    hk = (e[k + 1] - e[k]) / 2
    near = [i for i in range(N + 1) if e[k] - DNEAR <= e[i] <= e[k + 1] + DNEAR]
    nearset = set(near)
    Dn = [arb(0) for _ in range(7)]
    M7 = arb(0)
    hka = arb(hk)
    for i in range(N + 1):
        if i in nearset or D["C"][i] == 0:
            continue
        t = arb(e[i] - c)
        Ci = D["Csa"][i]
        L = abs(t).log()
        Hn = [t * t * (L / 2 - arb(3) / 4), t * (L - 1), L, 1 / t, -1 / (t * t), 2 / (t ** 3), -6 / (t ** 4)]
        for n in range(7):
            Dn[n] += Ci * Hn[n] * (1 if n % 2 == 0 else -1)
        M7 += abs(Ci) * 24 / (abs(t) - hka) ** 5
    return (c, hk, near, Dn, M7)


def poly_eval_taylor(Dn, start, u):
    tot = arb(0)
    for n in range(6, start - 1, -1):
        tot = tot * u + Dn[n] / FACT[n - start]
    return tot


def q_d2_lower(p, q):
    """Lower bound over [p, q] of sum_j w_j sum_rho Re 1/(x - rho)^2 (the second derivative of -sum w_j log|Q_j|)."""
    xb = acb(arb((p + q) / 2, (q - p) / 2))
    tot = arb(0)
    for j in D["W"]:
        wj = D["wa"][j]
        for (kind, r, lo, hi, mult) in D["roots"][j]:
            if kind == "r":
                dmax = qmax(qabs(p - lo), qabs(p - hi), qabs(q - lo), qabs(q - hi))
                tot += wj * mult / arb(dmax) ** 2
            else:
                val = (1 / ((xb - r) * (xb - r))).real
                if not val.is_finite():
                    return arb("nan")
                tot += wj * mult * val.lower()
    return tot


def quad_lb_positive(a0, a1, a2, w):
    hw = arb(w) / 2
    cands = [a0 - a1 * hw + a2 * hw * hw, a0 + a1 * hw + a2 * hw * hw]
    if a2 > 0:
        ustar = -a1 / (2 * a2)
        if not (ustar > hw or ustar < -hw):
            cands.append(a0 - a1 * a1 / (4 * a2))
    elif not (a2 <= 0):
        cands.append(a0 - abs(a1) * hw - abs(a2) * hw * hw)
    ok = all(cc > 0 for cc in cands)
    lows = [float(cc.lower()) if cc.is_finite() else float("-inf") for cc in cands]
    return ok, min(lows)


def lb_sub(p, q, model):
    """Rigorous test F > 0 (lower) / h > 0 (upper) on [p, q]; model = node-interval Taylor model or None (direct)."""
    e, Cs, Csa = D["e"], D["Cs"], D["Csa"]
    w = q - p
    m = (p + q) / 2
    four_w = 4 * w
    half_w = arb(w) / 2
    if model is not None:
        c, hk, near, Dn, M7 = model
        um = arb(m - c)
        uJ = arb(m - c, w / 2)
        aum = abs(um)
        Rm = poly_eval_taylor(Dn, 0, um) + arb(0, (M7 * aum ** 7 / 5040).abs_upper())
        Rpm = poly_eval_taylor(Dn, 1, um) + arb(0, (M7 * aum ** 6 / 720).abs_upper())
        uJmax = abs(uJ).abs_upper()
        F2 = poly_eval_taylor(Dn, 2, uJ) + arb(0, (M7 * uJmax ** 5 / 120).abs_upper())
        D2 = F2.lower() if F2.is_finite() else arb("nan")
        idx = near
    else:
        Rm, Rpm, D2 = arb(0), arb(0), arb(0)
        idx = range(len(e))
    Sp, Sq, Sm = arb(0), arb(0), arb(0)
    for i in idx:
        ei, Ci = e[i], Cs[i]
        if Ci == 0:
            continue
        if ei < p:
            dmin, dmax = p - ei, q - ei
        elif ei > q:
            dmin, dmax = ei - q, ei - p
        else:
            dmin, dmax = fmpq(0), qmax(ei - p, q - ei)
        if Ci > 0 and dmin <= four_w:
            if not dmax < fmpq(1, 2):
                return False, float("-inf"), False, arb("nan")
            Sp += Csa[i] * H_ball(arb(ei - p))
            Sq += Csa[i] * H_ball(arb(ei - q))
            Sm += Csa[i] * H_ball(arb(ei - m))
            continue
        t = arb(ei - m)
        L = abs(t).log()
        Rm += Csa[i] * (t * t * (L / 2 - arb(3) / 4))
        Rpm -= Csa[i] * (t * (L - 1))
        if dmin > four_w:
            D2 += Csa[i] * L - abs(Csa[i]) * half_w / arb(dmin)
        else:          # here Cs_i < 0
            D2 += Csa[i] * arb(dmax).log()
    for j in D["W"]:
        Qm = D["Q"][j](m)
        dQm = D["dQ"][j](m)
        Rm -= D["wa"][j] * arb(qabs(Qm)).log()
        Rpm -= D["wa"][j] * arb(dQm / Qm)
    D2 += q_d2_lower(p, q)
    if D["mode"] == "lower":
        ma = arb(m)
        if q <= ONE:
            Rm += (ma * ma - 1) / 2
            Rpm += ma
            D2 += 1
        elif p >= ONE:
            Rm += ma.log()
            Rpm += 1 / ma
            D2 += -1 / (arb(p) * arb(p))
        else:
            raise RuntimeError("interval crosses 1")
        Rm += D["const"]
    a0 = Rm + (Sp + Sq) / 2
    a1 = Rpm + (Sq - Sp) / arb(w)
    a2 = D2 / 2
    ok, lb = quad_lb_positive(a0, a1, a2, w)
    hm = Rm + Sm
    return ok, lb, bool(hm < 0), hm


def lb_bracket(p, q, j0, lo, hi):
    """F > 0 (h > 0) on [p, q] inside the (merged) root bracket [lo, hi]; j0 = list of the weighted Q_j with a root in it."""
    js = j0 if isinstance(j0, list) else [j0]
    xb = arb((p + q) / 2, (q - p) / 2)
    xr = arb((lo + hi) / 2, (hi - lo) / 2)
    val = D["sigma"] * U_real(xb)
    for j in D["W"]:
        if j in js:
            M1 = abs(D["dQ"][j](xr)).abs_upper()
            val -= D["wa"][j] * (M1 * arb(hi - lo)).log()
        else:
            val -= D["wa"][j] * abs(D["Q"][j](xb)).log()
    if D["mode"] == "lower":
        if q <= ONE:
            val += (xb * xb - 1) / 2
        elif p >= ONE:
            val += xb.log()
        else:                      # bracket around 1: g is increasing, so g >= g(p) = (p^2 - 1)/2 on [p, q]
            pa = arb(p)
            val += (pa * pa - 1) / 2
        val += D["const"]
    return bool(val > 0), float(val.lower()) if val.is_finite() else float("-inf")


def pieces_of(a, b):
    """Split [a, b] at 1 (lower mode) and at the weighted root brackets; returns (p, q, bracket-or-None)."""
    cuts = []
    if D["mode"] == "lower" and a < ONE < b:
        cuts.append((ONE, ONE, None))
    for (lo, hi, j) in D["brackets"]:
        if hi > a and lo < b:
            cuts.append((qmax(lo, a), qmin(hi, b), (j, lo, hi)))
    cuts.sort(key=lambda t: (t[0], t[1]))
    out, cur = [], a
    for (p1, q1, jb) in cuts:
        if p1 > cur:
            out.append((cur, p1, None))
        if jb is not None:
            out.append((p1, q1, jb))
            cur = q1
        else:
            cur = qmax(cur, p1)
    if cur < b:
        out.append((cur, b, None))
    return out


def run_interval(a, b, model, tag, wmin=fmpq(1, 10 ** 11)):
    n, minlb, where, fails, depth_max, covered = 0, float("inf"), None, [], 0, fmpq(0)
    for (p0, q0, jb) in pieces_of(a, b):
        if jb is not None:
            ok, lb = lb_bracket(p0, q0, jb[0], jb[1], jb[2])
            n += 1
            covered += q0 - p0
            if lb < minlb:
                minlb, where = lb, float(p0)
            if not ok:
                fails.append((float(p0), f"bracket of Q_{jb[0]}"))  # jb[0] = list of polynomials
            continue
        stack = [(p0, q0, 0)]
        while stack:
            p, q, dep = stack.pop()
            ok, lb, neg, hm = lb_sub(p, q, model)
            n += 1
            if neg:
                fails.append((float((p + q) / 2), f"value {hm.str(5)} < 0"))
                covered += q - p + sum((qq - pp for (pp, qq, dd) in stack), fmpq(0))
                stack = []
                break
            if n > MAXTESTS:
                fails.append((float((p + q) / 2), f"more than {MAXTESTS} tests in {tag}"))
                covered += q - p + sum((qq - pp for (pp, qq, dd) in stack), fmpq(0))
                stack = []
                break
            if ok:
                covered += q - p
                if lb < minlb:
                    minlb, where = lb, float((p + q) / 2)
                depth_max = max(depth_max, dep)
                continue
            if q - p < wmin:
                fails.append((float((p + q) / 2), f"lower {lb:.3e} at width {float(q - p):.1e}"))
                covered += q - p
                continue
            mid = (p + q) / 2
            stack.append((mid, q, dep + 1))
            stack.append((p, mid, dep + 1))
    return dict(tag=tag, n=n, minlb=minlb, where=where, fails=fails, depth=depth_max, covered=covered)


def task(spec):
    kind, a, b, k = spec
    model = far_model(k) if kind == "node" else None
    return run_interval(a, b, model, f"{kind} [{float(a):.6g}, {float(b):.6g}]")


def main():
    t0 = time.time()
    path = sys.argv[1]
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    setup(path)
    mode, e, Q, w, s, C, N, errs = D["mode"], D["e"], D["Q"], D["w"], D["s"], D["C"], D["N"], D["errs"]
    print(f"{path}: {D['d']['format']}, {len(Q)} polynomials ({len(D['W'])} weighted), {N + 1} nodes on [{float(e[0])}, "
          f"{float(e[-1])}], s = {float(s)!r}")
    ok_all = not errs
    m0 = sum(C, fmpq(0)); m1 = sum((C[i] * e[i] for i in range(N + 1)), fmpq(0)); m2 = sum((C[i] * e[i] ** 2 for i in range(N + 1)), fmpq(0))
    mom = m0 == 0 and m1 == 0 and m2 == 2 * s
    print(f"moments exact: {mom}; data errors: {errs if errs else 'none'}")
    ok_all = ok_all and mom
    degsum = sum((w[j] * Q[j].degree() for j in range(len(Q))), fmpq(0))
    if mode == "upper":
        u1 = degsum < s
        print(f"(U1) sum b deg = {float(degsum):.10f} < s = {float(s):.10f}: {u1}")
        ok_all = ok_all and u1
        ctx.prec = 256
        u2min, u2ok = (float("inf"), None), True
        for j in D["W"]:
            tot = arb(s) * arb(abs(int(Q[j].coeffs()[-1]))).log()
            for (kind, r, lo, hi, mult) in D["roots"][j]:
                tot += (U_complex(r) if kind == "c" else U_real(r.real)) * mult
            u2ok = u2ok and bool(tot > 0)
            if float(tot.mid()) < u2min[0]:
                u2min = (float(tot.mid()), j)
        print(f"(U2) all weighted: {u2ok}; smallest int log|Q| dmu = {u2min[0]:.4e} (Q_{u2min[1]})")
        exact, logpart = int_g_dmu()
        vex = (arb(exact) + logpart) / arb(s)
        vok = bool(arb(D["v"]) - vex > 0)
        ctx.prec = PREC
        print(f"(V) (1/s) int g dmu = {vex.str(15)} (exact rational part {float(exact):.12f}); claim v = {float(D['v'])!r}; "
              f"margin {float((arb(D['v']) - vex).mid()):.3e}: {vok}")
        ok_all = ok_all and u2ok and vok
        dom_a, dom_b = e[0], e[-1]
    else:
        kappa = 1 - degsum - 2 * s
        bud = kappa > 0
        print(f"budget sum c deg + 2s = {float(degsum + 2 * s):.12f} < 1: {bud}; kappa = {float(kappa):.10f}")
        ctx.prec = 192
        I = energy_lower()
        ctx.prec = PREC
        print(f"energy I = {I.str(12)}")
        Ilo = I.lower()
        D["const"] = -arb(D["m"]) + Ilo / arb(s)
        ok_all = ok_all and bud and I.is_finite()
        x_lo, x_hi = D["x_lo"], D["x_hi"]
        # left tail (0, x_lo]: U is decreasing on [0, a0], a0 = left end of the support of G (G = 0 on [e_0, a0])
        k0 = next(k for k in range(N + 1) if D["G"][k] > 0)
        a0 = e[k0 - 1]
        lt_ok = x_lo <= a0 and x_lo < ONE
        val = arb(-1) / 2 - arb(D["m"]) + Ilo / arb(s) - 2 * U_real(arb(0))
        for j in D["W"]:
            if j == D["jx"]:
                val -= D["wa"][j] * arb(x_lo).log()
                continue
            for (kind, r, lo, hi, mult) in D["roots"][j]:
                if kind == "r" and hi > 0 and lo < x_lo:
                    lt_ok = False
            mx = abs(D["Q"][j](arb(x_lo / 2, x_lo / 2))).abs_upper()
            val -= D["wa"][j] * mx.log()
        lt_ok = lt_ok and bool(val > 0)
        print(f"left tail (0, {x_lo}] (support of G starts at {float(a0)}): F >= {val.str(6)}: {lt_ok}")
        # right tail [x_hi, inf)
        rt_ok = x_hi >= ONE and x_hi >= e[-1]
        val = arb(kappa) * arb(x_hi).log() - arb(D["m"]) + Ilo / arb(s)
        for j in D["W"]:
            val -= D["wa"][j] * arb(sum(abs(int(a)) for a in D["polys"][j])).log()
            for (kind, r, lo, hi, mult) in D["roots"][j]:
                if kind == "r" and hi > x_hi:
                    rt_ok = False
        rt_ok = rt_ok and bool(val > 0) and bud
        print(f"right tail [{x_hi}, inf): F >= {val.str(6)} (increasing): {rt_ok}")
        ok_all = ok_all and lt_ok and rt_ok
        dom_a, dom_b = x_lo, x_hi
    # brackets of weighted real roots inside the domain
    br = []
    for j in D["W"]:
        for (kind, r, lo, hi, mult) in D["roots"][j]:
            if kind == "r" and hi > dom_a and lo < dom_b:
                if lo <= dom_a or hi >= dom_b:
                    print(f"root bracket of Q_{j} straddles the domain end"); ok_all = False
                br.append((lo, hi, j))
    br.sort()
    merged = []               # merge overlapping brackets (weighted polynomials sharing a root, e.g. reducible ones)
    for (lo, hi, j) in br:
        if merged and lo <= merged[-1][1]:
            mlo, mhi, js = merged[-1]
            merged[-1] = (mlo, qmax(mhi, hi), js + [j])
        else:
            merged.append((lo, hi, [j]))
    D["brackets"] = merged
    # tasks: node intervals (Taylor model) inside the support, direct intervals outside
    specs = []
    if dom_a < e[0]:
        specs.append(("direct", dom_a, e[0], None))
    for k in range(N):
        a_, b_ = qmax(e[k], dom_a), qmin(e[k + 1], dom_b)
        if b_ > a_:
            specs.append(("node", a_, b_, k))
    if dom_b > e[-1]:
        a = e[-1]
        while a < dom_b:
            b = qmin(dom_b, qmax(a * 2, a + 1))
            specs.append(("direct", a, b, None))
            a = b
    t1 = time.time()
    if nproc > 1:
        with Pool(nproc) as pool:
            res = pool.map(task, specs, chunksize=4)
    else:
        res = [task(sp) for sp in specs]
    nint = sum(r["n"] for r in res)
    fails = [f for r in res for f in r["fails"]]
    covered = sum((r["covered"] for r in res), fmpq(0))
    best = min(res, key=lambda r: r["minlb"])
    cov_ok = covered == dom_b - dom_a
    name = "(U3)" if mode == "upper" else "F >= 0 on the main range"
    print(f"{name} [{float(dom_a)}, {float(dom_b)}]: {nint} interval tests, coverage exact: {cov_ok}, failures: {len(fails)}, "
          f"brackets {len(br)}, smallest certified lower bound {best['minlb']:.4e} near {best['where']}  [{time.time() - t1:.0f}s]")
    for f in fails[:20]:
        print("   FAILED near", f)
    ok_all = ok_all and cov_ok and not fails
    if ok_all:
        ctx.prec = 200
        if mode == "upper":
            ev = arb(D["v"]).exp()
            print(f"VERIFIED: log lambda_A+ <= v = {D['v']}; lambda_A+ <= e^v = {ev.str(15)}")
        else:
            em = arb(D["m"]).exp()
            print(f"VERIFIED: log lambda_A+ >= m = {D['m']} = {arb(D['m']).str(15)}; lambda_A+ >= e^m = {em.str(15)}")
    else:
        print("NOT VERIFIED")
    print(f"[{time.time() - t0:.0f}s]")
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()
