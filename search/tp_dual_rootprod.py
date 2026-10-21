# Fixed-reference dual LP for a primal on arbitrary nodes (primalN_<tag>.npy = [val, len(e), e..., p...]).
# ROOT-PRODUCT version of tp_dual_fixedN_par.py: log|Q(x)| = log|lead| + sum log|x - zeta| over FLINT-isolated roots (float64
# Horner is wrong for the degree-32..38 Q_25..Q_27 near x ~ 1.95); the near-root grid also uses these roots. Output *_dualR*.
# Dual polynomial set: 'same' (the primal's subset) or a set name (e.g. F27). Numbers budget (H2') imposed unless 'integers'.
import os, numpy as np, sys, json, time
from scipy.optimize import linprog
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tpc import stable_hat_energy_matrix, stable_hat_potential
from tpsets import load, ROOT
from flint import fmpz_poly, acb
def flint_roots(q):
    """complex roots (with multiplicity) of the integer polynomial q (descending coefficients), to double precision."""
    out = []
    for r, mult in fmpz_poly([int(v) for v in q[::-1]]).complex_roots():
        out += [complex(float(r.real.mid()), float(r.imag.mid()))]*mult
    return np.array(out, dtype=complex)
def logabs_rootprod(q, x):
    z = flint_roots(q); v = np.full(len(x), np.log(abs(q[0])))
    for zz in z: v += np.log(np.abs(x - zz))
    return v
_G = {}
def _uchunk(x):
    return stable_hat_potential(_G['e'], x) @ _G['p']
def Unu(e, p, x, chunk=2000, nproc=None):
    """U_nu on many points, parallel over chunks (fork inherits e, p)."""
    import multiprocessing as mp, os
    _G['e'], _G['p'] = e, p
    nproc = nproc or int(os.environ.get('TP_NPROC', '32'))
    os.environ['OMP_NUM_THREADS'] = '1'
    pieces = [x[i0:i0 + chunk] for i0 in range(0, len(x), chunk)]
    with mp.get_context('fork').Pool(nproc) as pool:
        res = pool.map(_uchunk, pieces)
    return np.concatenate(res)
if __name__ == '__main__':
    prim, dualset = sys.argv[1], sys.argv[2]; numbers = not (len(sys.argv) > 3 and sys.argv[3] == 'integers')
    XF = 80.0
    a = np.load(prim); ne = int(a[1]); e = a[2:2 + ne]; p = a[2 + ne:]
    if dualset == 'same':
        idx = json.load(open(prim.replace('.npy', '_idx.json'))); P = load('F27'); polys = [P[j] for j in idx]
    else:
        polys = load(dualset); idx = list(range(len(polys)))
    degs = np.array([len(q) - 1 for q in polys], float)
    t0 = time.time(); K = stable_hat_energy_matrix(e); E = p @ K @ p
    roots = np.concatenate([flint_roots(q) for q in polys]); rr = roots[np.abs(roots.imag) < 1e-12].real; rr = rr[rr > 0]
    Xs = e[-1]
    x = np.unique(np.concatenate([np.geomspace(1e-7, 0.02, 3000), np.arange(0.02, Xs + 0.5, 0.0005), e[1:-1],
                                  0.5*(e[1:] + e[:-1]), np.geomspace(Xs + 0.5, XF, 3000)]))
    near = [np.concatenate([r + np.geomspace(1e-7, 0.02, 400), r - np.geomspace(1e-7, 0.02, 400)]) for r in rr if 0 < r < Xs + 0.5]
    sub = np.linspace(0, 1, 9)[1:-1]
    nodal = (e[:-1, None] + (e[1:] - e[:-1])[:, None]*sub[None, :]).ravel()
    x = np.unique(np.concatenate([x] + near + [nodal]))
    x = x[(x > 0)]
    x = x[np.min(np.abs(x[:, None] - rr[None, :]), axis=1) > 1e-9] if len(rr) else x
    U = Unu(e, p, x)
    L = np.array([logabs_rootprod(q, x) for q in polys]).T
    f = np.log(np.maximum(1.0, x)); ok = np.all(np.isfinite(L), axis=1); x, U, L, f = x[ok], U[ok], L[ok], f[ok]
    print(f"  E(nu) = {E:.4e}; {len(x)} grid points ({time.time()-t0:.0f}s)", flush=True)
    J = len(polys); sel = np.arange(0, len(x), 20)
    rows_extra = [np.concatenate([[0.0], degs, [1.0]])]; rhs_extra = [1.0 - 1e-6]
    if numbers: rows_extra.append(np.concatenate([[0.0], degs, [2.0]])); rhs_extra.append(1.0 - 1e-6)
    for it in range(60):
        A = np.hstack([np.ones((len(sel), 1)), L[sel], (2*U[sel] - E)[:, None]])
        A = np.vstack([A] + [r[None, :] for r in rows_extra])
        res = linprog(np.concatenate([[-1.0], np.zeros(J + 1)]), A_ub=A, b_ub=np.concatenate([f[sel], rhs_extra]),
                      bounds=[(None, None)] + [(0, None)]*(J + 1), method='highs',
                      options={'primal_feasibility_tolerance': 1e-10, 'dual_feasibility_tolerance': 1e-10})
        assert res.status == 0, res.message
        m, c, lam = res.x[0], res.x[1:J + 1], res.x[J + 1]
        F = f - m - L @ c - lam*(2*U - E); Fmin = F.min()
        print(f"  it {it}: m = {m:.8f} (Omega {np.exp(m):.7f}) lam = {lam:.5f} min F = {Fmin:.2e} at x = {x[np.argmin(F)]:.5f}", flush=True)
        if Fmin > -1e-10: break
        viol = np.nonzero(F < -1e-10)[0]
        lm = [i for i in viol if (i == 0 or F[i] <= F[i-1]) and (i == len(F)-1 or F[i] <= F[i+1])]
        sel = np.unique(np.concatenate([sel, lm]))
    print(f"RESULT m = {m:.8f} Omega = {np.exp(m):.7f}; sum c deg + lam = {degs @ c + lam:.4f}; sum c deg + 2 lam = {degs @ c + 2*lam:.4f}")
    # the 5 worst local minima of F (diagnostic for node refinement)
    lmins = [i for i in range(1, len(F) - 1) if F[i] <= F[i-1] and F[i] <= F[i+1]]
    worst = sorted(lmins, key=lambda i: F[i])[:12]
    print("lowest local minima of F:", [(round(float(x[i]), 5), f"{F[i]:.2e}") for i in sorted(worst)])
    tag = 'N' if numbers else 'I'
    np.save(prim.replace('.npy', f'_dualR{tag}_{dualset}.npy'), np.concatenate([[m, lam, E], c]))
    json.dump(idx, open(prim.replace('.npy', f'_dualR{tag}_{dualset}_idx.json'), 'w'))
