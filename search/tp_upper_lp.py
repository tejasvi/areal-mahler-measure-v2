# UPPER bound for lambda+ (totally positive Mahler measure) via Smith's criterion (arXiv:2111.12660, Thm 1.5), in the checkable form
# of Talebizadeh Sardari-Orloski (arXiv:2304.10021): a probability measure mu on [0, X] (X/4 > 1) and weights b_Q >= 0 with
#   (i) U_mu(x) - sum_Q b_Q log|Q(x)| >= 0 on [0, X]  (=> on C by the minimum principle),  (ii) sum_Q b_Q deg Q < 1,
#   (iii) int log|Q| dmu > 0 for every Q in A  =>  int log|P| dmu >= 0 for all P in Z[x]  =>  lambda+ <= exp(int log+ dmu).
# Float prototype: hat densities on nodes e, LP over (p, b), constraints on a grid; log|Q| by root products (FLINT roots).
import os, numpy as np, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scipy.optimize import linprog
from tpc import stable_hat_potential, _hat_quad_nodes
from tp_primal_hat import logplus_hat
from tp_dual_rootprod import flint_roots, logabs_rootprod
from tpsets import load
def logQ_hat_flint(q, e):
    """int log|Q| phi_k for every hat k: exact for real roots (hat potentials), quadrature for complex roots; FLINT roots."""
    z = flint_roots(q); out = np.full(len(e) - 2, np.log(abs(q[0])))
    real = z[np.abs(z.imag) < 1e-12].real; cpx = z[np.abs(z.imag) >= 1e-12]
    if len(real): out += stable_hat_potential(e, real).sum(0)
    if len(cpx):
        Xq, Wq = _hat_quad_nodes(e)
        for rho in cpx: out += (np.log(np.abs(Xq - rho))*Wq).sum(1)
    return out
if __name__ == '__main__':
    setname, X, h, gridh, tag = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
    polys = load(setname); t0 = time.time()
    e = np.round(np.arange(0.0, X + 1e-12, h), 10); n = len(e) - 2
    roots = np.concatenate([flint_roots(q) for q in polys]); rr = np.sort(roots[np.abs(roots.imag) < 1e-12].real); rr = rr[(rr > 0) & (rr < X)]
    x = np.arange(gridh/2, X, gridh)
    near = [np.concatenate([r + np.geomspace(1e-6, 0.01, 30), r - np.geomspace(1e-6, 0.01, 30)]) for r in rr]
    x = np.unique(np.concatenate([x, e[1:-1], 0.5*(e[1:] + e[:-1])] + near)); x = x[(x > 0) & (x < X)]
    x = x[np.min(np.abs(x[:, None] - rr[None, :]), axis=1) > 1e-9]
    U = stable_hat_potential(e, x)                                   # (len(x), n): potential of each unit-mass hat at x
    L = np.array([logabs_rootprod(q, x) for q in polys]).T           # (len(x), J)
    LQ = np.array([logQ_hat_flint(q, e) for q in polys])             # (J, n): int log|Q| phi_k
    f = logplus_hat(e); degs = np.array([len(q) - 1 for q in polys], float); J = len(polys)
    print(f"{setname}: {J} polys, {n} hats (h = {h}, X = {X}), {len(x)} grid points ({time.time()-t0:.0f}s)", flush=True)
    # variables: [p (n), b (J)];  minimize f.p
    c = np.concatenate([f, np.zeros(J)])
    A1 = np.hstack([-U, L])                                          # -(U p - L b) <= 0
    A2 = np.concatenate([np.zeros(n), degs])[None, :]                # sum b deg <= 1 - 1e-6
    A3 = np.hstack([-LQ, np.zeros((J, J))])                          # -(LQ p) <= -1e-7
    Aub = np.vstack([A1, A2, A3]); bub = np.concatenate([np.zeros(len(x)), [1 - 1e-6], -1e-7*np.ones(J)])
    Aeq = np.concatenate([np.ones(n), np.zeros(J)])[None, :]
    res = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=[1.0], bounds=[(0, None)]*(n + J), method='highs',
                  options={'primal_feasibility_tolerance': 1e-10, 'dual_feasibility_tolerance': 1e-10})
    print("status", res.status, res.message, flush=True)
    if res.status == 0:
        p, b = res.x[:n], res.x[n:]; val = res.fun
        sup = e[1:-1][p > 1e-9*p.max()]
        hmin = (U @ p - L @ b).min()
        print(f"UPPER LP value {val:.7f} (Omega {np.exp(val):.7f}); support [{sup.min():.3f}, {sup.max():.3f}]; min h on grid {hmin:.2e}; "
              f"sum b deg {degs @ b:.6f}; int log|Q| dmu min {(LQ @ p).min():.2e} ({time.time()-t0:.0f}s)", flush=True)
        print("   b:", [(j, round(float(v), 7)) for j, v in enumerate(b) if v > 1e-9], flush=True)
        np.save(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'upper_{tag}.npy'), np.concatenate([[val, len(e)], e, p, b]))
