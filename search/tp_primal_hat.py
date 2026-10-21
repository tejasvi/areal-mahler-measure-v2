# Primal SOCP over hat densities on uniform nodes e_k = k*h (k = 0..N), hats k = 1..N-1 (unit mass):
#   minimize int log+ x dmu  s.t.  int log|Q_j| dmu >= 0 (all j),  I(mu) >= 0.
import os, numpy as np, cvxpy as cp, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tpc import stable_hat_energy_matrix, stable_hat_potential, _hat_quad_nodes
from tpsets import load, ROOT
def logQ_hat(p, e):
    """int log|Q(x)| phi_k(x) dx for every hat k (exact for real roots via hat potentials, quadrature for complex roots)."""
    r = np.roots(np.array(p, float)); lead = abs(p[0])
    out = np.full(len(e) - 2, np.log(lead))
    real = r[np.abs(r.imag) < 1e-9].real; cpx = r[np.abs(r.imag) >= 1e-9]
    if len(real): out += stable_hat_potential(e, real).sum(0)
    if len(cpx):
        Xq, Wq = _hat_quad_nodes(e)
        for rho in cpx: out += (np.log(np.abs(Xq - rho))*Wq).sum(1)
    return out
def logplus_hat(e):
    Xq, Wq = _hat_quad_nodes(e)
    return (np.log(np.maximum(1.0, Xq))*Wq).sum(1)
def solve(e, polys, energy=True, C=None):
    n = len(e) - 2
    f = logplus_hat(e); L = np.array([logQ_hat(p, e) for p in polys])
    p = cp.Variable(n, nonneg=True); cons = [cp.sum(p) == 1, L @ p >= 0]
    if energy:
        K = stable_hat_energy_matrix(e)
        if C is None: C = max(1.0, 2*np.abs(K).max())
        P = C*np.ones((n, n)) - K; ev, V = np.linalg.eigh(0.5*(P + P.T)); R = (V*np.sqrt(np.maximum(ev, 0))).T
        cons.append(cp.sum_squares(R @ p) <= C)
    prob = cp.Problem(cp.Minimize(f @ p), cons); prob.solve(solver='CLARABEL')
    return prob.value, np.maximum(p.value, 0), [c.dual_value for c in cons], prob.status
if __name__ == '__main__':
    setname, X, N, tag = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    polys = load(setname)
    e = np.linspace(0.0, X, N + 1)
    t0 = time.time()
    for energy in (False, True):
        val, pv, du, st = solve(e, polys, energy)
        sup = e[1:-1][pv > 1e-7*pv.max()]
        print(f"{setname} X={X} N={N} energy={energy}: primal {val:.7f} (Omega {np.exp(val):.7f}) [{st}] support [{sup.min():.3f}, {sup.max():.3f}] ({time.time()-t0:.0f}s)", flush=True)
        if energy:
            print("   duals c_j:", [(j, round(float(c), 7)) for j, c in enumerate(du[1]) if c > 1e-8], " lambda_E:", du[2], flush=True)
        np.save(f'{ROOT}/scratch/primal_{tag}_{setname}_N{N}_E{int(energy)}.npy', np.concatenate([[val, X, N], pv]))
