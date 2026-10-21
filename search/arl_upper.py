# Upper-bound LP for lambda+ (Smith criterion; see tp_upper_lp.py) with cutting planes and parallel potentials.
# v3: the cutting-plane check evaluates U_mu by the closed form sum_i C_i H(e_i - x) (fast, vectorised) on a FINE grid over the
#     support (spacing FINEH) plus 1e-4 elsewhere; LP rows still use the per-hat stable potentials.
#   minimize int log+ dmu  over hat weights p >= 0 (sum p = 1) and b_Q >= 0, subject to
#   h(x) = U_mu(x) - sum_Q b_Q log|Q(x)| >= delta on the LP grid (grid refined by cutting planes from a fine check grid),
#   sum_Q b_Q deg Q <= 1 - 1e-6,  int log|Q| dmu >= eps_Q for every Q (Q must be irreducible: use an *irr set).
# Usage: tp_upper_lp3.py SET NODES(uniform:X:h | file.npy) GRIDH TAG [DELTA] [FINEH] [NPROC] [INIT full|sparse] [METHOD highs|highs-ds|highs-ipm]
import numpy as np, sys, json, time, os
import multiprocessing as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scipy.optimize import linprog
from tpc import stable_hat_potential
from arl_common import g_hat as logplus_hat   # areal objective g (name kept for a minimal diff)
from tp_dual_rootprod import flint_roots, logabs_rootprod
from tp_upper_lp import logQ_hat_flint
from tpsets import load
OUT = os.path.dirname(os.path.abspath(__file__))
_G = {}
def _urows(xc): return stable_hat_potential(_G['e'], xc)
def _uvec(xc): return stable_hat_potential(_G['e'], xc) @ _G['p']
def pool_map(fn, pieces, nproc):
    os.environ['OMP_NUM_THREADS'] = '1'
    with mp.get_context('fork').Pool(nproc) as pool:
        return pool.map(fn, pieces)
def Umat(e, x, nproc, chunk=400):
    _G['e'] = e
    return np.vstack(pool_map(_urows, [x[i:i + chunk] for i in range(0, len(x), chunk)], nproc))
def Uvec(e, p, x, nproc, chunk=2000):
    _G['e'], _G['p'] = e, p
    return np.concatenate(pool_map(_uvec, [x[i:i + chunk] for i in range(0, len(x), chunk)], nproc))
def Lmat(polys, x): return np.array([logabs_rootprod(q, x) for q in polys]).T
def U_closed(e, p, x, chunk=20000):
    """U_mu(x) for mu = sum_k p_k phi_k (unit-mass hats on nodes e): G_k = 2 p_k/(e_{k+1}-e_{k-1}), slope jumps C_i,
    U(x) = sum_i C_i H(e_i - x), H(t) = t^2/2 log|t| - 3 t^2/4."""
    G = np.zeros(len(e)); G[1:-1] = 2*p/(e[2:] - e[:-2])
    sl = np.diff(G)/np.diff(e); C = np.zeros(len(e)); C[0] = sl[0]; C[-1] = -sl[-1]; C[1:-1] = sl[1:] - sl[:-1]
    nz = np.nonzero(C)[0]; Cn, en = C[nz], e[nz]; out = np.empty(len(x))
    for i in range(0, len(x), chunk):
        t = en[None, :] - x[i:i + chunk, None]; at = np.abs(t)
        with np.errstate(divide='ignore', invalid='ignore'):
            H = np.where(at > 0, t*t/2*np.log(np.where(at > 0, at, 1.0)) - 0.75*t*t, 0.0)
        out[i:i + chunk] = H @ Cn
    return out
if __name__ == '__main__':
    setname, nodespec, gridh, tag = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4]
    delta = float(sys.argv[5]) if len(sys.argv) > 5 else 1e-6
    fineh = float(sys.argv[6]) if len(sys.argv) > 6 else 2e-5
    nproc = int(sys.argv[7]) if len(sys.argv) > 7 else 32
    init = sys.argv[8] if len(sys.argv) > 8 else 'full'
    method = sys.argv[9] if len(sys.argv) > 9 else 'highs'
    t0 = time.time(); polys = load(setname); J = len(polys); degs = np.array([len(q) - 1 for q in polys], float)
    if nodespec.endswith('.npy'): e = np.load(nodespec)
    else:
        _, X, h = nodespec.split(':'); e = np.round(np.arange(0.0, float(X) + 1e-12, float(h)), 10)
    n = len(e) - 2; f = logplus_hat(e); LQ = np.array([logQ_hat_flint(q, e) for q in polys])
    roots = np.concatenate([flint_roots(q) for q in polys]); rr = np.sort(roots[np.abs(roots.imag) < 1e-12].real)
    rr = rr[(rr > 0) & (rr < e[-1])]
    def clean(x):
        x = np.unique(x); x = x[(x > 0) & (x < e[-1])]
        return x[np.min(np.abs(x[:, None] - rr[None, :]), axis=1) > 1e-9] if len(rr) else x
    x = clean(np.concatenate([np.arange(gridh/2, e[-1], gridh), e[1:-1], 0.5*(e[1:] + e[:-1])]) if init == 'full'
              else np.arange(gridh/2, e[-1], gridh))
    U = Umat(e, x, nproc); L = Lmat(polys, x)
    xc = clean(np.concatenate([np.arange(1e-4/2, e[-1], 1e-4)] +
                              [np.concatenate([r + np.geomspace(1e-7, 0.003, 25), r - np.geomspace(1e-7, 0.003, 25)]) for r in rr]))
    Lc = Lmat(polys, xc)
    print(f"{setname}: J = {J}, {n} hats on [0, {e[-1]}], LP grid {len(x)}, check grid {len(xc)}, delta {delta} ({time.time()-t0:.0f}s)", flush=True)
    for it in range(40):
        c = np.concatenate([f, np.zeros(J)])
        Aub = np.vstack([np.hstack([-U, L]), np.concatenate([np.zeros(n), degs])[None, :], np.hstack([-LQ, np.zeros((J, J))])])
        bub = np.concatenate([-delta*np.ones(len(x)), [1 - 1e-6], -1e-7*np.ones(J)])
        res = linprog(c, A_ub=Aub, b_ub=bub, A_eq=np.concatenate([np.ones(n), np.zeros(J)])[None, :], b_eq=[1.0],
                      bounds=[(0, None)]*(n + J), method=method,
                      options={'primal_feasibility_tolerance': 1e-10, 'dual_feasibility_tolerance': 1e-10})
        assert res.status == 0, res.message
        p, b = res.x[:n], res.x[n:]; val = res.fun
        if it == 0:                                                    # fine check grid over the current support (+-0.02)
            sup0 = e[1:-1][p > 1e-12*p.max()]
            xf = clean(np.arange(max(sup0.min() - 0.02, fineh), min(sup0.max() + 0.02, e[-1] - fineh), fineh))
            xc = np.unique(np.concatenate([xc, xf])); Lc = Lmat(polys, xc)
            if it == 0:
                ref = Uvec(e, p, xc[::997], nproc); cl = U_closed(e, p, xc[::997])
                print(f'  closed-form U vs stable potentials: max |diff| {np.abs(ref - cl).max():.2e} on {len(cl)} points', flush=True)
        hc = U_closed(e, p, xc) - Lc @ b
        bad = np.nonzero(hc < 0.5*delta)[0]
        print(f"  it {it}: value {val:.8f} (Omega {np.exp(val):.7f}), rows {len(x)}, min h(check) {hc.min():.2e} at {xc[np.argmin(hc)]:.5f}, "
              f"{len(bad)} below delta/2 ({time.time()-t0:.0f}s)", flush=True)
        if len(bad) == 0: break
        lm = [i for i in bad if (i == 0 or hc[i] <= hc[i - 1]) and (i == len(hc) - 1 or hc[i] <= hc[i + 1])]
        if init == 'sparse' and len(bad) > 20*len(lm):          # broad violations: also add a thinned set of violated points
            lm = sorted(set(lm) | set(bad[::max(1, len(bad)//4000)].tolist()))
        xn = np.setdiff1d(xc[lm], x)
        if len(xn) == 0: xn = np.setdiff1d(xc[bad[:200]], x)
        U = np.vstack([U, Umat(e, xn, nproc)]); L = np.vstack([L, Lmat(polys, xn)]); x = np.concatenate([x, xn])
    sup = e[1:-1][p > 1e-9*p.max()]
    print(f"RESULT {tag}: value {val:.8f} Omega {np.exp(val):.7f}; support [{sup.min():.4f}, {sup.max():.4f}]; sum b deg {degs @ b:.6f}; "
          f"min int log|Q| dmu {(LQ @ p).min():.2e}; min h(check) {hc.min():.2e}", flush=True)
    print("   b:", [(j, float(v)) for j, v in enumerate(b) if v > 0], flush=True)
    np.save(f'{OUT}/arlup_{tag}.npy', np.concatenate([[val, len(e)], e, p, b]))
    json.dump({'set': setname, 'nodes': nodespec, 'delta': delta, 'fineh': fineh, 'J': J}, open(f'{OUT}/arlup_{tag}.json', 'w'))
