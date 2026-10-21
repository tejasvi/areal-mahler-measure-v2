# Areal analog, lower-bound side: hat primal SOCP  min int g dmu  s.t.  int log|Q_j| dmu >= 0, I(mu) >= 0 (energy), mu prob. on [0, X].
# Usage: arl_primal.py SET NODES(uniform:X:h | file.npy) TAG
import os, numpy as np, cvxpy as cp, sys, time, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tpc import stable_hat_energy_matrix
from tp_upper_lp import logQ_hat_flint
from tpsets import load
from arl_common import g_hat
setname, nodespec, tag = sys.argv[1], sys.argv[2], sys.argv[3]
polys = load(setname)
if nodespec.endswith('.npy'): e = np.load(nodespec)
else:
    _, X, h = nodespec.split(':'); e = np.round(np.arange(0.0, float(X) + 1e-12, float(h)), 10)
t0 = time.time(); n = len(e) - 2
f = g_hat(e); L = np.array([logQ_hat_flint(q, e) for q in polys])
for energy in (False, True):
    p = cp.Variable(n, nonneg=True); cons = [cp.sum(p) == 1, L @ p >= 0]
    if energy:
        K = stable_hat_energy_matrix(e); C = max(1.0, 2*np.abs(K).max())
        P = C*np.ones((n, n)) - K; ev, V = np.linalg.eigh(0.5*(P + P.T)); R = (V*np.sqrt(np.maximum(ev, 0))).T
        cons.append(cp.sum_squares(R @ p) <= C)
    prob = cp.Problem(cp.Minimize(f @ p), cons); prob.solve(solver='CLARABEL')
    pv = np.maximum(p.value, 0); sup = e[1:-1][pv > 1e-7*pv.max()]
    print(f"{setname} {nodespec} energy={energy}: primal {prob.value:.7f} (N0 {np.exp(prob.value):.7f}) [{prob.status}] support [{sup.min():.3f}, {sup.max():.3f}] ({time.time()-t0:.0f}s)", flush=True)
    if energy:
        print("   duals c:", [(j, round(float(cc), 7)) for j, cc in enumerate(cons[1].dual_value) if cc > 1e-8], " lambda_E:", cons[2].dual_value, flush=True)
        np.save(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'primal_{tag}.npy'), np.concatenate([[prob.value, len(e)], e, pv]))
