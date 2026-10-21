# Common closed forms for step densities (uniform on cells) and integer polynomials with real roots.
import numpy as np
def G(t):   # antiderivative of log|t|
    t = np.asarray(t, dtype=float); out = np.zeros_like(t); m = t != 0
    out[m] = t[m]*np.log(np.abs(t[m])) - t[m]; return out
def H(t):   # antiderivative of G
    t = np.asarray(t, dtype=float); out = np.zeros_like(t); m = t != 0
    out[m] = 0.5*t[m]**2*np.log(np.abs(t[m])) - 0.75*t[m]**2; return out
def cells(xmin, xmax, n):
    e = np.linspace(xmin, xmax, n+1); return e[:-1], e[1:]
def energy_matrix(a, b):
    # K_kl = (1/(w_k w_l)) int_{a_k}^{b_k} int_{a_l}^{b_l} log|x-y| dy dx
    A = a[:, None]; B = b[:, None]; C = a[None, :]; D = b[None, :]
    w = (b - a)
    num = H(D - A) - H(D - B) - H(C - A) + H(C - B)
    return num/(w[:, None]*w[None, :])
def potential(a, b, x):
    # U_k(x) = (1/w_k) int_{a_k}^{b_k} log|x-y| dy, shape (len(x), n)
    x = np.asarray(x, dtype=float)[:, None]; w = (b - a)[None, :]
    return (G(b[None, :] - x) - G(a[None, :] - x))/w
def poly_roots(coeffs):
    # coeffs highest-first integer list; returns lead, real roots (assert real)
    r = np.roots(coeffs); assert np.max(np.abs(r.imag)) < 1e-9, coeffs
    return coeffs[0], np.sort(r.real)
def logQ_cell_avg(coeffs, a, b):
    lead, roots = poly_roots(coeffs); w = b - a
    out = np.full_like(a, np.log(abs(lead)), dtype=float)
    for rho in roots: out += (G(b - rho) - G(a - rho))/w
    return out
def logQ(coeffs, x):
    lead, roots = poly_roots(coeffs); x = np.asarray(x, dtype=float)
    return np.log(abs(lead)) + sum(np.log(np.abs(x - rho)) for rho in roots)
A7 = [[1,0],[1,-1],[1,-2],[1,-3,1],[1,-5,6,-1],[1,-7,13,-7,1],[1,-7,14,-8,1]]
A11_PIN = [[1,0],[-1,1],[1,-3,1],[1,-5,6,-1],[1,-7,14,-8,1],[1,-7,13,-7,1],[1,-9,28,-35,15,-1],
           [1,-9,26,-29,11,-1],[1,-9,27,-32,13,-1],[1,-6,9,-1],[1,-9,27,-31,12,-1]]
# ---- piecewise-linear (hat) densities on uniform nodes x_k = x0 + k w, k = 1..N-1 (hats normalized to mass 1) ----
def J4(t):  # J'' = H, J(t) = t^4/24 log|t| - 25/288 t^4
    t = np.asarray(t, dtype=float); out = np.zeros_like(t); m = t != 0
    out[m] = t[m]**4/24.0*np.log(np.abs(t[m])) - 25.0/288.0*t[m]**4; return out
def hat_nodes(xmin, xmax, N):
    e = np.linspace(xmin, xmax, N+1); return e, e[1]-e[0]   # hats centred at e[1..N-1]
def hat_potential(e, w, x):
    # U_k(x) for hats centred at e[k], k=1..N-1: (1/w^2)[H(e_{k+1}-x) - 2H(e_k - x) + H(e_{k-1}-x)]
    x = np.asarray(x, dtype=float)[:, None]; c = e[1:-1][None, :]
    return (H(c + w - x) - 2*H(c - x) + H(c - w - x))/w**2
def hat_energy_matrix(e, w):
    k = np.arange(1, len(e)-1); D = (k[:, None] - k[None, :]).astype(float)
    out = np.zeros_like(D)
    for m, t in zip(range(-2, 3), [1, -4, 6, -4, 1]): out += t*J4((D + m)*w)
    return out/w**4
def hat_logQ(coeffs, e, w):
    lead, roots = poly_roots(coeffs)
    out = np.full(len(e)-2, np.log(abs(lead)))
    for rho in roots: out += hat_potential(e, w, np.array([rho]))[0]
    return out
# ---- non-uniform hats: nodes e[0..N]; hat k (k=1..N-1) on [e[k-1], e[k+1]], mass 1 ----
def hat_coefs(e):
    # second-derivative coefficients c (3 per hat) at positions e[k-1], e[k], e[k+1]
    wl = e[1:-1] - e[:-2]; wr = e[2:] - e[1:-1]; s = 2.0/(wl + wr)
    return np.stack([s/wl, -s*(1/wl + 1/wr), s/wr], axis=1)        # (N-1, 3)
def hatnu_potential(e, x):
    # U_k(x) = sum_i c_ki H(pos_ki - x)
    c = hat_coefs(e); x = np.asarray(x, dtype=float)[:, None]
    return c[None, :, 0]*H(e[None, :-2] - x) + c[None, :, 1]*H(e[None, 1:-1] - x) + c[None, :, 2]*H(e[None, 2:] - x)
def hatnu_energy_matrix(e):
    c = hat_coefs(e); pos = np.stack([e[:-2], e[1:-1], e[2:]], axis=1)   # (n,3)
    n = len(c); K = np.zeros((n, n))
    for i in range(3):
        for j in range(3):
            K += np.outer(c[:, i], c[:, j]) * J4(pos[:, i][:, None] - pos[:, j][None, :])
    return K
def hatnu_logQ(coeffs, e):
    lead, roots = poly_roots(coeffs)
    out = np.full(len(e)-2, np.log(abs(lead)))
    for rho in roots: out += hatnu_potential(e, np.array([rho]))[0]
    return out
# ---- numerically stable hat matrices: closed form for near pairs, Gauss-Legendre for far pairs ----
_GLX, _GLW = np.polynomial.legendre.leggauss(8)
def _hat_quad_nodes(e):
    # quadrature nodes/weights for each hat k: two segments [e_{k-1}, e_k], [e_k, e_{k+1}] with linear density
    l, m, r = e[:-2], e[1:-1], e[2:]; peak = 2.0/(r - l)
    X = []; W = []
    for a, b, rising in [(l, m, True), (m, r, False)]:
        xs = 0.5*(a + b)[:, None] + 0.5*(b - a)[:, None]*_GLX[None, :]
        dens = peak[:, None]*((xs - a[:, None])/(b - a)[:, None] if rising else (b[:, None] - xs)/(b - a)[:, None])
        X.append(xs); W.append(0.5*(b - a)[:, None]*_GLW[None, :]*dens)
    return np.hstack(X), np.hstack(W)        # (n, 16)
def stable_hat_potential(e, x, far=4.0):
    x = np.asarray(x, dtype=float)
    U = hatnu_potential(e, x)                                  # closed form (fine for near points)
    c = e[1:-1]; wmax = np.maximum(e[1:-1] - e[:-2], e[2:] - e[1:-1])
    Xq, Wq = _hat_quad_nodes(e)
    D = np.abs(x[:, None] - c[None, :]); farm = D > far*wmax[None, :]
    rows, cols = np.nonzero(farm)
    # quadrature: U_k(x) = sum_q W_kq log|x - X_kq|
    for k in np.unique(cols):
        rk = rows[cols == k]
        U[rk, k] = np.log(np.abs(x[rk][:, None] - Xq[k][None, :])) @ Wq[k]
    return U
def stable_hat_energy_matrix(e, far=4.0):
    K = hatnu_energy_matrix(e)
    c = e[1:-1]; wmax = np.maximum(e[1:-1] - e[:-2], e[2:] - e[1:-1]); n = len(c)
    Xq, Wq = _hat_quad_nodes(e)
    D = np.abs(c[:, None] - c[None, :]); farm = D > far*np.maximum(wmax[:, None], wmax[None, :])
    for k in range(n):
        ls = np.nonzero(farm[k])[0]
        if len(ls) == 0: continue
        # K_kl = sum_q sum_r W_kq W_lr log|X_kq - X_lr|
        L = np.log(np.abs(Xq[k][None, :, None] - Xq[ls][:, None, :]))     # (len(ls), 16, 16)
        K[k, ls] = np.einsum('q,lqr,lr->l', Wq[k], L, Wq[ls])
    return 0.5*(K + K.T)
def ld_hat_potential(e, x, chunk=20000):
    # closed-form hat potentials evaluated in extended precision (np.longdouble) to tame cancellation
    eL = np.asarray(e, dtype=np.longdouble); c = hat_coefs(eL)
    out = np.empty((len(x), len(e)-2))
    def HL(t):
        o = np.zeros_like(t); m = t != 0; o[m] = t[m]**2/2*np.log(np.abs(t[m])) - np.longdouble(3)/4*t[m]**2; return o
    for i0 in range(0, len(x), chunk):
        xx = np.asarray(x[i0:i0+chunk], dtype=np.longdouble)[:, None]
        U = c[None, :, 0]*HL(eL[None, :-2] - xx) + c[None, :, 1]*HL(eL[None, 1:-1] - xx) + c[None, :, 2]*HL(eL[None, 2:] - xx)
        out[i0:i0+chunk] = U.astype(float)
    return out
