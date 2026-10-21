# Certificate format tpmahler-upper-v1 (upper bound for Smyth's limit point of the Mahler measure of totally positive integers)

A certificate claims log lambda+ <= v (and log E+ <= v). lambda+ is the smallest limit point of
{Omega(alpha) = M(alpha)^{1/deg alpha} : alpha a totally positive algebraic integer}, M the Mahler measure, and E+ is the
essential minimum. The mathematics is Lemma U1 in `../PROGRESS.md` (section "UPPER BOUND"): Smith's theorem (arXiv:2111.12660,
Thm 1.5) in the checkable form of Talebizadeh Sardari-Orloski (arXiv:2304.10021). This file says how to rebuild every object from
the JSON.

## Fields (JSON object)
| field | type | meaning |
|---|---|---|
| `format` | string | must be `"tpmahler-upper-v1"` |
| `definition` | string | informal reminder; not used |
| `polys` | list of lists of JSON integers | Q_j in Z[x], coefficients in DESCENDING order (`[a_d, ..., a_0]`); each must be primitive and irreducible over Q; pairwise distinct |
| `b` | list of strings `"p/q"` | exact rationals b_j >= 0, one per polynomial |
| `nodes` | list of strings `"p/q"` | exact rationals 0 <= e_0 < e_1 < ... < e_N |
| `G` | list of strings `"p/q"` | exact rationals G_0..G_N >= 0 with G_0 = G_N = 0 |
| `X` | string `"p/q"` | rational X > 4 with e_N <= X (Sigma = [0, X], capacity X/4 > 1) |
| `v` | string `"p/q"` | the claimed bound for log lambda+ |
| `lp_value`, `mass_float` | numbers | informational only |

## Objects to rebuild
* G: [e_0, e_N] -> [0, inf) is the continuous piecewise-linear function with G(e_k) = G_k, extended by 0 outside [e_0, e_N].
* mu = G(x) dx, mass s = mu(R) = sum_k (e_{k+1} - e_k)(G_k + G_{k+1})/2 (exact rational).
* U(z) = int log|z - y| G(y) dy (the potential of mu, of mass s; "+log" convention), for z in C.
* h(x) = U(x) - sum_j b_j log|Q_j(x)| for real x (h = +inf at real roots of the Q_j with b_j > 0).

## What VERIFIED means
1. Data: as in the table (integer coefficients; primitive; irreducible over Q; distinct; b_j >= 0; nodes strictly increasing with
   e_0 >= 0; G_k >= 0, G_0 = G_N = 0; s > 0; X > 4 and e_N <= X).
2. (U1) sum_j b_j deg Q_j < s (strict, exact rational comparison).
3. (U2) int log|Q_j| dmu > 0 for every j with b_j > 0 (">= 0" suffices for the lemma).
4. (U3) h(x) >= 0 for every x in [e_0, e_N].
5. (V) (1/s) int log+ x dmu(x) <= v.

## Why this proves the claim
Lemma U1: (U3) extends to h >= 0 on C (h is superharmonic off [e_0, e_N] and tends to +inf at infinity by (U1)); then for every
irreducible P in Z[x] not among the weighted Q_j, int log|P| dmu = s log|lead P| + sum_rho U(rho) >= (s - sum b_j deg Q_j) log|lead P|
>= 0 (resultants are nonzero integers), and (U2) covers the weighted Q_j. So mu/s satisfies Smith's condition (1) on [0, X]; it is a
weak* limit of root distributions of totally positive algebraic integers alpha_k, and log Omega(alpha_k) -> (1/s) int log+ dmu <= v.
Hence E+ <= e^v, and lambda+ <= e^v (convexity of the set of arithmetic measures gives a continuum of limit values).

## Notes for an independent check
* G'' = sum_i C_i delta_{e_i} with exact slope jumps C_i (sum C_i = 0, sum C_i e_i = 0). Closed forms:
  U(x) = sum_i C_i H(e_i - x), H(t) = t^2/2 log|t| - 3 t^2/4 (real x); for complex z, U(z) = Re sum_i C_i Hc(e_i - z),
  Hc(t) = t^2/2 Log t - 3 t^2/4 (one branch for all i: the branch terms cancel because sum C_i e_i^2 is real).
  U''(x) = sum_i C_i log|e_i - x| (log singularities at the nodes; integrable).
* (U2): int log|Q| dmu = s log|a| + sum over the roots rho of Q of U(rho) (exact identity); some weighted Q_j have complex roots
  (pairs on the unit circle, from Flammang's construction).
* int log+ dmu: on each piece G(x) = alpha + beta x, int log x (alpha + beta x) dx with antiderivatives x log x - x and
  x^2/2 log x - x^2/4, over [max(e_k, 1), e_{k+1}].
* Only [e_0, e_N] needs checking (no tails): the minimum principle does the rest.

## Relation to the lower-bound format tpmahler-cert-v1 (for formalization)
Same encoding conventions: integer polynomials in descending order, exact rationals as "p/q" strings, a continuous piecewise-linear
density G on rational nodes with G_0 = G_N = 0, slope jumps C_i and the closed form U = sum_i C_i H(e_i - x).
Field correspondence: polys -> polys; c (weights) -> b; m (claim, lower) -> v (claim, upper); x_lo/x_hi (tails) -> X (capacity;
no tails are needed, the check runs on [e_0, e_N] only). The inequality is reversed: lower bound F = log+ - m - sum c log|Q| - 2U +
I/s >= 0 on (0, inf); upper bound h = U - sum b log|Q| >= 0 on [e_0, e_N], plus (U1), (U2), (V).
The interval-wide positivity argument (root windows with Sturm-matched brackets, chord bound with the Green's function, exact
endpoint values, closed-form knot singularities, bisection) is the same as in check_tp.py and is spelled out in the docstring of
check_up.py.
