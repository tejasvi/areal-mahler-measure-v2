# Certificate format `tpmahler-cert-v1`

The file `tp_certificate.json` is a JSON object. It proves log ℓ₊ ≥ m, where ℓ₊ is the smallest limit point of the values
Ω(α) = M(α)^{1/deg α} for totally positive algebraic integers α, and M(α) is the Mahler measure (leading coefficient
included). The mathematics is in Sections 2 and 3 of the paper. This file says how to rebuild every object from the JSON.

## Fields

| field | type | meaning |
|---|---|---|
| `format` | string | `"tpmahler-cert-v1"` |
| `definition` | string | a reminder of the definition of F; not used by the checkers |
| `polys` | list of lists of integers | the polynomials Q_j, coefficients from the leading one down to the constant term; deg Q_j = length − 1; must contain `[1, 0]` (the polynomial x) |
| `c` | list of strings `"p/q"` | the weights c_j ≥ 0, one per polynomial, in the same order; the weight of x must be positive |
| `nodes` | list of strings `"p/q"` | the nodes e_0 < e_1 < … < e_N, with e_0 ≥ 0 |
| `G` | list of strings `"p/q"` | the values G_0, …, G_N ≥ 0 at the nodes, with G_0 = G_N = 0 |
| `m` | string `"p/q"` | the bound: the claim is log ℓ₊ ≥ m |
| `x_lo`, `x_hi` | strings `"p/q"` | cut points for the two tails, used by `check_certificate.py` only; they do not affect the claim |
| `m_float_dual`, `margin` | numbers | information only: the value of the linear program and the margin subtracted from it |

## Objects

* G is the continuous function on ℝ that is linear on each [e_k, e_{k+1}], has G(e_k) = G_k, and vanishes outside
  [e_0, e_N].
* ρ = G(x) dx, and s = ρ(ℝ) = Σ_k (e_{k+1} − e_k)(G_k + G_{k+1})/2.
* U(x) = ∫ log|x − y| G(y) dy is the potential of ρ, and I = ∫∫ log|x − y| G(x) G(y) dx dy is its logarithmic energy.
* F(x) = log⁺ x − m − Σ_j c_j log|Q_j(x)| − 2 U(x) + I/s for x > 0, where log⁺ x = max(0, log x).

## What the checkers verify

1. The data: integer coefficients, nonzero leading coefficients, distinct polynomials; c_j ≥ 0, with a positive weight
   for x; strictly increasing nodes with e_0 ≥ 0; G_k ≥ 0 and G_0 = G_N = 0; s > 0.
2. Σ_j c_j deg Q_j + 2s < 1, as an exact comparison of rationals.
3. F(x) ≥ 0 for every x > 0, where F = +∞ at the positive roots of the Q_j with c_j > 0.

## Why this proves the claim

Put ν = ρ/s and λ = s. Then λ(2U_ν − I(ν)) = 2U − I/s, so F is the function of Lemma 2.2 of the paper, and ν satisfies
ν(D(y, r)) ≤ K r with K = 2 max G/s. Lemma 2.2 and Corollary 2.3 give: for every ε > 0, all but finitely many totally
positive algebraic numbers α satisfy log Ω(α) ≥ m − ε. Hence ℓ₊ ≥ e^m, and, through β ↦ β², all but finitely many totally
real algebraic numbers β have Weil height h(β) ≥ m/2 − ε.

## Formulas for a check

* G'' = Σ_i C_i δ_{e_i}, where C_i is the jump of the slope of G at e_i; Σ_i C_i = Σ_i C_i e_i = 0.
* U(x) = Σ_i C_i H(e_i − x) with H(t) = t²/2 · log|t| − 3t²/4, and I = Σ_{i,k} C_i C_k L(e_i − e_k) with
  L(t) = t⁴/24 · log|t| − 25t⁴/288 (Lemma 3.1 of the paper). Here H'' = log|t| and L'' = H.
* log⁺ x has a corner at x = 1, which is a root of Q_1 = (x − 1)², so F = +∞ there.
* Tails as in `check_certificate.py`. On (0, x_lo], with x_lo at most the left end of the support of G,
  F ≥ −m + I/s − c_x log x_lo − Σ_{Q_j ≠ x} c_j log max_{[0, x_lo]} |Q_j| − 2 U(0), because U decreases to the left of the
  support. On [x_hi, ∞), with x_hi ≥ max(1, e_N), F ≥ κ log x − b with κ = 1 − Σ_j c_j deg Q_j − 2s and
  b = m − I/s + Σ_j c_j log ‖Q_j‖₁, where ‖Q‖₁ is the sum of the absolute values of the coefficients of Q.

## Reference values for `tp_certificate.json`

* SHA-256 `9f1e51157fbe30ce8ab436c8552d7d6a9b5aaaa5483449de00ec68716eac3864`.
* m = 19152077587359311629/35184372088832000000 = 0.5443347841…, so e^m = 1.72346152….
* 19 polynomials: x and eighteen polynomials of Table 1 of V. Flammang, J. Number Theory 151 (2015), rebuilt as
  x^d U_j(x + 1/x − 2) for j = 1, …, 8, 10, 11, 13, 14, 16, 18, 20, 21, 24, 25 (Table 1 of the paper).
* 5473 nodes in [0, 56/5]; G vanishes outside [0.08525, 11.2] and on 401 node intervals inside it; s = 0.0528854151…;
  Σ_j c_j deg Q_j + 2s = 0.8067717903… < 1; κ = 0.1932282096…; max G/s = 2.4358996620… (computed by `checker/check_aux.py`).
* I = 2.8781306519508358577886530862…·10⁻¹¹.
