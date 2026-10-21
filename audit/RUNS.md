# Second checker: programs, runs and logs

This directory holds a second checker for the two certificates of the areal bounds (formats tpmahler-areal-lower-v1 and
tpmahler-areal-upper-v1). It was written independently of the first checkers (check_arl_lo_59a1c54c.py and
check_arl_up_a17af685.py, shipped with the certificates), from the two certificate formats with log⁺ replaced by g and
from the two criteria of the paper, before the first checkers were read. It shares with them only the file formats. For
the upper format its method is that of the second checker of the companion paper on the totally positive case.

Here g(x) = log⁺ x + (min(1, x)² − 1)/2, that is (x² − 1)/2 on [0, 1] and log x on [1, ∞). In the programs and logs,
(U1), (U2), (U3) are the conditions of the upper-bound criterion and (V) is the claim (1/s) ∫ g dμ ≤ v; for the lower
format, "budget", "left tail", "right tail" and "main range" are the parts of the condition F ≥ 0 of the lower-bound
criterion.

## Programs

- `check_arl_second.py CERT NPROC`: the checker, for both formats. Exit status 0 and the last line `VERIFIED` only if
  every check passes.
  - Data checks: integer coefficients, nonzero leading coefficient, degree ≥ 1, no two polynomials equal up to sign,
    weights ≥ 0; nodes strictly increasing with e_0 ≥ 0; G_k ≥ 0 with G_0 = G_N = 0; s > 0. Upper format: every
    polynomial primitive and irreducible over Q (FLINT factorization), X > 4 and e_N ≤ X. Lower format: the polynomial x
    has positive weight. The slope jumps C_i of G are exact rationals, and the program checks Σ C_i = 0,
    Σ C_i e_i = 0, Σ C_i e_i² = 2s, which make the closed form U(x) = Σ C_i H(e_i − x), H(t) = t²/2 log|t| − 3t²/4,
    equal to ∫ log|x − y| G(y) dy.
  - Upper format: (U1) is exact. (U2) is s log|a| + Σ_ρ U(ρ) over the roots ρ as certified isolating balls, with the
    complex form of U at non-real roots. (V) is the exact rational integral of the linear pieces of G against
    (x² − 1)/2 on [e_k, min(e_{k+1}, 1)], plus ball arithmetic for the pieces above 1. (U3) h = U − Σ b_j log|Q_j| ≥ 0 is
    proved on every point of [e_0, e_N] as in the companion paper's second checker: a Taylor model of order 6 for the
    nodes more than 0.05 away (remainder from |H⁽⁷⁾(t)| = 24/|t|⁵), chords for the concave near node terms, quadratic
    minorants R(m) + R′(m)u + inf R″ · u²/2 for all other terms, bisection, real roots of the weighted polynomials cut
    out by brackets of half-width 10⁻⁴⁰ (overlapping brackets merged) with h ≥ U − Σ_other b log|Q| − b_j log(sup|Q_j′| ·
    width), a stop at the first certainly negative value, and an exact coverage check.
  - Lower format: the budget Σ c_j deg Q_j + 2s < 1 is exact. I = Σ C_i C_j J(e_i − e_j), J(t) = t⁴/24 log|t| −
    25t⁴/288 (so J″ = H and J⁗ = log|t|), is a ball, and its lower end is used. F(x) = g(x) − m − Σ c_j log|Q_j(x)| −
    2U(x) + I/s ≥ 0 is proved on (0, ∞):
    - left tail (0, x_lo], with x_lo at most the left end of the support of G: g ≥ −1/2; U(x) ≤ U(0), since U is
      decreasing there; the term of the polynomial x is at least −c_x log x_lo; the other terms are bounded by
      max |Q_j| on [0, x_lo] in ball arithmetic; and no weighted real root lies in (0, x_lo];
    - right tail [x_hi, ∞), with x_hi ≥ max(1, e_N) and no weighted real root above x_hi:
      F ≥ (1 − Σ c_j deg Q_j − 2s) log x − m − Σ c_j log Σ_k |q_jk| + I/s, which is increasing;
    - main range [x_lo, x_hi]: the same minorant machinery as for (U3), with the node terms multiplied by −2, so that the
      node terms with C_i < 0 are the concave ones bounded by chords; g is handled exactly on each side of 1 (g″ = 1 on
      (0, 1), −1/x² on (1, ∞)), and no interval crosses 1.
- `validate_arl.py CERT NSAMP SEED`: tests the minorant itself. The interval test is run on NSAMP random node intervals;
  on up to 1,200 accepted sub-intervals the function (F or h) is evaluated by the full closed form (all nodes, no Taylor
  model) at 4 random points each and must not be below the minorant.
- `scan_contacts.py CERT`: evaluates F or h at node-interval midpoints and finer points and lists the smallest values
  (contact points). It was used to choose the sites of the hidden controls.
- `make_hidden_notch.py CERT X0_NUM X0_DEN A_NUM A_DEN WT_NUM WT_DEN OUT [grid]` (lower format): removes a
  flat-bottomed notch (flat part of half-width a, ramps of width wt) from G at x0. m and everything else are unchanged,
  and the exact change of I/s is included. The notch depth is tuned by bisection so that F < 0 at x0, while F > 0 at
  every node within 0.01 of x0, at the 4 inserted nodes and, with `grid`, at the 11 points of the grid of spacing
  2·10⁻⁵ around x0.
- `make_hidden_dips.py CERT XNUM XDEN NODE OUTDIR` (upper format): adds a flat-top bump of density at two sites, dip A
  with centre XNUM/XDEN + 1/(3·10⁷) and dip B with its left foot 10⁻⁷ to the right of the node with index NODE. Each
  bump's mass is tuned by bisection so that h < 0 at its centre, while h > 0 at every node within 0.01, at the 4 inserted
  nodes and at the 11 grid points around the bump. v is raised by 10⁻⁶, so that only (U3) can fail.

Requirements: Python 3.12 or later with python-flint 0.9.0 (FLINT 3.6.0); no other third-party package.
`check_arl_second.py` uses multiprocessing with fork and refuses to run under `python -O`.

## Inputs (SHA-256)

| item | file | SHA-256 |
| --- | --- | --- |
| lower certificate (λ_A⁺ ≥ 1.3823468) | cert_arllo_L1_4c5758a0.json | 4c5758a01151acbd31d0193517137a230b152539e43ce1337c869eb1f5c150ef |
| upper certificate (λ_A⁺ ≤ 1.3842933) | cert_arlup_U2_cb2c7960.json | cb2c7960028f60e3bcce2c9f76f08cadfe00eaa131fda60b6b18ffe7591080a3 |
| first checker, lower format | check_arl_lo_59a1c54c.py | 59a1c54c6735e0f9880c823c1feb4b9d851b3cc1c8ed5a57096e1136b43c8cfd |
| first checker, upper format | check_arl_up_a17af685.py | a17af685ce81ecc17fa22ff9f0943ffe6326a4922c7531615e16cda3e2baf6e4 |

The 34 controls of the certificate release (16 for the lower certificate, 18 for the upper one, including the narrow
notch and the narrow dip) were checked unchanged. Their SHA-256 values are in the SHA256SUMS file of that release.

## Results

All logs are in `logs/` and were produced by the commands in `logs/COMMANDS.txt`, with paths relative to a working
directory that contains the certificates and controls of the release.

| run | result |
| --- | --- |
| lower certificate (`logs/check_second_L1.log`) | VERIFIED. Budget Σ c_j deg Q_j + 2s = 0.902471165815 < 1; I = 5.50777014626·10⁻¹¹; left tail (0, 1/1000]: F ≥ 1.54691; right tail [4897, ∞): F ≥ 0.195062; main range [1/1000, 4897]: 10,220 interval tests, 189 root brackets, exact coverage, no failure, smallest certified lower bound 7.474·10⁻⁸. e^m = 1.38234689517067 ± 1.5·10⁻¹⁵ |
| upper certificate (`logs/check_second_U2.log`) | VERIFIED. (U1) Σ b_j deg Q_j = 0.9395242426 < s; (U2) smallest ∫ log\|Q_j\| dμ = 1.0·10⁻⁷; (V) (1/s) ∫ g dμ = 0.325189754087824 ± 4.7·10⁻¹⁶, margin 9.892·10⁻¹⁰ to v; (U3) 9,839 interval tests, 138 root brackets, exact coverage, no failure, smallest certified lower bound 3.06·10⁻¹⁰. e^v = 1.38429329774151 ± 4.3·10⁻¹⁵ |
| 16 lower controls of the certificate release (`logs/controls_L1/`) | 14 are NOT VERIFIED for their intended reasons: the data checks (G_N ≠ 0, a negative value of G, swapped nodes), the budget, the left tail (x_lo inside the support of G), or a certainly negative value of F (G × 0.99 and × 1.01, a changed coefficient, c_x × 1.001, dropped polynomials, m + 4·10⁻⁶, m + 10⁻⁵). m + 2·10⁻⁶ is VERIFIED (the margin of m is 3·10⁻⁶). The narrow notch is NOT VERIFIED with F = −7.8915·10⁻⁹ at x = 0.232750375 |
| 18 upper controls of the certificate release (`logs/controls_U2/`) | all 12 MUSTFAIL and 4 OPEN controls are NOT VERIFIED (data checks, (U1), (V), or a certainly negative value of h). POSITIVE (v + 10⁻⁶) is VERIFIED. The narrow dip is NOT VERIFIED with h = −5.0667·10⁻⁹ at x = 12.507110375 |
| hidden notches (`controls/lower_hidden/`, `logs/controls_L1/own_*.log`) | both NOT VERIFIED at the notch centre (F = −2.0120·10⁻⁸ at x = 0.22251033, −4.6649·10⁻⁸ at x = 3.1625003), while the budget and both tails pass. The first checker rejects both at the same points with the same values (`logs/first_checker/`) |
| hidden dips (`controls/upper_hidden/`, `logs/controls_U2/own_*.log`) | both NOT VERIFIED at the bump centre (h = −3.1788·10⁻⁸ at x = 7.16875003, −2.1806·10⁻⁸ at x = 2.2900011), while (U1), (U2) and (V) pass. The first checker rejects both at the same points with the same values |
| generators (`logs/make_hidden_*.log`, `logs/regen_cmp.txt`) | regenerate the four controls byte for byte. At the check points around each control, F or h is at least 5.7655·10⁻⁸ (notch A, 17 points), 5.9301·10⁻⁸ (notch B, 6 points), 1.2967·10⁻⁸ (dip A) and 3.0508·10⁻⁸ (dip B) |
| minorant validation (`logs/validate_arl_L1.log`, `logs/validate_arl_U2.log`) | 1,536 points in 384 accepted sub-intervals (lower) and 2,488 points in 622 (upper): the full closed form is never below the minorant; the smallest gaps are 7.1·10⁻¹⁵ and 2.2·10⁻¹⁴ |
| contact scans (`logs/scan_contacts_L1.log`, `logs/scan_contacts_U2.log`) | the smallest values of F and h at node-interval midpoints and finer points, from which the sites of the hidden controls were chosen |
