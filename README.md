# Bounds for the smallest limit point of the areal Mahler measure of totally positive algebraic integers

Certificates and programs for the paper
*Bounds for the smallest limit point of the areal Mahler measure of totally positive algebraic integers*
by Tejasvi Singh Tomar (Zenodo, [doi:10.5281/zenodo.23537179](https://doi.org/10.5281/zenodo.23537179)). This repository
is archived in the supplement record [doi:10.5281/zenodo.23537180](https://doi.org/10.5281/zenodo.23537180).

## Result

For an algebraic number `α` of degree `n` whose minimal polynomial `P` in `Z[x]` is primitive with leading coefficient
`a ≥ 1`, let `N₀(α) = ‖P‖₀^(1/n)`, where `‖P‖₀ = exp((1/π) ∫_D log|P| dA)` is Pritsker's areal analog of the Mahler measure
(`D` the unit disc). For totally positive `α` with conjugates `α_i`,

    n log N₀(α) = log a + Σ_i g(α_i),     g(x) = log⁺x + (min(1, x)² − 1)/2.

Let `λ_A⁺` be the smallest limit point of the values of `N₀` on totally positive algebraic integers and `E_A⁺` their
essential minimum. Flammang (J. Number Theory 151 (2015), Theorem 2) proved `N₀(α) ≥ 1.380047` with six explicit
exceptions. The paper proves

    1.3830134 < E_A⁺ ≤ λ_A⁺ < 1.3842933,

that all but finitely many totally positive algebraic numbers satisfy `N₀(α) > 1.3830134`, and that infinitely many totally
positive algebraic integers satisfy `N₀(α) < 1.3842933`.

- **Lower bound** (Lemma 2.1 of the paper, the method of auxiliary functions with a logarithmic energy term).
  `certificate/lower/cert_arllo_L3_5295faaa.json` gives 34 integer polynomials `Q_j` with weights `c_j > 0`, a continuous
  piecewise linear function `G` on 4,660 rational nodes in `[0, 24]`, with mass `s`, and a rational number `m`, such that
  `Σ c_j deg Q_j + 2s < 1` and

      F(x) = g(x) − m − Σ c_j log|Q_j(x)| − 2U(x) + I/s ≥ 0   for all x > 0 that are not roots of the Q_j,

  where `U(x) = ∫ log|x − y| G(y) dy` and `I = ∫ U G`. Then `log E_A⁺ ≥ m = 45636208945290293391/140737488355328000000`, and
  `e^m = 1.38301343…`.
- **Upper bound** (Lemma 3.1 of the paper, based on Smith's criterion, Ann. of Math. 200 (2024), Theorem 1.5 and
  Proposition 5.15). `certificate/upper/cert_arlup_U2_cb2c7960.json` gives a measure `μ = G(x) dx` on `[31/200, 16]` with
  mass `s` and 23 irreducible primitive polynomials with weights `b_j > 0`, such that (C1) `Σ b_j deg Q_j < s`, (C2)
  `∫ log|Q_j| dμ ≥ 0` for every `j`, (C3) `U(x) − Σ b_j log|Q_j(x)| ≥ 0` on `[31/200, 16]`, and
  `(1/s) ∫ g dμ ≤ v = 325189755077/10^12`. Then `λ_A⁺ ≤ e^v = 1.38429329…`.

All polynomials in both certificates are polynomials of Flammang's Table 2, and the certificates themselves list their
coefficients, so the checks do not depend on a transcription of the table.

## Files

| path | content |
|---|---|
| `certificate/lower/`, `certificate/upper/` | the two certificates (formats `tpmahler-areal-lower-v1` and `tpmahler-areal-upper-v1`) |
| `certificate/formats/` | the format `tpmahler-cert-v1` of the lower-bound companion paper (doi:10.5281/zenodo.23524856) and the format `tpmahler-upper-v1` of the upper-bound companion paper (doi:10.5281/zenodo.23529776); the areal formats have the same fields (see Notes) |
| `checker/lower/check_arl_lo_59a1c54c.py`, `checker/upper/check_arl_up_a17af685.py` | the checkers, in ball arithmetic (python-flint, 200 bits); their docstrings explain the method |
| `checker/*/diff_*.txt` | their differences from the companion papers' checkers, from which they were derived |
| `logs/` | checker runs on the certificates, certificate summaries by row number of Flammang's table, `COMMANDS.txt`, and a run-time diagnostic for the lower checker (`logs/lower/diag_L3.py` and its log) |
| `controls/lower/`, `controls/upper/` | modified certificates with the checkers' logs, `summary.txt`, the programs that made them, and `controls/REGEN.txt` |
| `search/` | the floating-point search that produced the certificates (not needed for checking; see Search) |
| `table2/` | two transcriptions of Flammang's Table 2 and the reproduction of her constant |
| `smyth/` | the values of `N₀(β_n²)` along Smyth's sequence |
| `audit/` | a second checker, `check_arl_second.py`, written independently of the checkers in `checker/`, with its runs: `audit/L3/` for the lower certificate and the rest of `audit/` for the upper one (see the Notes) |

## Checking the certificates

Requirements: Python 3 (we used 3.12.13) with python-flint (0.9.0). Both checkers refuse to run under `python -O`. The last
argument is the number of processes.

    python3 -I checker/lower/check_arl_lo_59a1c54c.py certificate/lower/cert_arllo_L3_5295faaa.json 24
    python3 -I checker/upper/check_arl_up_a17af685.py certificate/upper/cert_arlup_U2_cb2c7960.json 24

Each checker prints its intermediate results and ends with `VERIFIED` and exit status 0 only if every condition holds.
The lower checker checks the data and `Σ c_j deg Q_j + 2s < 1` exactly, bounds `F` from below on `(0, 1/1000]` and on
`[4357, ∞)`, and covers `[1/1000, 4357]` by windows around the roots and intervals between them, on which it uses a chord
bound with a rigorous bound for `∫ (F'')₊`. The upper checker checks (C1) exactly, (C2) from the closed form of `U` at
the roots of the `Q_j`, (C3) on the whole interval `[31/200, 16]`, and the bound on `(1/s) ∫ g dμ` from closed forms. The
logged runs took 3,665 s for the lower checker (48 processes, on a server with 96 threads) and 55 s for the upper checker
(40 processes). About 43,000 of the 69,298 interval checks of the lower checker fall on the single interval `[24, 4357]`,
far from the support of `G`, and one process does all of them, so more processes do not make the lower check much faster
(see `logs/COMMANDS.txt`). Runs on other servers printed the same lines apart from the times: 3,231 s for the lower check
with 24 processes (`audit/L3/logs/first_checker/check_first_L3.log`) and 28 s for the upper check with 24 processes.

To test the checkers, `controls/` holds modified certificates, made by the programs in the same directories;
`controls/REGEN.txt` records that all control files regenerate byte for byte. The checkers reject all of them except two,
both of which still state true bounds: `neg_m_plus_2e-6.json` raises `m` by `2·10^-6`, less than the margin of `3·10^-6`
that the search left, and `neg_POSITIVE_v_plus_1e-6.json` raises the upper bound `v` by `10^-6`. Two of the modified
certificates change the density by a narrow tent of mass about `10^-7` and width `4·10^-6`: `neg_NARROWNOTCH.json`
removes it at `x = 6129/25000`, inside `(0, 1)` where `g'' = 1`, and `neg_NARROWDIP.json` adds it at `x = 12.50711`,
between two points of the grid on which the search checked (C3). Although the inequality still holds at the nodes of the
tent and at nearby points, the checkers reject both certificates at the tent.

The second checker in `audit/` belongs to an audit of the checkers, certificates, logs and control files. It was written
from the certificate formats and the lemmas of the paper, without reference to the checkers in `checker/`. On
`[1/1000, 4357]` and on `[31/200, 16]` it bounds `F` and the function `U − Σ b_j log|Q_j|` of (C3) from below by
quadratic minorants from Taylor's formula with a rigorous bound for the second derivative, and it handles both formats:

    python3 -I audit/check_arl_second.py certificate/lower/cert_arllo_L3_5295faaa.json 16
    python3 -I audit/check_arl_second.py certificate/upper/cert_arlup_U2_cb2c7960.json 16

It ends with `VERIFIED` on both certificates and rejects the same modified certificates as the checkers in `checker/`.
Four further modified certificates, two in `audit/L3/controls/` and two in `audit/controls/upper_hidden/`, are narrow
notches and bumps at points where the inequalities are nearly tight, and all the checkers reject them. `audit/L3/RUNS.md`
lists the runs for the lower certificate and `audit/RUNS.md` those for the upper one.

## Search

The programs in `search/` are not needed to check the certificates. For the lower bound, `arl_nodes.py` chooses the nodes
from the solution of a preliminary run with the nodes `k/100` of `[0, 20]`, `arl_primal.py` minimizes `∫ g dθ` over
measures with piecewise linear densities on these nodes subject to `∫ log|Q| dθ ≥ 0` for the polynomials of Flammang's
Table 2 and to a nonnegative logarithmic energy, `arl_dual.py` solves the linear program for the auxiliary function with
this measure as the reference measure, and `export_arl_lo.py` writes the exact certificate. For the upper bound,
`arl_upper.py` solves the linear program for the measure and the weights by cutting planes, and `export_arl_up.py` writes
the exact certificate. `search/COMMANDS.txt` gives the command lines, and the solutions and logs of the runs are in
`search/lp/`. The export commands at the end of `search/COMMANDS.txt` reproduce both certificates, and the preliminary
run's lower certificate, byte for byte from the stored solutions (`search/REPRO_EXPORT.txt`). The search also
needs numpy, scipy, cvxpy and clarabel.

## Flammang's table

`table2/` holds two transcriptions of Table 2 of Flammang's paper: one from the HAL version's text layer
(`flammang2015_table2_pdftotext_79c01c8c.json`, made by `parse_table2_67d9bc59.py`) and one from images of its pages
(`flammang2015_table2_image_transcription_8720412f.json`). They agree except for row 58, to which the text layer appends
fragments of the page footer; the image transcription has the printed row. `repro_table2_19cdc09e.py` evaluates her
auxiliary function with the printed weights and re-optimizes the weights. With the image transcription, the printed weights
give `N₀ ≥ 1.3800359` and the re-optimized weights give `N₀ ≥ 1.3800500`. Its outputs for both transcriptions are the two
`.out` files. Most of the search used the text transcription (`search/data/A58.json`); the last linear program and the
export of the lower certificate used the image transcription (`search/data/A58v.json`). Row 58 has weight 0 in both
certificates and every polynomial in the certificates is a row of the table, so the wrong row does not affect them.

## Notes on the shipped files

- Some of the wording in the audited files predates the paper: Lemma 2.1 of the paper is called "Lemma TA1" there and
  Lemma 3.1 "Lemma UA1"; the lemmas of the companion papers from which they are adapted are called "Lemma TP1" (for the
  lower bound) and "Lemma U1" (for the upper bound). The conditions (C1)–(C3) are "(U1)–(U3)", and the bound on
  `(1/s) ∫ g dμ` is "(V)". (The docstrings also refer to `PROGRESS.md`, a file of working notes that is not part of this
  repository.)
- The formats `tpmahler-areal-lower-v1` and `tpmahler-areal-upper-v1` have exactly the fields of the formats in
  `certificate/formats/`. Only the left side of `F` changes, from `log⁺x` to `g(x)`, and the claim of the upper certificate
  bounds `(1/s) ∫ g dμ` instead of `(1/s) ∫ log⁺x dμ`.
- The runs in `audit/` outside `audit/L3/` that concern the lower bound (`audit/controls/lower_hidden/` and the files of
  `audit/logs/` with `L1` in their names) are for the preliminary run's lower certificate, `cert_arllo_L1_4c5758a0.json`,
  which gives the smaller bound `1.3823468`. It is archived at doi:10.5281/zenodo.23532710 and is not included here; the
  first command at the end of `search/COMMANDS.txt` reproduces it.
- The audited files are shipped unchanged. The shipped copies of the search programs, of `controls/*/make_dip_arl_*.py`
  and of `logs/lower/diag_L3.py` use paths relative to the program's directory, and `tp_primal_hat.py`, which
  `tp_upper_lp.py` imports, is copied from the release of the upper-bound companion paper (doi:10.5281/zenodo.23529778).

## Licence

Apache-2.0 for the programs (see `LICENSE`); the certificates, data and logs are released under CC BY 4.0 with the paper.
