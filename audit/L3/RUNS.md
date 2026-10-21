# Second checker on the sharper lower certificate (L3): runs and logs

This directory adds the runs of the second checker on the sharper lower certificate cert_arllo_L3_5295faaa.json. The
programs are byte-identical copies of those in the parent directory, and ../RUNS.md describes them:
check_arl_second.py, make_hidden_notch.py, validate_arl.py and scan_contacts.py (same SHA-256 as there). The certificate
format (tpmahler-areal-lower-v1) and the first checker check_arl_lo_59a1c54c.py are unchanged.

Claim checked: log λ_A⁺ ≥ m = 45636208945290293391/140737488355328000000, so λ_A⁺ ≥ e^m = 1.38301343570895… ≥ 1.3830134.

Requirements: Python 3.12 or later with python-flint 0.9.0 (FLINT 3.6.0); no other third-party package.

## Inputs (SHA-256)

| item | file | SHA-256 |
| --- | --- | --- |
| lower certificate L3 (λ_A⁺ ≥ 1.3830134) | cert_arllo_L3_5295faaa.json | 5295faaacff1340587ef77ae593b69432cecaf793ed7c5f599ae38bc48cac5e8 |
| first checker, lower format | check_arl_lo_59a1c54c.py | 59a1c54c6735e0f9880c823c1feb4b9d851b3cc1c8ed5a57096e1136b43c8cfd |

The 16 controls of the L3 certificate release were checked unchanged. Their SHA-256 values are in the SHA256SUMS file of
that release.

## Results

All logs are in `logs/` and were produced by the commands in `logs/COMMANDS.txt`, with paths relative to a working
directory that contains the L3 certificate and its controls.

| run | result |
| --- | --- |
| lower certificate L3 (`logs/check_second_L3.log`) | VERIFIED. Budget Σ c_j deg Q_j + 2s = 0.903938467019 < 1; I = 1.09650468663·10⁻¹⁰; left tail (0, 1/1000]: F ≥ 1.53243; right tail [4357, ∞): F ≥ 0.192145; main range [1/1000, 4357]: 39,833 interval tests, 230 root brackets, exact coverage, no failure, smallest certified lower bound 3.4862·10⁻⁸. e^m = 1.38301343570895 ± 3.5·10⁻¹⁵ |
| first checker on L3 (`logs/first_checker/check_first_L3.log`) | VERIFIED. Apart from the timings, the log equals the one in the certificate release: 230 root windows, 4,790 chord intervals, 69,298 interval checks, smallest endpoint lower bound 2.9632·10⁻⁶ |
| 16 controls of the certificate release (`logs/controls_L3/`) | 14 are NOT VERIFIED for their intended reasons: the data checks (G_N ≠ 0, a negative value of G, swapped nodes), the budget (1.403938 > 1), the left tail (x_lo inside the support of G), or a certainly negative value of F (G × 0.99 and × 1.01, a changed coefficient, c_x × 1.001, dropped polynomials, m + 4·10⁻⁶, m + 10⁻⁵). m + 2·10⁻⁶ is VERIFIED (the margin of m is 3·10⁻⁶). The narrow notch is NOT VERIFIED with F = −7.9155·10⁻⁹ at x = 0.245160375 |
| hidden notches (`controls/lower_hidden/`, `logs/controls_L3/own_*.log`) | both NOT VERIFIED at the notch centre (F = −7.4669·10⁻⁸ at x = 0.68820533, −8.5172·10⁻⁸ at x = 3.4833503), while the budget and both tails pass. The first checker also rejects both, each with one failed interval inside the notch and the same values, −7.467·10⁻⁸ and −8.517·10⁻⁸ (`logs/first_checker/`) |
| generators (`logs/make_hidden_notch_*.log`, `logs/regen_cmp.txt`) | regenerate both hidden notches byte for byte. At the check points around each notch, F is at least 9.5122·10⁻⁹ (notch A, 62 points) and 1.7981·10⁻⁸ (notch B, 27 points) |
| minorant validation (`logs/validate_arl_L3.log`) | 1,628 points in 407 accepted sub-intervals: the full closed form is never below the minorant; the smallest gap is 3.1·10⁻¹⁵ |
| contact scan (`logs/scan_contacts_L3.log`) | the smallest values of F at node-interval midpoints and finer points, from which the sites of the hidden notches were chosen |
