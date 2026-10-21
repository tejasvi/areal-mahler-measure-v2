# Negative (and two positive) controls for a tpmahler-upper-v1 certificate: perturbed copies; expected verdicts in the tag
# (MUSTFAIL / OPEN / POSITIVE). Usage: make_neg_up.py CERT.json OUTDIR
import json, copy, sys
from fractions import Fraction as Fr
src, outdir = sys.argv[1], sys.argv[2]
base = json.load(open(src)); v = Fr(base['v'])
def save(tag, c): json.dump(c, open(f"{outdir}/neg_{tag}.json", 'w'))
w = sorted(((Fr(x), i) for i, x in enumerate(base['b'])), reverse=True)
h1, h2 = w[0][1], w[1][1]
for tag, dv in [('MUSTFAIL_v_minus_2e-9', -Fr(2, 10**9)), ('MUSTFAIL_v_minus_1e-6', -Fr(1, 10**6)), ('POSITIVE_v_plus_1e-6', Fr(1, 10**6))]:
    c = copy.deepcopy(base); c['v'] = str(v + dv); save(tag, c)
for tag, fac in [('OPEN_b_heaviest_x1.01', Fr(101, 100)), ('OPEN_b_heaviest_x0.99', Fr(99, 100))]:
    c = copy.deepcopy(base); c['b'][h1] = str(Fr(base['b'][h1])*fac); save(tag, c)
for tag, i in [('MUSTFAIL_drop_heaviest', h1), ('MUSTFAIL_drop_second_heaviest', h2)]:
    c = copy.deepcopy(base); c['b'][i] = '0'; save(tag, c)
for tag, fac in [('MUSTFAIL_G_x0.99', Fr(99, 100)), ('OPEN_G_x1.01', Fr(101, 100)), ('MUSTFAIL_G_x0.5', Fr(1, 2))]:
    c = copy.deepcopy(base); c['G'] = [str(Fr(g)*fac) for g in base['G']]; save(tag, c)
c = copy.deepcopy(base); c['polys'].append([1] + [0]*99 + [-3]); c['b'].append('1/100'); save('MUSTFAIL_budget_violation', c)
c = copy.deepcopy(base); q = base['polys'][h1]; c['polys'][h1] = [a - b for a, b in zip(q + [0], [0] + q)]   # Q*(x - 1): reducible
save('MUSTFAIL_reducible_weighted', c)
c = copy.deepcopy(base); c['polys'][h1] = [2*a for a in base['polys'][h1]]; save('MUSTFAIL_nonprimitive', c)
c = copy.deepcopy(base); k = len(c['nodes'])//2; c['nodes'][k], c['nodes'][k + 1] = c['nodes'][k + 1], c['nodes'][k]; save('MUSTFAIL_nodes_swapped', c)
c = copy.deepcopy(base); c['G'][-1] = '1/1000000'; save('MUSTFAIL_G_N_nonzero', c)
c = copy.deepcopy(base); c['X'] = '4'; save('MUSTFAIL_X_capacity_4', c)
i2 = next((i for i, q in enumerate(base['polys']) if q == [1, -3, 1]), None)
if i2 is not None:
    c = copy.deepcopy(base); c['polys'][i2] = [1, -3, -1]; save('OPEN_coeff_x2m3xp1_to_x2m3xm1', c)
print('written to', outdir)
