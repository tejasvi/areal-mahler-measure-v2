# Negative controls for a tpmahler certificate: perturbed copies that must NOT verify.
import json, copy, sys
from fractions import Fraction as Fr
src, outdir = sys.argv[1], sys.argv[2]
base = json.load(open(src)); m = Fr(base['m'])
def save(tag, c): json.dump(c, open(f"{outdir}/neg_{tag}.json", 'w'))
for tag, dm in [('m_plus_1e-5', Fr(1, 10**5)), ('m_plus_4e-6', Fr(4, 10**6)), ('m_plus_2e-6', Fr(2, 10**6))]:
    c = copy.deepcopy(base); c['m'] = str(m + dm); save(tag, c)
# drop the heaviest non-x polynomial and the polynomial (x-1)^2 / x-1
w = [(Fr(x), i) for i, x in enumerate(base['c'])]
heavy = max((wi for wi in w if base['polys'][wi[1]] != [1, 0]))[1]
c = copy.deepcopy(base); c['c'][heavy] = '0'; save('drop_heaviest', c)
for i, q in enumerate(base['polys']):
    if q in ([1, -2, 1], [1, -1]):
        c = copy.deepcopy(base); c['c'][i] = '0'; save('drop_x_minus_1', c)
c = copy.deepcopy(base); c['G'] = [str(Fr(g)*Fr(101, 100)) for g in base['G']]; save('G_times_1.01', c)
c = copy.deepcopy(base); c['G'] = [str(Fr(g)*Fr(99, 100)) for g in base['G']]; save('G_times_0.99', c)
ix = [i for i, q in enumerate(base['polys']) if q == [1, 0]][0]
c = copy.deepcopy(base); c['c'][ix] = str(Fr(base['c'][ix])*Fr(1001, 1000)); save('cx_times_1.001', c)
c = copy.deepcopy(base); c['polys'].append([1] + [0]*99 + [-3]); c['c'].append('1/200'); save('budget_violation', c)
c = copy.deepcopy(base); k = next(i for i, g in enumerate(base['G']) if Fr(g) > 0); c['G'][k] = str(-Fr(base['G'][k])); save('negative_G', c)
print('written to', outdir)
