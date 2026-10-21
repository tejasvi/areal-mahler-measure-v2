# Extra negative controls (distinct from make_neg.py; for this certificate the heaviest non-x polynomial IS (x-1)^2, so
# make_neg.py's 'drop_heaviest' and 'drop_x_minus_1' coincide).
import json, copy, sys
from fractions import Fraction as Fr
src, outdir = sys.argv[1], sys.argv[2]
base = json.load(open(src))
def save(tag, c): json.dump(c, open(f"{outdir}/neg_{tag}.json", 'w'))
w = sorted(((Fr(x), i) for i, x in enumerate(base['c']) if base['polys'][i] != [1, 0]), reverse=True)
c = copy.deepcopy(base); c['c'][w[1][1]] = '0'; save('drop_second_heaviest', c)
i3 = base['polys'].index([1, -3, 1])
c = copy.deepcopy(base); c['polys'][i3] = [1, -3, 2]; save('coeff_x2m3xp1_to_x2m3xp2', c)
first = min(Fr(e) for e, g in zip(base['nodes'], base['G']) if Fr(g) > 0)
c = copy.deepcopy(base); c['x_lo'] = str(first + Fr(1, 100)); save('x_lo_inside_support', c)
c = copy.deepcopy(base); c['G'][-1] = '1/1000000'; save('G_N_nonzero', c)
c = copy.deepcopy(base); k = len(c['nodes'])//2; c['nodes'][k], c['nodes'][k + 1] = c['nodes'][k + 1], c['nodes'][k]; save('nodes_swapped', c)
print('written to', outdir)
