# Parse Flammang 2015 Table 2 (polynomials of her Theorem 2, Pritsker measure) from the pdftotext of her paper.
# Usage: parse_table2.py NORME_TXT OUT_JSON
import re, sys, json
from fractions import Fraction as Fr
src, out = sys.argv[1], sys.argv[2]
lines = open(src, encoding='utf-8').read().split('\n')
i0 = next(i for i, l in enumerate(lines) if 'Table 2:' in l)
entries, cur = [], None
for l in lines[i0 + 1:]:
    if l.strip().startswith('References') or l.strip().startswith('[A'): break
    m = re.match(r'^\s*(\d+)\s+(0\.\d+)\s+(.*)$', l)
    if m and (cur is None or int(m.group(1)) == cur['j'] + 1):
        if cur: entries.append(cur)
        cur = {'j': int(m.group(1)), 'c': m.group(2), 'txt': m.group(3).strip()}
    elif cur and re.match(r'^\s+[+−\-]\s*\d*x', l):
        cur['txt'] += ' ' + l.strip()
    elif cur and re.match(r'^\s+[+−\-]', l):
        cur['txt'] += ' ' + l.strip()
if cur: entries.append(cur)
def parse_poly(t):
    t = t.replace('−', '-').replace(' ', '')
    terms = re.findall(r'([+-]?)(\d*)(x(\d*))?', t)
    coef = {}
    for sgn, num, xs, ex in terms:
        if not num and not xs: continue
        c = int(num) if num else 1
        if sgn == '-': c = -c
        d = (int(ex) if ex else 1) if xs else 0
        coef[d] = coef.get(d, 0) + c
    deg = max(coef)
    return [coef.get(k, 0) for k in range(deg, -1, -1)]
res = []
for e in entries:
    q = parse_poly(e['txt'])
    res.append({'j': e['j'], 'c': e['c'], 'desc': q, 'txt': e['txt']})
json.dump(res, open(out, 'w'), indent=0)
print(len(res), 'entries; degrees:', [len(r['desc']) - 1 for r in res])
