# Polynomial sets for the totally positive Mahler-measure problem (descending integer coefficients).
import os, json
ROOT = os.path.dirname(os.path.abspath(__file__))
def flammang27():
    D = json.load(open(f'{ROOT}/data/flammang2015_tp_Q.json'))
    return [[1, 0]] + [e['desc'] for e in D]          # x plus her 27 reciprocal Q_j(x) = x^d U_j(x + 1/x - 2)
def load(name):
    if name == 'F27': return flammang27()
    return json.load(open(f'{ROOT}/data/{name}.json'))
