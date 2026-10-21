# Areal (Pritsker) analog: objective g(x) = log+ x + (min(1, x)^2 - 1)/2 (log N0(alpha) = (1/n)(log a + sum g(alpha_i))).
import os, numpy as np, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tpc import _hat_quad_nodes
def g(x):
    x = np.asarray(x, float)
    return np.where(x >= 1, np.log(np.maximum(x, 1.0)), (np.minimum(x, 1.0)**2 - 1)/2)
def g_hat(e):
    Xq, Wq = _hat_quad_nodes(e)
    return (g(Xq)*Wq).sum(1)
