# Context computation (not part of any proof): the absolute Pritsker measure along Smyth's sequence.
# beta_0 = 1, beta_n - 1/beta_n = beta_{n-1}; the 2^n conjugates of beta_n are the branch values (y +- sqrt(y^2 + 4))/2 over the
# conjugates y of beta_{n-1}; beta_n^2 is totally positive with conjugates c^2, so log N0(beta_n^2) = mean of g(c^2),
# g(x) = log+ x + (min(1, x)^2 - 1)/2.  Also prints Omega(beta_n) (its square is Smyth's ell^2 for the Mahler measure).
import numpy as np
g = lambda x: np.where(x >= 1, np.log(np.maximum(x, 1.0)), (np.minimum(x, 1.0)**2 - 1)/2)
c = np.array([1.0])
for n in range(1, 23):
    c = np.concatenate([(c + np.sqrt(c*c + 4))/2, (c - np.sqrt(c*c + 4))/2])
    print(f"n = {n:2d}: N0(beta_n^2) = {np.exp(np.mean(g(c*c))):.13f}   Omega(beta_n) = {np.exp(np.mean(np.log(np.maximum(1.0, np.abs(c))))):.13f}", flush=True)
