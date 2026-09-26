"""Discrete-grid coalescent likelihoods adapted from the supplied
real100_DIS_*.py source scripts (f = log Ne on Ngrid cells)."""
import numpy as np
from math import comb

def a_coeffs_kmax(kmax):
    if kmax in (1, 2):
        return np.array([1.0])
    a_prev = np.array([1.0])
    for k in range(3, kmax + 1):
        Mk, Mk_prev = comb(k - 1, 2), comb(k - 2, 2)
        a_cur = np.zeros(Mk + 1)
        for i in range(Mk + 1):
            if i == 0:
                a_cur[i] = 1.0
                continue
            a_im1_k = a_cur[i - 1]
            a_i_km1 = a_prev[i] if i <= Mk_prev else 0.0
            a_cur[i] = ((comb(k - 1, 2) - i + 1) * a_im1_k
                        + comb(k, 2) * a_i_km1) / (comb(k, 2) - i)
        a_prev = a_cur
    return a_prev

class DiscreteCoal:
    """Ngrid piecewise-constant cells on [0,T]; all ntip tips sampled at t=0."""
    def __init__(self, coal_times, ntip, T, Ngrid=100):
        self.coal = np.sort(np.asarray(coal_times, float))
        self.ntip, self.T, self.Ngrid = ntip, float(T), Ngrid
        self.grid = np.linspace(0.0, self.T, Ngrid + 1)
        self.loc = self.grid[1:]
        self.L = Ngrid
        self.a = a_coeffs_kmax(ntip)
        self.log_a = np.log(self.a)
        n_vec = np.arange(ntip, 1, -1)
        self.initC = n_vec * (n_vec - 1) / 2.0
        self.coal_idx = self._iv(self.coal)
        t_new = np.concatenate([[0.0], self.coal])
        self.interval_info = []
        for lo, hi in zip(t_new[:-1], t_new[1:]):
            interior = self.grid[(self.grid > lo) & (self.grid < hi)]
            pts = np.unique(np.concatenate([[lo], interior, [hi]]))
            left, right = pts[:-1], pts[1:]
            self.interval_info.append((self._iv((left + right) / 2.0), right - left))
        # exposure weights: sum_i initC_i * (time in cell j during interval i)
        w = np.zeros(self.L)
        for i, (gi, ln) in enumerate(self.interval_info):
            np.add.at(w, gi, self.initC[i] * ln)
        self.exposure = w

    def _iv(self, x):
        return np.clip(np.searchsorted(self.grid, x, side="left") - 1, 0, self.L - 1)

    def log_rhs(self, Lam):
        """log sum_i a_i x^i with x=exp(-Lam), computed in log space (all a_i>0)."""
        i = np.arange(len(self.a))
        return np.logaddexp.reduce(self.log_a - i * Lam)

    def loglik_sc(self, f):
        f = np.asarray(f, float)
        return -np.sum(f[self.coal_idx]) - np.sum(self.exposure * np.exp(-f))

    def log_bound_prob(self, f):
        Lam = np.sum(np.exp(-np.asarray(f, float))) * self.T / self.Ngrid
        if Lam <= 0:
            return -np.inf, Lam
        return (self.ntip - 1) * np.log(-np.expm1(-Lam)) + self.log_rhs(Lam), Lam

    def loglik_bc(self, f):
        lb, Lam = self.log_bound_prob(f)
        if not np.isfinite(lb):
            return -np.inf
        return self.loglik_sc(f) - lb
