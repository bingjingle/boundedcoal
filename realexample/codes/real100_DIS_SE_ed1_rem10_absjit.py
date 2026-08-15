#!/usr/bin/env python

import argparse
import math
import time
from pathlib import Path

import numpy as np
import pyreadr
from scipy.linalg import cho_solve
from scipy.special import comb
from scipy.stats import norm


# -----------------------------------------------------------------------------
# Arguments
# -----------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Bounded-coalescent discretized SE inference for the UPGMA100 real dataset."
    )
    p.add_argument("--proposal-sigma-l", type=float, required=True)
    p.add_argument("--nsim1", type=int, required=True,
                   help="Burn-in iterations.")
    p.add_argument("--nsim2", type=int, required=True,
                   help="Post-burn-in iterations.")
    p.add_argument("--jitter", type=float, default=1e-14)
    p.add_argument("--rem", type=int, default=10)
    p.add_argument(
        "--input-file",
        default="/home/groups/juliapr/Bingjing/code/boundedcoal/syndata1/UGPMA100new.rda",
    )
    p.add_argument(
        "--output-dir",
        default="/scratch/groups/juliapr/output_Bingjing/boundedcoaldata/real/",
    )
    return p.parse_args()


args = parse_args()

proposal_sigma_l = args.proposal_sigma_l
nsim1 = args.nsim1
nsim2 = args.nsim2
jitter = args.jitter
rem = args.rem

if proposal_sigma_l <= 0:
    raise ValueError("proposal_sigma_l must be positive")
if nsim1 <= 0 or nsim2 <= 0:
    raise ValueError("nsim1 and nsim2 must be positive")
if jitter <= 0:
    raise ValueError("jitter must be positive")
if rem <= 0:
    raise ValueError("rem must be positive")

np.random.seed(123)

# -----------------------------------------------------------------------------
# UPGMA100 real data
# -----------------------------------------------------------------------------
data = pyreadr.read_r(args.input_file)
points_inhomo = np.asarray(data["bt_adj"], dtype=float).squeeze()

T = 25.0
ntip = 100
Ngrid = 100

grid = np.linspace(0.0, T, Ngrid + 1)
loc = grid[1:]
L = len(loc)

t_vec = points_inhomo
t_vec_new = np.concatenate([[0.0], t_vec])

alpha = 0.1
beta = 0.1

# Same lineage-count coefficients as the uploaded discretized R code.
n_vec = np.arange(ntip, 1, -1)
initC = n_vec * (n_vec - 1) / 2.0


# -----------------------------------------------------------------------------
# Bounded-coalescent normalizing probability
# -----------------------------------------------------------------------------
def a_coeffs_kmax(kmax):
    if kmax == 1 or kmax == 2:
        return np.array([1.0])

    a_prev = np.array([1.0])
    for k in range(3, kmax + 1):
        Mk = int(comb(k - 1, 2))
        Mk_prev = int(comb(k - 2, 2))
        a_cur = np.zeros(Mk + 1, dtype=float)

        for i in range(Mk + 1):
            if i == 0:
                a_cur[i] = 1.0
                continue

            a_im1_k = a_cur[i - 1]
            a_i_km1 = a_prev[i] if i <= Mk_prev else 0.0
            numerator = (
                (comb(k - 1, 2) - i + 1) * a_im1_k
                + comb(k, 2) * a_i_km1
            )
            denominator = comb(k, 2) - i
            a_cur[i] = numerator / denominator

        a_prev = a_cur

    return a_prev


a = a_coeffs_kmax(ntip)


def rhs_value(x):
    poly = 0.0
    for coeff in reversed(a):
        poly = poly * x + coeff
    return poly


# R equivalent: findInterval(x, grid, left.open=TRUE, rightmost.closed=TRUE)
def interval_index(x):
    idx = np.searchsorted(grid, x, side="left") - 1
    return np.clip(idx, 0, L - 1)


coal_idx = interval_index(t_vec)

# Precompute the grid pieces inside each coalescent interval.
interval_info = []
for left_end, right_end in zip(t_vec_new[:-1], t_vec_new[1:]):
    interior = grid[(grid > left_end) & (grid < right_end)]
    points = np.unique(np.concatenate([[left_end], interior, [right_end]]))
    left = points[:-1]
    right = points[1:]
    mid = (left + right) / 2.0
    grid_idx = interval_index(mid)
    interval_info.append((grid_idx, right - left))


def log_lik(f):
    sum_f = np.sum(f[coal_idx])

    I = np.empty(len(interval_info), dtype=float)
    for i, (idx_i, len_i) in enumerate(interval_info):
        I[i] = np.sum(np.exp(-f[idx_i]) * len_i)

    Lambda = np.sum(np.exp(-f)) * T / Ngrid
    x = np.exp(-Lambda)
    val = rhs_value(x)

    if (not np.isfinite(val)) or val <= 0:
        return -np.inf

    logboundprob = (ntip - 1) * np.log(-np.expm1(-Lambda)) + np.log(val)
    return -sum_f - np.sum(initC * I) - logboundprob


# -----------------------------------------------------------------------------
# Sampling + GP covariance, translated from the uploaded DIS R implementation
# -----------------------------------------------------------------------------
def elliptical_slice(initial_theta, prior, lnpdf,
                     cur_lnpdf=None, angle_range=None):
    if cur_lnpdf is None:
        cur_lnpdf = lnpdf(initial_theta)

    nu = prior
    hh = np.log(np.random.uniform()) + cur_lnpdf

    if angle_range is None or angle_range == 0:
        phi = np.random.uniform() * 2.0 * np.pi
        phi_min = phi - 2.0 * np.pi
        phi_max = phi
    else:
        phi_min = -angle_range * np.random.uniform()
        phi_max = phi_min + angle_range
        phi = np.random.uniform() * (phi_max - phi_min) + phi_min

    while True:
        xx_prop = initial_theta * np.cos(phi) + nu * np.sin(phi)
        cur_lnpdf = lnpdf(xx_prop)

        if cur_lnpdf > hh:
            return xx_prop, cur_lnpdf

        if phi > 0:
            phi_max = phi
        elif phi < 0:
            phi_min = phi
        else:
            raise RuntimeError("BUG DETECTED")

        phi = np.random.uniform() * (phi_max - phi_min) + phi_min


def build_cov(l_value, jitter_start):
    theta1_value = np.exp(-2.0 * l_value)
    distance = loc[:, None] - loc[None, :]
    cov_K = np.exp(-theta1_value / 2.0 * distance ** 2)

    jitter_try = jitter_start
    while jitter_try <= 1.0:
        try:
            cov_K_chol = np.linalg.cholesky(
                cov_K + jitter_try * np.eye(L)
            )
            logdet_K = 2.0 * np.sum(np.log(np.diag(cov_K_chol)))
            return cov_K_chol, logdet_K
        except np.linalg.LinAlgError:
            jitter_try *= 10.0

    raise np.linalg.LinAlgError("Cholesky failed")


def quad_from_chol(f, cov_K_chol):
    return float(f @ cho_solve((cov_K_chol, True), f))


def log_posterior_l(l_value, f, theta, mu_l, scale_l, jitter_start):
    cov_K_chol, logdet_K = build_cov(l_value, jitter_start)
    quad = quad_from_chol(f, cov_K_chol)
    partial = norm.logpdf(l_value, loc=mu_l, scale=scale_l) - 0.5 * logdet_K
    return partial - 0.5 * theta * quad, partial, cov_K_chol, logdet_K


def metropolis_hastings_update_l(
        l_current, proposal_sigma_l, mu_l, scale_l, f, theta,
        cov_K_chol_current, logdet_current, jitter_start,
        partial_current=None):

    if partial_current is None:
        partial_current = (
            norm.logpdf(l_current, loc=mu_l, scale=scale_l)
            - 0.5 * logdet_current
        )

    quad_current = quad_from_chol(f, cov_K_chol_current)
    log_post_current = partial_current - 0.5 * theta * quad_current

    l_prop = l_current + np.random.normal(0.0, proposal_sigma_l)
    proposed = log_posterior_l(
        l_prop, f, theta, mu_l, scale_l, jitter_start
    )
    log_post_prop, partial_prop, chol_prop, logdet_prop = proposed

    if np.log(np.random.uniform()) < log_post_prop - log_post_current:
        return l_prop, partial_prop, chol_prop, logdet_prop, True

    return (
        l_current, partial_current,
        cov_K_chol_current, logdet_current, False,
    )


# ED1, matching the uploaded R method.
theta1 = 10.0
l = -0.5 * np.log(theta1)
mu_l = l
scale_l = 1e4

g_mk = np.full(L, 3.0)
theta = alpha / beta
cov_K_chol, logdet_K = build_cov(l, jitter)
partialloglik = None
curloglike = None
acc_count = 0

# Full latent chain, as in the uploaded R implementation.
g_mk_list = np.empty((nsim1 + nsim2, L), dtype=float)
theta_list = np.empty(nsim1 + nsim2, dtype=float)
theta1_list = np.empty(nsim1 + nsim2, dtype=float)
l_list = np.empty(nsim1 + nsim2, dtype=float)

start_time = time.time()

for ite in range(nsim1 + nsim2):
    if ite % 1000 == 0:
        print(ite, flush=True)

    prior = (cov_K_chol / np.sqrt(theta)) @ np.random.standard_normal(L)
    g_mk, curloglike = elliptical_slice(
        g_mk, prior, log_lik, curloglike
    )
    g_mk_list[ite, :] = g_mk

    if ite % rem == 0:
        alpha_pos = alpha + L / 2.0
        beta_pos = beta + 0.5 * quad_from_chol(g_mk, cov_K_chol)
        theta = float(np.random.gamma(alpha_pos, 1.0 / beta_pos))

        (
            l, partialloglik, cov_K_chol, logdet_K, acc
        ) = metropolis_hastings_update_l(
            l, proposal_sigma_l, mu_l, scale_l,
            g_mk, theta, cov_K_chol, logdet_K,
            jitter, partialloglik,
        )

        theta1 = np.exp(-2.0 * l)
        if acc:
            acc_count += 1

    theta_list[ite] = theta
    theta1_list[ite] = theta1
    l_list[ite] = l


timerun = time.time() - start_time
timerun_10000 = timerun / (nsim1 + nsim2) * 10000.0
n_hyper_updates = sum(
    1 for ite in range(nsim1 + nsim2) if ite % rem == 0
)
acc_rate = acc_count / n_hyper_updates

# -----------------------------------------------------------------------------
# Posterior summaries on N_e(t) = exp(g(t)); no synthetic truth for real data.
# -----------------------------------------------------------------------------
g_post = np.exp(g_mk_list[nsim1:, :])
med = np.quantile(g_post, 0.5, axis=0)


def hpd_interval(x, alpha=0.05):
    x = np.sort(np.asarray(x, dtype=float))
    n = len(x)
    k = int(np.floor((1.0 - alpha) * n))
    widths = x[k:] - x[:n-k]
    j = np.argmin(widths)
    return x[j], x[j + k]


low = np.zeros(L)
high = np.zeros(L)
for j in range(L):
    low[j], high[j] = hpd_interval(g_post[:, j], alpha=0.05)


def num_tag(x):
    s = f"{x:.12g}"
    s = s.replace("e-0", "e-").replace("e+0", "e+")
    return s


output_dir = Path(args.output_dir)
output_dir.mkdir(parents=True, exist_ok=True)

step_tag = num_tag(proposal_sigma_l)
output_file = output_dir / (
    f"real100_SE_DIS_stepsize_{step_tag}"
    f"_nsim1_{nsim1}_nsim2_{nsim2}.npz"
)

np.savez(
    output_file,
    g_mk_list=g_mk_list,
    theta_list=theta_list,
    theta1_list=theta1_list,
    l_list=l_list,
    c=points_inhomo,
    loc=loc,
    med=med,
    low=low,
    high=high,
    acc_rate=acc_rate,
    timerun_10000=timerun_10000,
    proposal_sigma_l=proposal_sigma_l,
    nsim1=nsim1,
    nsim2=nsim2,
    jitter=jitter,
    rem=rem,
    T=T,
    ntip=ntip,
)

print("Saved:", output_file)
print("MH acceptance rate:", acc_rate)
print("Runtime per 10,000 iterations:", timerun_10000)
