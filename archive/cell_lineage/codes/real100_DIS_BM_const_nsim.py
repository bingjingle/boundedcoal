#!/usr/bin/env python
# Locate the archived inputs independently of the working directory.
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _paths import data_file, output_file, OUTPUT_DIR


import argparse
import math
import os
import time
from math import comb

import numpy as np
import pyreadr

np.random.seed(123)


def parse_args():
    p = argparse.ArgumentParser(
        description="UPGMA100 bounded-coalescent discretized inference with Brownian-motion prior."
    )
    p.add_argument("--const", type=float, required=True,
                   help="BM intercept-scale multiplier: sigma0=sqrt(alpha/beta)*const.")
    p.add_argument("--nsim1", type=int, required=True,
                   help="Burn-in iterations.")
    p.add_argument("--nsim2", type=int, required=True,
                   help="Post-burn-in iterations.")
    p.add_argument(
        "--output-dir",
        default=str(OUTPUT_DIR)
    )
    return p.parse_args()


def num_tag(x):
    s = f"{x:.12g}"
    s = s.replace("e-0", "e-").replace("e+0", "e+")
    return s.replace(".", "p").replace("+", "")


args = parse_args()
const = args.const
nsim1 = args.nsim1
nsim2 = args.nsim2

if const <= 0:
    raise ValueError("const must be positive.")
if nsim1 <= 0 or nsim2 <= 0:
    raise ValueError("nsim1 and nsim2 must be positive.")

# Fixed hyperparameters from the supplied DIS-BM code.
alpha = 0.1
beta = 0.1
rem = 1

# UPGMA100 real data.
input_file = data_file("UGPMA100new.rda")
data = pyreadr.read_r(input_file)
points_inhomo = np.asarray(data["bt_adj"]).squeeze().astype(float)

T = 25.0
ntip = 100
Ngrid = 100

if points_inhomo.size != ntip - 1:
    raise ValueError(
        f"Expected {ntip - 1} UPGMA coalescent times in bt_adj, found {points_inhomo.size}."
    )

grid = np.linspace(0.0, T, Ngrid + 1)
loc = grid[1:]
L = len(loc)

# R findInterval(..., left.open=TRUE, rightmost.closed=TRUE), converted to 0-based indices.
idx = np.searchsorted(grid, points_inhomo, side="left") - 1
idx = np.clip(idx, 0, L - 1)


def a_coeffs_kmax(kmax):
    if kmax in (1, 2):
        return np.array([1.0])

    a_prev = np.array([1.0])
    for k in range(3, kmax + 1):
        Mk = comb(k - 1, 2)
        Mk_prev = comb(k - 2, 2)
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


# choose(ntip,2), choose(ntip-1,2), ..., choose(2,2)
lineages = np.arange(ntip, 1, -1)
initC = lineages * (lineages - 1) / 2.0

# Precompute the exact piecewise-grid exposure used by the supplied R likelihood.
# This is algebraically the same as summing initC[i] * I[i] at each likelihood call.
weighted_exposure = np.zeros(L, dtype=float)
t_vec_new = np.concatenate([[0.0], points_inhomo])
for i in range(len(points_inhomo)):
    left_interval = t_vec_new[i]
    right_interval = t_vec_new[i + 1]
    if right_interval <= left_interval:
        continue

    # Add overlap length of (left_interval,right_interval] with each grid cell.
    first = max(0, np.searchsorted(grid, left_interval, side="right") - 1)
    last = min(L - 1, np.searchsorted(grid, right_interval, side="left") - 1)
    for j in range(first, last + 1):
        cell_left = grid[j]
        cell_right = grid[j + 1]
        overlap = max(0.0, min(right_interval, cell_right) - max(left_interval, cell_left))
        if overlap > 0:
            weighted_exposure[j] += initC[i] * overlap


def log_lik(f):
    f = np.asarray(f, dtype=float)
    sum_f = np.sum(f[idx])
    no_coal = np.sum(weighted_exposure * np.exp(-f))

    Lambda = np.sum(np.exp(-f)) * T / Ngrid
    x = np.exp(-Lambda)
    val = rhs_value(x)
    if val <= 0 or Lambda <= 0:
        return -np.inf

    log_one_minus_x = np.log(-np.expm1(-Lambda))
    logboundprob = (ntip - 1) * log_one_minus_x + np.log(val)
    return -sum_f - no_coal - logboundprob


def chol_sample(mean, cov_chol):
    return mean + cov_chol @ np.random.standard_normal(mean.size)


def elliptical_slice(initial_theta, prior, lnpdf,
                     cur_lnpdf=None, angle_range=None):
    initial_theta = np.asarray(initial_theta, dtype=float)
    D = len(initial_theta)

    if cur_lnpdf is None:
        cur_lnpdf = lnpdf(initial_theta)

    if len(prior.shape) == 1:
        nu = prior
    else:
        if prior.shape[0] != D or prior.shape[1] != D:
            raise IOError("Prior must be a D-element sample or DxD chol(Sigma)")
        nu = prior @ np.random.normal(size=D)

    hh = math.log(np.random.uniform()) + cur_lnpdf

    if angle_range is None or angle_range == 0.0:
        phi = np.random.uniform() * 2.0 * math.pi
        phi_min = phi - 2.0 * math.pi
        phi_max = phi
    else:
        phi_min = -angle_range * np.random.uniform()
        phi_max = phi_min + angle_range
        phi = np.random.uniform() * (phi_max - phi_min) + phi_min

    while True:
        xx_prop = initial_theta * math.cos(phi) + nu * math.sin(phi)
        proposed_lnpdf = lnpdf(xx_prop)
        if proposed_lnpdf > hh:
            return xx_prop, proposed_lnpdf

        if phi > 0:
            phi_max = phi
        elif phi < 0:
            phi_min = phi
        else:
            raise RuntimeError("BUG DETECTED: shrunk to current position.")
        phi = np.random.uniform() * (phi_max - phi_min) + phi_min


# Brownian covariance on the 100 discretization locations.
cov_K = np.minimum.outer(loc, loc)
ones_vec = np.ones((L, 1))
sigma0 = np.sqrt(alpha / beta) * const
cov_K_corrected = cov_K + (ones_vec @ ones_vec.T) * sigma0 ** 2
cov_K_mod_chol = np.linalg.cholesky(cov_K_corrected)
cov_K_mod_inv = np.linalg.inv(cov_K_corrected)

print("const =", const, "alpha =", alpha, "beta =", beta, "sigma0 =", sigma0)

g_mk = np.full(L, 3.0)
theta = alpha / beta
n_total = nsim1 + nsim2

# The supplied DIS-BM implementation retains every iteration (rem=1).
g_mk_list = np.empty((n_total, L), dtype=np.float64)
theta_list = np.empty(n_total, dtype=np.float64)
curloglike = None

start_time = time.time()
for ite in range(n_total):
    if ite % 10000 == 0:
        print(f"iteration {ite}/{n_total}", flush=True)

    cov_K_chol_final = cov_K_mod_chol / np.sqrt(theta)
    prior = chol_sample(np.zeros(L), cov_K_chol_final)

    g_mk, curloglike = elliptical_slice(
        g_mk,
        prior,
        log_lik,
        cur_lnpdf=curloglike,
        angle_range=None,
    )

    g_mk_list[ite] = g_mk

    alpha_pos = alpha + L / 2
    beta_pos = beta + 0.5 * g_mk @ cov_K_mod_inv @ g_mk
    theta = float(np.random.gamma(alpha_pos, 1.0 / beta_pos))
    theta_list[ite] = theta

runtime = time.time() - start_time
timerun_10000 = runtime / n_total * 10000

# g_mk is log N_e for the discretized BM method, exactly as in the supplied R code.
g_post = np.exp(g_mk_list[nsim1:])
med = np.quantile(g_post, 0.5, axis=0)


def hpd_interval(x, alpha=0.05):
    x = np.sort(np.asarray(x, dtype=float))
    n = len(x)
    k = int(np.floor((1 - alpha) * n))
    widths = x[k:] - x[:n - k]
    j = np.argmin(widths)
    return x[j], x[j + k]


low = np.zeros(L)
high = np.zeros(L)
for j in range(L):
    low[j], high[j] = hpd_interval(g_post[:, j], alpha=0.05)

output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(
    output_dir,
    f"real100_BM_DIS_const_{num_tag(const)}_nsim1_{nsim1}_nsim2_{nsim2}.npz",
)

np.savez(
    output_file,
    g_mk_list=g_mk_list,
    theta_list=theta_list,
    c=points_inhomo,
    d=loc,
    p=timerun_10000,
    med=med,
    low=low,
    high=high,
    const=const,
    sigma0=sigma0,
    nsim1=nsim1,
    nsim2=nsim2,
    rem=rem,
    alpha=alpha,
    beta=beta,
    T=T,
)

print("Saved:", output_file)
print("const =", const, "sigma0 =", sigma0)
print("nsim1 =", nsim1, "nsim2 =", nsim2)
print("Runtime per 10,000 iterations =", timerun_10000)
