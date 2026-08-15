#!/usr/bin/env python

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
        description="UPGMA100 bounded-coalescent RI inference with Brownian-motion prior."
    )
    p.add_argument("--const", type=float, required=True,
                   help="BM intercept-scale multiplier: sigma0=sqrt(alpha/beta)*const.")
    p.add_argument("--nsim1", type=int, required=True,
                   help="Burn-in iterations.")
    p.add_argument("--nsim2", type=int, required=True,
                   help="Post-burn-in iterations.")
    p.add_argument(
        "--output-dir",
        default="/scratch/groups/juliapr/output_Bingjing/boundedcoaldata/real/"
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

# Fixed hyperparameters
alpha = 0.1
beta = 0.1
rem = 10

# UPGMA100 real data
input_file = "/home/groups/juliapr/Bingjing/code/boundedcoal/syndata1/UGPMA100new.rda"
data = pyreadr.read_r(input_file)
points_inhomo = np.asarray(data["bt_adj"]).squeeze().astype(float)

T = 25.0
ntip = 100
N = ntip - 1
K = ntip
Ngrid = 100

if points_inhomo.size != N:
    raise ValueError(
        f"Expected {N} UPGMA coalescent times in bt_adj, found {points_inhomo.size}."
    )

x4 = np.linspace(0.0, T, Ngrid + 1)[1:]
xxx = np.concatenate([x4, points_inhomo])
days = np.concatenate([[0.0], points_inhomo, [T]])
diff = np.diff(days)

Nfinal = Ngrid + N + K
g_mk = np.full(Nfinal, 3.0)
g_mk[Ngrid + N:Ngrid + N + K] = 3.0 * T * diff / np.sum(diff)


# Brownian-motion covariance functions.
def expo_quad_kernel(xn, xm):
    return min(float(xn), float(xm))


def expo_quad_kernel2(xn, t1, t2):
    xn, t1, t2 = float(xn), float(t1), float(t2)
    h = t2 - t1
    if xn <= t1:
        return xn * h
    elif xn < t2:
        u = xn - t1
        return t1 * h + u * h - 0.5 * u ** 2
    else:
        return h * (t1 + 0.5 * h)


def expo_quad_kernel3(t1, t2, t3, t4):
    t1, t2, t3, t4 = map(float, (t1, t2, t3, t4))
    h1 = t2 - t1
    h2 = t4 - t3
    if t1 == t3 and t2 == t4:
        return t1 * h1 ** 2 + h1 ** 3 / 3.0
    # In the construction below i <= j, so the intervals are ordered.
    return h1 * h2 * (t1 + 0.5 * h1)


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


n_vec = np.arange(ntip, 1, -1)
com_vec = n_vec * (n_vec - 1) / 2


def chol_sample(mean, cov_chol):
    return mean + cov_chol @ np.random.standard_normal(mean.size)


def log_lik(f, ns):
    if np.any(f <= 0):
        return -np.inf

    term1 = np.sum(np.log(f[:, Ngrid:Ngrid + ns]))
    term2 = -np.sum(com_vec * f[:, Ngrid + ns:Nfinal - 1])

    Lambda = np.sum(f[:, Ngrid + ns:Nfinal])
    x = np.exp(-Lambda)
    val = rhs_value(x)
    if val <= 0 or Lambda <= 0:
        return -np.inf

    # Stable version of log(1-exp(-Lambda)).
    log_one_minus_x = np.log(-np.expm1(-Lambda))
    logboundprob = (ntip - 1) * log_one_minus_x + np.log(val)
    return term1 + term2 - logboundprob


def elliptical_slice(initial_theta, prior, lnpdf, pdf_params=(),
                     cur_lnpdf=None, angle_range=None):
    D = len(initial_theta)
    if cur_lnpdf is None:
        cur_lnpdf = lnpdf(initial_theta, *pdf_params)

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
        cur_lnpdf = lnpdf(xx_prop, *pdf_params)
        if cur_lnpdf > hh:
            break
        if phi > 0:
            phi_max = phi
        elif phi < 0:
            phi_min = phi
        else:
            raise RuntimeError("BUG DETECTED: shrunk to current position.")
        phi = np.random.uniform() * (phi_max - phi_min) + phi_min

    return xx_prop, cur_lnpdf


# Build RI Brownian covariance matrix.
cov_K = np.zeros((Nfinal, Nfinal), dtype=np.float64)
for i in range(Nfinal):
    for j in range(i, Nfinal):
        if j < Ngrid + N and i < Ngrid + N:
            kij = expo_quad_kernel(xxx[i], xxx[j])
        elif j >= Ngrid + N and i < Ngrid + N:
            kij = expo_quad_kernel2(
                xxx[i],
                days[j - Ngrid - N],
                days[j - Ngrid - N + 1],
            )
        else:
            kij = expo_quad_kernel3(
                days[i - Ngrid - N],
                days[i - Ngrid - N + 1],
                days[j - Ngrid - N],
                days[j - Ngrid - N + 1],
            )
        cov_K[i, j] = kij
        cov_K[j, i] = kij

ones_vec = np.ones(Nfinal)
ones_vec[Ngrid + N:] = diff
ones_vec = ones_vec.reshape(-1, 1)

# Same BM correction used by the supplied synthetic-data code.
sigma0 = np.sqrt(alpha / beta) * const
cov_K_corrected = cov_K + (ones_vec @ ones_vec.T) * sigma0 ** 2
cov_K_mod_chol = np.linalg.cholesky(cov_K_corrected)
cov_K_mod_inv = np.linalg.inv(cov_K_corrected)

theta = alpha / beta
g_mk_list2 = []
g_mk_list3 = []
theta_list = []
curloglike = None

start_time = time.time()
for ite in range(nsim1 + nsim2):
    if ite % 10000 == 0:
        print(f"iteration {ite}/{nsim1 + nsim2}", flush=True)

    cov_K_chol_final = cov_K_mod_chol / np.sqrt(theta)
    prior = chol_sample(np.zeros(Nfinal), cov_K_chol_final)

    g_mk, curloglike = elliptical_slice(
        g_mk.reshape(1, -1),
        prior,
        log_lik,
        pdf_params=[N],
        cur_lnpdf=curloglike,
        angle_range=None,
    )

    if ite % rem == 0:
        g_mk_list2.append(1.0 / g_mk[0][0:Ngrid])
        g_mk_list3.append(1.0 / g_mk[0][Ngrid:Ngrid + N])

        alpha_pos = alpha + Nfinal / 2
        beta_pos = beta + 0.5 * g_mk[0] @ cov_K_mod_inv @ g_mk[0]
        theta = float(np.random.gamma(alpha_pos, 1.0 / beta_pos))
        theta_list.append(theta)

runtime = time.time() - start_time
timerun_10000 = runtime / (nsim1 + nsim2) * 10000

# Burn-in is measured in original MCMC iterations; the saved RI chain is thinned by rem.
g_mk_list2 = np.asarray(g_mk_list2)
g_mk_list3 = np.asarray(g_mk_list3)
theta_list = np.asarray(theta_list)
nn = int(nsim1 / rem)
samples = g_mk_list2[nn:]
med = np.quantile(samples, 0.5, axis=0)


def hpd_interval(x, alpha=0.05):
    x = np.sort(np.asarray(x, dtype=float))
    n = len(x)
    k = int(np.floor((1 - alpha) * n))
    widths = x[k:] - x[:n - k]
    j = np.argmin(widths)
    return x[j], x[j + k]


low = np.zeros(Ngrid)
high = np.zeros(Ngrid)
for j in range(Ngrid):
    low[j], high[j] = hpd_interval(samples[:, j], alpha=0.05)

output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(
    output_dir,
    f"real100_BM_RI_const_{num_tag(const)}_nsim1_{nsim1}_nsim2_{nsim2}.npz",
)

np.savez(
    output_file,
    aaa=g_mk_list3,
    aa=g_mk_list2,
    c=points_inhomo,
    d=x4,
    p=timerun_10000,
    s=theta_list,
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
