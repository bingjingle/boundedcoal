#!/usr/bin/env python

import argparse
import math
import time
from pathlib import Path

import numpy as np
from scipy.linalg import cho_solve
from scipy.special import comb
from scipy.stats import norm
from statsmodels.sandbox.distributions.extras import mvstdnormcdf


# -----------------------------------------------------------------------------
# Arguments
# -----------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Random-integral inference with a squared-exponential kernel."
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
        required=True,
        help="CSV with coal_times column, or an R data file (requires pyreadr).",
    )
    p.add_argument("--output-dir", default="./out")
    p.add_argument("--ntip", type=int, default=100)
    p.add_argument("--bound", type=float, default=25.0,
                   help="tau; for the standard coalescent it only sets the grid extent.")
    p.add_argument("--theta1-init", type=float, default=10.0)
    p.add_argument("--standard", action="store_true",
                   help="standard-coalescent likelihood instead of bounded.")
    p.add_argument("--seed", type=int, default=123)
    p.add_argument("--tag", default=None)
    p.add_argument("--hyper-every", type=int, default=1,
                   help="run the lengthscale MH update on every k-th theta update "
                        "(1 = upstream behaviour). Each update costs two 305-dimensional "
                        "Gaussian orthant probabilities, ~0.5 s apiece.")
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

np.random.seed(args.seed)
BOUNDED = not args.standard

# -----------------------------------------------------------------------------
# UPGMA100 real data: same setup as the pasted real-data Python code
# -----------------------------------------------------------------------------
if args.input_file.endswith(".rda"):
    import pyreadr
    data = pyreadr.read_r(args.input_file)
    points_inhomo = np.asarray(data["bt_adj"], dtype=float).squeeze()
else:
    import pandas as pd
    points_inhomo = np.asarray(pd.read_csv(args.input_file)["coal_times"], dtype=float)
points_inhomo = np.sort(points_inhomo)

T = args.bound
ntip = args.ntip
N = ntip - 1
Ngrid = 100
K = ntip if BOUNDED else ntip - 1

x4 = np.linspace(0.0, T, Ngrid + 1)[1:]
xxx = np.concatenate([x4, points_inhomo])

days = (np.concatenate([[0.0], points_inhomo, [T]]) if BOUNDED
        else np.concatenate([[0.0], points_inhomo]))
diff = np.diff(days)
Nfinal = Ngrid + N + K

g_mk = 3.0 + np.zeros(Nfinal)
g_mk[Ngrid + N:Ngrid + N + K] = 3.0 * T * diff / np.sum(diff)

alpha = 0.1
beta = 0.1


# -----------------------------------------------------------------------------
# Bounded-coalescent likelihood ingredients
# -----------------------------------------------------------------------------
def expo_quad_kernel(theta1, xn, xm):
    return np.exp(-theta1 / 2.0 * np.sum((xn - xm) ** 2))


def expo_quad_kernel2(theta1, xn, t1, t2):
    rt = np.sqrt(theta1 / 2.0)
    return np.sqrt(np.pi / 2.0 / theta1) * (
        math.erf(rt * (t2 - xn)) - math.erf(rt * (t1 - xn))
    )


def expo_quad_kernel3(theta1, t1, t2, t3, t4):
    rt = np.sqrt(theta1 / 2.0)
    term = (
        (t2 - t3) * math.erf(rt * (t2 - t3))
        - (t2 - t4) * math.erf(rt * (t2 - t4))
        - (t1 - t3) * math.erf(rt * (t1 - t3))
        + (t1 - t4) * math.erf(rt * (t1 - t4))
    )
    exp_part = (
        np.exp(-theta1 / 2.0 * (t2 - t3) ** 2)
        - np.exp(-theta1 / 2.0 * (t1 - t3) ** 2)
        - np.exp(-theta1 / 2.0 * (t2 - t4) ** 2)
        + np.exp(-theta1 / 2.0 * (t1 - t4) ** 2)
    )
    return np.sqrt(np.pi / 2.0 / theta1) * term + exp_part / theta1


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


n_vec = np.arange(ntip, 1, -1)
com_vec = n_vec * (n_vec - 1) / 2.0


def chol_sample(mean, cov_chol):
    return mean + cov_chol @ np.random.standard_normal(mean.size)


def log_lik(f, ns):
    if np.any(f <= 0):
        return -np.inf

    term1 = np.sum(np.log(f[:, Ngrid:Ngrid + ns]))
    if not BOUNDED:
        return term1 - np.sum(com_vec * f[:, Ngrid + ns:Nfinal])
    term2 = -np.sum(com_vec * f[:, Ngrid + ns:Nfinal - 1])

    Lambda = np.sum(f[:, Ngrid + ns:Nfinal])
    x = np.exp(-Lambda)
    val = rhs_value(x)

    if (not np.isfinite(val)) or val <= 0:
        return -np.inf

    log_one_minus_x = np.log(-np.expm1(-Lambda))
    logboundprob = (ntip - 1) * log_one_minus_x + np.log(val)
    return term1 + term2 - logboundprob


def elliptical_slice(initial_theta, prior, lnpdf, pdf_params=(),
                     cur_lnpdf=None, angle_range=None):
    if cur_lnpdf is None:
        cur_lnpdf = lnpdf(initial_theta, *pdf_params)

    if len(prior.shape) == 1:
        nu = prior
    else:
        nu = prior @ np.random.normal(size=initial_theta.size)

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
            raise RuntimeError("BUG DETECTED")

        phi = np.random.uniform() * (phi_max - phi_min) + phi_min

    return xx_prop, cur_lnpdf


def mvn_prob_from_chol_fast(lower, upper, mean, L,
                            maxpts=2_000_000, abseps=1e-3, releps=1e-3):
    lower = np.asarray(lower, float).ravel()
    upper = np.asarray(upper, float).ravel()
    mean = np.asarray(mean, float).ravel()

    s = np.linalg.norm(L, axis=1)
    s = np.maximum(s, 1e-12)

    lower_z = (lower - mean) / s
    upper_z = (upper - mean) / s
    L_tilde = L / s[:, None]
    R = L_tilde @ L_tilde.T
    R = 0.5 * (R + R.T)

    p = mvstdnormcdf(
        lower_z, upper_z, R,
        maxpts=int(maxpts), abseps=abseps, releps=releps,
    )
    return float(max(p, 1e-300))


def build_chol(theta1, jitter_start):
    Kmat = np.zeros((Nfinal, Nfinal), dtype=np.float64)

    for i in range(Nfinal):
        for j in range(i, Nfinal):
            if j < Ngrid + N and i < Ngrid + N:
                kij = expo_quad_kernel(theta1, xxx[i], xxx[j])
            elif j >= Ngrid + N and i < Ngrid + N:
                kij = expo_quad_kernel2(
                    theta1, xxx[i],
                    days[j - Ngrid - N],
                    days[j - Ngrid - N + 1],
                )
            else:
                kij = expo_quad_kernel3(
                    theta1,
                    days[i - Ngrid - N],
                    days[i - Ngrid - N + 1],
                    days[j - Ngrid - N],
                    days[j - Ngrid - N + 1],
                )

            Kmat[i, j] = kij
            Kmat[j, i] = kij

    jitter_try = jitter_start
    while jitter_try <= 1.0:
        try:
            return np.linalg.cholesky(
                Kmat + jitter_try * np.eye(Nfinal)
            )
        except np.linalg.LinAlgError:
            jitter_try *= 10.0

    raise np.linalg.LinAlgError("Cholesky failed")


def log_posterior(l_value, mu_l, scale_l, f, theta0, jitter_start,
                  cov_K_chol=None,
                  mvn_tol=(2_000_000, 1e-3, 1e-3)):
    theta11 = np.exp(-2.0 * l_value)

    built = False
    if cov_K_chol is None:
        cov_K_chol = build_chol(theta11, jitter_start)
        built = True

    comp2 = -0.5 * float(theta0) * (
        f @ cho_solve((cov_K_chol, True), f)
    )
    comp1 = norm.logpdf(l_value, loc=mu_l, scale=scale_l)

    d = len(f)
    mean = np.zeros(d)
    lower = np.zeros(d)
    upper = np.full(d, np.inf)

    prob = mvn_prob_from_chol_fast(
        lower, upper, mean, cov_K_chol,
        maxpts=mvn_tol[0], abseps=mvn_tol[1], releps=mvn_tol[2],
    )

    comp3 = np.log(prob) + np.sum(np.log(np.diag(cov_K_chol)))
    subtotal = comp1 - comp3
    total = subtotal + comp2

    if built:
        return total, subtotal, cov_K_chol
    return total, subtotal


def metropolis_hastings_update(
        l_current, proposal_sigma_l, mu_l, scale_l,
        f, theta0, cov_K_chol_current, jitter_start,
        log_post_current_partial=None,
        mvn_tol=(2_000_000, 1e-3, 1e-3)):

    if log_post_current_partial is None:
        log_post_current, log_post_current_partial = log_posterior(
            l_current, mu_l, scale_l, f, theta0,
            jitter_start, cov_K_chol_current, mvn_tol,
        )
    else:
        log_post_current = (
            log_post_current_partial
            - 0.5 * float(theta0)
            * (f @ cho_solve((cov_K_chol_current, True), f))
        )

    l_prop = l_current + np.random.normal(0.0, proposal_sigma_l)
    lp_prop, log_post_prop_partial, L_prop = log_posterior(
        l_prop, mu_l, scale_l, f, theta0,
        jitter_start=jitter_start,
        cov_K_chol=None,
        mvn_tol=mvn_tol,
    )

    if np.log(np.random.uniform()) < lp_prop - log_post_current:
        return l_prop, log_post_prop_partial, L_prop, True

    return l_current, log_post_current_partial, cov_K_chol_current, False


# -----------------------------------------------------------------------------
# ED1 + MCMC
# -----------------------------------------------------------------------------
theta1 = args.theta1_init
l = -0.5 * np.log(theta1)
mu_l = l
scale_l = 1e4

cov_K_chol = build_chol(theta1, jitter)
theta = alpha / beta

g_mk_list2 = []
g_mk_list3 = []
theta_list = []
theta1_list = []
l_list = []

curloglike = None
partialloglik = None
acc_count = 0
mvn_tol = (2_000_000, 1e-3, 1e-3)

start_time = time.time()

for ite in range(nsim1 + nsim2):
    if ite % 1000 == 0:
        print(ite, flush=True)

    cov_K_chol_final = cov_K_chol / np.sqrt(theta)
    prior = chol_sample(np.zeros(Nfinal), cov_K_chol_final)

    g_mk, curloglike = elliptical_slice(
        g_mk.reshape(1, -1),
        prior,
        log_lik,
        pdf_params=[N],
        cur_lnpdf=curloglike,
    )

    g_mk_list2.append(1.0 / g_mk[0][0:Ngrid])
    g_mk_list3.append(1.0 / g_mk[0][Ngrid:Ngrid + N])

    if ite % rem == 0:
        alpha_pos = alpha + Nfinal / 2.0
        u = cho_solve((cov_K_chol, True), g_mk[0])
        beta_pos = beta + 0.5 * (g_mk[0] @ u)

        theta = float(np.random.gamma(alpha_pos, 1.0 / beta_pos))
        theta_list.append(theta)

        if (ite // rem) % args.hyper_every == 0:
            l, partialloglik, cov_K_chol, acc = metropolis_hastings_update(
                l, proposal_sigma_l, mu_l, scale_l,
                g_mk[0], theta, cov_K_chol,
                jitter, partialloglik, mvn_tol,
            )
        else:
            acc = False

        theta1 = np.exp(-2.0 * l)
        theta1_list.append(theta1)
        l_list.append(l)

        if acc:
            acc_count += 1


timerun = time.time() - start_time
timerun_10000 = timerun / (nsim1 + nsim2) * 10000.0
n_hyper_updates = len(theta_list)
acc_rate = acc_count / n_hyper_updates

# -----------------------------------------------------------------------------
# Posterior summaries: real data, so no synthetic truth/SSE/coverage is used.
# -----------------------------------------------------------------------------
all_Ne_grid = np.asarray(g_mk_list2)
samples = all_Ne_grid[nsim1:, :]
med = np.quantile(samples, 0.5, axis=0)


def hpd_interval(x, alpha=0.05):
    x = np.sort(np.asarray(x, dtype=float))
    n = len(x)
    k = int(np.floor((1.0 - alpha) * n))
    widths = x[k:] - x[:n-k]
    j = np.argmin(widths)
    return x[j], x[j + k]


low = np.zeros(Ngrid)
high = np.zeros(Ngrid)
for j in range(Ngrid):
    low[j], high[j] = hpd_interval(samples[:, j], alpha=0.05)


def num_tag(x):
    s = f"{x:.12g}"
    s = s.replace("e-0", "e-").replace("e+0", "e+")
    return s


output_dir = Path(args.output_dir)
output_dir.mkdir(parents=True, exist_ok=True)

step_tag = num_tag(proposal_sigma_l)
output_file = output_dir / (
    (args.tag + ".npz") if args.tag else
    f"real100_SE_RI_stepsize_{step_tag}_nsim1_{nsim1}_nsim2_{nsim2}.npz"
)

np.savez(
    output_file,
    aaa=np.asarray(g_mk_list3),
    aa=all_Ne_grid,
    c=points_inhomo,
    d=x4,
    s=np.asarray(theta_list),
    z=np.asarray(theta1_list),
    l=np.asarray(l_list),
    t=acc_rate,
    p=timerun_10000,
    med=med,
    low=low,
    high=high,
    proposal_sigma_l=proposal_sigma_l,
    nsim1=nsim1,
    nsim2=nsim2,
    jitter=jitter,
    rem=rem,
    T=T,
    ntip=ntip,
    grid=x4,
    bounded=BOUNDED,
    theta1_init=args.theta1_init,
    coal=points_inhomo,
)
print(f"DONE {args.tag}  med[0]={med[0]:.4f} med[50]={med[50]:.4f} "
      f"med[-1]={med[-1]:.4f}  hi_max={high.max():.4g}", flush=True)

print("Saved:", output_file)
print("MH acceptance rate:", acc_rate)
print("Runtime per 10,000 iterations:", timerun_10000)
