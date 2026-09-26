#!/usr/bin/env python
# Resolve data and results independently of the working directory.
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from _paths import data_file, output_file

from scipy.stats import expon
from scipy.stats import uniform
from scipy.stats import norm
from scipy.stats import multivariate_normal
from numpy.random import multinomial
from numpy.random import uniform
import numpy as np
import matplotlib.pyplot as plt
import scipy 
from scipy.stats import truncnorm
import math
from scipy.stats import mvn
import pyreadr
import time
import sys
from math import comb
from numpy.linalg import LinAlgError
import random
from scipy.linalg import cho_solve, cho_factor
from scipy.stats import mvn
from statsmodels.sandbox.distributions.extras import mvstdnormcdf
np.set_printoptions(formatter={'float': '{:0.17f}'.format})
np.random.seed(123)


# -------------------------- data & setup -------------------------------------
data = pyreadr.read_r(data_file('syn3_std_ntip50.rda'))
points_inhomo = np.array(data["coal_time"]).squeeze()
points_inhomo=points_inhomo[int(sys.argv[1])]#args[1] range from 0 to 
T = np.max(points_inhomo)
ntip = 50

N = ntip - 1
Ngrid = 100
K = N
x4 = np.linspace(0, T, Ngrid + 1)[1:]
xxx = np.concatenate([x4, points_inhomo])

days = np.concatenate([[0], points_inhomo])
diff = np.diff(days)
Nfinal = Ngrid + N + K

# Construct g_mk
com_vec = np.zeros(ntip - 1)
for mmm in range(ntip - 1):
    com_vec[mmm] = comb(ntip - mmm, 2)
g_mk = 3 + np.zeros(Nfinal)
g_mk[Ngrid + N:Ngrid + N + K] = 3 * days[N] * diff / sum(diff)

nsim1 = 100000
nsim2 = 100000

def inten2(t):
    return 25*np.exp(-5*t)

# -------------------------- kernels ------------------------------------------
def expo_quad_kernel(theta1, xn, xm):
    return np.exp(-theta1/2*np.sum((xn - xm)**2))

def expo_quad_kernel2(theta1, xn, t1, t2):
    rt = np.sqrt(theta1/2)
    return np.sqrt(np.pi/2/theta1)*(math.erf(rt*(t2-xn)) - math.erf(rt*(t1-xn)))

def expo_quad_kernel3(theta1, t1, t2, t3, t4):
    rt = np.sqrt(theta1/2)
    term = (t2-t3)*math.erf(rt*(t2-t3)) - (t2-t4)*math.erf(rt*(t2-t4)) \
         - (t1-t3)*math.erf(rt*(t1-t3)) + (t1-t4)*math.erf(rt*(t1-t4))
    exp_part = np.exp(-theta1/2*(t2-t3)**2) - np.exp(-theta1/2*(t1-t3)**2) \
             - np.exp(-theta1/2*(t2-t4)**2) + np.exp(-theta1/2*(t1-t4)**2)
    return np.sqrt(np.pi/2/theta1)*term + (1/theta1)*exp_part

# -------------------------- sampling utils -----------------------------------
def chol_sample(mean, cov_chol):
    return mean + cov_chol @ np.random.standard_normal(mean.size)

def log_lik(f, ns):
    if np.prod(f > 0) == 0:
        return float('-inf')
    else:
        return np.sum(np.log(f[:, Ngrid:Ngrid+ns])) - np.sum(com_vec * (f[:, (Ngrid+ns):Nfinal]))

def elliptical_slice(initial_theta, prior, lnpdf, pdf_params=(),
                     cur_lnpdf=None, angle_range=None):
    D = len(initial_theta)
    if cur_lnpdf is None:
        cur_lnpdf = lnpdf(initial_theta, *pdf_params)

    if len(prior.shape) == 1:
        nu = prior
    else:
        if not (prior.shape[0] == D and prior.shape[1] == D):
            raise IOError("Prior must be a D-element sample or DxD chol(Sigma)")
        nu = prior @ np.random.normal(size=D)

    hh = math.log(np.random.uniform()) + cur_lnpdf

    if angle_range is None or angle_range == 0.:
        phi = np.random.uniform()*2.*math.pi
        phi_min = phi - 2.*math.pi
        phi_max = phi
    else:
        phi_min = -angle_range*np.random.uniform()
        phi_max = phi_min + angle_range
        phi = np.random.uniform()*(phi_max - phi_min) + phi_min

    while True:
        xx_prop = initial_theta*math.cos(phi) + nu*math.sin(phi)
        cur_lnpdf = lnpdf(xx_prop, *pdf_params)
        if cur_lnpdf > hh:
            break
        if phi > 0:
            phi_max = phi
        elif phi < 0:
            phi_min = phi
        else:
            raise RuntimeError('BUG DETECTED: Shrunk to current position and still not acceptable.')
        phi = np.random.uniform()*(phi_max - phi_min) + phi_min
    return (xx_prop, cur_lnpdf)

# ---------- fast MVN prob using Cholesky only (no full Sigma) ----------------
def mvn_prob_from_chol_fast(lower, upper, mean, L,
                            maxpts=2_000_000, abseps=1e-3, releps=1e-3):
    lower = np.asarray(lower, float).ravel()
    upper = np.asarray(upper, float).ravel()
    mean  = np.asarray(mean,  float).ravel()

    # std devs (sqrt of diag(Sigma)) = row norms of L
    s = np.linalg.norm(L, axis=1)
    s = np.maximum(s, 1e-12)

    lower_z = (lower - mean) / s
    upper_z = (upper - mean) / s

    L_tilde = L / s[:, None]
    R = L_tilde @ L_tilde.T
    R = 0.5 * (R + R.T)

    p = mvstdnormcdf(lower_z, upper_z, R,
                     maxpts=int(maxpts), abseps=abseps, releps=releps)
    return float(max(p, 1e-300))

# ------------------------ build Cholesky (discard K) -------------------------
def build_chol(theta1, jitter):
    """Build kernel K (locally), return only its jittered Cholesky; discard K."""
    K = np.zeros((Nfinal, Nfinal), dtype=np.float64)
    for i in range(Nfinal):
        for j in range(i, Nfinal):
            if j < Ngrid + N and i < Ngrid + N:
                kij = expo_quad_kernel(theta1, xxx[i], xxx[j])
            elif j >= Ngrid + N and i < Ngrid + N:
                kij = expo_quad_kernel2(theta1, xxx[i],
                                        days[j - Ngrid - N], days[j - Ngrid - N + 1])
            else:
                kij = expo_quad_kernel3(theta1,
                                        days[i - Ngrid - N], days[i - Ngrid - N + 1],
                                        days[j - Ngrid - N], days[j - Ngrid - N + 1])
            K[i, j] = kij
            if j != i:
                K[j, i] = kij


    cur_jitter = jitter
    while True:
        try:
            L = np.linalg.cholesky(K + np.eye(Nfinal) * cur_jitter)
            return L
        except LinAlgError as e:
            print(f"Caught LinAlgError: {e}. Increasing jitter to {cur_jitter*5}")
            cur_jitter *= 5

# ------------------------ posterior & MH (chol-only) -------------------------
def log_posterior(l, mu_l, scale_l, f, theta0, jitter, cov_K_chol=None,mvn_tol=(2_000_000, 1e-3, 1e-3)):
    theta11 = np.exp(-2*l)

    built = False
    if cov_K_chol is None:
        cov_K_chol = build_chol(theta11, jitter)
        built = True

    theta0_scalar = np.asarray(theta0).item()
    comp2 = -0.5 * theta0_scalar * (f @ cho_solve((cov_K_chol, True), f))
    comp1 = norm.logpdf(l, loc=mu_l, scale=scale_l)

    d = len(f)
    mean = np.zeros(d)
    lower = np.zeros(d)
    upper = np.full(d, np.inf)
    maxpts, abseps, releps = mvn_tol
    
    prob = mvn_prob_from_chol_fast(lower, upper, mean, cov_K_chol,
                                   maxpts=maxpts, abseps=abseps, releps=releps)
    
    comp3 = np.log(prob) + np.sum(np.log(np.diag(cov_K_chol)))
    subtotal=comp1-comp3
    total = subtotal + comp2 
    if built:
        return total,subtotal, cov_K_chol
    return total,subtotal

def metropolis_hastings_update(l_current, proposal_sigma_l, mu_l, scale_l,
                               f, theta0, cov_K_chol_current, 
                               jitter, log_post_current_partial=None,mvn_tol=(2_000_000, 1e-3, 1e-3)):
    if log_post_current_partial is None:
        log_post_current,log_post_current_partial = log_posterior(l_current, mu_l, scale_l, f, theta0,jitter, 
                                         cov_K_chol_current, mvn_tol)
    else:
        log_post_current = log_post_current_partial-0.5 * np.asarray(theta0).item() * (f @ cho_solve((cov_K_chol_current, True), f))
    

    l_prop = l_current + np.random.normal(0, proposal_sigma_l)
    lp_prop,log_post_prop_partial, L_prop = log_posterior(l_prop, mu_l, scale_l, f, theta0,jitter=jitter,
                                                 cov_K_chol=None,  mvn_tol=mvn_tol)

    log_alpha = lp_prop - log_post_current 
   
    if np.log(np.random.uniform()) < log_alpha:
        return l_prop,log_post_prop_partial, L_prop,  True
    else:
        return l_current,log_post_current_partial, cov_K_chol_current, False

# ---------------------------- init & run -------------------------------------
theta1 = np.float64(10.0)
l = np.float64(-np.log(theta1) / 2)
proposal_sigma_l = 0.01
mu_l = l
scale_l = 1e4
jitter = 1e-12
cov_K_chol = build_chol(theta1, jitter)
alpha = 0.1
beta  = 0.1
theta = 100.0
rem = 10
g_mk_list2, g_mk_list3, theta_list, theta1_list = [], [], [], []
curloglike = None
partialloglik=None
mvn_tol = (2_000_000, 1e-3, 1e-3)
acc_count=0
start_time = time.time()
for ite in range(nsim1 + nsim2):
    cov_K_chol_final = cov_K_chol / np.sqrt(theta)
    prior = chol_sample(mean=np.zeros(Nfinal), cov_chol=cov_K_chol_final)
    g_mk, curloglike = elliptical_slice(g_mk.reshape(1, -1), prior, log_lik,
                                        pdf_params=[N], cur_lnpdf=curloglike, angle_range=None)

    if (ite % rem) == 0:
        g_mk_list2.append(1/g_mk[0][0:Ngrid])
        g_mk_list3.append(1/g_mk[0][Ngrid:Ngrid+N])
        alpha_pos = alpha + Nfinal/2
        u = cho_solve((cov_K_chol, True), g_mk[0])
        beta_pos = beta + 0.5 * (g_mk[0] @ u)
        theta = float(np.random.gamma(alpha_pos, 1/beta_pos))
        theta_list.append(theta)
        l, partialloglik,cov_K_chol,acc = metropolis_hastings_update(
            l, proposal_sigma_l, mu_l, scale_l, g_mk[0], theta,
            cov_K_chol, jitter, partialloglik,mvn_tol)
        theta1 = np.exp(-2*l)
        theta1_list.append(theta1)
        if acc:
            acc_count+=1

timerun = time.time() - start_time
timerun_10000 = timerun / (nsim2 + nsim1) * 10000
acc_rate=acc_count/int((nsim2+nsim1)/rem)
nn = int(nsim1/rem)
low  = np.quantile(np.array(g_mk_list2)[nn:,], 0.025, axis=0)
high = np.quantile(np.array(g_mk_list2)[nn:,], 0.975, axis=0)
med  = np.quantile(np.array(g_mk_list2)[nn:,], 0.5, axis=0)
truth = inten2(x4)
l2_dist1 = sum((np.array(med).squeeze() - truth)**2)
coverage1 = np.sum((truth >= low.squeeze()) * (truth <= high.squeeze())) / len(x4)
width1 = sum(high - low) / Ngrid
np.savez(output_file('stdcoaldata/syn3/SE50/syn3_tips50_data'+sys.argv[1]+'fullbayes_step001_jitter1e-12.npz'), aaa=g_mk_list3,aa=g_mk_list2,c=points_inhomo,d=x4,e=truth,f=coverage1,i=jitter,j=l2_dist1,p=timerun_10000,q=width1,s=theta_list,t=acc_rate,z=theta1_list)

