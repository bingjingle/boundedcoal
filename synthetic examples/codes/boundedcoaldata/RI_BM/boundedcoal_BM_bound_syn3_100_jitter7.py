#!/usr/bin/env python
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
np.set_printoptions(suppress=True)
import pyreadr
import time
import sys
from math import comb
from numpy.linalg import LinAlgError
import random
from scipy.linalg import cho_solve, cho_factor
import math
from math import lgamma, exp
from scipy.special import comb 
np.random.seed(123)



data = pyreadr.read_r('/home/groups/juliapr/Bingjing/code/boundedcoal/syndatafinal/syn3_bounded_ntip100.rda')
points_inhomo = np.array(data["coal_time"]).squeeze()
points_inhomo=points_inhomo[int(sys.argv[1])]#args[1] range from 0 to 29
T=0.71
ntip=100
N = ntip - 1
Ngrid = 100
K= ntip  # counts of integrals
x4 = np.linspace(0, T, Ngrid + 1)[1:]
xxx=np.concatenate([x4, points_inhomo])

days = np.concatenate([[0], points_inhomo,[T]])
diff = np.diff(days)
Nfinal = Ngrid + N + K
g_mk = 3 + np.zeros(Nfinal)
g_mk[ Ngrid+N:Ngrid+N+K ]=3*T*diff/sum(diff)
nsim1=1000000
nsim2=1000000

def inten2(t):
    return 25*np.exp(-5*t)


# Build covariance matrix with stable kernel
def expo_quad_kernel(xn, xm):
    return min(xn, xm)
def expo_quad_kernel2(xn, t1, t2):
    # Ensure type consistency
    xn, t1, t2 = float(xn), float(t1), float(t2)
    if xn <= t1:
        return xn * (t2 - t1)
    elif t1 < xn < t2:
        return -0.5 * (xn**2) - 0.5 * (t1**2) + xn * t2
    else:  # xn >= t2
        return 0.5 * (t2 - t1) * (t2 + t1)  # Better than (t2**2 - t1**2)/2

def expo_quad_kernel3(t1, t2, t3, t4):
    # Use stable equality comparison
    if (t1==t3 and t2==t4 ):
        return (t2**3)/3 + 2*(t1**3)/3 - 0.5*(t1**2)*t2 - 0.5*t2*(t1**2)
    else:
        return 0.5 * (t4 - t3) * (t2 + t1) * (t2 - t1)  # Rewritten for stability

dtype = np.longdouble
LD = dtype  # alias

def r_stable(k: int, j: int) -> np.longdouble:
    # All constants and ops in longdouble
    if j == 1:
        return LD(1)

    sign  = LD(-1) if ((j - 1) % 2) else LD(1)
    term1 = LD(2*j - 1)

    # log ∏_{m=1}^{j-1} (k-m)/(k+m) = ∑_{m=1}^{j-1} [log(k-m) - log(k+m)]
    # do it directly in longdouble instead of math.lgamma
    m = np.arange(1, j, dtype=LD)
    kLD = LD(k)
    log_ratio = np.sum(np.log(kLD - m) - np.log(kLD + m), dtype=LD)

    return LD(sign * term1 * np.exp(log_ratio))

# ensure integer
ntip = int(ntip)

# comb(n,2) = n*(n-1)/2 done in longdouble (avoid scipy.special.comb exact=False)
com_vec = np.fromiter(
    (LD(ntip - m) * LD(ntip - m - 1) / LD(2) for m in range(ntip-1)),
    dtype=LD
)

# r_stable values as longdouble (iterator yields longdouble)
r_vec = np.fromiter(
    (r_stable(ntip, ntip - m) for m in range(ntip-1)),
    dtype=LD
)


def chol_sample(mean, cov_chol):
    return mean + cov_chol @ np.random.standard_normal(mean.size)


def log_lik(f, ns):
    """
    Numerically stable version of your log likelihood.
    - Stabilizes sum r_vec*exp(-com_vec * sum(f[...] ))
    - Uses log1p
    - Accumulates in longdouble if available
    """
    dtype = np.longdouble
    f = np.asarray(f, dtype=dtype)

    # 1) domain check only where log is taken
    #    (avoid rejecting just because some unused entries of f are <= 0)
    if np.any(f <= 0):
        return -np.inf
    # ----- term 1: sum log f -----
    term1 = np.sum(np.log(f[:, Ngrid:Ngrid+ns], dtype=dtype), dtype=dtype)

    # ----- term 2: - sum(com_vec * f_slice) -----
    # Your original does: -np.sum(com_vec * f[:, Ngrid+ns : Nfinal-1])
    # This relies on broadcasting; be explicit about accumulation dtype.
    term2 = -np.sum(com_vec * f[:, Ngrid+ns:Nfinal-1], dtype=dtype)

    # ----- term 3: -log(1 + sum(r_vec * exp(-com_vec * S))) -----
    # where S = sum(f[:, Ngrid+ns:Nfinal])  (your code used a scalar S)
    S = np.sum(f[:, Ngrid+ns:Nfinal], dtype=dtype)       # scalar
    s = -com_vec * S                                     # vector of exponents

    # Max-shift to keep exp() well-scaled, even if s is very negative
    m = np.max(s)                                        # scalar
    # sum_r_exp = sum_j r_j * exp(s_j - m)
    sum_r_exp = np.sum(r_vec * np.exp(s - m, dtype=dtype), dtype=dtype)

    # We need -log(1 + exp(m) * sum_r_exp) computed stably.
    y = np.exp(m, dtype=dtype) * sum_r_exp               # may be tiny or huge
    denom = 1 + y                                        # use longdouble here

    # Guard: if alternating signs in r_vec make denom <= 0, likelihood is invalid
    if not np.isfinite(denom) or denom <= 0:
        return -np.inf

    term3 = -np.log(denom, dtype=dtype)                  # log1p equivalent already done
    return term1 + term2 + term3
    
def elliptical_slice(initial_theta,prior,lnpdf,pdf_params=(),
                     cur_lnpdf=None,angle_range=None):
    """
    NAME:
       elliptical_slice
    PURPOSE:
       Markov chain update for a distribution with a Gaussian "prior" factored out
    INPUT:
       initial_theta - initial vector
       prior - cholesky decomposition of the covariance matrix 
               (like what numpy.linalg.cholesky returns), 
               or a sample from the prior
       lnpdf - function evaluating the log of the pdf to be sampled
       pdf_params= parameters to pass to the pdf
       cur_lnpdf= value of lnpdf at initial_theta (optional)
       angle_range= Default 0: explore whole ellipse with break point at
                    first rejection. Set in (0,2*pi] to explore a bracket of
                    the specified width centred uniformly at random.
    OUTPUT:
       new_theta, new_lnpdf
    HISTORY:
       Originally written in matlab by Iain Murray (http://homepages.inf.ed.ac.uk/imurray2/pub/10ess/elliptical_slice.m)
       2012-02-24 - Written - Bovy (IAS)
    """
    D= len(initial_theta)
    if cur_lnpdf is None:
        cur_lnpdf= lnpdf(initial_theta,*pdf_params)

    # Set up the ellipse and the slice threshold
    if len(prior.shape) == 1: #prior = prior sample
        nu= prior
    else: #prior = cholesky decomp
        if not prior.shape[0] == D or not prior.shape[1] == D:
            raise IOError("Prior must be given by a D-element sample or DxD chol(Sigma)")
        nu= np.dot(prior,np.random.normal(size=D))
    hh = math.log(np.random.uniform()) + cur_lnpdf

    # Set up a bracket of angles and pick a first proposal.
    # "phi = (theta'-theta)" is a change in angle.
    if angle_range is None or angle_range == 0.:
        # Bracket whole ellipse with both edges at first proposed point
        phi= np.random.uniform()*2.*math.pi
        phi_min= phi-2.*math.pi
        phi_max= phi
    else:
        # Randomly center bracket on current point
        phi_min= -angle_range*np.random.uniform()
        phi_max= phi_min + angle_range
        phi= np.random.uniform()*(phi_max-phi_min)+phi_min

    # Slice sampling loop
    while True:
        # Compute xx for proposed angle difference and check if it's on the slice
        xx_prop = initial_theta*math.cos(phi) + nu*math.sin(phi)
        cur_lnpdf = lnpdf(xx_prop,*pdf_params)
        if cur_lnpdf > hh:
            # New point is on slice, ** EXIT LOOP **
            break
        # Shrink slice to rejected point
        if phi > 0:
            phi_max = phi
        elif phi < 0:
            phi_min = phi
        else:
            raise RuntimeError('BUG DETECTED: Shrunk to current position and still not acceptable.')
        # Propose new angle difference
        phi = np.random.uniform()*(phi_max - phi_min) + phi_min
    return (xx_prop,cur_lnpdf)


cov_K = np.zeros((Nfinal, Nfinal), dtype=np.float64)
theta=100
rem=10
alpha=.1
beta=.1
g_mk_list2=[]
g_mk_list3=[]
theta_list=[]

for i in range(Nfinal):
        for j in range(i,Nfinal):
            if (j<Ngrid+N) and (i<Ngrid+N):
                cov_K[i][j]=expo_quad_kernel(xxx[i],xxx[j])
            if (j>(Ngrid+N-1)) and (i<(Ngrid+N)):
                cov_K[i][j]=expo_quad_kernel2(xxx[i],   days[j-Ngrid-N], days[j-Ngrid-N+1])
            if  (i>(Ngrid+N-1)):
                cov_K[i][j]=expo_quad_kernel3(days[i-Ngrid-N],days[i-Ngrid-N+1],  days[j-Ngrid-N],days[j-Ngrid-N+1 ]          )  
            if j!=i:
                cov_K[j][i]=cov_K[i][j]
# Construct ones_vec
ones_vec = np.ones(Nfinal, dtype=np.float64)
ones_vec[Ngrid + N:] = diff
ones_vec = ones_vec.reshape(-1, 1)

jitter =1e-7
    
cov_K_jittered = cov_K + np.eye(Nfinal) * jitter

# Cholesky solve
cov_K_chol = np.linalg.cholesky(cov_K_jittered)
w = cho_solve((cov_K_chol, True), ones_vec)

# Accuracy diagnostics
A_star_chol = cov_K_chol @ cov_K_chol.T
chol_error = np.linalg.norm(A_star_chol - cov_K, ord='fro')
solve_error = np.linalg.norm(cov_K @ w - ones_vec, ord='fro')
d = ones_vec.T @ w

print("Cholesky approx error:", chol_error)
print("Linear solve error:", solve_error)
print("Min eigenvalue:", np.min(np.linalg.eigvalsh(cov_K)))
print("Condition number:", np.linalg.cond(cov_K))
print("dtype:", cov_K.dtype)
print("scalar d:", d.item())

noise_var = 1e-15
K = cov_K * noise_var + np.eye(Nfinal)
K = K.astype(LD)
L, lower = cho_factor(K, lower=True)

# Helper to avoid forming full inverse
def apply_Kinv(B):
    return cho_solve((L, lower), B)

A_inv_eye = cov_K @ apply_Kinv(np.eye(Nfinal))
A_inv_w = cov_K @ apply_Kinv(w)
w_T_A_inv_w = np.float64((w.T @ A_inv_w).item())
cov_K_mod = A_inv_eye + (A_inv_w @ A_inv_w.T) / (np.float64(d.item())-w_T_A_inv_w)

cov_K_mod_chol = np.linalg.cholesky(cov_K_mod + jitter * np.eye(Nfinal))
error_chol = np.linalg.norm(cov_K_mod_chol @ cov_K_mod_chol.T - cov_K_mod, ord='fro')
print("Cholesky reconstruction error:", error_chol)



curloglike = None
start_time=time.time()
for ite in range(nsim1+nsim2):
    cov_K_chol_final=cov_K_mod_chol/np.sqrt(theta)
    prior=chol_sample(mean=np.zeros(Nfinal), cov_chol=cov_K_chol_final)#nu
    g_mk,curloglike=elliptical_slice(g_mk.reshape(1,-1),prior,log_lik,pdf_params=[N],cur_lnpdf=curloglike,angle_range=None)

    if (ite%rem==0):
        g_mk_list2.append(1/g_mk[0][0:Ngrid])
        g_mk_list3.append(1/g_mk[0][Ngrid:Ngrid+N])
        alpha_pos=alpha+Nfinal/2
        u = cho_solve((cov_K_mod_chol, True), g_mk[0])
        beta_pos=beta+1/2*g_mk[0]@u
        theta=np.random.gamma(alpha_pos, 1/beta_pos,1)
        theta_list.append(theta)
timerun=time.time()-start_time
timerun_10000=timerun/(nsim2+nsim1)*10000

       

cond_number = np.linalg.cond(cov_K)
print(cond_number)
cond_number = np.linalg.cond(cov_K_mod)
print(cond_number)
cond_number = np.linalg.cond(cov_K_mod+np.eye(Nfinal)*jitter)
print(cond_number)
cond_number = np.linalg.cond(cov_K_mod_chol @ cov_K_mod_chol.T )
print(cond_number)


nn=int(nsim1/rem)       
#grids
low=np.quantile(np.array(g_mk_list2)[nn:,], 0.025, axis=0)
high=np.quantile(np.array(g_mk_list2)[nn:,], 0.975, axis=0)
med=np.quantile(np.array(g_mk_list2)[nn:,], 0.5, axis=0)
truth=inten2(x4)
l2_dist1=sum((np.array(med).squeeze()-truth)**2)
coverage1=np.sum((truth>=low.squeeze()) * (truth<=high.squeeze()))/len(x4)
width1=sum(high-low)/Ngrid
np.savez('/scratch/groups/juliapr/output_Bingjing/boundedcoaldata/bounded/syn3/RI_BM/syn3_tips100_data_'+sys.argv[1]+'jitter7.npz', aaa=g_mk_list3,aa=g_mk_list2,c=points_inhomo,d=x4,e=truth,f=coverage1,i=jitter,j=l2_dist1,p=timerun_10000,q=width1,s=theta_list)
