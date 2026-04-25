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

data = pyreadr.read_r('/home/groups/juliapr/Bingjing/code/boundedcoal/syndata1/UGPMA100new.rda')
points_inhomo = np.array(data["bt_adj"]).squeeze()
T=25
print(points_inhomo)

ntip=100  
N = ntip - 1
K= ntip  # counts of integrals

xxx=points_inhomo

days = np.concatenate([[0], points_inhomo,[T]])
diff = np.diff(days)
Nfinal =  N + K
g_mk = 3 + np.zeros(Nfinal)
g_mk[ N:N+K ]=3*T*diff/sum(diff)
nsim1=1500000
nsim2=1500000

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

def r_stable(k: int, j: int) -> float:
    if j == 1:
        return 1.0

    sign  = -1.0 if ((j - 1) % 2) else 1.0
    term1 = float(2*j - 1)

    m = np.arange(1, j, dtype=np.float64)
    kf = float(k)

    log_ratio = np.sum(
        np.log(kf - m) - np.log(kf + m),
        dtype=np.float64
    )

    return sign * term1 * np.exp(log_ratio)

# ensure integer
ntip = int(ntip)

m = np.arange(ntip, 1, -1, dtype=np.float64)
com_vec = m * (m - 1.0) / 2.0


r_vec = np.array(
    [r_stable(ntip, ntip - m) for m in range(ntip - 1)],
    dtype=np.float64
)

# -------------------------- sampling utils -----------------------------------
def chol_sample(mean, cov_chol):
    return mean + cov_chol @ np.random.standard_normal(mean.size)


def log_lik(f, ns):
    """
    Numerically stable version of your log likelihood.
    Uses double precision (float64).
    """
    dtype = np.float64
    f = np.asarray(f, dtype=dtype)

    # domain check only where log is taken
    if np.any(f <= 0):
        return -np.inf

    # ----- term 1: sum log f -----
    term1 = np.sum(np.log(f[:, 0:ns], dtype=dtype), dtype=dtype)

    # ----- term 2: - sum(com_vec * f_slice) -----
    term2 = -np.sum(com_vec * f[:, ns:Nfinal-1], dtype=dtype)

    # ----- term 3: -log(1 + sum(r_vec * exp(-com_vec * S))) -----
    S = np.sum(f[:, ns:Nfinal], dtype=dtype)
    s = -com_vec * S

    m = np.max(s)
    sum_r_exp = np.sum(r_vec * np.exp(s - m, dtype=dtype), dtype=dtype)

    y = np.exp(m, dtype=dtype) * sum_r_exp

    denom = 1 + y                                        # use longdouble here

    # Guard: if alternating signs in r_vec make denom <= 0, likelihood is invalid
    if not np.isfinite(denom) or denom <= 0:
        return -np.inf

    term3 = -np.log(denom, dtype=dtype)             
    
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
theta=50
rem=10
alpha=.1
beta=.1
g_mk_list2=[]
g_mk_list3=[]
theta_list=[]

for i in range(Nfinal):
        for j in range(i,Nfinal):
            if (j<N) and (i<N):
                cov_K[i][j]=expo_quad_kernel(xxx[i],xxx[j])
            if (j>(N-1)) and (i<(N)):
                cov_K[i][j]=expo_quad_kernel2(xxx[i],   days[j-N], days[j-N+1])
            if  (i>(N-1)):
                cov_K[i][j]=expo_quad_kernel3(days[i-N],days[i-N+1],  days[j-N],days[j-N+1 ]          )  
            if j!=i:
                cov_K[j][i]=cov_K[i][j]
# Construct ones_vec
ones_vec = np.ones(Nfinal, dtype=np.float64)
ones_vec[N:] = diff
ones_vec = ones_vec.reshape(-1, 1)

jitter =1e-9
    
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
print("Condition number:", np.linalg.cond(cov_K_jittered))
print("dtype:", cov_K.dtype)
print("scalar d:", d.item())

noise_var = 1e-16
KK = cov_K * noise_var + np.eye(Nfinal)
L, lower = cho_factor(KK, lower=True)

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

import math
from scipy import special
from scipy import optimize

EPS = 10e-15


class TruncatedMVN:
    """
    Create a normal distribution :math:`X  \sim N ({\mu}, {\Sigma})` subject to linear inequality constraints
    :math:`lb < X < ub` and sample from it using minimax tilting. Based on the MATLAB implemention by the authors
    (reference below).

    :param np.ndarray mu: (size D) mean of the normal distribution :math:`\mathbf {\mu}`.
    :param np.ndarray cov: (size D x D) covariance of the normal distribution :math:`\mathbf {\Sigma}`.
    :param np.ndarray lb: (size D) lower bound constrain of the multivariate normal distribution :math:`\mathbf lb`.
    :param np.ndarray ub: (size D) upper bound constrain of the multivariate normal distribution :math:`\mathbf ub`.
    :param Union[int, None] seed: a random seed.

    Note that the algorithm may not work if 'cov' is close to being rank deficient.

    Reference:
    Botev, Z. I., (2016), The normal law under linear restrictions: simulation and estimation via minimax tilting,
    Journal of the Royal Statistical Society Series B, 79, issue 1, p. 125-148,

    Example:
        >>> d = 10  # dimensions
        >>>
        >>> # random mu and cov
        >>> mu = np.random.rand(d)
        >>> cov = 0.5 - np.random.rand(d ** 2).reshape((d, d))
        >>> cov = np.triu(cov)
        >>> cov += cov.T - np.diag(cov.diagonal())
        >>> cov = np.dot(cov, cov)
        >>>
        >>> # constraints
        >>> lb = np.zeros_like(mu) - 2
        >>> ub = np.ones_like(mu) * np.inf
        >>>
        >>> # create truncated normal and sample from it
        >>> n_samples = 100000
        >>> samples = TruncatedMVN(mu, cov, lb, ub).sample(n_samples)

    Reimplementation by Paul Brunzema
    """

    def __init__(self, mu, cov, lb, ub, seed=None):
        self.dim = len(mu)
        if not cov.shape[0] == cov.shape[1]:
            raise RuntimeError("Covariance matrix must be of shape DxD!")
        if not (self.dim == cov.shape[0] and self.dim == len(lb) and self.dim == len(ub)):
            raise RuntimeError("Dimensions D of mean (mu), covariance matric (cov), lower bound (lb) "
                               "and upper bound (ub) must be the same!")

        self.cov = cov
        self.orig_mu = mu
        self.orig_lb = lb
        self.orig_ub = ub
        
        # permutated
        self.lb = lb - mu  # move distr./bounds to have zero mean
        self.ub = ub - mu  # move distr./bounds to have zero mean
        if np.any(self.ub <= self.lb):
            raise RuntimeError("Upper bound (ub) must be strictly greater than lower bound (lb) for all D dimensions!")

        # scaled Cholesky with zero diagonal, permutated
        self.L = np.empty_like(cov)
        self.unscaled_L = np.empty_like(cov)

        # placeholder for optimization
        self.perm = None
        self.x = None
        self.mu = None
        self.psistar = None

        # for numerics
        self.eps = EPS

        # a random state
        self.random_state = np.random.RandomState(seed)

    def sample(self, n):
        """
        Create n samples from the truncated normal distribution.

        :param int n: Number of samples to create.
        :return: D x n array with the samples.
        :rtype: np.ndarray
        """
        if not isinstance(n, int):
            raise RuntimeError("Number of samples must be an integer!")

        # factors (Cholesky, etc.) only need to be computed once!
        if self.psistar is None:
            self.compute_factors()

        # start acceptance rejection sampling
        rv = np.array([], dtype=np.float64).reshape(self.dim, 0)
        accept, iteration = 0, 0
        while accept < n:
            logpr, Z = self.mvnrnd(n, self.mu)  # simulate n proposals
            idx = -np.log(self.random_state.rand(n)) > (self.psistar - logpr)  # acceptance tests
            rv = np.concatenate((rv, Z[:, idx]), axis=1)  # accumulate accepted
            accept = rv.shape[1]  # keep track of # of accepted
            iteration += 1
            if iteration == 10 ** 3:
                print('Warning: Acceptance prob. smaller than 0.001.')
            elif iteration > 10 ** 4:
                accept = n
                rv = np.concatenate((rv, Z), axis=1)
                print('Warning: Sample is only approximately distributed.')

        # finish sampling and postprocess the samples!
        order = self.perm.argsort(axis=0)
        rv = rv[:, :n]
        rv = self.unscaled_L @ rv
        rv = rv[order, :]

        # retransfer to original mean
        rv += np.tile(self.orig_mu.reshape(self.dim, 1), (1, rv.shape[-1]))  # Z = X + mu
        return rv
    
    def compute_factors(self):
        # compute permutated Cholesky factor and solve optimization

        # Cholesky decomposition of matrix with permuation
        self.unscaled_L, self.perm = self.colperm()
        D = np.diag(self.unscaled_L)
        if np.any(D < self.eps):
            print('Warning: Method might fail as covariance matrix is singular!')

        # rescale
        scaled_L = self.unscaled_L / np.tile(D.reshape(self.dim, 1), (1, self.dim))
        self.lb = self.lb / D
        self.ub = self.ub / D

        # remove diagonal
        self.L = scaled_L - np.eye(self.dim)

        # get gradient/Jacobian function
        gradpsi = self.get_gradient_function()
        x0 = np.zeros(2 * (self.dim - 1))

        # find optimal tilting parameter non-linear equation solver
        sol = optimize.root(gradpsi, x0, args=(self.L, self.lb, self.ub), method='hybr', jac=True)
        if not sol.success:
            print('Warning: Method may fail as covariance matrix is close to singular!')
        self.x = sol.x[:self.dim - 1]
        self.mu = sol.x[self.dim - 1:]

        # compute psi star
        self.psistar = self.psy(self.x, self.mu)
        
    def reset(self):
        # reset factors -> when sampling, optimization for optimal tilting parameters is performed again

        # permutated
        self.lb = self.orig_lb - self.orig_mu  # move distr./bounds to have zero mean
        self.ub = self.orig_ub - self.orig_mu

        # scaled Cholesky with zero diagonal, permutated
        self.L = np.empty_like(self.cov)
        self.unscaled_L = np.empty_like(self.cov)

        # placeholder for optimization
        self.perm = None
        self.x = None
        self.mu = None
        self.psistar = None

    def mvnrnd(self, n, mu):
        # generates the proposals from the exponentially tilted sequential importance sampling pdf
        # output:     logpr, log-likelihood of sample
        #             Z, random sample
        mu = np.append(mu, [0.])
        Z = np.zeros((self.dim, n))
        logpr = 0
        for k in range(self.dim):
            # compute matrix multiplication L @ Z
            col = self.L[k, :k] @ Z[:k, :]
            # compute limits of truncation
            tl = self.lb[k] - mu[k] - col
            tu = self.ub[k] - mu[k] - col
            # simulate N(mu,1) conditional on [tl,tu]
            Z[k, :] = mu[k] + self.trandn(tl, tu)
            # update likelihood ratio
            logpr += lnNormalProb(tl, tu) + .5 * mu[k] ** 2 - mu[k] * Z[k, :]
        return logpr, Z

    def trandn(self, lb, ub):
        """
        Sample generator for the truncated standard multivariate normal distribution :math:`X \sim N(0,I)` s.t.
        :math:`lb<X<ub`.

        If you wish to simulate a random variable 'Z' from the non-standard Gaussian :math:`N(m,s^2)`
        conditional on :math:`lb<Z<ub`, then first simulate x=TruncatedMVNSampler.trandn((l-m)/s,(u-m)/s) and set
        Z=m+s*x.
        Infinite values for 'ub' and 'lb' are accepted.

        :param np.ndarray lb: (size D) lower bound constrain of the normal distribution :math:`\mathbf lb`.
        :param np.ndarray ub: (size D) upper bound constrain of the normal distribution :math:`\mathbf lb`.

        :return: D samples if the truncated normal distribition x ~ N(0, I) subject to lb < x < ub.
        :rtype: np.ndarray
        """
        if not len(lb) == len(ub):
            raise RuntimeError("Lower bound (lb) and upper bound (ub) must be of the same length!")

        x = np.empty_like(lb)
        a = 0.66  # threshold used in MATLAB implementation
        # three cases to consider
        # case 1: a<lb<ub
        I = lb > a
        if np.any(I):
            tl = lb[I]
            tu = ub[I]
            x[I] = self.ntail(tl, tu)
        # case 2: lb<ub<-a
        J = ub < -a
        if np.any(J):
            tl = -ub[J]
            tu = -lb[J]
            x[J] = - self.ntail(tl, tu)
        # case 3: otherwise use inverse transform or accept-reject
        I = ~(I | J)
        if np.any(I):
            tl = lb[I]
            tu = ub[I]
            x[I] = self.tn(tl, tu)
        return x

    def tn(self, lb, ub, tol=2):
        # samples a column vector of length=len(lb)=len(ub) from the standard multivariate normal distribution
        # truncated over the region [lb,ub], where -a<lb<ub<a for some 'a' and lb and ub are column vectors
        # uses acceptance rejection and inverse-transform method

        sw = tol  # controls switch between methods, threshold can be tuned for maximum speed for each platform
        x = np.empty_like(lb)
        # case 1: abs(ub-lb)>tol, uses accept-reject from randn
        I = abs(ub - lb) > sw
        if np.any(I):
            tl = lb[I]
            tu = ub[I]
            x[I] = self.trnd(tl, tu)

        # case 2: abs(u-l)<tol, uses inverse-transform
        I = ~I
        if np.any(I):
            tl = lb[I]
            tu = ub[I]
            pl = special.erfc(tl / np.sqrt(2)) / 2
            pu = special.erfc(tu / np.sqrt(2)) / 2
            x[I] = np.sqrt(2) * special.erfcinv(2 * (pl - (pl - pu) * self.random_state.rand(len(tl))))
        return x

    def trnd(self, lb, ub):
        # uses acceptance rejection to simulate from truncated normal
        x = self.random_state.randn(len(lb))  # sample normal
        test = (x < lb) | (x > ub)
        I = np.where(test)[0]
        d = len(I)
        while d > 0:  # while there are rejections
            ly = lb[I]
            uy = ub[I]
            y = self.random_state.randn(len(uy))  # resample
            idx = (y > ly) & (y < uy)  # accepted
            x[I[idx]] = y[idx]
            I = I[~idx]
            d = len(I)
        return x

    def ntail(self, lb, ub):
        # samples a column vector of length=len(lb)=len(ub) from the standard multivariate normal distribution
        # truncated over the region [lb,ub], where lb>0 and lb and ub are column vectors
        # uses acceptance-rejection from Rayleigh distr. similar to Marsaglia (1964)
        if not len(lb) == len(ub):
            raise RuntimeError("Lower bound (lb) and upper bound (ub) must be of the same length!")
        c = (lb ** 2) / 2
        n = len(lb)
        f = np.expm1(c - ub ** 2 / 2)
        x = c - np.log(1 + self.random_state.rand(n) * f)  # sample using Rayleigh
        # keep list of rejected
        I = np.where(self.random_state.rand(n) ** 2 * x > c)[0]
        d = len(I)
        while d > 0:  # while there are rejections
            cy = c[I]
            y = cy - np.log(1 + self.random_state.rand(d) * f[I])
            idx = (self.random_state.rand(d) ** 2 * y) < cy  # accepted
            x[I[idx]] = y[idx]  # store the accepted
            I = I[~idx]  # remove accepted from the list
            d = len(I)
        return np.sqrt(2 * x)  # this Rayleigh transform can be delayed till the end

    def psy(self, x, mu):
        # implements psi(x,mu); assumes scaled 'L' without diagonal
        x = np.append(x, [0.])
        mu = np.append(mu, [0.])
        c = self.L @ x
        lt = self.lb - mu - c
        ut = self.ub - mu - c
        p = np.sum(lnNormalProb(lt, ut) + 0.5 * mu ** 2 - x * mu)
        return p

    def get_gradient_function(self):
        # wrapper to avoid dependancy on self

        def gradpsi(y, L, l, u):
            # implements gradient of psi(x) to find optimal exponential twisting, returns also the Jacobian
            # NOTE: assumes scaled 'L' with zero diagonal
            d = len(u)
            c = np.zeros(d)
            mu, x = c.copy(), c.copy()
            x[0:d - 1] = y[0:d - 1]
            mu[0:d - 1] = y[d - 1:]

            # compute now ~l and ~u
            c[1:d] = L[1:d, :] @ x
            lt = l - mu - c
            ut = u - mu - c

            # compute gradients avoiding catastrophic cancellation
            w = lnNormalProb(lt, ut)
            pl = np.exp(-0.5 * lt ** 2 - w) / np.sqrt(2 * math.pi)
            pu = np.exp(-0.5 * ut ** 2 - w) / np.sqrt(2 * math.pi)
            P = pl - pu

            # output the gradient
            dfdx = - mu[0:d - 1] + (P.T @ L[:, 0:d - 1]).T
            dfdm = mu - x + P
            grad = np.concatenate((dfdx, dfdm[:-1]), axis=0)

            # construct jacobian
            lt[np.isinf(lt)] = 0
            ut[np.isinf(ut)] = 0

            dP = - P ** 2 + lt * pl - ut * pu
            DL = np.tile(dP.reshape(d, 1), (1, d)) * L
            mx = DL - np.eye(d)
            xx = L.T @ DL
            mx = mx[:-1, :-1]
            xx = xx[:-1, :-1]
            J = np.block([[xx, mx.T],
                          [mx, np.diag(1 + dP[:-1])]])
            return (grad, J)

        return gradpsi

    def colperm(self):
        perm = np.arange(self.dim)
        L = np.zeros_like(self.cov)
        z = np.zeros_like(self.orig_mu)

        for j in perm.copy():
            pr = np.ones_like(z) * np.inf  # compute marginal prob.
            I = np.arange(j, self.dim)  # search remaining dimensions
            D = np.diag(self.cov)
            s = D[I] - np.sum(L[I, 0:j] ** 2, axis=1)
            s[s < 0] = self.eps
            s = np.sqrt(s)
            tl = (self.lb[I] - L[I, 0:j] @ z[0:j]) / s
            tu = (self.ub[I] - L[I, 0:j] @ z[0:j]) / s
            pr[I] = lnNormalProb(tl, tu)
            # find smallest marginal dimension
            k = np.argmin(pr)

            # flip dimensions k-->j
            jk = [j, k]
            kj = [k, j]
            self.cov[jk, :] = self.cov[kj, :]  # update rows of cov
            self.cov[:, jk] = self.cov[:, kj]  # update cols of cov
            L[jk, :] = L[kj, :]  # update only rows of L
            self.lb[jk] = self.lb[kj]  # update integration limits
            self.ub[jk] = self.ub[kj]  # update integration limits
            perm[jk] = perm[kj]  # keep track of permutation

            # construct L sequentially via Cholesky computation
            s = self.cov[j, j] - np.sum(L[j, 0:j] ** 2, axis=0)
            if s < -0.01:
                raise RuntimeError("Sigma is not positive semi-definite")
            elif s < 0:
                s = self.eps
            L[j, j] = np.sqrt(s)
            new_L = self.cov[j + 1:self.dim, j] - L[j + 1:self.dim, 0:j] @ L[j, 0:j].T
            L[j + 1:self.dim, j] = new_L / L[j, j]

            # find mean value, z(j), of truncated normal
            tl = (self.lb[j] - L[j, 0:j - 1] @ z[0:j - 1]) / L[j, j]
            tu = (self.ub[j] - L[j, 0:j - 1] @ z[0:j - 1]) / L[j, j]
            w = lnNormalProb(tl, tu)  # aids in computing expected value of trunc. normal
            z[j] = (np.exp(-.5 * tl ** 2 - w) - np.exp(-.5 * tu ** 2 - w)) / np.sqrt(2 * math.pi)
        return L, perm


def lnNormalProb(a, b):
    # computes ln(P(a<Z<b)) where Z~N(0,1) very accurately for any 'a', 'b'
    p = np.zeros_like(a)
    # case b>a>0
    I = a > 0
    if np.any(I):
        pa = lnPhi(a[I])
        pb = lnPhi(b[I])
        p[I] = pa + np.log1p(-np.exp(pb - pa))
    # case a<b<0
    idx = b < 0
    if np.any(idx):
        pa = lnPhi(-a[idx])  # log of lower tail
        pb = lnPhi(-b[idx])
        p[idx] = pb + np.log1p(-np.exp(pa - pb))
    # case a < 0 < b
    I = (~I) & (~idx)
    if np.any(I):
        pa = special.erfc(-a[I] / np.sqrt(2)) / 2  # lower tail
        pb = special.erfc(b[I] / np.sqrt(2)) / 2  # upper tail
        p[I] = np.log1p(-pa - pb)
    return p


def lnPhi(x):
    # computes logarithm of  tail of Z~N(0,1) mitigating numerical roundoff errors
    out = -0.5 * x ** 2 - np.log(2) + np.log(special.erfcx(x / np.sqrt(2)) + EPS)  # divide by zeros error -> add eps
    return out




#GP_regression_ess(points_inhomo,g_mk,theta0,theta1,noise_var,rang,num_points,cov_K_noise)
def GP_regression_ess(xi,rang,num_points,cov_chol):
    NN=N+K
   
     #cov is sigma_22
    
    x1=np.linspace(0,rang,num_points+1)[1:]      # prediction points, integer is to make it easy
    M=len(x1)
    mean=np.zeros((1,M))[0]
    posterior_cov=np.zeros((M,M)) #final covariance returned
    k_matrix=np.zeros((M,NN)) #Sigma_12

    k_matrix_pre=np.zeros((M,M)) #Sigma_11
    for i in range(M):
        for j in range(NN):
            if j<(ntip-1):
                k_matrix[i][j]=expo_quad_kernel(x1[i],xi[j])
            else:
                k_matrix[i][j]=expo_quad_kernel2(x1[i],days[j-ntip+1], days[j-ntip+2])
              
    #k_C=cho_solve((cov_chol, True), yi)
    k_C2=cho_solve((cov_chol, True), k_matrix.T)
   
    for i in range(1,M):
        for j in range(i,M):
            k_matrix_pre[i][j]=expo_quad_kernel(x1[i],x1[j])
            k_matrix_pre[j][i]=k_matrix_pre[i][j]
    #posterior_cov=(k_matrix_pre-np.dot(k_matrix,k_C2))/theta+np.eye(M)*noise_var
    posterior_cov=k_matrix_pre-np.dot(k_matrix,k_C2)
    
    
    return x1,k_matrix, posterior_cov


curloglike = None
num_points=100
jitter2=1e-3
x1,k_matrix, posterior_cov=GP_regression_ess(np.array(points_inhomo),T,num_points,cov_K_mod_chol)
for ite in range(nsim1+nsim2):
    cov_K_chol_final=cov_K_mod_chol/np.sqrt(theta)
    prior=chol_sample(mean=np.zeros(Nfinal), cov_chol=cov_K_chol_final)#nu
    g_mk,curloglike=elliptical_slice(g_mk.reshape(1,-1),prior,log_lik,pdf_params=[N],cur_lnpdf=curloglike,angle_range=None)
    if (ite%rem==0):
        g_mk_list3.append(1/g_mk[0][0:N])
        k_C=cho_solve((cov_K_mod_chol, True), g_mk[0])
        mean=np.dot(k_matrix,k_C)
        cov=posterior_cov/theta+np.eye(num_points)*jitter2
        lb =np.zeros_like(mean)
        ub =np.ones_like(mean) * np.inf
        line = TruncatedMVN(mean, cov, lb, ub).sample(1)
        g_mk_list2.append(1/line)
        alpha_pos=alpha+Nfinal/2
        u = cho_solve((cov_K_mod_chol, True), g_mk[0])
        beta_pos=beta+1/2*g_mk[0]@u
        theta=np.random.gamma(alpha_pos, 1/beta_pos,1)
        theta_list.append(theta)


np.savez('/scratch/groups/juliapr/output_Bingjing/boundedcoaldata/real_100_BM_jitter9.npz', aaa=g_mk_list3,aa=g_mk_list2,c=points_inhomo,d=x1,i=jitter,s=theta_list)    