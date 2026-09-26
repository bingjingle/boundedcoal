#!/usr/bin/env Rscript
# Discretised bounded / standard coalescent via JuliaPalacios/phylodyn.
# The precision prior is explicit. Both equal-tailed and minimal-width 95%
# intervals are saved; the submitted paper uses equal-tailed intervals.
# This historical driver preserves its original iteration settings.
suppressMessages({library(ape); library(phylodyn)})
a <- commandArgs(trailingOnly = TRUE)
which <- a[1]; prec <- as.numeric(a[2]); tag <- a[3]
csv <- a[4]; ntip <- as.integer(a[5]); tau <- as.numeric(a[6]); outdir <- a[7]
NS <- 200000; NB <- 100000; Ng <- 100
# BC_SMOKE=1 shrinks the chain so the pipeline can be tested end-to-end quickly.
if (nzchar(Sys.getenv("BC_SMOKE"))) { NS <- 400; NB <- 200 }

ct <- sort(read.csv(csv)$coal_times)
stopifnot(length(ct) == ntip - 1)
data <- list(coal_times = ct, lineages = ntip:2,
             intercoal_times = c(ct[1], ct[-1] - ct[-length(ct)]),
             samp_times = 0, n_sampled = ntip)
x <- seq(0, tau, length.out = Ng + 1)[-1]

# Seed is overridable.  phylodyn's bounded sampler can break down numerically on a
# particular chain realisation -- `bound_prob` goes NA inside coal_loglik_bounded and
# the run aborts -- and which datasets that hits depends on the precision prior:
# at seed 123, syn1 fails at prec 0.1 and syn4 fails at prec 0.01.  Re-drawing the
# chain with another seed avoids that realisation; it does not fix the underlying
# fragility, which is worth reporting upstream.
set.seed(as.integer(Sys.getenv("BC_SEED", "123")))
zz <- file("/dev/null", "w"); sink(zz)
t0 <- Sys.time()
r <- if (which == "bc")
  phylodyn:::mcmc_sampling(data, alg = "bound_ESS", nsamp = NS, nburnin = NB, ngrid = Ng,
                           bound = tau, printevery = 1e9,
                           prec_alpha = prec, prec_beta = prec) else
  phylodyn:::mcmc_sampling(data, alg = "ESS", nsamp = NS, nburnin = NB, ngrid = Ng,
                           printevery = 1e9, prec_alpha = prec, prec_beta = prec)
el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
sink(); close(zz)

# ---- credible bands, both conventions -------------------------------------
# Ne draws on the ngrid-1 interior cells; columns follow phylodyn's logfmat.
NeS <- exp(r$samp[, 1:(Ng - 1), drop = FALSE])

hpd <- function(v, alpha = 0.05) {             # minimal-width interval
  v <- sort(v); n <- length(v); k <- floor((1 - alpha) * n)
  j <- which.min(v[(k + 1):n] - v[1:(n - k)])
  c(v[j], v[j + k])
}
hp  <- apply(NeS, 2, hpd)
stepf <- function(z) stats::stepfun(r$grid, c(z[1], z, utils::tail(z, 1)))
hpd_lo_fun <- stepf(hp[1, ]); hpd_hi_fun <- stepf(hp[2, ])

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
out <- data.frame(x = x,
                  med    = r$med_fun(x),
                  low    = r$low_fun(x),    # equal-tailed 2.5%  (phylodyn default)
                  hi     = r$hi_fun(x),     # equal-tailed 97.5%
                  hpd_lo = hpd_lo_fun(x),   # minimal-width 95%
                  hpd_hi = hpd_hi_fun(x))
write.csv(out, file.path(outdir, paste0(tag, ".csv")), row.names = FALSE)
write.csv(data.frame(mid = exp(r$samp[, 50])),
          file.path(outdir, paste0(tag, "_trace.csv")), row.names = FALSE)
saveRDS(list(grid = r$grid, x = x, NeS = NeS, secs = el, prec = prec,
             nsamp = NS, nburnin = NB, ngrid = Ng, tau = tau, ntip = ntip,
             which = which, coal_times = ct),
        file.path(outdir, paste0(tag, "_chain.rds")), compress = "xz")
cat("DONE", tag, round(el), "s  prec", prec, " per10k", round(el / (NS + NB) * 10000, 2),
    " med[1]/[50]/[100]:", signif(out$med[1], 4), signif(out$med[50], 4),
    signif(out$med[100], 4), " hi_max:", signif(max(out$hi), 4), "\n")
