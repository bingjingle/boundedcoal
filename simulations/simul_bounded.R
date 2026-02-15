##Simulation from bounded coalescent

#' Simulate from inhomogeneous, heterochronous bounded coalescent
#' 
#' @param samp_times numeric vector of sampling times.
#' @param n_sampled numeric vector of samples taken per sampling time.
#' @param traj function that returns effective population size at time t.
#' @param bound upper bound on the TMRCA
#'   
#' @return A list containing vectors of coalescent times \code{coal_times}, 
#'   intercoalescent times \code{intercoal_times}, and number of active lineages
#'   \code{lineages}, as well as passing along \code{samp_times} and
#'   \code{n_sampled}.
#' @export
#' 
#' @examples
#' coalsim_bounded(0:2, 3:1, unif_traj, bound=1)
coalsim_bounded <- function(samp_times, n_sampled, traj, bound, ...)
{
  #This function is not functional for heterochronous sampling
  #the r_j factors need to be updated
  ##hazard and inv. hazard for dominating intensity
  if (stats::is.stepfun(traj)) {
    knots = knots(traj)
    midpts = c(min(knots) - 1, knots[-1] - diff(knots)/2, 
               max(knots) + 1)
    traj_inv <- stats::stepfun(x = knots, y = 1/traj(midpts))
    hazard_simple<-function(t,start,target) integrate_step_fun(traj_inv,start,start+t) - target
    hazard_simple_bound<-hazard_simple(bound,0,0)
    #dom_rate<-function(t) traj_inv(t)*(1/(1-exp(hazard_simple(t,0,0)-hazard_simple_bound)))
    #hazard_dom <- function(t, lins, start, target) 0.5 * lins * 
    #  (lins - 1) * integrate_step_fun(dom_rate, start, 
    #                                  start + t) - target
    is_stepfun = TRUE
  }else {
    traj_inv <- function(t) 1/traj(t, ...)
  #  traj_inv <- function(t) 1/traj(t)
    hazard_simple<-function(t,start,target) stats::integrate(traj_inv,start,start+t)$value - target
    hazard_simple_bound<-hazard_simple(bound,0,0)
    #dom_rate<-function(t) traj_inv(t)*(1/(1-exp(hazard_simple(t,0,0)-hazard_simple_bound)))
    #hazard_dom <- function(t, lins, start, target) 0.5 * lins * 
    #  (lins - 1) * stats::integrate(dom_rate, start, start + t)$value - target
    is_stepfun = FALSE
  }
  val_upper<-2*traj_inv(bound)
  
  ##hazard target
  r_func <- function(k,j) {
    if(j == 1) return(1)
    prod <- 1
    for(m in 1:(j-1)) {
      prod <- prod * ((2*m + 1)/(2*m - 1)) * ((k - m)/(k + m))
    }
    return((-1)^(j-1) * prod)
  }
  ntip<-sum(n_sampled)
  result_list <- lapply(seq_len(ntip), function(k) {
    sapply(seq_len(k), function(i) r_func(k, i))
  })
  com_vec <- choose(1:ntip, 2)
  
  coal_times = NULL
  lineages = NULL
  curr = 1
  active_lineages = n_sampled[curr]
  time = samp_times[curr]
  while (time <= max(samp_times) || active_lineages > 1) {
    if (active_lineages == 1) {
      curr <- curr + 1
      active_lineages <- active_lineages + n_sampled[curr]
      time <- samp_times[curr]
    }
    w <- stats::rexp(1) / (0.5*active_lineages * (active_lineages-1))
    target<-hazard_simple_bound-log(1+exp(-w)*(exp(hazard_simple_bound-hazard_simple(time,0,0))-1))
    if (is_stepfun) {
      y <- hazard_uniroot_stepfun(traj_inv_stepfun = hazard_simple, 
                                   start = 0, target = target)
    }else {
      y <- stats::uniroot(hazard_simple, 
                          start = 0, target = target, lower=0, upper = val_upper, 
                          extendInt = "upX")$root
    }
    while (curr < length(samp_times) && y >= samp_times[curr + 
                                                               1]) {
     # target <- -hazard_dom(t = samp_times[curr + 1] - time, 
    #                    lins = active_lineages, start = time, target = target)
      curr <- curr + 1
      active_lineages <- active_lineages + n_sampled[curr]
      time <- samp_times[curr]
      w <- stats::rexp(1) / (0.5*active_lineages * (active_lineages-1))
      target<-hazard_simple_bound-log(1+exp(-w)*(exp(hazard_simple_bound-hazard_simple(time,0,0))-1))
      if (is_stepfun) {
        y <- hazard_uniroot_stepfun(traj_inv_stepfun = hazard_simple, 
                                    start = 0, target = target)
      }else {
        y <- stats::uniroot(hazard_simple, 
                            start = 0, target = target, lower=0, upper = val_upper, 
                            extendInt = "upX")$root
      }
    }
    time <- y
    ##now thinning here
    con<-hazard_simple(time,0,0)-hazard_simple_bound
    denom=sum(result_list[[active_lineages]]*exp(com_vec[1:(active_lineages)]*con))
    numer=sum(result_list[[active_lineages-1]]*exp(com_vec[1:(active_lineages-1)]*con))
    if(stats::runif(1) <= (numer/denom)*(1-exp(con))) {
      coal_times = c(coal_times, time)
      lineages = c(lineages, active_lineages)
      active_lineages = active_lineages - 1
    }
    
    
  }
  return(list(coal_times = coal_times, lineages = lineages, 
              intercoal_times = c(coal_times[1], diff(coal_times)), 
              samp_times = samp_times, n_sampled = n_sampled))
  
}


####MLE estimation of piece-wise constant under bounded coal

  
# Returns an objective (fn) and gradient (gr) suitable for optim()
bound_coal_loglik <- function(init) {
  
  # precompute terms that do not depend on par
  ntip <- sum(init$ns)
  
  r_func <- function(k, j) {
    if (j == 1) return(1)
    prod <- 1
    for (m in 1:(j - 1)) {
      prod <- prod * ((2*m + 1)/(2*m - 1)) * ((k - m)/(k + m))
    }
    (-1)^(j - 1) * prod
  }
  
  r_ntip <- sapply(seq_len(ntip), function(i) r_func(ntip, i))
  com_vec <- choose(seq_len(ntip), 2)
  
  fn <- function(par) {
    if (length(par) != init$ng)
      stop(sprintf("Incorrect length for par; should be %d", init$ng))
    
    f <- rep(par, init$gridrep)
    
    llnocoal  <- init$D * init$C * exp(-f)
    sllnocoal <- init$D * exp(-f)
    
    Lambda <- sum(sllnocoal)
    bound_prob <- sum(r_ntip * exp(-com_vec * Lambda))
    
    # your original scalar: ll = sum(-y*f - llnocoal - log(bound_prob))
    ll_vec <- -init$y * f - llnocoal - log(bound_prob)
    ll <- sum(ll_vec[!is.nan(ll_vec)])
    
    # optim minimizes, so return negative log-likelihood (nll)
    -ll
  }
  
  gr <- function(par) {
    if (length(par) != init$ng)
      stop(sprintf("Incorrect length for par; should be %d", init$ng))
    
    f <- rep(par, init$gridrep)
    
    llnocoal  <- init$D * init$C * exp(-f)
    sllnocoal <- init$D * exp(-f)
    
    Lambda <- sum(sllnocoal)
    bound_prob <- sum(r_ntip * exp(-com_vec * Lambda))
    grad_bound <- sum(r_ntip * com_vec * exp(-com_vec * Lambda))
    
    # this is your original dll (gradient of ll wrt par)
    dll <- apply(init$rep_idx, 1, function(idx) {
      sum(-init$y[idx[1]:idx[2]] + llnocoal[idx[1]:idx[2]])
    }) - (grad_bound / bound_prob) * apply(init$rep_idx, 1, function(idx) {
      sum(sllnocoal[idx[1]:idx[2]])
    })
    
    # gradient of nll = -ll is -dll
    -as.numeric(dll)
  }
  
  list(fn = fn, gr = gr)
}


#' MLE estimation of Ne from bounded coalescent.
#' 
#' @param data \code{phylo} object or list containing vectors of coalescent 
#'   times \code{coal_times}, sampling times \code{samp_times}, and number 
#'   sampled per sampling time \code{n_sampled}.
#' @param lengthout numeric specifying number of grid points.
#' @param bound proivde uper bound on TMRCA
#'   
#' @return Phylodynamic reconstruction of effective population size at grid 
#'   points. 
#' 
#' @examples
#' data<-coalsim_bounded(c(0),c(10),constant,bound=0.5)
#' bounded_skyline(data, lengthout=10, bound=0.5)
bounded_skyline <- function(data, lengthout = 5, bound = 1) {
  if (inherits(data, "phylo")) {
    phy <- summarize_phylo(data)
  } else if (all(c("coal_times", "samp_times", "n_sampled") %in% names(data))) {
    phy <- with(data, list(
      samp_times = samp_times,
      coal_times = coal_times,
      n_sampled  = n_sampled
    ))
  } else stop("data must be a phylo or a list/data.frame with coal_times, samp_times, n_sampled")
  
  grid <- seq(0, bound + 1e-4, length.out = lengthout)
  
  lik_init <- phylodyn:::coal_lik_init(
    samp_times = phy$samp_times,
    n_sampled  = phy$n_sampled,
    coal_times = phy$coal_times,
    grid       = grid
  )
  
  obj <- bound_coal_loglik(lik_init)
  par0 <- rep(0, lik_init$ng)
  
  fit <- optim(
    par     = par0,
    fn      = obj$fn,
    gr      = obj$gr,
    method  = "BFGS",
    control = list(maxit = 1000, reltol = 1e-10)
  )
  
  Ne=exp(fit$par)
  return(list(Ne=Ne,grid=grid))
}

#Example

exp_traj = function(t, scale=1000, rate=1)
{
  return(scale * exp(-t*rate))
}

set.seed(123)
nsamp<-20
bound<-0.5
data<-coalsim_bounded(c(0),c(nsamp),exp_traj,bound=bound,scale=25,rate=5)

res1<-bounded_skyline(data, lengthout = 7, bound =bound)
plot(res1$grid,c(res1$Ne[1],res1$Ne),type="S",ylab="Ne(t)",xlab="Time",ylim=c(0,25))
points(seq(0,max(res1$grid),by=.05),exp_traj(seq(0,max(res1$grid),by=.05),scale=25,rate=5),type="l",col="red")

tree<-generate_newick(data)$newick
skylineplot(tree)
points(seq(0,max(res1$grid),by=.05),exp_traj(seq(0,max(res1$grid),by=.05),scale=25,rate=5),type="l",col="red")

bnpr<-BNPR(tree)
plot_BNPR(bnpr)
points(seq(0,max(res1$grid),by=.05),exp_traj(seq(0,max(res1$grid),by=.05),scale=25,rate=5),type="l",col="red")


