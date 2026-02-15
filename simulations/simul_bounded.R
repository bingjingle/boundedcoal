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






