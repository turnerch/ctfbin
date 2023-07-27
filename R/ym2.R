#' DALEC2 Ecosystem With Real Data
#'
#' todo
#' 
#' @author Turner Haugen
#'
#' @param filter object of class pomp, contains particle filter used to model data
#' @param pars list of parameter values
#' @param n.part number of particles used in particle filter
#'
#' @return obj: list, contains estimate of marginal log-likelihood obtained from particle filter (obj$logLik) and trajectory (obj$traj)
#'
#' @examples todo
#'
#' @import pomp
#'
#' @export

ym2 <- function(filter, pars, n.part = 200) {
  
  if(!("pomp" %in% class(filter))) stop("Custom filter must be a pomp object")
  
  if(!is.vector(pars)) stop("Custom parameters must be a vector")
  
  if(!is.numeric(n.part)) stop("n.part parameter must be numeric!")

  pars['LAI_0'] = 3.788
  pars['lma'] = 75
  pars['Clab_0'] = 199.8883
  pars['omega_obs_f'] = .25
  obj <- pfilter(filter, Np = n.part, filter.traj = TRUE, params = pars[c('Clab_0', 'lma', 'LAI_0', 'omega_obs_f',
                                                             'd_fall', 'c_rfall', 'c_lf', 'd_onset', 'c_ronset', 'c_eff',
                                                             'omega_f', 'omega_lab', 'f_f', 'f_lab')])
  
  return(list(logLik = logLik(obj), traj = filter_traj(obj)))
}