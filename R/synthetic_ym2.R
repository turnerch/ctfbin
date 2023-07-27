#' Synthetic DALEC2 ecosystem simulation
#'
#' This script is designed to generate two years of synthetic dataset from a
#' state space formulation of the DALEC2 ecosystem model. Meteorological driver data
#' from the University of Notre Dame Environmental Research Center East (UNDE) were obtained
#' from NEON and are used during the generation of the sythetic data
#' 
#' @author Turner Haugen
#'
#' @param filter object of class pomp, particle filter model used to generate estimate of marginal log-likelihood
#' @param pars named list of parameter values
#' @param n.part integer, number of particles (will be 200 by default)
#' @param scaled boolean, true if parameters are scaled to a [0,1]^10 hypercube (will be true by default)
#' @param verbose boolean, true if verbose warning about setting scaled to default value is needed
#'
#' @return obj: list, contains estimate of marginal log-likelihood obtained from particle filter (obj$logLik) and trajectory (obj$traj)
#'
#' @examples
#' synthetic_ym2(synthetic_ym2_pf(), synthetic_ym2_pars(), verbose = FALSE)
#' synthetic_ym2(synthetic_ym2_pf(seed = 13), synthetic_ym2_pars(FALSE), n.part = 150, scaled = FALSE)
#'
#' @import pomp
#'
#' @export

synthetic_ym2 <- function(filter, pars, n.part = 200, scaled = TRUE, verbose = TRUE) {
  
  # Checks all parameters provided to make sure they are the correct types

  if(!("pomp" %in% class(filter))) stop("Custom filter must be a pomp object")

  if(!is.vector(pars)) stop("Custom parameters must be a vector")
  
  if(!is.numeric(n.part)) stop("n.part parameter must be numeric!")
  
  if(!is.logical(scaled)) stop("scaled parameter must be logical!")
  
  if(!is.logical(verbose)) stop("verbose parameter must be logical!")
  
  # Verbose warning can be turned off by switching verbose parameter to FALSE
  
  if(verbose & scaled) print("Parameters currently being passed as scaled values")

  ## list of true parameter values, scaled to live on a [0,1]^10 
  ## hypercube
  
  ## set values for constant parameters
  pars['LAI_0'] = 2.78
  pars['lma'] = 95
  pars['Clab_0'] = 129.88
  pars['omega_obs_f'] = .2
  
  ## convert parameters from [0,1] to range expected by model if scaled is true
  if (scaled){
    pars['d_fall'] = pars['d_fall']*364 + 1
    pars['d_onset'] = pars['d_onset']*364 + 1
    pars['c_lf'] = pars['c_lf'] * .875 + .125
    pars['c_rfall'] = pars['c_rfall'] * 130 + 20
    pars['c_ronset'] = pars['c_ronset'] * 90 + 10
    pars['c_eff'] = pars['c_eff'] * 90 + 10
    pars['f_f'] = pars['f_f'] * .49 + .01
    pars['f_lab'] = pars['f_lab'] * .49 + .01
  }
  
  ## note that omega_f and omega_lab live on [0,1] natively
  ## so no conversion is required
  
  ## objective function evaluation from particle filter. objective
  ## function is a stochastic approximation of the marginal log-likelihood
  ## of a non-linear non-gaussian state space model
  obj <- pfilter(filter, Np = n.part, filter.traj = TRUE, 
                 params = pars[c('d_fall', 'd_onset', 'c_lf', 'c_rfall',
                                                              'c_ronset', 'c_eff', 'f_f', 'f_lab', 'omega_f', 'omega_lab', 
                                                             'LAI_0', 'Clab_0', 'lma', 'omega_obs_f')])
  ## return objective function
  return(list(logLik = logLik(obj), traj = filter_traj(obj)))
}