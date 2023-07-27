#' YM2 default parameters
#'
#' This script defines a default set of parameters for ym2.  Since this uses real data, these are not the TRUE parameters!  They provide the highest marginal log-likelihood of all 50000 particle MCMC iterations.
#' 
#' @author Turner Haugen
#'
#' @param N/A
#'
#' @return A list of parameter values
#'
#' @examples
#' ym2_pars()
#'
#' @import pomp
#'
#' @export

ym2_pars <- function() {
  return(c(
    Clab_0 = 199.88826735,
    lma = 75,
    LAI_0 = 3.788,
    omega_obs_f = 0.05,
    d_fall = 273.31828977,
    c_rfall = 48.98778074,
    c_lf = 0.92769078,
    d_onset = 95.14796422,
    c_ronset = 44.59112732,
    c_eff = 81.24987320,
    omega_f = 0.01445637,
    omega_lab = 0.054312313,
    f_f = 0.01353197,
    f_lab = 0.45956056))
}