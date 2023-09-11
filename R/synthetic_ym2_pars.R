#' Synthetic YM2 default parameters
#'
#' This script defines the parameters for synthetic_ym2, scaled to live on a [0,1]^10 hypercube.
#' 
#' @author Turner Haugen
#'
#' @param scaled logical that decides if parameters are scaled (True) or unscaled (False)
#'
#' @return A list of parameter values
#'
#' @examples
#' synthetic_ym2_pars()
#' synthetic_ym2_pars(F)
#'
#' @import pomp
#'
#' @export

synthetic_ym2_pars <- function(scaled = TRUE) {
  ## obtaining default parameters uses similar method as the beginning of ym2_pf()
  
  ## check if value for variable scaled is valid
  if(!is.logical(scaled)) stop("scaled parameter must be logical!")
  
  ## use unde drivers
  drivers_unde <- loglai::drivers_unde
  
  ## take subset
  synth_drivers <- drivers_unde[1074:1954,]
  
  ## initialize state vectors
  Cfol <- Clab <- rep(NA, nrow(synth_drivers))
  
  Cfol[1] <- 264.1
  Clab[1] <- 129.88
  
  ## set process parameters
  ## note: process parameters chosen to give a FOREST, and not a fucking PARKING LOT
  
  ## day leaves start fallong
  d_f <- 275
  ## day leaves start growing again
  d_o <- 110
  ## total proportion of leaves that fall
  c_lf <- .89
  ## number of days of leaf fall period
  c_rf <- 60
  ## number of days of labile regrowth
  c_ro <- 45
  ## measure of nitrogen use efficiency of forest canopy
  c_eff <- 20
  ## proportion of GPP allocated to foliage
  f_f <- .03
  ## proportion of GPP allocated to labile carbon
  f_lab <- .12
  
  ## set variance governing parameter
  omega_f <- omega_lab <-  .02
  
  # If scaled values are needed, returns scaled true parameter values
  # Otherwise returns unscaled values
  if(scaled) {
    ## list of true parameter values, scaled to live on a [0,1]^10 
    ## hypercube
    
    return(c(d_fall = (d_f- 1) / 364,
             d_onset = (d_o - 1) / 364,
             c_lf = (c_lf - .125) / .875,
             c_rfall = (c_rf - 20) / 130,
             c_ronset = (c_ro - 10) / 90,
             c_eff = (c_eff - 10) / 90,
             f_f = (f_f - .01) / .49,
             f_lab = (f_lab - .01) / .49,
             omega_f = omega_f,
             omega_lab = omega_lab))
  }
  else(
    ## list of true parameter values, unscaled
    
    return(c(d_fall = d_f,
             d_onset = d_o,
             c_lf = c_lf,
             c_rfall = c_rf,
             c_ronset = c_ro,
             c_eff = c_eff,
             f_f = f_f,
             f_lab = f_lab,
             omega_f = omega_f,
             omega_lab = omega_lab))
  )
}