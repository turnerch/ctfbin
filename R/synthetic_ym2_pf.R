#' Synthetic YM2 particle filter
#'
#' This script defines a particle filter for use in the synthetic_ym2() function, which was designed
#' by John W. Smith Jr.  It embeds a DALEC2 process model for carbon evolution.
#' 
#' @author Turner Haugen
#'
#' @param seed changes randomization seed, will make particle filter give different results (is 3141 by default)
#'
#' @return A Csnippet of the particle filter
#'
#' @examples
#' filter <- synthetic_ym2_pf()
#' filter2 <- synthetic_ym2_pf(seed = 16766)
#'
#' @import pomp
#' @importFrom lamW lambertW0
#'
#' @export
 
synthetic_ym2_pf <- function(seed = 3141) {
  ## seed servo
  
  if(!is.numeric(seed)) stop("seed parameter must be numeric!")
  
  set.seed(seed)
  
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
  
  ## set variance governing parameters
  omega_f <- omega_lab <-  .02
  omega_f_obs <- .2
  omega_lab_obs <- .09
  
  ## fixed parameters
  lma <- 95
  s <- 365.25 / pi
  
  phi_f <- rep(NA, length(Cfol))
  
  for (i in 2:length(Cfol)){
    ## simulate gross photosynthetic production for day i
    G <- loglai::GPP(lai = Cfol[i-1] / lma, c_eff = c_eff, maxt = synth_drivers$MaxT[i],
             mint = synth_drivers$MinT[i], Ca = synth_drivers$ca[i], lat = 46.23,
             yearday = synth_drivers$yearday[i], rad = synth_drivers$rad[i])
    
    ## psi_f computation using lamW package
    psi_f <- - c_rf * sqrt(lambertW0(1 / (2*pi*(log(1-c_lf)^2)))) / sqrt(2)
    
    ## get phi_f and phi_o values
    phi_f[i] <- sqrt(2 / pi) * ((-log(1 - c_lf))/c_rf) 
    phi_f[i] <- phi_f[i] * exp(- (sin((synth_drivers$yearday[i] - d_f + psi_f)/s) * sqrt(2)*s / c_rf)^2)
    
    phi_o <- sqrt(2 / pi) * (6.9088 / c_ro)
    phi_o <- phi_o * exp(- (sin((synth_drivers$yearday[i] - d_o - c_ro * .6245)/s) * sqrt(2)*s / c_ro)^2)
    
    ## deterministic process evolution DALEC2 equations
    cf_evo <- (1 - phi_f[i]) * Cfol[i-1] + f_f * G + phi_o * Clab[i-1]
    clab_evo <- (1 - phi_o) * Clab[i-1] + f_lab * G
    
    ## mean function for lognormal embeddings
    cf_mean <- log(cf_evo) - .5*log(1 + omega_f^2)
    clab_mean <- log(clab_evo) - .5*log(1 + omega_lab^2)
    
    ## sd function for lognormal embeddings
    cf_sd <- sqrt(log(1 + omega_f^2))
    clab_sd <- sqrt(log(1 + omega_lab^2))
    
    ## simulate state evolution
    Cfol[i] <-  rlnorm(1, meanlog = cf_mean, sdlog = cf_sd)
    Clab[i] <- rlnorm(1, meanlog = clab_mean, sdlog = clab_sd)
  }
  
  ## generate observed data
  LAI_obs <- rlnorm(length(Cfol), log(Cfol) - .5*log(1 + omega_f_obs^2), sqrt(log(1 + omega_f_obs^2))) / lma 
  Clab_obs <- rlnorm(length(Clab), log(Clab) - .5*log(1 + omega_lab_obs^2), sqrt(log(1 + omega_lab_obs^2)))
  
  ## measurement component for particle filter, written 
  ## using Csnippet functionality from pomp package
  dmeas_D2_mm <- loglai::dmeas_mm()
  
  Clab_obs[!(1:880 %% 30 == 1)] <- NA
  
  ## create some synthetic data
  LAI_obs[seq(from = 1, to = 882, by = 4)] <- NA
  
  ## coerce to data frame
  pomp_df = data.frame(LAI = LAI_obs, 
                       Clab_obs = Clab_obs)
  ## no leaf-fall for now
  pomp_df$Lf = NA
  
  ## drivers in format accepted by pomp
  Driver_df = as.data.frame(cbind(time = 0:880, 
                                  MaxT = synth_drivers$MaxT, 
                                  MinT = synth_drivers$MinT,
                                  yearday = synth_drivers$yearday, 
                                  ca = synth_drivers$ca,
                                  rad = synth_drivers$rad))
  
  ## prior evaluation. not necessary for BO per se but useful
  ## to define them now to get other algo benchmarks a la pMCMC
  dprior_D2 = Csnippet("
                      double lik3;
                      double lik4; 
                      double lik5;
                      double lik6;
                      double lik7; 
                      double lik8;
                      double lik9;
                      double lik10; 
                      double lik11;
                      double lik12;
                      lik3 = dunif(d_onset, 1.0, 365.0, give_log);
                      lik4 = dunif(d_fall, 1.0, 365.0, give_log);
                      lik5 = dunif(c_eff, 10.0, 100.0, give_log);
                      lik6 = dunif(c_lf, .125, 1.0, give_log);
                      lik7 = dunif(c_ronset, 10.0, 100.0, give_log);
                      lik8 = dunif(c_rfall, 20.0, 150.0, give_log);
                      lik9 = dunif(omega_lab, 0.0, 1.0, give_log);
                      lik10 = dunif(omega_f, 0.0, 1.0, give_log);
                      lik11 = dunif(f_lab, .01, .5, give_log);
                      lik12 = dunif(f_f, .01, .5, give_log);
                      lik = (give_log) ? lik3 + lik4 + lik5 + lik6 + lik7 + lik8 + lik9 + lik10 + lik11 + lik12 : lik3 * lik4 * lik5 * lik6 * lik7 * lik8 * lik9 * lik10 * lik11 * lik12;
                     ")
  
  ## transform function to map bounded interval onto real line
  ## useful
  Csnip_transform = Csnippet("
  T_d_onset = logit((d_onset - 1.0) / (365.0 - 1.0));
  T_d_fall = logit((d_fall - 1.0) / (365.0 - 1.0));
  T_c_eff = logit((c_eff - 10.0) / (100.0 - 10.0));
  T_c_ronset = logit((c_ronset - 10.0) / (100.0 - 10.0));
  T_c_rfall = logit((c_rfall - 20.0) / (150.0 - 20.0));
  T_c_lf = logit((c_lf - .125) / (.875));
  T_omega_lab = log(omega_lab);
  T_omega_f = log(omega_f);
  T_f_f = logit((f_f - .01)/.49);
  T_f_lab = logit((f_lab - .01)/.49);
  ")
  
  ## inverse transform function to map real line back into bounded interval
  Csnip_inv = Csnippet("
  d_onset = 364 * expit(T_d_onset) + 1.0;
  d_fall = 364 * expit(T_d_fall) + 1.0;
  c_eff = 90 * expit(T_c_eff) + 10.0;
  c_ronset = 90 * expit(T_c_ronset) + 10.0;
  c_rfall = 130 * expit(T_c_rfall) + 20.0;
  c_lf = .875 * expit(T_c_lf) + .125;
  omega_lab = exp(T_omega_lab);
  omega_f = exp(T_omega_f);
  f_f = (.49) * expit(T_f_f) + 0.1;
  f_lab = (.49) * expit(T_f_lab) + 0.1;
  ")
  
  ## particle filter function for synthetic data, to be used later
  ## in objective function. embeds DALEC2 process model for carbon
  ## evolution
  
  filter_sym <- pomp(data = pomp_df, times = 0:880, t0 = 0,
                     rprocess = discrete_time(loglai::synthetic_ym2_evo(), delta.t = 1),
                     statenames = c('Clab', 'Cfol'),
                     covar = covariate_table(Driver_df, times = 0:880),
                     tcovar = 'time',
                     covarnames = c('MaxT', 'MinT', 'yearday', 'ca', 'rad'),
                     rinit = function(Clab_0, LAI_0, lma, ...){
                       c(Clab = rnorm(1, Clab_0, .01*Clab_0),
                         Cfol = rnorm(1, LAI_0*lma, .01*LAI_0*lma))
                     },
                     dmeasure = dmeas_D2_mm,
                     dprior = dprior_D2,
                     params = c(d_fall = d_f,
                                d_onset = d_o,
                                c_lf = c_lf,
                                c_rfall = c_rf,
                                c_ronset = c_ro,
                                c_eff = c_eff,
                                f_f = f_f,
                                f_lab = f_lab,
                                omega_f = omega_f,
                                omega_lab =  omega_lab,
                                LAI_0 = Cfol[1] / lma,
                                Clab_0 = Clab[1],
                                lma = lma,
                                omega_obs_f = omega_f_obs),
                     partrans = parameter_trans(toEst = Csnip_transform, fromEst = Csnip_inv),
                     paramnames = c('d_fall', 'd_onset', 'c_lf', 'c_rfall',
                                    'c_ronset', 'c_eff', 'f_f', 'f_lab', 'omega_f', 'omega_lab', 
                                    'LAI_0', 'Clab_0', 'lma', 'omega_obs_f'))
  
  ## Returns pomp object filter
  return(filter_sym)
}
