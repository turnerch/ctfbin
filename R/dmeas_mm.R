#' Dmeas
#'
#' measurement component for particle filter, written using Csnippet functionality from pomp package
#' 
#' @author John W. Smith Jr
#'
#' @param N/A
#'
#' @return Csnippet of measurements for particle filter
#'
#' @examples
#' dmeas_mm(270.1)
#' 
#' @import pomp
#'
#' @export

dmeas_mm <- function() {
dmeas_D2_mm <- Csnippet("
                  double lik1;
                  double lik2;
                  double lik3; 
                  if (ISNA(LAI)) {
                    lik1 = (give_log) ? 0 : 1;
                  } else if (!ISNA(LAI)){
                    double mean_fol_l = log(Cfol) - .5*log(1 + pow(.2,2));
                    double sd_fol_l = pow(log(1 + pow(.2,2)), .5);
                    lik1 = dlnorm(lma*LAI, mean_fol_l, sd_fol_l, give_log);
                  }
                  if (ISNA(Lf)) {
                    lik2 = (give_log) ? 0 : 1;
                  } else if (!ISNA(Lf)){
                    double pi = atan(1) * 4;
                    double s = 365.25/pi;
                    double psi_f = -1.3588480 + 4.5994549 * c_lf - 9.5964931 * pow(c_lf, 2.0) + 12.1567793 * pow(c_lf, 3.0) -6.8903864 * pow(c_lf, 4.0) -0.3576296 * pow(c_lf, 5.0) + 1.3941622 * pow(c_lf, 6.0);
                    double fall = ((pow(2/pi, .5) * (-log(1 - c_lf)/c_rfall) * exp(-pow((sin((yearday - d_fall + psi_f)/s) * pow(2, .5) * s / c_rfall), 2))));
                    lik2 = dnorm(Lf, Cfol*fall, (5 + 1*Cfol*fall), give_log);
                  }
                  if (ISNA(Clab_obs)) {
                    lik1 = (give_log) ? 0 : 1;
                  } else if (!ISNA(Clab_obs)){
                    double mean_lab_l = log(Clab) - .5*log(1 + pow(.09,2));
                    double sd_lab_l = pow(log(1 + pow(.09,2)), .5);
                    lik3 = dlnorm(Clab_obs, mean_lab_l, sd_lab_l, give_log);
                  }
                  lik = (give_log) ? lik1 + lik2 + lik3 : lik1 * lik2 * lik3 ;
                ")
  return(dmeas_D2_mm)
}