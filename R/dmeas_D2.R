#' todo
#'
#' todo
#' 
#' @author Turner Haugen
#'
#' @param todo
#'
#' @return Csnippet of measurements
#'
#' @examples todo
#'
#' @import pomp
#'
#' @export

dmeas_D2 <- function() {
  return(Csnippet("
                  double lik1;
                  double lik2;
                  if (ISNA(LAI)) {
                    lik1 = (give_log) ? 0 : 1;
                  } else if (!ISNA(LAI)){
                    lik1 = dnorm(lma*LAI, Cfol, Cfol*.25, give_log);
                  }
                  if (ISNA(Lf)) {
                    lik2 = (give_log) ? 0 : 1;
                  } else if (!ISNA(Lf)){
                    double pi = atan(1) * 4;
                    double s = 365.25/pi;
                    double psi_f = -1.3588480 + 4.5994549 * c_lf - 9.5964931 * pow(c_lf, 2.0) + 12.1567793 * pow(c_lf, 3.0) -6.8903864 * pow(c_lf, 4.0) -0.3576296 * pow(c_lf, 5.0) + 1.3941622 * pow(c_lf, 6.0);
                    double fall = ((pow(2/pi, .5) * (-log(1 - c_lf)/c_rfall) * exp(-pow((sin((yearday - d_fall + c_rfall*psi_f)/s) * pow(2, .5) * s / c_rfall), 2))));
                    lik2 = dnorm(Lf, Cfol*fall, (5 + 1*Cfol*fall), give_log);
                  }
                  lik = (give_log) ? lik1 + lik2 : lik1 * lik2 ;
                "))
}