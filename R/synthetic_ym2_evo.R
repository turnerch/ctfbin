#' Synthetic YM2 Evolution Function
#'
#' Simulates one step forwards in time in a time series process.
#' 
#' @author John W. Smith Jr
#'
#' @param N/A
#'
#' @return Csnippet of evolution function
#'
#' @examples
#' synthetic_ym2_evo()
#' 
#' @import pomp
#'
#' @export

synthetic_ym2_evo <- function() {
  return(Csnippet("
                  double psid = -2;
                  double rtot = 1;
                  double trange;
                  double gs;
                  double pp;
                  double qq = -204.6453;
                  double ci;
                  double e0;
                  double dec;
                  double mult;
                  double dayl;
                  double cps;
                  double lat_e = 46.23391;
                  double a0 = c_eff; 
                  double a1 = 0.0156935; 
                  double a2 = 4.22273;
                  double a3 = 208.868;
                  double a4 = 0.0453194; 
                  double a5 = 0.37836;
                  double a6 = 7.19298;
                  double a7 = 0.011136;
                  double a8 = 2.1001; 
                  double a9 = 0.789798;
                  double pi = 3.141593;
                  double s = 365.25/pi;
                  double G;
                  double psi_f = -1.3588480 + 4.5994549 * c_lf - 9.5964931 * pow(c_lf, 2.0) + 12.1567793 * pow(c_lf, 3.0) -6.8903864 * pow(c_lf, 4.0) -0.3576296 * pow(c_lf, 5.0) + 1.3941622 * pow(c_lf, 6.0);
                  double onset = pow(2/pi, .5) * (6.9088/c_ronset) * exp(-pow(sin((yearday - d_onset - .6245*c_ronset)/s) * pow(2, .5) * s / c_ronset, 2));
                  double fall = ((pow(2/pi, .5) * (-log(1 - c_lf)/c_rfall) * exp(-pow((sin((yearday - d_fall + c_rfall * psi_f)/s) * pow(2, .5) * s / c_rfall), 2))));
                  trange = .5*(MaxT - MinT);
                  gs = pow(fabs(psid),(0.789798)) / (0.37836*rtot + trange);
                  pp = (Cfol / lma) /gs*a0*exp(a7*(MaxT));
                  ci = .5*(ca + qq - pp + pow(pow(ca+qq-pp,2) - 4*(ca*qq - pp*a2), .5 ));
                  e0 = (a6*pow(Cfol/lma,2)) / (pow(Cfol/lma,2) + a8);
                  dec = -23.4*cos((360*(yearday + 10)/365)*pi/180)*pi/180;
                  mult = tan(lat_e*pi/180)*tan(dec);
                  if (mult >= 1){
                    dayl = 24;
                  } else if(mult <= -1){
                    dayl = 0;
                  } else{
                    dayl = 24*acos(-mult)/pi;
                  }
                  cps = e0*rad*gs*(ca - ci) / (e0*rad + gs*(ca - ci));
                  G = cps*(a1*dayl + a4);
                  double eps_omega = .000001;
                  double mu_lab = (1 - onset) * Clab + f_lab * G;
                  double mean_lab = log(mu_lab) - .5*log(1 + pow(omega_lab,2));
                  double sd_lab = pow(log(1 + pow(omega_lab,2)), .5) + eps_omega;
                  double mu_fol = (1 - fall)*Cfol + onset*Clab + f_f * G;
                  double mean_fol = log(mu_fol) - .5*log(1 + pow(omega_f,2) + eps_omega);
                  double sd_fol = pow(log(1 + pow(omega_f,2)), .5);
                  
                  Cfol = rlnorm(mean_fol, sd_fol + eps_omega);
                  Clab = exp(rnorm(mean_lab, sd_lab + eps_omega));
                  "))
}