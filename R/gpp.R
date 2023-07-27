#' Modified GPP Function
#'
#' Creates GPP function in order to compute the GPP for a day based on values found in the drivers_unde data.
#' GPP function from Aggregated Canopy Model, modified for the reparameterization
#' of the DALEC2 model.  For use in synthetic_ym2 default particle filter.
#' 
#' @author John W. Smith Jr
#'
#' @param lai Leaf Area Index, is unitless and is the total surface of leaves per area of ground
#' @param c_eff Canopy Efficiency, is unitless but describes how efficiently a tree's canopy uses light
#' @param maxt Maximum temperature, in celsius
#' @param mint Minimum temperature, in celsius
#' @param Ca Atmospheric carbon dioxide, in parts per million
#' @param lat global latitude of location of observation
#' @param yearday which day of the year observation was recorded, as an integer
#' @param rad daily shortwave radiation, in gCm ^ -1
#'
#' @return GPP (gross primary productivity) from input provided
#'
#' @examples
#' GPP(lai = 2.641, c_eff = 20, maxt = 12.83247, mint = 6.867838, Ca = 400.95, lat = 46.23, yearday = 268, rad = 11.26568)
#'
#' @export

GPP = function(lai, c_eff, maxt, mint, Ca, lat, yearday, rad){
  ## Check if function parameters are legal
  if(!is.numeric(lai)) stop("lai parameter must be numeric!")
  if(!is.numeric(c_eff)) stop("c_eff parameter must be numeric!")
  if(!is.numeric(maxt)) stop("maxt parameter must be numeric!")
  if(!is.numeric(mint)) stop("mint parameter must be numeric!")
  if(!is.numeric(Ca)) stop("Ca parameter must be numeric!")
  if(!is.numeric(lat)) stop("lat parameter must be numeric!")
  if(!is.numeric(yearday)) stop("yearday parameter must be numeric!")
  if(!is.numeric(rad)) stop("rad parameter must be numeric!")
  
  ## Sets parameters
  psid = -2
  rtot = 1
  trange = .5*(maxt - mint)
  a = c(c_eff, 
        0.0156935, 
        4.22273,
        208.868,
        0.0453194, 
        0.37836,
        7.19298,
        0.011136,
        2.1001, 
        0.789798)
  gs = abs((psid))^a[10] / (a[6]*rtot + trange)
  pp = lai/gs*a[1]*exp(a[8]*maxt)
  qq = a[3] - a[4]
  ci = .5*(Ca + qq - pp + ((Ca+qq-pp)^2 - 
                             4*(Ca*qq - pp*a[3]))^(.5) )
  e0 = (a[7]*lai**2 / (lai**2 + a[9]) )
  dec = -23.4*cos((360*(yearday + 10)/365)*pi/180)*pi/180
  mult = tan(lat*pi/180)*tan(dec)
  
  ## Changes dayl variable depending on value of mult
  if (mult >= 1){
    dayl = 24
  } else if(mult <= -1){
    dayl = 0
  } else{
    dayl = 24*acos(-mult)/pi
  }
  cps = e0*rad*gs*(Ca - ci) / (e0*rad + gs*(Ca - ci))
  gpp = cps*(a[2]*dayl + a[5])
  return(gpp)
}