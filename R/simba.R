#' Synthetic YM2 default parameters
#'
#' todo
#' 
#' @author Munir Winkel
#'
#' @param X: todo
#'
#' @return todo
#'
#' @examples
#' todo
#'
#' @import
#'
#' @export

simba <- function(X)
    #### library(fields)
    #### mm = 40; X = makepts( mm = mm , p = 6)
    x1 = X[,1]; x2 = X[,2]; x3 = X[,3];
    x4 = X[,4]; x5 = X[,5]; x6 = X[,6];
    # par(mfrow=c(3,2))
    # for( v in seq(0,0.3, length.out = 6)){
    #   x3  =   v
    # x4 = x5 = x6 = 0.1
    #upper bumps, involving all variables
    Y = 3.14749 + sin(  2*pi*(x1^2 - 2*x2*(1 + x3) ) )*
      ( pnorm( 30 *(x2 -.3))  + pnorm(30*(.8 - x2)) - 1 )*
      2*sin( 4*pi*x1 + 3*pi*(1+x3) + 2*pi*(x4 + x5) + 3*pi*(1 + x6))  +
      
      ### Simba lion king when x3 near 0.2
      (4 + 6*(x1) )*
      ((pnorm( 30*(x2 - .0)) + pnorm(30*(.2 - x2))  ) - 1  )*
      ((pnorm( 30*(x1 - .0)) + pnorm(30*(.6 - x1))  ) - 1  )*
      (pnorm( 10*(.2 - x3))) +
      
      ### other area below simba lion king
      ### Simba lion king when x3 near 0.2
      (  1 - 8*(x1 + x2 - x4 - x5 - x6)^2 )*
      ((  pnorm( 40*(x2 - .0)) + pnorm(40*(.2 - x2))  ) - 1  )*
      ((  pnorm( 40*(x1 - .6)) + pnorm(40*(1 - x1))  ) - 1  )*
      (   pnorm( 10*(.2 - x3))  ) +
      
      #### unassuming, locally active region when x3 > .2
      .5*(  1 - sin( 8*pi*x1 + 7*pi*x2*x3 - 4*pi*x4*x5*x6) )*
      ((  pnorm( 30*(x2 - .0)) + pnorm(30*(.3 - x2))  ) - 1  )*
      (   pnorm( 8*(x3 - .3))  ) +
      
      #### deceptive maximum, when x2 > .8, x3 > .2, involving
      #### all variables
      (5*cos( 2*(x2 + .5)*(-x4 + .5)*(-x5+.5)^2 )*(-x6 - .5) - 
           .02*((1-x2)^2 + 
                  (1-x1)^2 + 
                  (1-x3  - .3*x4)^2  + 
                  (1-x5  + .5*x4)^2 + 
                  (.8-x6 - .4*x4)^2 )   )*
      (pnorm(5*(x2 - 1) + pnorm(10*(.5 - x3))))
    # # 
    #     best = which.max(Y); X[best,]
    #     image.plot(matrix(Y,mm), main = paste("Simba, X3 = ", mean(x3), sep = "",
    #                 ",  f(x*) = ", round(max(Y),2)),
    #                xlab = "X1", ylab = "X2")
    # 
    #     points(x1[best], x2[best] , pch = "*", cex = 2, col = 'white')
    # }
    
    #     return(Y)
    
    
    
  }# end 
