d_rsize_Weibull <- function(x,TRpar,r){
  shape <- TRpar[1]
  scale <- TRpar[2]
  x^r*dweibull(x, shape=shape, scale = scale, log = FALSE)/r_moment_gamma_Weib(TRpar,r, "weib")
}
