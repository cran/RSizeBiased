r_moment_gamma_Weib <- function(TRpar,r, dist){
  shape <- TRpar[1]
  scale <- TRpar[2]
  rm <- switch(dist,
               "weib"=  scale**r*gamma((r + shape)/shape),
               "gamma"= scale**r*gamma(r + shape)/gamma(shape)
  )
  return(rm)
}
