log_Lik_Weib_gamma_weighted <- function(TRpar, datain, r, dist){
  shape <- TRpar[1]
  scale <- TRpar[2]
  x <- datain
  loglik <- switch(dist,
           "weib"= sum(-(x**shape/scale**shape) - r*log(scale) +
                   log(shape) + (-1 + r)*log(x) + shape*(-log(scale) + log(x)) - log(gamma((r + shape)/shape))),
            "gamma" = sum(log(dgamma(x,shape = shape+r, scale = scale, log = FALSE)))
  )
  return(-loglik)
}


