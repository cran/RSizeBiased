r_rsize_Weibull <- function(n, TRpar,r) {
  u<-rgamma(n, shape = r/TRpar[1]+1, scale =1)
  y<-TRpar[2]*u^(1/TRpar[1])
  return (y)  
}