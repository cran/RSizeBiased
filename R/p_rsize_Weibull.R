p_rsize_Weibull <- function(q,TRpar,r){
  shape <- TRpar[1]
  scale <- TRpar[2]
  1 - gammainc((q/scale)**shape,(r + shape)/shape)[2]/gamma((r + shape)/shape)
}
