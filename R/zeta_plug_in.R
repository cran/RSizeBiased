zeta_plug_in <- function(null_value, datain,r,EST_par,type, dist){
  n<-length(datain)
  if(type==1){
    z <- (T1T2.Mean.Var(datain,r, type)- null_value)/sqrt(s11.s22(EST_par,r, "s11", dist)/n)
  } else {
    z <- (T1T2.Mean.Var(datain,r, type)- null_value)/sqrt(s11.s22(EST_par,r, "s22", dist)/n)
  }
  return(z)
}

