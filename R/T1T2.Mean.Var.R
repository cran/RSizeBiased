T1T2.Mean.Var <- function(datain,r, type)
  {
   test.stat<- ifelse(type==1, sum(datain^(1-r))/sum(datain^(-r)),
                        sum(datain^(2-r))/sum(datain^(-r))-(sum(datain^(1-r))/sum(datain^(-r)))^2)
   return(test.stat)
}


