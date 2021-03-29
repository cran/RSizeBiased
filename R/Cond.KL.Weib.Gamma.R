Cond.KL.Weib.Gamma <- function(par,nullvalue,hata,hatb,type,dist)
  {
  #K:shape
  #L:scale
  #type==1 = mean
  #(hatk-1)*psi(0, hatk)-log(halL)-hatk-log(Gamma(hatk))
  K=par
  L = switch(dist,
         "weib"= ifelse(type==1, nullvalue/gamma((1+K)/K), sqrt(nullvalue/(gamma((2+K)/K)-(gamma((1+K)/K))^2))),
         "gamma" = ifelse(type==1, nullvalue/K, sqrt(nullvalue/K))
         )

  result = switch(dist,
     "weib"=   -log(K/L^K)+(hata-K)*(log(hatb)+digamma(1)/hata)+(hatb/L)^K*gamma(K/hata+1),
    "gamma" =  log(gamma(K))+K*log(L)-(K-1)*(psi(0, hata)+log(hatb))+hatb*hata/L
 )

  return(result)
 }


