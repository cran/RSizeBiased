Size.BiasedMV.Tests<-function(datain_r,r,nullMEAN,nullVAR, start_par,nboot,alpha,prior_sel, distr){

  if(!is.vector(datain_r)){
    stop("'datain_r' should be a vector")
  }

  if(length(datain_r)<2){
    stop("'datain_r' should have at least 2 observations")
  }

  if(r<0){
    stop("'r' should be non-negative")
  }

  if(nullMEAN<=0){
    stop("'nullMEAN' should be positive")
  }
  if(nullVAR<=0){
    stop("'nullVAR' should be positive")
  }

  if(distr=="weib"){
    # func<-log_Lik_weib_gamma_weighted
    # plugin<-zeta_plug_in
    # KLmean<-Cond.KL.Weib.Gamma
    # KLVAR<-Cond.KL.Weib.Gamma
  }else if(distr=="gamma"){
    # func<-log_Lik_weib_gamma_weighted
    # plugin<-zeta_plug_in
    # KLmean<-Cond.KL.Weib.Gamma
    # KLVAR<-Cond.KL.Weib.Gamma
  }else {
    stop("'distr' not recognized")
  }

  if(!is.vector(start_par)){
    stop("'start_par' should be a vector")
  }

  if(length(start_par)!=2){
    stop("'start_par' should be a two elements vector containing the starting values for the MLE for the two parameter distribution (weibull or gamma)")
  }

  if(nboot<0){
    stop("'nboot' should be a non-negative integer")
  }

  if(nboot!=round(nboot)){
    stop("'nboot' should be non-negative integer")
  }

  if(any(alpha<=0)|any(alpha>=1)){
    stop("All elements of 'alpha' should be in (0,1)")
  }

  if(!(prior_sel=="normal"|prior_sel=="gamma")){
    stop("'prior_sel' not recognized")
  }

  alpha=sort(alpha)
  boot_T1_mean<-rep(NA,nboot)
  boot_T2_var<-rep(NA,nboot)
  pvaluebT1<-NA
  pvaluebT2<-NA


  est_as_r_biased <- optim(start_par,log_Lik_Weib_gamma_weighted,datain=datain_r,r=r, dist=distr, hessian = TRUE)
  EST_par <- est_as_r_biased$par
  CovMatrix<-inv(est_as_r_biased$hessian)
  n=length(datain_r)

  asymptotic<-0
  if(EST_par[1]>r){
    asymptotic<-1
  }
  #Step 2. compute and store the observed values of the test statistics
  T1value<-T1T2.Mean.Var(datain_r,r, 1)
  T2value<-T1T2.Mean.Var(datain_r,r, 2)
  Tivalues<-c(T1value,T2value)

  if(asymptotic==1){
    #Step 2. compute and store the observed values of the these statistics
    Zeta_i<-c(zeta_plug_in(nullMEAN, datain_r,r,EST_par,1, distr),  zeta_plug_in(nullVAR, datain_r,r,EST_par,2, distr))
  }else{
    Zeta_i<-c(NA,NA)
  }

  #Step 3. determine the closest to the estimated distribution, member of the class
  #that satisfies the equation for the tested null mean or variance value,
  shape_c_mean<-optimize(Cond.KL.Weib.Gamma,c(0.1*EST_par[1], 10*EST_par[1]),nullvalue=nullMEAN,hata=EST_par[1],hatb=EST_par[2], type=1, dist=distr)$minimum
  if(distr=="weib"){
    scale_c_mean<-nullMEAN/gamma((1+shape_c_mean)/shape_c_mean)
  }else if(distr=="gamma"){
    scale_c_mean<-nullMEAN/shape_c_mean
  }

  shape_c_var<-optimize(Cond.KL.Weib.Gamma,c(0.1*EST_par[1], 10*EST_par[1]),nullvalue=nullVAR,hata=EST_par[1],hatb=EST_par[2], type=2, dist=distr)$minimum
  if(distr=="weib"){
    scale_c_var<-sqrt(nullVAR/(gamma((2+shape_c_var)/shape_c_var)-(gamma((1+shape_c_var)/shape_c_var))^2))
  }else if(distr=="gamma"){
    scale_c_var<-sqrt(nullVAR/shape_c_var)
  }

  if(nboot>0){
    for (ib in 1:nboot){
      #step 4 and 5_mean. assign normal prior distribution to the distribution parameters and generate values from the prior distributions
        rscale_mean=-99
        while(rscale_mean<10^(-7)|is.na(rscale_mean)){
          if(prior_sel=="normal"){
            rshape_mean=-99
            while(rshape_mean<10^(-7)){
              rshape_mean <- rnorm(1, mean = shape_c_mean, sd = sqrt(CovMatrix[1,1]))
            }
          }else if(prior_sel=="gamma"){
            rscale_mean=-99
            agamma=shape_c_mean^2/(CovMatrix[1,1])
            bgamma=(CovMatrix[1,1])/shape_c_mean
            rshape_mean=-99
            while(rshape_mean<10^(-7)){
              rshape_mean <- rgamma(1, shape = agamma, scale  = bgamma)
            }
          }
          if(distr=="weib"){
            rscale_mean <- nullMEAN/(gamma((1+rshape_mean)/rshape_mean))
          }else if(distr=="gamma"){
            rscale_mean <- nullMEAN/rshape_mean
          }
      }

      trueboot_mean<-c(rshape_mean,rscale_mean)

      #step 6_mean. generate an r-size biased bootstrap sample of size n from the proposed distribution of Step 5
      if(distr=="weib"){
        boot_data_r_mean <-  r_rsize_Weibull(n,c(rshape_mean,rscale_mean),r)
      }else if(distr=="gamma"){
        boot_data_r_mean <- rgamma(n, shape = rshape_mean+r, scale =rscale_mean)
      }
      #Step 8_mean. compute and store the observed values of the test statistic
      boot_T1_mean[ib]<-T1T2.Mean.Var(boot_data_r_mean,r, type=1)

      #step 4 and 5_var. assign normal prior distribution to the distribution parameters and generate values from the prior distributions
      rscale_var=-99
        while(rscale_var<10^(-7)|is.na(rscale_var)){
          if(prior_sel=="normal"){
            rshape_var=-99
            while(rshape_var<10^(-7)){
              rshape_var <- rnorm(1, mean = shape_c_var, sd = sqrt(CovMatrix[1,1]))
            }
          }else if(prior_sel=="gamma"){
            rscale_var=-99
            agamma=shape_c_var^2/(CovMatrix[1,1])
            bgamma=(CovMatrix[1,1])/shape_c_var
            rshape_var=-99
            while(rshape_var<10^(-7)){
              rshape_var <- rgamma(1, shape = agamma, scale  = bgamma)
            }
          }
          if(distr=="weib"){
            rscale_var<-sqrt(nullVAR/(gamma((2+rshape_var)/rshape_var)-(gamma((1+rshape_var)/rshape_var))^2))
          }else if(distr=="gamma"){
            rscale_var<-sqrt(nullVAR/rshape_var)
          }
        }

      trueboot_var<-c(rshape_var,rscale_var)
      if(distr=="weib"){
        boot_data_r_var <-r_rsize_Weibull(n,c(rshape_var,rscale_var),r)
      }else if(distr=="gamma"){
        boot_data_r_var <- rgamma(n, shape = rshape_var+r, scale =rscale_var)
      }
      boot_T2_var[ib]<-T1T2.Mean.Var(boot_data_r_var,r, type=2)
    }
    T1qunat<-c(quantile(boot_T1_mean,alpha/2),rev(quantile(boot_T1_mean,1-alpha/2)))
    T2qunat<-c(quantile(boot_T2_var,alpha/2),rev(quantile(boot_T2_var,1-alpha/2)))
    pvaluebT1<-2*min(length(which(Tivalues[1]>boot_T1_mean)),length(which(Tivalues[1]<boot_T1_mean)))/nboot
    pvaluebT2<-2*min(length(which(Tivalues[2]>boot_T2_var)),length(which(Tivalues[2]<boot_T2_var)))/nboot

  }
  nalpha=length(alpha)
  decision <- data.frame(matrix(nrow =2, ncol = nalpha))
  decisionasympt <- data.frame(matrix(nrow =2, ncol = nalpha))

  rownames(decision)<-c("m0","var0")
  rownames(decisionasympt)<-c("m0","var0")
  colnames(decision)<-alpha
  colnames(decisionasympt)<-alpha


  for(iii in 1:nalpha){
    decision[1,iii] <- ifelse(Tivalues[1]<T1qunat[iii] | Tivalues[1]>T1qunat[2*nalpha-iii+1],1,0)
    decision[2,iii] <- ifelse(Tivalues[2]<T2qunat[iii] | Tivalues[2]>T2qunat[2*nalpha-iii+1],1,0)
  }

  if(asymptotic==1){
    asymptotic_p_mean<-2*(1-pnorm(abs(Zeta_i[1]),0,1))
    asymptotic_p_var<-2*(1-pnorm(abs(Zeta_i[2]),0,1))
  }else{
    asymptotic_p_mean<-NA
    asymptotic_p_var<-NA
  }

  if(asymptotic==1){
    for(iii in 1:nalpha){
      decisionasympt[1,iii] <- ifelse(asymptotic_p_mean<alpha[iii],1,0)
      decisionasympt[2,iii] <- ifelse(asymptotic_p_var<alpha[iii],1,0)
    }
  }else{
    decisionasympt[1,iii] <- NA
    decisionasympt[2,iii] <- NA
  }
  return(list(par = EST_par, loglik=-est_as_r_biased$value, CovMatrix = CovMatrix,Zeta_i=Zeta_i,Tivalues=Tivalues
              ,T1_bootstrap_quan=T1qunat,T2_bootstrap_quan=T2qunat,NullValues=c(nullMEAN,nullVAR),distribution=distr,alpha=alpha,bootstrap_p_mean=pvaluebT1,
              bootstrap_p_var=pvaluebT2,decision=decision,asymptotic_p_mean=asymptotic_p_mean,asymptotic_p_var=asymptotic_p_var,
              decisionasympt=decisionasympt,prior_sel=prior_sel))
}
