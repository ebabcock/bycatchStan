printStanCode<-function() {
  
  #write out binomial stan model for survival/mortality. No predictions
write("data{
 int N;
 int Ncoef;
 array[N] int Y;
 matrix[N,Ncoef] xMatrix;
}
parameters{
 vector[Ncoef] b;
}
transformed parameters{
  vector[N] logitmu;
  logitmu = xMatrix*b;
}
model{
  b~normal(0,10);
  Y~bernoulli_logit(logitmu);
}
generated quantities {
  array[N] real LL;
  array[N] real Yrep;
  for(i in 1:N) {
   Yrep[i] = bernoulli_logit_rng(logitmu[i]);
   LL[i] = bernoulli_logit_lpmf(Y[i]|logitmu[i]);
  }
}
", file="stan/binomial.stan")
  
  #Write out binomial stan model with new data to predict
  write(
    "data{
 int N;
 int Nall;
 int Ncoef;
 array[N] int Y;
 matrix[N,Ncoef] xMatrix;
 matrix[Nall,Ncoef] xMatrixAll;
}
parameters{
 vector[Ncoef] b;
}
transformed parameters{
  vector[N]  logitmu;
  logitmu = xMatrix*b;
}
model{
  b~normal(0,10);
  Y~bernoulli_logit(logitmu);
}
generated quantities {
  array[N] real LL;
  array[N] real Yrep;
  array[Nall] real logitMuAll;
  array[Nall] real strataProb;
  for(i in 1:N) {
   Yrep[i] = bernoulli_logit_rng(logitmu[i]);
   LL[i] = bernoulli_logit_lpmf(Y[i]|logitmu[i]);
  }
  for(i in 1:Nall) {
    logitMuAll[i] = xMatrixAll[i,]*b;
    strataProb[i] = inv_logit(logitMuAll[i]);
  }
}

", file="stan/binomialP.stan")
  
  
#Write out the stan file for negative binomial with
#no simulation of effort, includes priors
write("data{
 int N;
 int Ncoef;
 array[N] int Y;
 real interceptSD;
 real coefficientSD;
 int phiType;
 real phiPar;
 vector[N] Effort;
 matrix[N,Ncoef] xMatrix;
}
parameters{
 real  b0;
 vector[Ncoef] b;
 real<lower=0.00001,upper=100> phi;
}
transformed parameters{
  vector[N] logmu;
  vector[N] mu;
  logmu = b0+xMatrix*b;
  for(i in 1:N) {
   mu[i] = exp(logmu[i])*Effort[i];
  }
}
model{
  b0~normal(0,interceptSD);
  b~normal(0,coefficientSD);
  if(phiType==1) {
    phi~exponential(phiPar);
  }  else  {
    phi~normal(0,phiPar);
  }
  Y~neg_binomial_2(mu,phi);
}
generated quantities {
  array[N] real LL;
  array[N]real Yrep;
  for(i in 1:N) {
   Yrep[i] = neg_binomial_2_rng(mu[i],phi);
   LL[i] = neg_binomial_2_lpmf(Y[i]|mu[i],phi);
  }
}

",file="stan/NB2matrixNoEffort.stan")
  
#No effort intercept only
  write(
    "data{
 int N;
 int Ncoef;
 real interceptSD;
 int phiType;
 real phiPar;
 array[N] int Y;
 vector[N] Effort;
}
parameters{
 real  b0;
 real<lower=0.00001,upper=100> phi;
}
transformed parameters{
  vector[N] logmu;
  vector[N] mu;
  for(i in 1:N) {
   logmu[i] = b0;
   mu[i] = exp(logmu[i])*Effort[i];
  }
}
model{
  b0~normal(0,interceptSD);
  if(phiType==1) {
    phi~exponential(phiPar);
  }  else  {
    phi~normal(0,phiPar);
  }
  Y~neg_binomial_2(mu,phi);
}
generated quantities {
  array[N] real LL;
  array[N] real Yrep;
  for(i in 1:N) {
   Yrep[i] = neg_binomial_2_rng(mu[i],phi);
   LL[i] = neg_binomial_2_lpmf(Y[i]|mu[i],phi);
  }
}

",file="stan/NB2matrixNoEffort1.stan")
  
  
}