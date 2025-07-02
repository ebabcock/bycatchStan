data{
 int N;
 int Ncoef;
 real interceptSD;
 int phiType;
 real phiPar;
 int Y[N];
 vector[N] offset;
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
   mu[i] = exp(logmu[i])*offset[i];
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
  real LL[N];
  real Yrep[N];
  for(i in 1:N) {
   Yrep[i] = neg_binomial_2_rng(mu[i],phi);
   LL[i] = neg_binomial_2_lpmf(Y[i]|mu[i],phi);
  }
}


