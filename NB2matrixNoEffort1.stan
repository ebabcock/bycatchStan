data{
 int N;
 int Ncoef;
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
  b0~normal(0,10);
  phi~normal(0,1);
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


