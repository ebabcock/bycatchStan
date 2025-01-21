data{
 int N;
 int Nall;
 int Ncoef;
 int Y[N];
 vector[N] offset;
 matrix[N,Ncoef] xMatrix;
 matrix[Nall,Ncoef] xMatrixAll;
 vector[Nall] EffortAll;
}
parameters{
 vector[Ncoef] b;
 real<lower=0.00001,upper=100> phi;
}
transformed parameters{
  vector[N] logmu;
  vector[N] mu;
  logmu = xMatrix*b;
  for(i in 1:N) {
   mu[i] = exp(logmu[i])*offset[i];
  }
}
model{
  b~normal(0,10);
  phi~normal(0,1);
  Y~neg_binomial_2(mu,phi);
}
generated quantities {
  real LL[N];
  real Yrep[N];
  real muAll[Nall];
  real StrataBycatch[Nall];
  for(i in 1:N) {
   Yrep[i] = neg_binomial_2_rng(mu[i],phi);
   LL[i] = neg_binomial_2_lpmf(Y[i]|mu[i],phi);
  }
  for(i in 1:Nall) {
    muAll[i] = exp(xMatrixAll[i,]*b)*EffortAll[i];
    StrataBycatch[i] = neg_binomial_2_rng(muAll[i],phi);
  }
}



