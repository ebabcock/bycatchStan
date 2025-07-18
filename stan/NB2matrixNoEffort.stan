data{
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


