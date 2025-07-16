data{
 int N;
 int Ncoef;
 array[N] real Y;
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

