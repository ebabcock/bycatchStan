data{
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


