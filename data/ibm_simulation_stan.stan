data{
  int N;
  int k;
  vector [N] sizet;
  vector [N] growth_std;
  int type[N];
  real sizemat[N,N];
  real distmat[N,N];
  
  int nn;
  real Itype[nn,k];
  real sigma_Itype[nn];
  real a0[nn];
  real c[nn,k];
  real c0[nn];
  real sigma[nn];
  real a1[nn,k];
  real a01[nn];
  real a3[nn,k];
  real a03[nn];
  real a2[nn];
  real sigma_a1[nn];
  real sigma_a3[nn];
  real sigma_c[nn];
}
transformed data{
  vector[N] sizet_std;
  sizet_std=(sizet-mean(sizet))/(2*sd(sizet));

}
transformed parameters{
  real mu_dd[nn,N];
  real mu_base[nn,N];
  
  {
  vector[N] cf;
  vector[N] cf_std;
  vector[N] smat[N];
  vector[N] dmat[N];
  
  for(n in 1:nn){
  
    for(i in 1:N){
      for(j in 1:N){
        smat[i,j] = sizemat[i,j]^a1[nn,type[i]];
        dmat[i,j] = distmat[i,j]^2*a2[nn];
      }
      cf[i] = sum(smat[i] ./ exp(dmat[i]));
    }
    
    cf_std = (cf - mean(cf)) / (2*sd(cf));
    
    for(l in 1:N){
      mu_dd[n,l]=Itype[n,type[l]] 
        + (c[n,type[l]]/(2*sd(sizet)))*mean(sizet) 
        + (c[n,type[l]]/(2*sd(sizet)))*sizet[l]
        + (a3[n,type[l]]/(2*sd(cf)))*mean(cf) 
        + (a3[n,type[l]]/(2*sd(cf)))*cf[l];
      mu_base[n,l]=Itype[n,type[l]]
        + (c[n,type[l]]/(2*sd(sizet)))*mean(sizet) 
        + (c[n,type[l]]/(2*sd(sizet)))*sizet[l];
    }
    
  }
  }
}
generated quantities{
  real pred_dd[nn,N];
  real pred_base[nn,N];
  for(n in 1:nn){
    for(i in 1:N){
      pred_dd[n,i] = normal_rng(mu_dd[n,i], sigma[n]);
      pred_base[n,i] = normal_rng(mu_base[n,i], sigma[n]);
  }}
}
