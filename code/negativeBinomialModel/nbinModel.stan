data {
  // Number of RNAs transcipts searched for
  int<lower = 1> NRNA;     
  
  int <lower = 1> NReplicates;
  
  // Array of integers for the counts data
  // for each transcript for each sample
  int<lower = 0> counts[NRNA,NReplicates,2];
}
transformed data{
  int K = (2 * NReplicates);
  
}

parameters {
  vector[NRNA] muWild;
  
  vector[NRNA] muMut;
  
  
  // scale factor accounting for differing total RNA counts
  vector[(2 * NReplicates - 1)] scaleFactorsRaw;
  
  // dispersion parameter for counts
  vector<lower = 0>[NRNA] phiWild;
  vector<lower = 0>[NRNA] phiMut;
  
  
  // must explicitly state support of priors, i.e. upper and lower bounds, otherwise samplier does not act as expected!!!

}

transformed parameters {
      // centre scale factors so that they are all compariable!
      vector[K] scaleFactors;  // centered
      for (k in 1:(K-1)) {
        scaleFactors[k] = scaleFactorsRaw[k];
      }
      scaleFactors[K] = -sum(scaleFactorsRaw);
}

model{
  // mean coefficient for dispersion regression
  real beta1;
  
  // intercept for dispersion regression
  real beta0;
  
  // variance for dispersion regression
  real<lower = 0> sigma;
  
   // sample the means from a log uniform distribution
   muWild ~ normal(4,4);
   
   muMut ~ normal(4,4);
   
   sigma ~ exponential(1);
   
   beta1 ~ normal(-1,1);
   
   beta0 ~ normal(8,1);
   
    // Linearly regress phi over mean
    log2(phiWild) ~ normal( muWild * beta1 + beta0, sigma);
    target += - log(fabs(phiWild) * log(2));
    
    log2(phiMut) ~ normal(muMut * beta1 + beta0, sigma);
    target += - log(fabs(phiMut) * log(2));
   
   // scale factor prior; one for each condition
   scaleFactorsRaw ~ double_exponential(0,0.1);
  
  for(i in 1:NRNA) for(j in 1:NReplicates){ 
    counts[i,j,1] ~ neg_binomial_2(2^(scaleFactors[j]+muWild[i]),phiWild[i]);
  }
  
  for(i in 1:NRNA) for(j in 1:NReplicates){
    counts[i,j,2] ~ neg_binomial_2(2^(scaleFactors[(NReplicates + j)]+muMut[i]),phiMut[i]);
  }
}
