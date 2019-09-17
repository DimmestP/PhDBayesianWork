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
  real<lower = 0> phi;
  
  
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
  
   // sample the means from a log uniform distribution
   muWild ~ normal(8,2);
   
   muMut ~ normal(8,2);
   
   // scale factor prior; one for each condition
   scaleFactorsRaw ~ double_exponential(0,0.1);
    
  // Cauchy prior for phi parameter;
   phi ~ cauchy(0,3);
  
  for(i in 1:NRNA) for(j in 1:NReplicates){ 
    counts[i,j,1] ~ neg_binomial_2(2^(scaleFactors[j]+muWild[i]),phi);
  }
  
  for(i in 1:NRNA) for(j in 1:NReplicates){
    counts[i,j,2] ~ neg_binomial_2(2^(scaleFactors[(NReplicates + j)]+muMut[i]),phi);
  }
}
