data {
  // Number of RNAs transcipts searched for
  int<lower = 1> NRNA;     
  
  int <lower = 1> NReplicates;
  
  // Array of integers for the counts data
  // for each transcript for each sample
  int<lower = 0> counts[NRNA,NReplicates,2];
}
parameters {
  vector[NRNA] muWild;
  
  vector[NRNA] muMut;

  // scale factor accounting for differing total RNA counts
  vector[NReplicates] aMut;
  
  vector[NReplicates] aWild;
  
  // dispersion parameter for counts
  real<lower = 0> phi;
  
  
  // must explicitly state support of priors, i.e. upper and lower bounds, otherwise samplier does not act as expected!!!

}

model{
  
   // sample the means from a log uniform distribution
   muWild ~ normal(8,2);
   
   muMut ~ normal(8,2);
   
   // scale factor prior; one for each condition
   aMut ~ normal(0,0.01);
   
   aWild ~ normal(0,0.01);
    
  // Cauchy prior for phi parameter;
   phi ~ cauchy(0,3);
  
  for(i in 1:NRNA) for(j in 1:NReplicates){ 
    counts[i,j,1] ~ neg_binomial_2(2^(aWild[j]+muWild[i]),phi);
  }
  
  for(i in 1:NRNA) for(j in 1:NReplicates){
    counts[i,j,2] ~ neg_binomial_2_log(2^(aMut[j]+muMut[i]),phi);
  }
}
