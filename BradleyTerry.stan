data{
int<lower=0> N_dyads;       // number of dyads
int<lower=0> N_ids;         // number of individuals
int focal_id;               // the focal individual whose dominance value is fixed
int<lower=0> ind1[N_dyads];   // first individual in the dyad
int<lower=0> ind2[N_dyads];   // second individual in the dyad
int<lower=0> win1[N_dyads];   // number of wins by the first individual
int<lower=0> win2[N_dyads];   // number of wins by the second individual
}

transformed data{
int<lower=0> n_ints[N_dyads];     // total interactions between each dyad
for(n in 1:N_dyads) {
n_ints[n] = win1[n] + win2[n];
  }
}

parameters{
vector[N_ids] d;            // latent dominance value parameter
}

transformed parameters{
vector[N_ids] d_fix;                  // transform d into d_fix, to fix the focal individual's d value
d_fix[focal_id] = 0.0;
for( n in 1:(focal_id-1) ) d_fix[n] = d[n];
for( n in (focal_id+1):N_ids ) d_fix[n] = d[n];
}

model{
vector[N_dyads] p;        // local parameter

d ~ normal( 0 , 1);        // prior on d

for( n in 1:N_dyads) {
    // probability of winning = the difference in dominance values between individuals
    // uses the logistic/inverse logit function to transform to probability
  p[n] = inv_logit( d_fix[ ind1[n] ] - d_fix[ ind2[n] ] );
 }

  win1 ~ binomial( n_ints , p );
}
