data{
int<lower=0> N_dyads;       // number of dyads
int<lower=0> N_ids;         // number of individuals
int<lower=0> P;              // number of covariates
int<lower=0> ind1[N_dyads];   // first individual in the dyad
int<lower=0> ind2[N_dyads];   // second individual in the dyad
int<lower=0> win1[N_dyads];   // number of wins by the first individual
int<lower=0> win2[N_dyads];   // number of wins by the second individual
row_vector[P] X[N_ids];       // covariates; 0 must be meaningful!
}

transformed data{
int<lower=0> n_ints[N_dyads];     // total interactions between each dyad
for(n in 1:N_dyads) {
  n_ints[n] = win1[n] + win2[n];
  }
}

parameters{
vector[N_ids] d_raw;                  // latent dominance value parameter
vector[P] Beta;                   // beta coefficients
real<lower=0> sigma;        // standard deviation of dominance values
}

transformed parameters{
vector[N_ids] d;
for(j in 1:N_ids ) {
  d[j] = (X[j] * Beta) + d_raw[j] * sigma;     // non-centered parameterisation
  }
}

model{
vector[N_dyads] p;        // local parameter

Beta ~ normal(0, 1);
sigma ~ cauchy(0, 2);

d_raw ~ normal( 0, 1 );

for( n in 1:N_dyads) {

    // probability of winning = the difference in dominance values between individuals
    // uses the logistic/inverse logit function to transform to probability

  p[n] = inv_logit( d[ ind1[n] ] - d[ ind2[n] ] );
 }

  win1 ~ binomial( n_ints , p );    // vectorised
}
