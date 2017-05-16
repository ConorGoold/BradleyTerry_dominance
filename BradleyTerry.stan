data{
int<lower=0> N_dyads;
int<lower=0> N_ids;
int focal_id;
int<lower=0> ind1[N_dyads];
int<lower=0> ind2[N_dyads];
int<lower=0> win1[N_dyads];
int<lower=0> win2[N_dyads];
}

transformed data{
int<lower=0> n_ints[N_dyads];
for(n in 1:N_dyads) {
n_ints[n] = win1[n] + win2[n];
  }
}

parameters{
vector[N_ids] d;
}

transformed parameters{
vector[N_ids] d_fix;
d_fix[focal_id] = 0.0;
for( n in 1:(focal_id-1) ) d_fix[n] = d[n];
for( n in (focal_id+1):N_ids ) d_fix[n] = d[n];
}

model{
vector[N_dyads] p;

d ~ normal( 0 , 1);

for( n in 1:N_dyads) {
  p[n] = inv_logit( d_fix[ ind1[n] ] - d_fix[ ind2[n] ] );

  win1[n] ~ binomial( n_ints[n] , p[n] );
 }
}
