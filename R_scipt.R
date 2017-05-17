library(rstsan)
library(rethinking)

setwd("~/Documents/PhD_NMBU/BradleyTerry_dominance")

d <- read.csv("Adams1.csv")

stan_data_list <- list(N_dyads = nrow(d), N_ids = 5,
                       ind1 = d$ind1, ind2 = d$ind2, 
                       win1 = d$win1, win2 = d$win2, 
                       focal_id = 3)

fit_stan <- stan(file = "BradleyTerry.stan", data = stan_data_list, 
                 chains = 4, cores = 4, iter = 2000, warmup = 1000
                )

print(fit_stan)

p_samples <- as.data.frame(fit_stan)

d_values <- p_samples[ , grep("d_fix", colnames(p_samples))]
colnames(d_values) <- 1:5
orderings <- as.vector(apply(d_values, 1, function(z) names(sort(z))))
orderings <- data.frame( chain_id = rep(1:nrow(p_samples), each=5), 
                            individual = orderings)
orderings_list <- lapply(split(orderings, orderings$chain_id), 
                         function(z) all(as.vector(z$individual) == 1:5) == TRUE ) 
possible_orders <- list( 1:5, c(1,2,3,5,4), c(2,1,3,4,5), c(1,2,4,3,5), c(2,1,3,5,4) ,
                      c(2,1,4,3,5), c(1,2,4,5,3), c(1,2,5,3,4))


