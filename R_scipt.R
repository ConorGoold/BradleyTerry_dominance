library(rstsan)
library(rethinking)
library(reshape2)

setwd("~/Documents/PhD_NMBU/BradleyTerry_dominance")

#==========================================================
# Fit the first Bradley-Terry model from Table 3
#==========================================================

d <- read.csv("Adams1.csv")

stan_data_list <- list(N_dyads = nrow(d), N_ids = 5,
                       ind1 = d$ind1, ind2 = d$ind2, 
                       win1 = d$win1, win2 = d$win2, 
                       focal_id = 3)

fit_stan <- stan(file = "BradleyTerry_1.stan", data = stan_data_list, 
                 chains = 4, cores = 4, iter = 5000, warmup = 2500
                )

print(fit_stan)

p_samples <- as.data.frame(fit_stan)    # get the samples

d_values <- p_samples[ , grep("d_fix", colnames(p_samples))]
  
colnames(d_values) <- LETTERS[1:5]       # re-name individuals

# get the order of individuals from each step in the MCMC chain
orderings <- as.vector(apply(d_values, 1, function(z) names(sort(z, decreasing = TRUE))))

# create a data frame from orderings
orderings_df <- data.frame( chain_id = rep(1:nrow(p_samples), each=5), 
                            individual = orderings)

# split the data frame and get the orderings as one vector
orderings_vecs <- lapply(split(orderings_df, orderings_df$chain_id) , 
                         function(z) as.vector(paste(z[,"individual"])) )

# what are the unique orders?
unique_orders <- unique(orderings_vecs)

# find the probabilities of each of the unique orders from the MCMC chains
ordering_probs <- rep(list(list()), length(unique_orders))
for(i in 1:length(unique_orders)) { 
  ordering_probs[[i]] <- sum(unlist( 
    lapply(split(orderings_df, orderings_df$chain_id), 
           function(z) sum(as.vector(z$individual) == unique_orders[[i]]) == 5 ) )
    )/nrow(p_samples)
}

# final data frame of dominance orders and their probability
dominance_orders <- data.frame(order = matrix(unlist(unique_orders), 
                                              nrow = length(unique_orders), byrow = TRUE),
                               probability = unlist(ordering_probs)
                              )
colnames(dominance_orders) <- c(1:5, "probability")
dominance_orders <- dominance_orders[order(dominance_orders$probability, decreasing = TRUE), ]                                
rownames(dominance_orders) <- 1:length(unique_orders)

## NB: The results replicate Adams (2005) partly, but simulation variance in the MCMC chain
## means that the lower probability dominance orders are likely to change. 
## This is not a fully Bayesian approach!
