library(rstsan)
library(rethinking)
library(ggplot2)

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
## This is not really a fully Bayesian approach. 

#==========================================================
# Fit the second Bradley-Terry model, using data from Table 2.

# This data set is of wins and losses between 20 female bighorn sheep, where the
# interest is in relating dominance to age.
#==========================================================

d2 <- read.csv("Adams2.csv")
d2_ages <- read.csv("Adams2_ages.csv")

d2$ind1 <- rep(1:length(unique(d2$ind1)),
                     sort(unlist(lapply(split(d2, d2$ind1), 
                         function(z) nrow(z)) ) , decreasing = TRUE))

ind2_new <- rep(list(list()), length(unique(d2$ind2)))
for(i in 1:19){
  ind2_new[[i]] <- c(unique(d2$ind1),20)[-c(1:i)]
}

d2$ind2 <- unlist(ind2_new)

d2_ages$ageZ <- (d2_ages$age - mean(d2_ages$age))/sd(d2_ages$age)
X_mat = as.matrix(d2_ages$ageZ)

stan_data_list <- list(N_dyads = nrow(d2), N_ids = length(unique(d2$ind1))+1,
                       ind1 = d2$ind1, ind2 = d2$ind2, 
                       win1 = d2$win1, win2 = d2$win2, 
                       X = X_mat, P = ncol(X_mat) ) 

fit_stan2 <- stan(file = "BradleyTerry_withCovariate.stan", 
                  data = stan_data_list, 
                  chains = 4, cores = 4, iter = 5000, warmup = 2500
)

print(fit_stan2)

# NB: Stan tells us that there are a number of divergent transitions, indicating problems
## in the convergence of the MCMC chains. Using the non-centered parameterisation doesn't get rid 
## of the divergent transitions. This may occur because of the large number of 'missing' dyads, where
## no interactions are observed. 

p_samples <- as.data.frame(fit_stan2)

#=========== Re-create figure 2 ===============================================

age_comb <- d2_ages$ageZ
pred <- sapply(age_comb, 
               function(z) p_samples$`Beta[1]` * z )
pred_mu <- apply(pred, 2, mean)
pred_hdi <- apply(pred, 2, function(z) HPDI(z,0.89) )

plotting_df <- data.frame(id = 1:20, 
                          age = d2_ages$age,
                          ageZ = d2_ages$ageZ,
                          mu = pred_mu, 
                          hdi_low = pred_hdi[1,], hdi_high = pred_hdi[2,],
                          dom_mu = as.vector(apply(p_samples[,grep("d",colnames(p_samples))][21:40],
                                                   2, mean ) ),
                          dom_low = as.vector(apply(p_samples[,grep("d",colnames(p_samples))][21:40],
                                                   2, function(z) HPDI(z) )[1,] ),
                          dom_high = as.vector(apply(p_samples[,grep("d",colnames(p_samples))][21:40],
                                                    2, function(z) HPDI(z) )[2,] )
                          )

ggplot(plotting_df, aes(age, mu)) + 
  geom_line() + 
  geom_point(aes(age, dom_mu, group=id), position = position_dodge(0.5)) + 
  geom_errorbar(aes(ymin=dom_low, ymax=dom_high, group=id), 
                width=0, position = position_dodge(0.5)) + 
  labs(x = "age", y="dominance value") + 
  theme_bw() + 
  theme(panel.grid = element_blank())


# NB: The dominance values and regression line seem to lie over a smaller range than 
## in Adams (2005). This is probably (not tested) a result of using regularising priors, rather than the vague priors used in BUGS. 

