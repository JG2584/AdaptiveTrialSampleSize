#BRAR_Example

library(rjags)
library(posterior)
library(dplyr)
#parameters they report
max_pt <- 80
min_pt <- 30
n_first_fr <- 20 #numer of fairly randomized patients
prob_minr <- 0.1 #minimum probability of randomization each arm

# Data initialization
nA <- n_first_fr / 2; nB <- n_first_fr / 2

#Scenario 1
pA <- 0.3; pB <- 0.12
#Scenario 2
pA <- 0.25; pB <- 0.12
#Scenario 3
pA <- 0.12; pB <- 0.12

yA <- sum(rbinom(nA, 1, pA))
yB <- sum(rbinom(nB, 1, pB))

# Formulate Dat 
te_dat <- list(
  yA = yA,
  yB = yB,
  nA = nA,
  nB = nB)

# Model definition
te_jagsmod <- "model {
# Binomial likelihood
  yA ~ dbin(pA,nA)
  yB ~ dbin(pB,nB)
  
  # Prior distributions
    pA ~ dunif(0,1)
    pB ~ dbeta(1.2,8.8)
  
  # Calculated probabilities
    p_supr <- step(pA - pB) #whether arm A is superior than B
    p_ef_gt_0.2 <- step(pA - 0.2) #whether prob_A is larger than 0.2
}"

#JAGS model
model <- jags.model(textConnection(te_jagsmod), 
                    data = te_dat,
                    n.chains = nChains)

update(model,nBurn)
te_post <- coda.samples(model, variable.names = te_params,
                        n.iter = nIter)
te_draws <- as_draws(te_post)
summary(te_draws)

# numbers of chains, burn-in iterations and iterations to keep
nChains <- 2
nBurn <- 200
nIter <- 1000
te_params <- c("pA","pB","p_supr","p_ef_gt_0.2")

# Simulation 
sim_procedure <- function() {
  # Data initialization
  nA <- nB <- n_first_fr / 2
  yA <- rbinom(nA,1,runif(1))
  yB <- rbinom(nB,1,rbeta(1,1.2,8.8))
  
  #JAGS model
  model <- jags.model(textConnection(te_jagsmod), 
                      data = te_dat,
                      n.chains = nChains)
  
  #Adaptive randomization phase
  for (n in (n_first_fr+1):max_pt) {
    update(model,nBurn)
    
    #Posterior Samples
    te_post <- coda.samples(model, variable.names = te_params,
                           n.iter = nIter)
    te_draws <- as_draws(te_post)
    summary(te_draws)
    
    te_pA_post <- as_tibble(as_draws_matrix(te_post), 
                               rownames = "Iteration") %>%
      select("Iteration",starts_with("pA"))  
    te_pB_post <- as_tibble(as_draws_matrix(te_post), 
                            rownames = "Iteration") %>%
      select("Iteration",starts_with("pB"))  
    
    #Check stopping rules
    futility_check <- sapply(te_pA_post, function(x) 
      (sum(x > 0.2) / (nIter*nChains)))
    efficacy_check <- sum(te_pA_post > te_pB_post)/ (nIter*nChains)
    
    if (futility_check < 0.05 ｜ efficacy_check > 0.95) {
      break
    }
    
    #Adaptive Randomization 
    p_randomize_A <- max(prob_min_randomization, 
                         min(1 - prob_min_randomization, 
                             mean(params[, "pA"])))
    
  }
  
}

