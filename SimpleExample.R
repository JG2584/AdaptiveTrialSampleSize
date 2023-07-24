# Basic Sample Size Calculation

# Data Generating Function
create_rct_data <- function(n, mu_0, mu_1, sigma_0, sigma_1) {
  # Randomization procedure 1
  group <- sample(rep(c(0,1),n)) 
  
  # Randomization procedure 2: Binomial
  # most random way to generate the sequence:note that this result in non-equal 
  # sample size per arm
  group <- (runif(2*n) > 0.5)*1
  
  # Randomization procedure 3: Truncated Binomial
  group <- c()
  for (i in 1:(2*n)) {
    a <- sum(group == 1)
    b <- sum(group == 0)
    if (max(a, b) >= n) {
      group <- c(group, rep(1, n - a), rep(0, n - b))
      break
    }
    group <- c(group, (runif(1) > 0.5)*1)
  }
  
  outcome <- (1-group) * rnorm(n=n, mean=mu_0, sd=sigma_0) +
    group * rnorm(n=n, mean=mu_1, sd=sigma_1)
  return(data.frame("group"=group, "outcome"=outcome))
}

# Test data-generating function
create_rct_data(n=3, mu_0=3, mu_1=4, sigma_0=0.1, sigma_1=0.1)

# Test whether to reject null-hypothesis
run_test <- function(data) {
  test_result <- t.test(outcome~group, data=data)
  return(as.integer(test_result$p.value<0.05))
}

# Simulation Study
set.seed(43)
kSim <- 1000
power_sim <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(power_sim) <- c("SampleSize","Power")
sample_size <- seq(10,150,20)
start_time <- Sys.time()
for (i in 1:length(sample_size)){
  df_temp <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(df_temp) <- c("nSim","sampleSize","test")
  #print(df_temp)
  for (k in 1:kSim){
    data <- create_rct_data(n = sample_size[i], mu_0=17, mu_1=18, sigma_0=2, sigma_1=2)
    result <- run_test(data)
    # Do not use rbind, as it changes the original column names of the dataframe
    df_temp[k,] <- c(nSim = k,sampleSize = sample_size[i],test = result)
    #print(df_temp)
  }
  power_sim[i,] <- c(sample_size[i],sum(df_temp$test)/kSim)
}
print(power_sim)
end_time <- Sys.time()
print(end_time - start_time)

# Check with the formula
power_formula <- sapply(sample_size, function (n){
  pnorm(sqrt((n*(17-18)^2)/(2^2+2^2)) - qnorm(0.025, lower.tail=F))
})

library(ggplot2)
ggplot(data.frame(
  n = rep(sample_size, 2),
  power = c(power_sim$Power, power_formula),
  which = rep(c("Simulation","Formula"), each=length(sample_size))
), aes(x=n, y=power, color=factor(which))) +
  geom_line() +
  labs(color="Method", y="Power", x="Sample size (per group)")

## Using SimEngine Package
library(SimEngine)
sim <- new_sim()
sim %<>% set_script(function() {
  data <- create_rct_data(n=L$n, mu_0=17, mu_1=18, sigma_0=2, sigma_1=2)
  reject <- run_test(data)
  return (list("reject"=reject))
})

sim %<>% set_levels(n=c(50,80,100,150))
sim %<>% set_config(num_sim=1000)
sim %<>% run()

power_sim <- sim %>% summarize(
  list(stat="mean", name="power", x="reject")
)
print(power_sim)

## Package the simulation study in a function

# Input: number of simulation study, mean of control group, mean of experimental group, 
# variance of control group and experimental group, minimized sample size, maximized sample size and 
# the threshold of power

# Output: List of two elements: Table of sample size with its power, the minimum
# sample size that is over the threshold
sim_cpmean <- function(kSim,mu_0,mu_1, sigma_0,sigma_1,min_sz,max_sz,threshold){
  set.seed(43)
  power_sim <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(power_sim) <- c("SampleSize","Power")
  # Consider change it to bisection method
  sample_size <- seq(min_sz,max_sz,10)
  start_time <- Sys.time()
  for (i in 1:length(sample_size)){
    df_temp <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(df_temp) <- c("nSim","sampleSize","test")
    #print(df_temp)
    for (k in 1:kSim){
      data <- create_rct_data(n = sample_size[i], mu_0=mu_0, mu_1=mu_1, sigma_0=sigma_0, sigma_1=sigma_1)
      result <- run_test(data)
      # Do not use rbind, as it changes the original column names of the dataframe
      df_temp[k,] <- c(nSim = k,sampleSize = sample_size[i],test = result)
      #print(df_temp)
    }
    power_sim[i,] <- c(sample_size[i],sum(df_temp$test)/kSim)
  }
  elgb_sz <- c()
  for (row in 1:dim(power_sim)[1]){
    if (power_sim[row,2] > threshold){
      elgb_sz <- c(elgb_sz,power_sim[row,1])
    }
  }
  end_time <- Sys.time()
  print(end_time - start_time)
  return(list(power_sim,min(elgb_sz)))
}
print(sim_cpmean(100,17,18,2,2,100,150,0.9))
