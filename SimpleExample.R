# Basic Sample Size Calculation

#### Data Generating Function ####
create_rct_data <- function(n, procedure, mu_0, mu_1, sigma_0, sigma_1) {
  # Randomization procedure 1
  if (procedure == 1){
    group <- sample(rep(c(0,1),n)) 
  }
  
  # Randomization procedure 2: Binomial
  # most random way to generate the sequence:note that this result in non-equal 
  # sample size per arm
  if (procedure == 2){
    group <- (runif(2*n) > 0.5)*1
  }
  
  
  # Randomization procedure 3: Truncated Binomial
  else{
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
  }
  
  
  #outcome <- (1-group) * rnorm(n=n, mean=mu_0, sd=sigma_0) +
  #  group * rnorm(n=n, mean=mu_1, sd=sigma_1)
  outcome <- c()
  for (i in 1: length(group)){
    if (group[i] == 0){
      outcome[i] <- rnorm(1,mean=mu_0, sd=sigma_0)
    }
    else{
      outcome[i] <- rnorm(1,mean=mu_1, sd=sigma_1)
    }
  }
  return(data.frame("group"=group, "outcome"=outcome))
}

# Test data-generating function
data1 <- create_rct_data(n=4, 3, mu_0=2, mu_1=4, sigma_0=0.1, sigma_1=0.1)

# Test whether to reject null-hypothesis using z-test
run_test <- function(data) {
  #test_result <- t.test(outcome~group, data=data)
  #return(as.integer(test_result$p.value<0.05))
  
  group_means <- tapply(data$outcome, data$group, mean)
  group_sds <- tapply(data$outcome, data$group, sd)
  n1 <- sum(data$group == levels(as.factor(data$group))[1])
  n2 <- sum(data$group == levels(as.factor(data$group))[2])
  
  z_score <- (group_means[1] - group_means[2]) / sqrt((group_sds[1]^2 / n1) + (group_sds[2]^2 / n2))
  p_value <- 2 * (1 - pnorm(abs(z_score)))
  
  return(as.integer(p_value < 0.05))
}

#### Simulation Study ####
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
    data <- create_rct_data(n = sample_size[i], 3, mu_0=17, mu_1=18, sigma_0=2, sigma_1=2)
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

#### Check with the formula ####
power_formula <- function (n){
  pnorm(sqrt((n*(17-18)^2)/(2^2+2^2)) - qnorm(0.025, lower.tail=F))
}
power_fm <- sapply(sample_size, power_formula)

n_formula <- function (pw){
  (qnorm(0.025, lower.tail=F)+qnorm((1-pw), lower.tail = F))^2*(2^2+2^2)/((17-18)^2)
}
n_fm <- sapply(power_sim$Power, n_formula)
n_fm

library(ggplot2)
# x-sample size y-power plot
ggplot(data.frame(
  n = rep(sample_size, 2),
  power = c(power_sim$Power, power_fm),
  which = rep(c("Simulation","Formula"), each=length(sample_size))
), aes(x=n, y=power, color=factor(which))) +
  geom_line() +
  labs(color="Method", y="Power", x="Sample size (per group)")

# x-power y-sample size plot
ggplot(data.frame(
  power = rep(power_sim$Power, 2),
  n = c(power_sim$SampleSize, n_fm),
  which = rep(c("Simulation","Formula"), each=length(sample_size))
), aes(y=n, x=power, color=factor(which))) +
  geom_line() +
  labs(color="Method", x="Power", y="Sample size (per group)")

#### Package the simulation study in a function ####

# Input: number of simulation study, mean of control group, mean of experimental group, 
# variance of control group and experimental group, minimized sample size, maximized sample size and 
# the threshold of power

# Output: List of two elements: Table of sample size with its power, the minimum
# sample size that is over the threshold
sim_cpmean <- function(kSim,rp,step,mu_0,mu_1, sigma_0,sigma_1,min_sz,max_sz,threshold){
  #set.seed(43)
  power_sim <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(power_sim) <- c("SampleSize","Power")
  # Consider change it to bisection method
  sample_size <- seq(min_sz,max_sz,step)
  start_time <- Sys.time()
  for (i in 1:length(sample_size)){
    df_temp <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(df_temp) <- c("nSim","sampleSize","test")
    #print(df_temp)
    for (k in 1:kSim){
      data <- create_rct_data(n = sample_size[i], rp,mu_0=mu_0, mu_1=mu_1, sigma_0=sigma_0, sigma_1=sigma_1)
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
  if (length(elgb_sz)==0){
    elgb_sz <- c(max_sz)
    print("No eligible sample size is found to exceed the threshold!")
  }
  end_time <- Sys.time()
  print(end_time - start_time)
  return(list(power_sim,min(elgb_sz)))
}
#print(sim_cpmean(100,10,17,18,2,2,100,150,0.9))

#### Metrics ####
# Variability of n
min_n <- c()
for (i in 1:100){
  min_n <- c(min_n,sim_cpmean(1000,2,1,17,18,2,2,95,105,0.9)[[2]])
}
print(min_n)
print(var(min_n))
print(sd(min_n))

##### Sensitivity Analysis #####
kEvl = 100;kSim = 1000
# Changing randomization procedure

# Function to calculate min_nr for a specific randomization procedure
calc_min_nr <- function(rp) {
  min_nr <- c()
  for (i in 1:kEvl) {
    min_nr <- c(min_nr, sim_cpmean(kSim = kSim, rp = rp, step = 1,
                                   17, 18, 2, 2,
                                   min_sz = 95, max_sz = 105, threshold = 0.9)[[2]])}
  return(min_nr)
}

min_nr1 <- calc_min_nr(rp = 1)
min_nr2 <- calc_min_nr(rp = 2)
min_nr3 <- calc_min_nr(rp = 3)
min_nr <- list(min_nr1,min_nr2,min_nr3)
var_nr <- c(var(min_nr1),var(min_nr2),var(min_nr3))

print(var_nr)
freq_table1 <- table(min_nr3)
max_frequency1 <- max(freq_table1)
most_frequent_elements <- names(freq_table1)[freq_table1== max_frequency1]

cat("Most frequently occurring element(s): ", most_frequent_elements, "\n")
print(c(abs(n_formula(0.9)-mean(min_nr1)),abs(n_formula(0.9)-mean(min_nr2)),
        abs(n_formula(0.9)-mean(min_nr3))))

# Histogram
par(mfrow = c(3,1))
for (i in 1:3){
  hist(min_nr[[i]], main = paste("Histogram of n for Randomization Procedure " , i),
       xlab = "Minimum Sample Size")
}

# Number of simulation study
calc_min_nk <- function(kSim) {
  min_nk <- c()
  for (i in 1:kEvl) {
    min_nk <- c(min_nk, sim_cpmean(kSim = kSim, rp = 3, step = 1,
                                   17, 18, 2, 2,
                                   min_sz = 95, max_sz = 105, threshold = 0.9)[[2]])}
  return(min_nk)
}

kSims <- c(100,500,1000)
min_nk <- list();var_nk <- c();abs_nk <- c()

for (i in 1:length(kSims)) {
   min_nk0 <- calc_min_nk(kSims[i])
   min_nk[[i]] <- min_nk0
   var_nk <- c(var_nk,var(min_nk0))
   abs_nk <- c(abs_nk,abs(n_formula(0.9)-mean(min_nk0)))
}

print(abs_nk)
par(mfrow = c(3,1))
for (i in 1:length(kSims)){
  hist(min_nk[[i]], main = paste("Histogram of n for" , kSims[i], "of simulation replicates"),
       xlab = "Minimum Sample Size")
}

par(mfrow = c(1,1))
ggplot(data.frame(
  kSim = kSims,
  var = var_nk), aes(y=var_nk, x=kSim)) +
    geom_line() +
    labs(x="Number of Simulation Study", y="Variability of n")

# Threshold
# Note that at first we set step larger to see a small range of sample size and
# then set it smaller as 1 to find the more accurate estimation and run the simulation
min_nt1 <- c();min_nt2 <- c();min_nt3 <- c()
thresholds <- c(0.8,0.9,0.95)
for (k in 1:kEvl){
    min_nt1 <- c(min_nt1,sim_cpmean(kSim = kSim,rp = 3,step = 1,
                                    17,18,2,2,
                                    min_sz = 70,max_sz = 80,threshold = thresholds[1])[[2]])
    min_nt2 <- c(min_nt2,sim_cpmean(kSim = kSim,rp = 3,step = 1,
                                    17,18,2,2,
                                    min_sz = 95,max_sz = 110,threshold = thresholds[2])[[2]])
    min_nt3 <- c(min_nt3,sim_cpmean(kSim = kSim,rp = 3,step = 1,
                                    17,18,2,2,
                                    min_sz = 120,max_sz = 130,threshold = thresholds[3])[[2]])
  }

# for (k in 1:kEvl){
#   min_nt3 <- c(min_nt3,sim_cpmean(kSim = kSim,rp = 3,step = 1,
#                                   17,18,2,2,
#                                   min_sz = 120,max_sz = 130,threshold = thresholds[3])[[2]])
# }
var_nt <- c(var(min_nt1),var(min_nt2),var(min_nt3))
print(c(abs(n_formula(0.8)-mean(min_nt1)),abs(n_formula(0.9)-mean(min_nt2)),
        abs(n_formula(0.95)-mean(min_nt3))))

par(mfrow = c(3,1))
min_nt <- list(min_nt1,min_nt2,min_nt3)
for (i in 1:3){
  hist(min_nt[[i]], main = paste("Histogram of Power Threshold ", thresholds[i]),
       xlab = "Minimum Sample Size")
}

par(mfrow = c(1,1))
ggplot(data.frame(
  threshold = c(0.8,0.9,0.95),
  var = var_nt), aes(y=var, x=threshold))+
    geom_line() +
    labs(x="Power Threshold", y="Variability of n")

save(min_nr,min_nk,min_nt, 
    file = "/Users/beep/Downloads/申研/Cambridge/MRC BSU Internship/results2.Rdata")

### Retetet
