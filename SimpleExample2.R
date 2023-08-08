# Basic Sample Size Calculation

## Randomization Procedures ##

proc1 <- function(n){
  group <- sample(rep(c(0,1),n)) 
  group
}

proc2 <- function(n){
  group <- (runif(2*n) > 0.5)*1
  group
}

proc3 <- function(n){
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
  group
}

#### Data Generating Function ####
create_rct_data <- function(n, procedure, kSim){
  if(procedure==1){
    group <- replicate(kSim, proc1(n)) # could use for multiple sims
    #group <- proc1(n)
  }else if(procedure==2){
    group <- replicate(kSim, proc2(n)) # could use for multiple sims
    #group <- proc2(n)
  }else{
    group <- replicate(kSim, proc3(n)) # could use for multiple sims
    #group <- proc3(n)
  }
  
  outcome <- apply(group,2,function(x) {
    new_col <- ifelse(x==0, rnorm(n=sum(x==0), mean=mu_0, sd=sigma_0),NA)
    new_col <- ifelse(is.na(new_col),rnorm(n=sum(x==1), mean=mu_1, sd=sigma_1),new_col)
    return(new_col)
  })
  return(list(group,outcome))
}

#### Test Function using z-test ####

counter <- 0 #have to set counter = 0 every time before running run-test2, which might be a problem
run_test <- function(group,outcome){
  apply(outcome,2,function(col){
    counter <<- counter+1
    # print(counter)
    grp0 <- group[,counter] == 0
    grp1 <- !grp0
    group0_mean <- mean(col[grp0])
    group1_mean <- mean(col[grp1])
    group0_sd <- sd(col[grp0])
    group1_sd <- sd(col[grp1])
    n1 <- sum(grp0)
    n2 <- sum(grp1)
    
    z_score <- (group0_mean - group1_mean) / sqrt((group0_sd^2 / n1) + (group1_sd^2 / n2))
    p_value <- 2 * (1 - pnorm(abs(z_score)))
    return(as.integer(p_value < 0.05))
  })
}

#### Simulation study: General trend of different sizes ####

set.seed(43)
#profvis({
kSim = 1000
power_sim <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(power_sim) <- c("SampleSize","Power", "MCSE")

sample_size <- seq(10,150,20)
start_time <- Sys.time()

for (i in 1:length(sample_size)){
  dl <- create_rct_data(sample_size[i],3,kSim)
  counter <- 0 #have to set counter = 0 every time before running run-test2
  result <- run_test(dl[[1]],dl[[2]])
  current.power <-  sum(result)/kSim
  mcse <- round(sqrt(current.power*(1-current.power)/kSim), 5)
  power_sim[i, ] <- c(sample_size[i],current.power, mcse)
}
print(power_sim)
end_time <- Sys.time()
print(end_time - start_time)
#})

##### Check with the formula #####

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
# x-power y-sample size plot
ggplot(data.frame(
  power = rep(power_sim$Power, 2),
  n = c(power_sim$SampleSize, n_fm),
  which = rep(c("Simulation","Formula"), each=length(sample_size))
), aes(y=n, x=power, color=factor(which))) +
  geom_line() +
  labs(color="Method", x="Power", y="Sample size (per group)")

#### Package the simulation study in a function ####

# Input: number of simulation study, randomization procedure, step between min_n 
# and max_n, min_n, max_n and the threshold of power

# Output: List of two elements: Table of sample size with its power and MC error, 
# the minimum sample size that is over the threshold

sim_cpmean <- function(kSim,rp,step,min_sz,max_sz,threshold){
  
  power_sim <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(power_sim) <- c("SampleSize","Power", "MCSE")
  sample_size <- seq(min_sz,max_sz,step)
  start_time <- Sys.time()
  
  for (i in 1:length(sample_size)){
    dt <- create_rct_data(sample_size[i],3,kSim)
    counter <- 0 #have to set counter = 0 every time before running run-test2
    result <- run_test(dt[[1]],dt[[2]])
    current.power <-  sum(result)/kSim
    mcse <- round(sqrt(current.power*(1-current.power)/kSim), 5)
    power_sim[i, ] <- c(sample_size[i],current.power, mcse)
  }
  
  print(power_sim)
  end_time <- Sys.time()
  print(end_time - start_time)
  
  #select eligible sample size attaining the threshold
  elgb_sz <- c()
  for (row in 1:dim(power_sim)[1]){
    if (power_sim[row,2] >= threshold){
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

sim_cpmean(1000,3,1,80,90,0.9)

#### Metrics ####
# Variability of n
min_n <- c()
for (i in 1:100){
  min_n <- c(min_n,sim_cpmean(1000,3,1,80,90,0.9)[[2]])
}
print(min_n)
print(var(min_n))
print(sd(min_n))

##### Sensitivity Analysis #####
kEvl = 100;kSim = 1000;mu_0=17; mu_1=18; sigma_0=2; sigma_1=2

###### Randomization procedure ######

# Function to calculate min_nr for a specific randomization procedure
calc_min_nr <- function(rp) {
  min_nr <- c()
  for (i in 1:kEvl) {
    min_nr <- c(min_nr, sim_cpmean(kSim = kSim, rp = rp, step = 1,
                                   min_sz = 80, max_sz = 90, threshold = 0.9)[[2]])}
  return(min_nr)
}

min_nr1 <- calc_min_nr(rp = 1)
min_nr2 <- calc_min_nr(rp = 2)
min_nr3 <- calc_min_nr(rp = 3)
min_nr <- list(min_nr1,min_nr2,min_nr3)
var_nr <- c(var(min_nr1),var(min_nr2),var(min_nr3))

print(var_nr)

# Histogram
par(mfrow = c(3,1))
for (i in 1:3){
  hist(min_nr[[i]], main = paste("Histogram of n for Randomization Procedure " , i),
       xlab = "Minimum Sample Size")
}

###### Number of simulation study ######
calc_min_nk <- function(kSim) {
  min_nk <- c()
  for (i in 1:kEvl) {
    min_nk <- c(min_nk, sim_cpmean(kSim = kSim, rp = 3, step = 1,
                                   min_sz = 80, max_sz = 90, threshold = 0.9)[[2]])}
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
  hist(min_nk[[i]], main = paste("Histogram of n for" , kSims[i], "simulation replicates"),
       xlab = "Minimum Sample Size")
}

par(mfrow = c(1,1))
ggplot(data.frame(
  kSim = kSims,
  var = var_nk), aes(y=var_nk, x=kSim)) +
  geom_line() +
  labs(x="Number of Simulation Study", y="Variability of n")

###### Power Threshold ######
# Note that at first we set step larger to see a small range of sample size and
# then set it smaller as 1 to find the more accurate estimation and run the simulation
min_nt1 <- c();min_nt2 <- c();min_nt3 <- c()
thresholds <- c(0.8,0.9,0.95)
for (k in 1:kEvl){
  min_nt1 <- c(min_nt1,sim_cpmean(kSim = kSim,rp = 3,step = 1,
                                  min_sz = 55,max_sz = 70,threshold = thresholds[1])[[2]])
  min_nt2 <- c(min_nt2,sim_cpmean(kSim = kSim,rp = 3,step = 1,
                                  min_sz = 80,max_sz = 90,threshold = thresholds[2])[[2]])
  min_nt3 <- c(min_nt3,sim_cpmean(kSim = kSim,rp = 3,step = 1,
                                  min_sz = 95,max_sz = 110,threshold = thresholds[3])[[2]])
}

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

#save(min_nr,min_nk,min_nt, 
#     file = "/Users/beep/Downloads/申研/Cambridge/MRC BSU Internship/results_corrected.Rdata")