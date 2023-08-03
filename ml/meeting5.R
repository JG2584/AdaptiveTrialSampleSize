
#### Simulation Study, speeding up df_temp ####
set.seed(0208)
kSim <- 10000

power_sim <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(power_sim) <- c("SampleSize","Power", "MCSE")
sample_size <- seq(80,90,10)
start_time <- Sys.time()

test <- rep(NA, kSim)
profvis({
  for (i in 1:length(sample_size)){
    for (k in 1:kSim){
      data <- create_rct_data(n = sample_size[i], 1, mu_0=17, mu_1=18, sigma_0=2, sigma_1=2)
      test[k] <- run_test(data)
    }
    df_temp <- data.frame(nSim=1:kSim, samplesize=sample_size[i], test=test)
    current.power <- sum(df_temp$test)/kSim
    mcse <- round(sqrt(current.power*(1-current.power)/kSim), 5)
    power_sim[i, ] <- c(sample_size[i],current.power, mcse)
  }
})

print(power_sim)
end_time <- Sys.time()
print(end_time - start_time)


#### Simulation Study, removing df_temp ####
set.seed(0208)
kSim <- 10000

sample_size <- seq(80,90,10)
power_sim <- data.frame(matrix(ncol = 3, nrow = length(sample_size)))
colnames(power_sim) <- c("SampleSize","Power", "MCSE")
start_time <- Sys.time()

test <- rep(NA, kSim)
profvis({
  for (i in 1:length(sample_size)){
    for (k in 1:kSim){
      data <- create_rct_data(n = sample_size[i], 1, mu_0=17, mu_1=18, sigma_0=2, sigma_1=2)
      test[k] <- run_test(data)
    }
    current.power <- sum(test)/kSim
    mcse <- round(sqrt(current.power*(1-current.power)/kSim), 5)
    power_sim[i, ] <- c(sample_size[i],current.power, mcse)
  }
})



#############################

####### Speed up rand proc 3 #####

proc1 <- function(n){
  group <- sample(rep(c(0,1),n)) 
  group
}

proc2 <- function(n){
  group <- (runif(2*n) > 0.5)*1
  group
}

proc3_2 <- function(n){
  group <- rep(NA, 2*n) # so its length doesn't change
  group[1:n] <- (runif(n) > 0.5)*1
  for (i in (n+1):(2*n)) {
    a <- sum(group == 1, na.rm = TRUE)
    b <- sum(group == 0, na.rm = TRUE)
    if (max(a, b) >= n) {
      group[i:(2*n)] <- c(rep(1, n - a), rep(0, n - b))
      break
    }
    group[i] <- (runif(1) > 0.5)*1
  }
  group
}

create_rct_data <- function(n=10, procedure=3, mu_0=1, mu_1=2, sigma_0=1, sigma_1=1, kSim=5){
  if(procedure==1){
    #group <- replicate(kSim, proc1(n)) # could use for multiple sims
    group <- proc1(n)
  }else if(procedure==2){
    #group <- replicate(kSim, proc2(n)) # could use for multiple sims
    group <- proc2(n)
  }else{
    #group <- replicate(kSim, proc3(n)) # could use for multiple sims
    group <- proc3(n)
  }
  
  df <- data.frame(group=group, outcome=NA)
  df$outcome[df$group==0] <-  rnorm(n=sum(df$group==0), mean=mu_0, sd=sigma_0)
  df$outcome[df$group==1] <-  rnorm(n=sum(df$group==1), mean=mu_1, sd=sigma_1)
  return(df)
}
