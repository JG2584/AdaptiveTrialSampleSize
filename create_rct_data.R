
create_rct_data <- function(n, procedure, kSim, mu_0, sigma_0, mu_1, sigma_1){
  if(procedure==1){
    #group <- replicate(kSim, proc1(n)) # could use for multiple sims
    group <- replicate(kSim, proc1(n), FALSE) # could use for multiple sims
    #group <- proc1(n)
  }else if(procedure==2){
    group <- replicate(kSim, proc2(n), FALSE) # could use for multiple sims
    #group <- proc2(n)
  }else{
    group <- replicate(kSim, proc3(n), FALSE) # could use for multiple sims
    #group <- proc3(n)
  }
  
outcome.list <- vector("list", kSim)
for(i in 1:kSim){
  outcome.vec <- rep(NA, 2*n)
  single.group <- group[[i]]
  outcome.vec[single.group==0] <-  rnorm(n=sum(single.group==0), mean=mu_0, sd=sigma_0)
  outcome.vec[single.group==1] <-  rnorm(n=sum(single.group==1), mean=mu_1, sd=sigma_1)
  outcome.list[[i]] <- outcome.vec
}

output <- Map(cbind, group, outcome.list)
output

}

