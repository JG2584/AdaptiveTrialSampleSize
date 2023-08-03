library(librarian)
librarian::shelf(profvis, microbenchmark)

profvis({
min_n <- c(min_n,sim_cpmean(1000,2,1,17,18,2,2,95,105,0.9)[[2]])
})


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
  
  outcome <- (1-group) * rnorm(n=n, mean=mu_0, sd=sigma_0) +
    group * rnorm(n=n, mean=mu_1, sd=sigma_1)
  return(data.frame("group"=group, "outcome"=outcome))
}

create_rct_data(n=10, procedure=3, mu_0=1, mu_1=2, sigma_0=1, sigma_1=1)




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


create_rct_data_2 <- function(n=10, procedure=3, mu_0=1, mu_1=2, sigma_0=1, sigma_1=1, kSim=5){
  if(procedure==1){
    #group <- replicate(kSim, proc1(n)) # use for multiple sims
    group <- proc1(n)
  }else if(procedure==2){
    #group <- replicate(kSim, proc2(n)) # use for multiple sims
    group <- proc2(n)
  }else{
    #group <- replicate(kSim, proc3(n)) # use for multiple sims
    group <- proc3(n)
  }
  
  df <- data.frame(group=group, outcome=NA)
  df$outcome[df$group==0] <-  rnorm(n=sum(df$group==0), mean=mu_0, sd=sigma_0)
  df$outcome[df$group==1] <-  rnorm(n=sum(df$group==1), mean=mu_1, sd=sigma_1)
  return(df)
}
create_rct_data_2()

### Interesting line:
outcome <- (1-group) * rnorm(n=n, mean=mu_0, sd=sigma_0) +
    group * rnorm(n=n, mean=mu_1, sd=sigma_1)
# ^^^ Issues:
# The first rnorm(n) gets recycled and used twice. This also happens with the
# second rnorm(n) call: there are 2n patients in the object "group".
# also
# There will be further recycling for whichever group has more than n patients.
# simple test:
p1 <- create_rct_data(n=10000, procedure=1, mu_0=0, mu_1=10, sigma_0=1, sigma_1=1)
tapply(p1$outcome, p1$group, summary)
p1X <- create_rct_data_2(n=100000, procedure=1, mu_0=0, mu_1=10, sigma_0=1, sigma_1=1)
tapply(p1X$outcome, p1X$group, summary)

p1$code <- "original"
p1X$code <- "modified"
p1.df <- rbind(p1, p1X)

ggplot(data=p1.df, mapping=aes(col=code))+
  geom_histogram(mapping=aes(x=outcome))+
  facet_wrap("group")


nsims <- 1000
orig.outp <- lapply(X=rep(25,nsims), FUN=create_rct_data, procedure=1, mu_0=0, mu_1=10, sigma_0=2, sigma_1=4)
mod.outp <-  lapply(X=rep(25,nsims), FUN=create_rct_data_2, procedure=1, mu_0=0, mu_1=10, sigma_0=2, sigma_1=4)

orig.var0 <- rep(NA, nsims)
mod.var0 <- rep(NA, nsims)
orig.var1 <- rep(NA, nsims)
mod.var1 <- rep(NA, nsims)

for(i in 1:nsims){
  orig.var0[i] <- var(orig.outp[[i]]$outcome[orig.outp[[i]]$group==0])
  mod.var0[i] <- var(mod.outp[[i]]$outcome[mod.outp[[i]]$group==0])
  orig.var1[i] <- var(orig.outp[[i]]$outcome[orig.outp[[i]]$group==1])
  mod.var1[i] <- var(mod.outp[[i]]$outcome[mod.outp[[i]]$group==1])
}

summary(orig.var0)
summary(mod.var0)
summary(orig.var1)
summary(mod.var1)

variances <- data.frame(orig0=orig.var0, orig1=orig.var1, mod0=mod.var0, mod1=mod.var1)

par(mfrow=c(1,2))
boxplot(variances$orig0, variances$mod0)
boxplot(variances$orig1, variances$mod1)



# consider:
n <- 5
group <- c(1,1,0,1,0,1,0,0,1,0) 
(1-group) * rnorm(n, mean=10) + group * rnorm(n, mean=20)

create_rct_data(n=10, procedure=3, mu_0=1, mu_1=2, sigma_0=1, sigma_1=1)
create_rct_data_2(n=10, procedure=2, mu_0=1, mu_1=2, sigma_0=1, sigma_1=1)


# Speed no different
microbenchmark::microbenchmark(create_rct_data(n=100, procedure=3, mu_0=1, mu_1=2, sigma_0=1, sigma_1=1),
                               create_rct_data_2(n=100, procedure=3, mu_0=1, mu_1=2, sigma_0=1, sigma_1=1),
                                times=1000)



######### NOTE: the code below will take forever to run, I think.
 outcome <- c()
  for (i in 1: length(group)){
    if (group[i] == 0){
      outcome[i] <- rnorm(1,mean=mu_0, sd=sigma_0)
    }
    else{
      outcome[i] <- rnorm(1,mean=mu_1, sd=sigma_1)
    }
  }
 
 