# run the Wenlock glyphis TMB code
# this has maternal and paternal skip spawning parameters now
#require(TMB)
library(tmbstan) # loads TMB and tmbstan without the f-ups... 
#compile("wlock_ssest.cpp") # ,"-O1 -g",DLLFLAGS="")
dyn.load("wlock_ssest.so")

###################
# set up the data #
###################

input.file <- "glyph_Wenlock_ckmr.dat" 

ckmr.df <- read.table(input.file,header=T)

# code for kin types: up = 0, hsp = 1, fsp = 2

ckmr.df$kcode <- rep(0)
ckmr.df$kcode[ckmr.df$kin.type == 'hsp'] <- 1
ckmr.df$kcode[ckmr.df$kin.type == 'fsp'] <- 2
data <- as.list(ckmr.df)

# tmax: number of cohorts including min(c1,c2)-1

crng <- c(min(c(data$c1,data$c2)),max((c(data$c1,data$c2))))
tmax <- length(crng[1]:crng[2])+1
tzero <- crng[1]-1
data$tzero <- tzero
data$tmax <- tmax
# data$nriv <- if( WENLOCK) 1 else 2 
data$nsex <- 2
data$p_hsp <- 0.92
data$kin.type <- NULL # this wigs out TMB on latest versions 
#data$psi <- c(0,0) # skip spawning factors for females and males

#########################
# set up the parameters #
#########################
pars <- list()


## lambda - log-scale growth rate
pars$lambda <- 0
  
## log(N0) for each river (total males and females) 
pars$lN0 <- log(600)

## iphi - logit-scale survival (s)
pars$iphi <- 2

## leta - log-scale river-specific litter effect
pars$lnu <- log(2)

## izeta - logit-scale sex ratio for each river
pars$izeta <- 0

## itheta - logit-scale multiple paternity parameter
pars$itheta <- -3

pars$lgamma <- 0

pars$ipsi <- log(0.7/(1-0.7))

########################
# create the AD object #
########################

map.obj <- list(
    iphi = factor(NA),
    lambda= factor(NA)
 )
    
obj <- MakeADFun(data=data,parameters=pars,map=map.obj,DLL='wlock_ssest')#, silent=TRUE)

fit0 <- do.call("optim", obj)
obftmp <- obj$fn()
reppo0 <- sdreport(obj)
reppo0.summary <- summary(reppo0)

# 

library(tmbstan)

if ( file.exists("FIT_POSTSTAN.Rdata")) {
    load(
        "FIT_POSTSTAN.Rdata"
    )
} else { 
    stan.hi <- c(log(6000), 3, 1, 5, 5, 5)
    stan.low <- c(0.001, -1, -7, -5, -2, -5)
    fit0.mcmc <-tmbstan( obj, 
        lower = stan.low, 
        upper = stan.hi, silent=TRUE)

    fit0.stan <- as.data.frame(fit0.mcmc)
    save(fit0.stan, file="FIT_POSTSTAN.Rdata")
}
#source("wlock_post_sample_v3.R")
#Nf.post <- post.Nf() 


