#+-------------------------------------------------
#
# This script demonstrates model congruence in
# SSE models for diversification
#
#+------------------------------------------------


#   ____________________________________________________________________________
#   Load libraries                                                          ####

# Install rphenoscate packages that contains functions to automatize hidden state expansions
# install.packages("devtools")
# require(devtools)
# install_github("phenoscape/rphenoscape")
# install_github("uyedaj/rphenoscate")
library(rphenoscate)

# source dependencies and install them if they are not
source('R/dependencies.R')

#-------- For parallel computations (to speed up)
library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)
#stopCluster(cl) #stop cluster
#--------

#   ____________________________________________________________________________
#   Simulate Data under BiSSE model                                         ####

set.seed(123)
NSIM=10 # N of trees to simulate


## Generate data
# This script uses much of code from HISSE package vignette
# PARAMETERS: lambda0, lambda1, mu0, mu1, q01, q10
# NOTE!!! the order of q01, q10 is different between BiSSE and HiSSE
pars <- c(0.1, 0.2, 0.03, 0.03, 0.02, 0.01)
phy<- list()

for (i in 1:NSIM){

  bisse <- NULL
  while (is.null(bisse)) {
    bisse <- tree.bisse(pars, max.taxa=100, x0=NA)
  }

  phy[[i]] <- bisse
}

# OR
#saveRDS(phy, file='data/phy_SSE.RDS')
phy <- readRDS(file='data/phy_SSE.RDS')


#   ____________________________________________________________________________
#   Likelihood under the Original BiSSE model                               ####

# constant mu, varibale extinction and trait evolution
Ln.bisse <- c()

for (i in 1:NSIM){
  phy.i <- phy[[i]]
  sim.dat <- data.frame(names(phy.i$tip.state), phy.i$tip.state)

  # PARAMETERS: turnover.0, turnover.1, eps.0, eps.1, q.10, q.01
  lik.hisse <- makeHiSSELikelihood(phy = phy.i, data = sim.dat, hidden.states = FALSE, ode.eps=0, root.p=rep(1/4, 4))
  likf <- lik.hisse$log.lik
  pars <- lik.hisse$pars

  ## Set the parameter values. Note that we have turnover and eps, not speciation and extinction!
  spec <- c(0.1, 0.2)
  extin <- rep(0.03, 2)
  turnover <- spec+extin
  eps <- extin/spec
  qs <- c(0.01, 0.02)
  pars <- setNames(c(turnover, eps, qs), names(pars))

  ## Compute the log-likelihood
  Ln.bisse[[i]] <- likf(pars)
}


#   ____________________________________________________________________________
#   Likelihood under Congruent EHE Model                                     ####

# Now we use a congruent EHE model with hidden states for the traits

# Set-up the model
Q0 <- initQ(c(0,1), c(0.01,0.01))
Qsmm <- amaSMM(Q0, Q0, controlling.state = NULL)

# make Q with equal rates (EHE). when it is lumped with the respect to the observable states it produces
# the original Bisse Q wit two rates
Qsmm[1,4] <- Qsmm[2,3] <- 0.01
Qsmm

# convert to HiSSE ordering of states
v=c(1,3, 2,4)
Qsmm[v,v]
Qsmm <- Qsmm[v,v]

# extract rates
qs <- c(Qsmm)
qs <- qs[qs>=0]

Ln.ehe <- c()
for (i in 1:NSIM){

  phy.i <- phy[[i]]
  sim.dat <- data.frame(names(phy.i$tip.state), phy.i$tip.state)

  # PARAMETERS: turnover.0, turnover.1, eps.0, eps.1, q.10, q.01
  lik.hisse <- makeHiSSELikelihood(phy = phy.i, data = sim.dat, hidden.states = 1, ode.eps=0, root.p=rep(1/4, 4))
  likf <- lik.hisse$log.lik
  pars <- lik.hisse$pars

  ## Set the parameter values. Note that we have turnover and eps, not speciation and extinction!
  spec <- c(0.1, 0.2, 0.1,  0.2)
  extin <- rep(0.03, 4)
  turnover <- spec+extin
  eps <- extin/spec
  pars <- setNames(c(turnover, eps, qs), names(pars))

  ## Compute the log-likelihood
  Ln.ehe[[i]] <- likf(pars)

}




#   ____________________________________________________________________________
#   Likelihood under non-Congruent Model                                     ####

# Now we slightly break the lumpubulity of the congruent EHE to make Nearly Lumpuble
# and to show that likelihoods are different

# Set-up the model
Q0 <- initQ(c(0,1), c(0.01,0.01))
Qsmm <- amaSMM(Q0, Q0, controlling.state = NULL)

# make Q with eaual rates. when it is lumped with the respect to the observable states it produces
# the original Bisse Q wit two rates
Qsmm[1,4] <- 0.01
Qsmm

# convert to HiSSE ordering of states
v=c(1,3, 2,4)
Qsmm[v,v]
Qsmm <- Qsmm[v,v]

# extract rates
qs <- c(Qsmm)
qs <- qs[qs>=0]

Ln.nlehe <- c()
for (i in 1:NSIM){

  phy.i <- phy[[i]]
  sim.dat <- data.frame(names(phy.i$tip.state), phy.i$tip.state)

  # PARAMETERS: turnover.0, turnover.1, eps.0, eps.1, q.10, q.01
  lik.hisse <- makeHiSSELikelihood(phy = phy.i, data = sim.dat, hidden.states = 1, ode.eps=1e-8, root.p=rep(1/4, 4))
  likf <- lik.hisse$log.lik
  pars <- lik.hisse$pars

  ## Set the parameter values. Note that we have turnover and eps, not speciation and extinction!
  spec <- c(0.1, 0.2, 0.1,  0.2)
  extin <- rep(0.03, 4)
  turnover <- spec+extin
  eps <- extin/spec
  pars <- setNames(c(turnover, eps, qs), names(pars))

  ## Compute the log-likelihood
  Ln.nlehe[[i]] <- likf(pars)

}


##  ............................................................................
##  Results: Likelihood with known parameters                               ####


print(Ln.ehe)
print(Ln.bisse)
print(Ln.nlehe)

SSE.ln <- cbind(
  'bisse'=Ln.bisse,
  'EHE'=Ln.ehe,
  'non_congruent'=Ln.nlehe
)

print(SSE.ln )
# the likelihoods are not perfectly the same due to numerical optimization
print(apply(SSE.ln , 2, function(x) x-SSE.ln [,1]))
# root mean squarred error
apply(SSE.ln, 2, function(x) (x-SSE.ln[,1])^2 %>% mean %>% sqrt() )



#   ____________________________________________________________________________
#   Maximum Likelihood Inference                                            ####

# We construct and run two models Q.IND (i.e., CID model), and character-dependent model (Q.COR). These models are
# congruent and have the same number of parameters

# Run yourself or load data
# load("data/congruent_SSE.RData")

##  ............................................................................
##  CID model                                                               ####

# Set-up the model
Qhis <- TransMatMakerHiSSE(hidden.traits=1)

Q0 <- initQ(c(0,1), c(0.01,0.02))
Q1 <- initQ(c(0,1), c(0.02,0.02))
Qsmm <- amaSMM(Q0, Q1, controlling.state = NULL)
Qsmm
Q.IND <- Q2model(Qsmm, diag.as = NA)
Q.IND[row(Q.IND)+col(Q.IND)==5] <- 0
Q.IND
v=c(1,3, 2,4)
Q.IND <- Q.IND[v,v]
colnames(Q.IND) <- rownames(Q.IND) <- colnames(Qhis)
print(Q.IND)
#       (0A) (1A) (0B) (1B)
# (0A)   NA    2    1    0
# (1A)    1   NA    0    1
# (0B)    1    0   NA    2
# (1B)    0    1    1   NA


### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Inference                                                               ####

# Use old HiSSE function, since matrix left diagonal is omitted in the new HiSSE
# Note matrix specification in old HiSSE is different from the new one.

# We use parallel computing to save time

turnover <- c(1,1,2,2)
extinction.fraction <- rep(1, 4)
f = c(1,1)
HiSSE.ind <- list()


HiSSE.ind <- foreach(i=1:NSIM, .packages=c("nloptr", "hisse") ) %dopar% {

  phy.i <- phy[[i]]
  sim.dat <- data.frame(names(phy.i$tip.state), phy.i$tip.state)
  HiSSE.ind[[i]] = hisse.old(phy.i, sim.dat, f=c(1,1), hidden.states=TRUE, turnover.anc=turnover,
                        eps.anc=extinction.fraction, trans.rate=Q.IND, root.p=rep(1/4, 4), sann=TRUE,
                        turnover.upper=1, eps.upper=1, trans.upper=1)
  # HiSSE.ind = hisse.old(phy.i, sim.dat, f=c(1,1), hidden.states=TRUE, turnover.anc=turnover,
  #                       eps.anc=extinction.fraction, trans.rate=Q.IND, root.p=rep(1/4, 4), sann=TRUE)

}


##  ............................................................................
##  Character-dependent model                                               ####

Q.COR <- Qsmm
Q.COR[1,4] <- Q.COR[2,3] <- 0.01
Q.COR
Q.COR <- Q2model(Q.COR, diag.as = NA)
Q.COR[3,2] <- Q.COR[4,1] <- 0
v=c(1,3, 2,4)
Q.COR <- Q.COR[v,v]
colnames(Q.COR) <- rownames(Q.COR) <- colnames(Qhis)
print(Q.COR)
#        (0A) (1A) (0B) (1B)
# (0A)   NA    2    1    2
# (1A)    1   NA    0    1
# (0B)    1    2   NA    2
# (1B)    0    1    1   NA

turnover <- c(1,1,2,2)
extinction.fraction <- rep(1, 4)
f = c(1,1)
HiSSE.cor <- list()

HiSSE.cor <- foreach(i=1:NSIM, .packages=c("nloptr", "hisse") ) %dopar% {

  phy.i <- phy[[i]]
  sim.dat <- data.frame(names(phy.i$tip.state), phy.i$tip.state)
  HiSSE.cor[[i]] = hisse.old(phy.i, sim.dat, f=c(1,1), hidden.states=TRUE, turnover.anc=turnover,
                             eps.anc=extinction.fraction, trans.rate=Q.COR, root.p=rep(1/4, 4), sann=TRUE)
  # HiSSE.cor = hisse.old(phy.i, sim.dat, f=c(1,1), hidden.states=TRUE, turnover.anc=turnover,
  #                       eps.anc=extinction.fraction, trans.rate=Q.COR, root.p=rep(1/4, 4), sann=TRUE)

}


stopCluster(cl)


##  ............................................................................
##  Results: Maximum Likelihood Inference                                   ####


SSE.Mln <- cbind(
'CID'=lapply(HiSSE.ind, function(x) x$loglik) %>% unlist(),
'COR'=lapply(HiSSE.cor, function(x) x$loglik) %>% unlist(),
'CID-COR'=lapply(HiSSE.ind, function(x) x$loglik) %>% unlist()-lapply(HiSSE.cor, function(x) x$loglik) %>% unlist()
)


# the likelihoods are not perfectly the same due to numerical optimization
# In the precooked data
# some are quite different, for example, in rows 1,4,9,10
print(SSE.Mln) %>% round(., 4)

# to prove that these differences are numerical errors
# Let us take the inferred parameters from CID model, convert them into COR model and estimate Ln again

Ln.cor.from.ind <- c()

for (i in 1:NSIM){

  phy.i <- phy[[i]]
  sim.dat <- data.frame(names(phy.i$tip.state), phy.i$tip.state)

  # PARAMETERS: turnover.0, turnover.1, eps.0, eps.1, q.10, q.01
  lik.hisse <- makeHiSSELikelihood(phy = phy.i, data = sim.dat, hidden.states = 1, ode.eps=0, root.p=rep(1/4, 4))
  likf <- lik.hisse$log.lik
  pars <- lik.hisse$pars

  ## Set the parameter values. Note that we have turnover and eps, not speciation and extinction!

  # v=c(1,3, 2,4)
  # Q.COR[v,v]
  # Q.IND[v,v]

  ind <- HiSSE.ind[[i]]
  sol <- ind$solution
  corr.rate <- sol['q0A1A']/2
  names(corr.rate) <- NULL

  sol['q0A1A'] <- corr.rate
  sol['q0A1B'] <- corr.rate
  sol['q0B1A'] <- corr.rate
  sol['q0B1B'] <- corr.rate

  Ln.cor.from.ind[[i]] <- likf(ind$solution)

}


SSE.Mln <- cbind(SSE.Mln, 'CORfromCID'=Ln.cor.from.ind, 'CID-COR2'=SSE.Mln[,1]-Ln.cor.from.ind)
# now the likelihhods between CID model and COR are the same
print(SSE.Mln)

#save.image("data/congruent_SSE.RData")
