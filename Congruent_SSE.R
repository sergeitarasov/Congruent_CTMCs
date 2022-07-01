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

# source dependecies and install them if they are not
source('R/dependencies.R')


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
phy
#length(phy)
#plot(phy[[2]])

saveRDS(phy, file='data/phy_SSE.RDS')
phy <- readRDS(file='data/phy_SSE.RDS')

#   ____________________________________________________________________________
#   Likelihood under BiSSE (Original model)                                  ####

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

Ln.bisse

