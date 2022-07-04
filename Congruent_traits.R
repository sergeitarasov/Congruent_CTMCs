#+-------------------------------------------------
#
# This script demonstrates model congruence in
# CTMCs for discrete characters on known phylogeny
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
# slightly modified version of rayDISC() function from corHMM package
# with disabled ACE reconstruction to save computational time
source('R/mod_ray_disc.R')


#   ____________________________________________________________________________
#   Simulate Data                                                           ####

# The analyses may take a while if you are not willing to wait proceed
# directly to Results at the bottom.

#Simulate a tree and ten replicates of a two-state character
set.seed(123)
tree<-pbtree(n=200, scale=100, b=1, d=0)
Q <- initQ(c(1,2), c(.3,.2))
print(Q)
hist <- sim.history(tree, Q, nsim=10)



#   ____________________________________________________________________________
#   ML inference under the original model                                   ####

Q.par <- Q2model(Q, diag.as = NA)
print(Q.par)
#   1  2
# 1 NA  2
# 2  1 NA

Q.out <- list()
for (i in 1:10){
  taxa <- cbind(hist[[i]]$tip.label, hist[[i]]$states)
  Q.out[[i]] <- mod_rayDISC(hist[[i]], taxa, rate.mat=Q.par, root.p=c(0.5, 0.5), verbose = TRUE)
}


#   ____________________________________________________________________________
#   ML inference under Equivalent Expansion Model                           ####

# We make a three-state equivalent model where state 1 has two hidden states
# This Equivalent Expansion Model has the same number of parameters

M.eq <- initQ(c(1:3), c(0))
M.eq[1,3] <- M.eq[2,3] <- 2
M.eq[3,1] <- M.eq[3,2] <- 1
M.eq <- Q2model(M.eq, diag.as = NA)
print(M.eq)
#    1  2  3
# 1 NA NA  2
# 2 NA NA  2
# 3  1  1 NA

M.eq.out <- list()
for (i in 1:10){
  states<- mapvalues(hist[[i]]$states, from = c("1", "2"), to=c("1&2", "3") )
  taxa<-cbind(names(states), states)
  root <- c(1/4, 1/4, 1/2)
  sum(root)
  M.eq.out[[i]] <- mod_rayDISC(hist[[i]], taxa, rate.mat=M.eq, root.p=root, verbose = TRUE)
}

#   ____________________________________________________________________________
#   ML inference under Equivalent nWTP  Model                              ####


M.eq.nwtp <- initQ(c(1:3), c(3))
M.eq.nwtp[1,3] <- M.eq.nwtp[2,3] <- 2
M.eq.nwtp[3,1] <- M.eq.nwtp[3,2] <- 1

M.eq.nwtp <- Q2model(M.eq.nwtp/.1, diag.as = NA)
print(M.eq.nwtp)
#    1  2  3
# 1 NA  1  3
# 2  1 NA  3
# 3  2  2 NA

M.eq.nwtp.out <- list()
for (i in 1:10){
  states<- mapvalues(hist[[i]]$states, from = c("1", "2"), to=c("1&2", "3") )
  taxa<-cbind(names(states), states)
  root <- c(1/4, 1/4, 1/2)
  sum(root)
  M.eq.nwtp.out[[i]] <- mod_rayDISC(hist[[i]], taxa, rate.mat=M.eq.nwtp, root.p=root, verbose = TRUE)
}


#   ____________________________________________________________________________
#   ML inference under Super-Expansion Model                                ####

# This super-extension model is similar to M.eq but has three parameters
M.s <- initQ(c(1:3), c(0))
M.s[1,3] <- M.s[2,3] <- 2
M.s[3,1] <- 1
M.s[3,2] <- 3
M.s <- Q2model(M.s/.1, diag.as = NA)
print(M.s)
#    1  2  3
# 1 NA NA  3
# 2 NA NA  3
# 3  1  2 NA

M.s.out <- list()
for (i in 1:10){
  states<- mapvalues(hist[[i]]$states, from = c("1", "2"), to=c("1&2", "3") )
  taxa<-cbind(names(states), states)
  root <- c(1/4, 1/4, 1/2)
  sum(root)
  M.s.out[[i]] <- mod_rayDISC(hist[[i]], taxa, rate.mat=M.s, root.p=root, verbose = TRUE)
}


#   ____________________________________________________________________________
#   Likelihoods under Sub-Expansions                                        ####

# Since general parametric form does not exist for sub-expansions
# we do not run ML inference here since sub-expansions for paericular
# parameter estimates may super-large number of hidden states.

# Instead we check congruence for given parameter values
# by specifying the intial parameters (can be anything) values to
# mod_rayDISC() through p=...

#   ____________________________________________________________________________
#   Original model                                                          ####

print(Q.par)

Q.hat.out <- list()
for (i in 1:10){
  taxa <- cbind(hist[[i]]$tip.label, hist[[i]]$states)
  Q.hat.out[[i]] <- mod_rayDISC(hist[[i]], taxa, rate.mat=Q.par, root.p=c(0.5, 0.5),  p=c(.2,.3), verbose = TRUE)
}


##  ............................................................................
##  EHE transformation                                                      ####

Q.ehe <- EHEtransform(Q)
Q.ehe <- Q2model(Q.ehe, diag.as = NA)
print(Q.ehe)
#      1.1 1.2 2.1 2.2 2.3
# 1.1  NA  NA   1   1   1
# 1.2  NA  NA   1   1   1
# 2.1   1   1  NA  NA  NA
# 2.2   1   1  NA  NA  NA
# 2.3   1   1  NA  NA  NA

Q.ehe.out <- list()
for (i in 1:10){
  states<- mapvalues(hist[[i]]$states, from = c("1", "2"), to=c("1&2", "3&4&5") )
  taxa<-cbind( names(states), states)
  root <- c(1/4, 1/4, rep(1/6,3))
  sum(root)
  Q.ehe.out[[i]] <- mod_rayDISC(hist[[i]], taxa, rate.mat=Q.ehe, root.p=root, p=c(.1), verbose = TRUE)
}


##  ............................................................................
##  EHE nWTP transformation                                                 ####

# off-diagonal paramters can be any values

Q.ehe.nwtp <- EHEtransform(Q)
Q.ehe.nwtp[Q.ehe.nwtp==0] <- 20
Q.ehe.nwtp <- Q2model(Q.ehe.nwtp, diag.as = NA)
Q.ehe.nwtp
#      1.1 1.2 2.1 2.2 2.3
# 1.1  NA   1   2   2   2
# 1.2   1  NA   2   2   2
# 2.1   2   2  NA   1   1
# 2.2   2   2   1  NA   1
# 2.3   2   2   1   1  NA

Q.ehe.nwtp.out1 <- list()
Q.ehe.nwtp.out2 <- list()
for (i in 1:10){
  states<- mapvalues(hist[[i]]$states, from = c("1", "2"), to=c("1&2", "3&4&5") )
  taxa<-cbind( names(states), states)
  root <- c(1/4, 1/4, rep(1/6,3))
  sum(root)
  Q.ehe.nwtp.out1[[i]] <- mod_rayDISC(hist[[i]], taxa, rate.mat=Q.ehe.nwtp, root.p=root, p=c(1, .1), verbose = TRUE)
  Q.ehe.nwtp.out2[[i]] <- mod_rayDISC(hist[[i]], taxa, rate.mat=Q.ehe.nwtp, root.p=root, p=c(20, .1), verbose = TRUE)
}




##  ............................................................................
##  CHE transformation                                                      ####


Q.che <- CHEtransform(Q)$Q
Q.che <- Q2model(Q.che, diag.as = NA)
print(Q.che)
     # 1.1 1.2 1.3 1.4 2.1 2.2 2.3 2.4
# 1.1  NA  NA  NA  NA   1   1   1  NA
# 1.2  NA  NA  NA  NA   1   1  NA   1
# 1.3  NA  NA  NA  NA   1  NA   1   1
# 1.4  NA  NA  NA  NA  NA   1   1   1
# 2.1  NA   1   1  NA  NA  NA  NA  NA
# 2.2   1  NA  NA   1  NA  NA  NA  NA
# 2.3   1  NA  NA   1  NA  NA  NA  NA
# 2.4  NA   1  NA   1  NA  NA  NA  NA

Q.che.out <- list()
for (i in 1:10){
  states<- mapvalues(hist[[i]]$states, from = c("1", "2"), to=c("1&2&3&4", "5&6&7&8") )
  taxa<-cbind( names(states), states)
  root <- rep(1/8, 8)
  sum(root)
  Q.che.out[[i]] <- mod_rayDISC(hist[[i]], taxa, rate.mat=Q.che, root.p=root, p=c(.1), verbose = TRUE)
}



#   ____________________________________________________________________________
#   Results                                                                 ####

# if you did not run analyses yourself, uncomment this to load my results
# load("data/congruent_traits.RData")

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Maximum Likelihoods are the same                                        ####

Mln <- cbind(
'original'=lapply(Q.out, function(x) x$loglik) %>% unlist(),
'M.eq'=lapply(M.eq.out, function(x) x$loglik) %>% unlist(),
'M.eq.nwtp'=lapply(M.eq.nwtp.out, function(x) x$loglik) %>% unlist(),
'M.s.out'=lapply(M.s.out, function(x) x$loglik) %>% unlist()
)

print(Mln)
# the likelihoods are not perfectly the same due to numerical optimization
print(apply(Mln, 2, function(x) x-Mln[,1]))
# root mean squarred error
apply(Mln, 2, function(x) (x-Mln[,1])^2 %>% mean %>% sqrt() )



### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Likelihoods with specified rates are the same too                       ####

Ln <- cbind(
'original'=lapply(Q.hat.out, function(x) x$loglik) %>% unlist(),
'ehe'=lapply(Q.ehe.out, function(x) x$loglik) %>% unlist(),
'ehe.nwtp1'=lapply(Q.ehe.nwtp.out1, function(x) x$loglik) %>% unlist(),
'ehe.nwtp2'=lapply(Q.ehe.nwtp.out2, function(x) x$loglik) %>% unlist(),
'che'=lapply(Q.che.out, function(x) x$loglik) %>% unlist()
)

print(Ln)
# the likelihoods are not perfectly the same due to numerical optimization
print(apply(Ln, 2, function(x) x-Ln[,1]))
# root mean squarred error
apply(Ln, 2, function(x) (x-Ln[,1])^2 %>% mean %>% sqrt() )


#save.image("data/congruent_traits.RData")

