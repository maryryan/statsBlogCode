###################
######################### CONFOUDNING BIAS SIMULATIONS
######################### FOR STATS N CATS
######################### MARY RYAN
######################### CREATED: 9.22.2018
###################

#### LOAD PACKAGES ####
library(tidyverse)
library(MASS)

#### CREATE SIMULATED DATASET UNDER CONFOUNDING MODEL ####
set.seed(1234)
## number of pop datapoints, sims, and sample size ##
N <- 1000
nsims <- 100
Ni <- 100
## relationship between X and W ##
#odds.ratio <- seq( 0, 3, by=0.2 )
corr.xw <- seq( 0, 0.9, by=0.1 )
## number of conditions ##
#K <- length( odds.ratio )
K <- length( corr.xw )
trueMeans <- c( 3, 7 )
var.x <- 1.2
var.w <- 3
## true betas for response model ##
beta.true <- c( 4, 2, 1.2 ) #intercept, x, w
## std dev for random noise ##
sig.ep <- 0.5

## pre-allocate w prob, w, D, and y lists ##
#w.p <- vector( mode="list", length=K )
XW <- vector( mode="list", length=K )
D <- vector( mode="list", length=K )
y <- vector( mode="list", length=K )

for( k in seq(K) ){
  
  #w.p[[k]] <- 1/( 1 + exp( 0.5 + odds.ratio[k]*x ) ) 
  #w[[k]] <- rbinom( N, 1, w.p[[k]] )
  cov.xw <- corr.xw[[k]] * sqrt(var.x) * sqrt(var.w)
  cov.mat <- matrix( c(var.x, cov.xw, cov.xw, var.w), ncol=2 )
  XW <- mvrnorm( N,  trueMeans, cov.mat )
  ## design matrix ##
  #D[[k]] <- cbind( rep(1, N), x, w[[k]] )
  D[[k]] <- cbind( rep(1, N), XW )
  ## response data ##
  y[[k]] <- D[[k]] %*% beta.true + rnorm( N, 0, sig.ep )
  
  
}

#### SIMULATIONS ####
## pre-allocate vectors ##
unadjust.beta1.est <- matrix( NA, ncol=nsims, nrow=K )
unadjust.beta1.var <- NULL
adjust.beta1.est <- matrix( NA, ncol=nsims, nrow=K )
adjust.beta1.var <- NULL

for( k in seq(K) ){
  
  for( i in seq(nsims) ){
    
    ## SAMPLE DATA ##
    index <- sample( 1:N, Ni, replace=T )
    x.obs <- D[[k]][index, 2]
    w.obs <- D[[k]][index, 3]
    y.obs <- y[[k]][index]
    
    ## ESTIMATE UNADJUSTED MODEL ##
    unadjust.model <- lm( y.obs ~ x.obs )
    unadjust.beta1.est[k,i] <- summary( unadjust.model )$coef[2,1]

    ## ESTIMATE ADJUSTED MODEL ##
    adjust.model <- lm( y.obs ~ x.obs + w.obs )
    adjust.beta1.est[k,i] <- summary( adjust.model )$coef[2,1]
    
  }
  unadjust.beta1.var[k] <- var( unadjust.beta1.est[k,] )
  adjust.beta1.var[k] <- var( adjust.beta1.est[k,] )
  
}

unadjust.beta1.est.m <- rowMeans( unadjust.beta1.est )
adjust.beta1.est.m <- rowMeans( adjust.beta1.est )

#### CREATE RESULTS TABLE ####
betaResults.unadjust <- cbind( round(unadjust.beta1.est.m, 5),
                      round(beta.true[2] - unadjust.beta1.est.m, 5) )
betaResults.adjust <- cbind( round(adjust.beta1.est.m, 5),
                               round(beta.true[2] - adjust.beta1.est.m, 5) )
varResults <- cbind( round(unadjust.beta1.var, 5), round(adjust.beta1.var, 5) )
bias.var <- round(varResults[,1] - varResults[,2], 5)

results.table.est <- cbind( rep(beta.true[2], length(corr.xw)),
                            betaResults.unadjust, betaResults.adjust )
colnames(results.table.est) <- c( "True Beta1","Unadjusted Beta1 Est.",
                                  "Unadjusted Beta1 Est. Bias",
                                  "Adjusted Beta1 Est.",
                                  "Adjusted Beta1 Est. Bias" )
#rownames(results.table) <- paste0(rep("Odds Ratio = ", K),
#                                  odds.ratio)
rho <- '\u03c1'
rownames(results.table.est) <- rep( paste0( rep(rho, length(corr.xw) ),
                                       rep("=", length(corr.xw)),
                                       corr.xw) )

results.table.var <- cbind( varResults, bias.var )
colnames(results.table.var) <- c( "Unadjusted Beta1 Var.",
                                  "Adjusted Beta1 Var.",
                                  "Unadjusted - Adjusted Var." )
rownames(results.table.var) <- rep( paste0( rep(rho, length(corr.xw) ),
                                            rep("=", length(corr.xw)),
                                            corr.xw) )

results.table.est.lim <- results.table.est[c(1,6),]
results.table.var.lim <- results.table.var[c(1,6),]
