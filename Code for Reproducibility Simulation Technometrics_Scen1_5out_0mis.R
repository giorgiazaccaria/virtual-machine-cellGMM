################################################################################
#        Code for replication - Technometrics Review article TCH-24-160        #
#       Simulation Study - Main Article - Table 1, Scenario 1 (well sep)       #
################################################################################
## n = 200, p = 5, G = 2, unbalanced clusters, non-spherical components
## 5% outliers, 0% missing
rm(list = ls())
setwd("/workspaces/virtual-machine-cellGMM")
## INSTALL PACKAGES (IF NECESSARY)
install.packages("cellWise")
install.packages("tclust")
remotes::install_version("snipEM", version = "1.0.1", upgrade = "never")
remotes::install_version("MixtureMissing", version = "2.0.0", upgrade = "never")

## LOAD THE PACKAGES
library(MASS)
library(pracma)
library(purrr)
library(mvnfast)
library(cellWise)
library(mclust)
library(tclust)
library(snipEM)
library(MixtureMissing)
library(doParallel)
library(parallel)
library(foreach)

##  UPLOAD THE cellGMM CODE (ATTACHED)
source("cellGMM.R", echo=TRUE)

set.seed(17387)

## PARAMETERS DEPENDING ON THE SCENARIO
nsample <- 10 
n.obs <- 200
G <- 2
p <- 5
maxdist <- 5
alpha.out <- 0.05
alpha.mis <- 0
out.val <- alpha.out*n.obs*p
NA.val <- alpha.mis*n.obs*p
## TUNING PARAMETER SETTING
alpha_tclust <- alpha.out*2
alpha_1 = alpha_2 <- alpha.out
alpha.A1 <- alpha.out
alpha.A2 <- alpha.out*2
nrep <- 40
nstart <- 10
niter <- 10
zero.tol <- 1e-16
tuning_param_init <- data.frame(alpha_tclust, alpha_1, alpha_2, alpha.A1, alpha.A2, 
                                nrep, nstart, niter)
rndstart <- 3
## PARAMETER SETTING
# Prior probabilities
pp <- c(0.3, 0.7)
# Labels
label.th <- NULL
for (g in 1:G) {
  label.th <- c(label.th, rep(g, n.obs*pp[g]))
}
Uth <- diag(G)
Uth <- Uth[label.th, , drop = FALSE]
# Component mean vectors
muval <- 10*rand(G-1, p-1)
mu <- matrix(0, G, p)
gg <- 1
for (g in 2:G) {
  mu[g, ] <- c((g + gg), muval[g-1, ])
}
D <- pdist(mu)
diag(D) <- NA
while (min(D, na.rm = TRUE) < maxdist) {
  muval <- 10*rand(G-1, p-1)
  mu <- matrix(0, G, p)
  gg <- 1
  for (g in 2:G) {
    mu[g, ] <- c((g + gg), muval[g-1, ])
  }
  D <- pdist(mu)
  diag(D) <- NA
}
# Component covariance matrices
Sigma <- vector (mode = "list", length = G)
Sigma[[1]] <- matrix(0, p, p)
for (j in 1:(p-1)) {
  for (jj in (j+1):p) {
    Sigma[[1]][j ,jj] <- (0.9)^(abs(j-jj))
  }
}
Sigma[[1]] <- Sigma[[1]] + t(Sigma[[1]]) + diag(p)
Sigma[[2]] <- Sigma[[1]]
Sigma.inv <- lapply(Sigma, solve)
maxfact.th <- max(unlist(lapply(Sigma, function(x){eigen(x)$values})))/min(unlist(lapply(Sigma, function(x){eigen(x)$values})))
maxfact <- ceiling(maxfact.th) + 5
## GENERATION PROCESS
X <- vector(mode = "list", length = nsample)
Xout <- vector(mode = "list", length = nsample)
Wth <- vector(mode = "list", length = nsample)
for (samp in 1:nsample) {
  # Data generation
  for (g in 1:G) {
    if (g == 1) {
      X[[samp]] <- mvrnorm(n.obs*pp[g], mu[g, ], Sigma[[g]], tol = .Machine$double.xmin)
    } else {
      X[[samp]] <- rbind(X[[samp]], mvrnorm(n.obs*pp[g], mu[g, ], Sigma[[g]], tol = .Machine$double.xmin))
    }
  }
  ## Contamination
  # Outlier generation
  repl <- sample(n.obs*p, out.val)
  Xout[[samp]] <- X[[samp]]
  Xout[[samp]][repl] <- runif(out.val, min = -10, max = 10)
  Wth[[samp]] <- matrix(1, nrow = n.obs, ncol = p)
  Wth[[samp]][repl] <- 0
  # Check outlyingness
  if (alpha.out > 0) {
    out.rows <- which(rowSums(Wth[[samp]]) < p)
    for (i in 1:length(out.rows)) {
      out <- which(Wth[[samp]][out.rows[i], ] == 0)
      D <- sweep(mu, 2, Xout[[samp]][out.rows[i], ], FUN = "-")
      while (any(!as.matrix(lapply(1:G, function(ii) {D[ii, , drop = FALSE] %*% Sigma.inv[[ii]] %*% t(D[ii, , drop = FALSE])})) > qchisq(p = 0.99, df = p))) {
        Xout[[samp]][out.rows[i], out] <- runif(length(out), min = -10, max = 10)
        D <- sweep(mu, 2, Xout[[samp]][out.rows[i], ], FUN = "-")
      }
    }
  }
  # Missing generation
  miss <- 1:(n.obs*p)
  miss <- miss[!(miss %in% repl)]
  replNA <- sample(miss, NA.val)
  Xout[[samp]][replNA] <- NA
  Wth[[samp]][replNA] <- 0
}

################################################################################
## MODELS' IMPLEMENTATION
## HETEROGENEOUS POPULATION METHODOLOGIES
Permut <- pracma:: perms(1:G)

## cellGMM
# No penalty (cellGMM.pen0)
internal.cellGMM <- function(X, G, Ug, Permut, tuningp, maxfact, rndstart){
  cellgmm <- cellGMM(X, G, tuningp, maxfact = maxfact, rndstart = rndstart, showprogress = FALSE)
  # Solve label switching  
  dif.pattern <- cellgmm$W[!duplicated(cellgmm$W), , drop = FALSE]
  pat.unit <- matrix(0, nrow(X), 1)
  for (t in 1:nrow(dif.pattern)) {
    pat.unit[apply(cellgmm$W, 1, identical, dif.pattern[t, ]), ] <- t
  }
  logl1 <- rep(0, dim(Permut)[1])
  for (l in 1:dim(Permut)[1]){
    permprior <- cellgmm$pp[Permut[l, ]]
    permmu <- cellgmm$mu[Permut[l, ], ]
    permSigma <- cellgmm$sigma[Permut[l, ]]
    w <- update_post(X, dif.pattern, pat.unit, permprior, permmu, permSigma)
    logl1[l] <- sum(sum(log(w)*Ug))
  }
  maxid <- which.max(logl1)
  labord <- Permut[maxid, ]
  # Order according to the label switching
  cellgmm.ord <- cellgmm
  cellgmm.ord$ordcomp <- labord
  cellgmm.ord$X.imputed <- cellgmm$X.imputed[, , labord]
  cellgmm.ord$pp <- cellgmm$pp[labord]
  cellgmm.ord$mu <- cellgmm$mu[labord, ]
  cellgmm.ord$sigma <- cellgmm$sigma[labord]
  cellgmm.ord$post <- cellgmm$post[, labord]
  return(cellgmm.ord)
}

result.cellGMM <- vector(mode = "list", length = nsample)
cl <- makeCluster(4)
registerDoParallel(cl)
result.cellGMM <- foreach (i = 1:nsample, .packages = c("MASS", "tclust", "purrr"), .errorhandling = "pass") %dopar% {
  internal.cellGMM(Xout[[i]], G, Uth, Permut, tuning_param_init, maxfact, rndstart)         
}
stopCluster(cl)

for (i in 1:nsample) {
  UU <- diag(G)
  UU <- UU[result.cellGMM[[i]]$label, ] 
  result.cellGMM[[i]]$U <- UU[, result.cellGMM[[i]]$ordcomp]
  abs.error.cellGMM <- abs((result.cellGMM[[i]]$X.imputed.best - X[[i]]))
  result.cellGMM[[i]]$error <- c(mean(abs.error.cellGMM), mean(abs.error.cellGMM / abs(X[[i]])), sqrt(mean(abs.error.cellGMM^2)))
}

# Penalty (cellGMM.penb)
internal.cellGMM.penopt <-function(X, G, penalty, init, Ug, Permut, tuningp, maxfact, rndstart){
  cellgmm <- cellGMM(X, G, tuningp, penalty = penalty, maxfact = maxfact, rndstart = rndstart, manual_initparam = init, showprogress = FALSE)
  # Solve label switching  
  dif.pattern <- cellgmm$W[!duplicated(cellgmm$W), , drop = FALSE]
  pat.unit <- matrix(0, nrow(X), 1)
  for (t in 1:nrow(dif.pattern)) {
    pat.unit[apply(cellgmm$W, 1, identical, dif.pattern[t, ]), ] <- t
  }
  logl1 <- rep(0, dim(Permut)[1])
  for (l in 1:dim(Permut)[1]){
    permprior <- cellgmm$pp[Permut[l, ]]
    permmu <- cellgmm$mu[Permut[l, ], ]
    permSigma <- cellgmm$sigma[Permut[l, ]]
    w <- update_post(X, dif.pattern, pat.unit, permprior, permmu, permSigma)
    logl1[l] <- sum(sum(log(w)*Ug, n.obs, G))
  }
  maxid <- which.max(logl1)
  labord <- Permut[maxid, ]
  # Order according to the label switching
  cellgmm.ord <- cellgmm
  cellgmm.ord$ordcomp <- labord
  cellgmm.ord$X.imputed <- cellgmm$X.imputed[, , labord]
  cellgmm.ord$pp <- cellgmm$pp[labord]
  cellgmm.ord$mu <- cellgmm$mu[labord, ]
  cellgmm.ord$sigma <- cellgmm$sigma[labord]
  cellgmm.ord$post <- cellgmm$post[, labord]
  return(cellgmm.ord)
}

penalty <- vector(mode = "list", length = nsample)
for (i in 1:nsample) {
  penalty[[i]] <- matrix(0, nrow = n.obs, ncol = p)
  for (ii in 1:n.obs) {
    logEst <- matrix(0, 1, p)
    for (g in 1:G) {
      logEst <- logEst + result.cellGMM[[i]]$post[ii, g]*log(diag(solve(result.cellGMM[[i]]$sigma[[g]])))
    }
    penalty[[i]][ii, ] <- 0.5 * (- logEst + qchisq(p = 0.99, df = 1))
  }
}

result.cellGMM.penopt <- vector(mode = "list", length = nsample)
cl <- makeCluster(4)
registerDoParallel(cl)
result.cellGMM.penopt <- foreach (i = 1:nsample, .packages = c("MASS", "tclust", "purrr"), .errorhandling = "pass") %dopar% {
  internal.cellGMM.penopt(Xout[[i]], G, penalty[[i]], result.cellGMM[[i]]$init, Uth, Permut, tuning_param_init, maxfact, rndstart)
}
stopCluster(cl)

for (i in 1:nsample) {
  UU <- diag(G)
  UU <- UU[result.cellGMM.penopt[[i]]$label, ] 
  result.cellGMM.penopt[[i]]$U <- UU[, result.cellGMM.penopt[[i]]$ordcomp]
  abs.error.cellGMM.penopt <- abs((result.cellGMM.penopt[[i]]$X.imputed.best - X[[i]]))
  result.cellGMM.penopt[[i]]$error <- c(mean(abs.error.cellGMM.penopt), mean(abs.error.cellGMM.penopt / abs(X[[i]])), sqrt(mean(abs.error.cellGMM.penopt^2)))
}

################################################################################
## TCLUST (Garcia et al, 2008) [R package: tclust]
post <- function(X,
                 pp,
                 mu,
                 Sigma) {
  G <- dim(mu)[1]
  w <- matrix(as.double(NA), dim(X)[1], G)
  for (g in 1:G) {
    w[,g] = log(pp[g]) + mclust::dmvnorm(X, mu[g,], as.matrix(Sigma[[g]]), log = TRUE)
  }
  wnorm <- apply(w, 1, max)
  w <- exp(sweep(w, 1, wnorm, "-"))
  w <- w / rowSums(w)
  w[which(w < sqrt(.Machine$double.eps), arr.ind = TRUE)] <- sqrt(.Machine$double.eps)
  return(w)
}

internal.tclust <- function(X, G, Ug, Permut, maxfact){
  model.tclust <- tclust(X, k = G, alpha = 0.25, iter.max = 500, restr.fact = maxfact)
  # Solve label switching  
  logl1 <- rep(0, dim(Permut)[1])
  for (l in 1:dim(Permut)[1]){
    permprior <- model.tclust$weights[Permut[l, ]]
    permmu <- t(model.tclust$centers)[Permut[l, ], ]
    permSigma <- array_tree(model.tclust$cov, 3)[Permut[l, ]]
    w <- post(X[model.tclust$cluster > 0, ], permprior, permmu, permSigma)
    logl1[l] <- sum(sum(log(w)*Ug[model.tclust$cluster > 0, ]))
  }
  maxid <- which.max(logl1)
  labord <- Permut[maxid, ]
  # Order according to the label switching
  tclust.ord <- model.tclust
  tclust.ord$ordcomp <- labord
  tclust.ord$weights <- model.tclust$weights[labord]
  tclust.ord$centers <- t(model.tclust$centers)[labord, ]
  tclust.ord$cov <- array_tree(model.tclust$cov, 3)[labord]
  tclust.ord$post <- post(X = X, pp = tclust.ord$weights, mu = tclust.ord$centers, Sigma = tclust.ord$cov)
  tclust.ord$label.compl <- apply(tclust.ord$post, 1, which.max)
  return(tclust.ord)
}

result.tclust <- vector(mode = "list", length = nsample)
cl <- makeCluster(4)
registerDoParallel(cl)
result.tclust <- foreach (i = 1:nsample, .packages = c("tclust", "purrr"), .errorhandling = "pass") %dopar% {
  internal.tclust(Xout[[i]], G, Uth, Permut, maxfact)
}
stopCluster(cl)

for (i in 1:nsample) {
  UU <- diag(G)
  UU <- UU[result.tclust[[i]]$label.compl, ] 
  result.tclust[[i]]$U <- UU
  result.tclust[[i]]$W <- matrix(1, n.obs, p)
  result.tclust[[i]]$W[result.tclust[[i]]$cluster == 0, ] <- rep(0, p)
}

################################################################################
## sclust (Farcomeni, 2014) [R package: snipEM]
# 25% of contamination out 
internal.sclust.25 <- function(X, G, Ug, Permut, maxfact){
  Vinit <- matrix(1, nrow(X), ncol(X))
  Vinit[which(X > quantile(X, 0.875) | X < quantile(X, 0.125))] <- 0
  Rinit <- kmeans(X, G)$clust
  model.sclust <- snipEM::sclust(X, k = G, V = Vinit, R = Rinit, restr.fact = maxfact, tol = 1e-6, maxiters = 500)
  model.sclust$sigma <- vector(mode = "list", length = G)
  for (g in 1:G){
    model.sclust$sigma[[g]] <- matrix(0, ncol(X), ncol(X))
    for (j in 1:ncol(X)) {
      model.sclust$sigma[[g]][j, ] <- model.sclust$S[g, j, ]
    }
  }
  # Solve label switching  
  dif.pattern <- model.sclust$V[!duplicated(model.sclust$V), , drop = FALSE]
  pat.unit <- matrix(0, nrow(X), 1)
  for (t in 1:nrow(dif.pattern)) {
    pat.unit[apply(model.sclust$V, 1, identical, dif.pattern[t, ]), ] <- t
  }
  logl1 <- rep(0, dim(Permut)[1])
  for (l in 1:dim(Permut)[1]){
    permprior <- model.sclust$pi[Permut[l, ]]
    permmu <- model.sclust$mu[Permut[l, ], ]
    permSigma <- model.sclust$sigma[Permut[l, ]]
    w <- update_post(X, dif.pattern, pat.unit, permprior, permmu, permSigma)
    logl1[l] <- sum(sum(log(w)*Ug))
  }
  maxid <- which.max(logl1)
  labord <- Permut[maxid, ]
  # Order according to the label switching
  sclust.ord <- model.sclust
  sclust.ord$ordcomp <- labord
  sclust.ord$pi <- model.sclust$pi[labord]
  sclust.ord$mu <- model.sclust$mu[labord, ]
  sclust.ord$sigma <- model.sclust$sigma[labord]
  sclust.ord$label.compl <- model.sclust$R
  sclust.ord$post <- update_post(X, dif.pattern, pat.unit, sclust.ord$pi, sclust.ord$mu, sclust.ord$sigma)
  if (any(model.sclust$R == 0)) {
    sclust.ord$post[model.sclust$R == 0, ] <- post(X = X[model.sclust$R == 0, , drop = FALSE], pp = sclust.ord$pi, mu = sclust.ord$mu, Sigma = sclust.ord$sigma)
    sclust.ord$label.compl[model.sclust$R == 0] <- apply(sclust.ord$post[model.sclust$R == 0, ], 1, which.max)
  }
  return(sclust.ord)
}

result.sclust.25 <- vector(mode = "list", length = nsample)
cl <- makeCluster(4)
registerDoParallel(cl)
result.sclust.25 <- foreach (i = 1:nsample, .packages = c("snipEM", "purrr", "tclust"), .errorhandling = "pass") %dopar% {
  internal.sclust.25(Xout[[i]], G, Uth, Permut, maxfact)       
}
stopCluster(cl)

for (i in 1:nsample) {
  UU <- diag(G)
  UU <- UU[result.sclust.25[[i]]$label.compl, ] 
  result.sclust.25[[i]]$U <- UU[, result.sclust.25[[i]]$ordcomp]
}

# 5% of contamination out 
internal.sclust <- function(X, G, Ug, Permut, maxfact){
  Vinit <- matrix(1, nrow(X), ncol(X))
  Vinit[which(X > quantile(X, 0.975) | X < quantile(X, 0.025))] <- 0 
  Rinit <- kmeans(X, G)$clust
  model.sclust <- snipEM::sclust(X, k = G, V = Vinit, R = Rinit, restr.fact = maxfact, tol = 1e-6, maxiters = 500)
  model.sclust$sigma <- vector(mode = "list", length = G)
  for (g in 1:G){
    model.sclust$sigma[[g]] <- matrix(0, ncol(X), ncol(X))
    for (j in 1:ncol(X)) {
      model.sclust$sigma[[g]][j, ] <- model.sclust$S[g, j, ]
    }
  }
  # Solve label switching  
  dif.pattern <- model.sclust$V[!duplicated(model.sclust$V), , drop = FALSE]
  pat.unit <- matrix(0, nrow(X), 1)
  for (t in 1:nrow(dif.pattern)) {
    pat.unit[apply(model.sclust$V, 1, identical, dif.pattern[t, ]), ] <- t
  }
  logl1 <- rep(0, dim(Permut)[1])
  for (l in 1:dim(Permut)[1]){
    permprior <- model.sclust$pi[Permut[l, ]]
    permmu <- model.sclust$mu[Permut[l, ], ]
    permSigma <- model.sclust$sigma[Permut[l, ]]
    w <- update_post(X, dif.pattern, pat.unit, permprior, permmu, permSigma)
    logl1[l] <- sum(sum(log(w)*Ug))
  }
  maxid <- which.max(logl1)
  labord <- Permut[maxid, ]
  # Order according to the label switching
  sclust.ord <- model.sclust
  sclust.ord$ordcomp <- labord
  sclust.ord$pi <- model.sclust$pi[labord]
  sclust.ord$mu <- model.sclust$mu[labord, ]
  sclust.ord$sigma <- model.sclust$sigma[labord]
  sclust.ord$label.compl <- model.sclust$R
  sclust.ord$post <- update_post(X, dif.pattern, pat.unit, sclust.ord$pi, sclust.ord$mu, sclust.ord$sigma)
  if (any(model.sclust$R == 0)) {
    sclust.ord$post[model.sclust$R == 0, ] <- post(X = X[model.sclust$R == 0, , drop = FALSE], pp = sclust.ord$pi, mu = sclust.ord$mu, Sigma = sclust.ord$sigma)
    sclust.ord$label.compl[model.sclust$R == 0] <- apply(sclust.ord$post[model.sclust$R == 0, ], 1, which.max)
  }
  return(sclust.ord)
}

result.sclust <- vector(mode = "list", length = nsample)
cl <- makeCluster(4)
registerDoParallel(cl)
result.sclust <- foreach (i = 1:nsample, .packages = c("snipEM", "purrr", "tclust"), .errorhandling = "pass") %dopar% {
  internal.sclust(Xout[[i]], G, Uth, Permut, maxfact)       
}
stopCluster(cl)

for (i in 1:nsample) {
  UU <- diag(G)
  UU <- UU[result.sclust[[i]]$label.compl, ] 
  result.sclust[[i]]$U <- UU[, result.sclust[[i]]$ordcomp]
}

################################################################################
## MNM (Multivariate Normal Mixture) [R package: MixtureMissing]
internal.MNM <- function(X, G, Ug, Permut){
  model.MNM <- MNM(X, G, max_iter = 500, epsilon = 1e-6)  
  # Solve label switching  
  logl1 <- rep(0, dim(Permut)[1])
  for (l in 1:dim(Permut)[1]){
    permprior <- model.MNM$pi[Permut[l, ]]
    permmu <- model.MNM$mu[Permut[l,], ]
    permSigma <- array_tree(model.MNM$sigma, 3)[Permut[l, ]]
    w <- post(X, permprior, permmu, permSigma)
    logl1[l] <- sum(sum(log(w)*Ug))
  }
  maxid <- which.max(logl1)
  labord <- Permut[maxid, ]
  # Order according to the label switching
  MNM.ord <- model.MNM
  MNM.ord$ordcomp <- labord
  MNM.ord$pi <- model.MNM$pi[labord]
  MNM.ord$mu <- model.MNM$mu[labord, ]
  MNM.ord$sigma <- array_tree(model.MNM$sigma, 3)[labord]
  MNM.ord$z_tilde <- model.MNM$z_tilde[, labord]
  return(MNM.ord)
}

result.MNM <- vector(mode = "list", length = nsample)
cl <- makeCluster(4)
registerDoParallel(cl)
result.MNM <- foreach (i = 1:nsample, .packages = c("MixtureMissing", "mclust", "purrr"), .errorhandling = "pass") %dopar% {
  internal.MNM(Xout[[i]], G, Uth, Permut)      
}
stopCluster(cl)

for (i in 1:nsample) {
  UU <- diag(G)
  UU <- UU[result.MNM[[i]]$clusters, ] 
  result.MNM[[i]]$U <- UU[, result.MNM[[i]]$ordcomp]
}

################################################################################
## MCNM (Tong and Tortora, 2022) [R package: MixtureMissing]
post_mcnm <- function(X,
                      pp,
                      mu,
                      Sigma, 
                      alpha,
                      eta) {
  G <- dim(mu)[1]
  w <- matrix(as.double(NA), dim(X)[1], G)
  for (g in 1:G) {
    w[,g] = log(pp[g]) + log(ContaminatedMixt::dCN(X, mu = mu[g,], Sigma = as.matrix(Sigma[[g]]), alpha = alpha[g], eta = eta[g]))
  }
  wnorm <- apply(w, 1, max)
  w <- exp(sweep(w, 1, wnorm, "-"))
  w <- w / rowSums(w)
  w[which(w < sqrt(.Machine$double.eps), arr.ind = TRUE)] <- sqrt(.Machine$double.eps)
  return(w)
}

internal.MCNM <- function(X, G, Ug, Permut){
  model.MCNM <- MCNM(X, G, max_iter = 500, epsilon = 1e-6)  
  # Solve label switching  
  logl1 <- rep(0, dim(Permut)[1])
  for (l in 1:dim(Permut)[1]){
    permprior <- model.MCNM$pi[Permut[l, ]]
    permmu <- model.MCNM$mu[Permut[l, ], ]
    permSigma <- array_tree(model.MCNM$sigma, 3)[Permut[l, ]]
    permalpha <- model.MCNM$alpha[Permut[l, ]]
    permeta <- model.MCNM$eta[Permut[l, ]]
    w <- post_mcnm(X, permprior, permmu, permSigma, permalpha, permeta)
    logl1[l] <- sum(sum(log(w)*Ug))
  }
  maxid <- which.max(logl1)
  labord <- Permut[maxid, ]
  # Order according to the label switching
  MCNM.ord <- model.MCNM
  MCNM.ord$ordcomp <- labord
  MCNM.ord$pi <- model.MCNM$pi[labord]
  MCNM.ord$mu <- model.MCNM$mu[labord, ]
  MCNM.ord$sigma <- array_tree(model.MCNM$sigma, 3)[labord]
  MCNM.ord$z_tilde <- model.MCNM$z_tilde[, labord]
  MCNM.ord$alpha <- model.MCNM$alpha[labord]
  MCNM.ord$eta <- model.MCNM$eta[labord]
  return(MCNM.ord)
}

result.MCNM <- vector(mode = "list", length = nsample)
cl <- makeCluster(4)
registerDoParallel(cl)
result.MCNM <- foreach (i = 1:nsample, .packages = c("MixtureMissing", "mclust", "purrr"), .errorhandling = "pass") %dopar% {
  internal.MCNM(Xout[[i]], G, Uth, Permut)      
}
stopCluster(cl)

for (i in 1:nsample) {
  UU <- diag(G)
  UU <- UU[result.MCNM[[i]]$clusters, ] 
  result.MCNM[[i]]$U <- UU[, result.MCNM[[i]]$ordcomp]
  result.MCNM[[i]]$W <- matrix(1, n.obs, p)
  result.MCNM[[i]]$W[result.MCNM[[i]]$outliers, ] <- rep(0, p)
}

################################################################################
## MtM (Wang et al, 2004) [R package: MixtureMissing]
post_t <- function(X,
                   pp,
                   mu,
                   Sigma, 
                   df) {
  G <- dim(mu)[1]
  w <- matrix(as.double(NA), dim(X)[1], G)
  for (g in 1:G) {
    w[,g] = log(pp[g]) + mvnfast::dmvt(X, mu[g,], as.matrix(Sigma[[g]]), df[g], log = TRUE)
  }
  wnorm <- apply(w, 1, max)
  w <- exp(sweep(w, 1, wnorm, "-"))
  w <- w / rowSums(w)
  w[which(w < sqrt(.Machine$double.eps), arr.ind = TRUE)] <- sqrt(.Machine$double.eps)
  return(w)
}

internal.MtM <- function(X, G, Ug, Permut){
  model.MtM <- MtM(X, G, max_iter = 500, epsilon = 1e-6, outlier_cutoff = 0.75)  
  # Solve label switching  
  logl1 <- rep(0, dim(Permut)[1])
  for (l in 1:dim(Permut)[1]){
    permprior <- model.MtM$pi[Permut[l, ]]
    permmu <- model.MtM$mu[Permut[l, ], ]
    permSigma <- array_tree(model.MtM$sigma, 3)[Permut[l, ]]
    permdf <- model.MtM$df[Permut[l, ]]
    w <- post_t(X, permprior, permmu, permSigma, permdf)
    logl1[l] <- sum(sum(log(w)*Ug))
  }
  maxid <- which.max(logl1)
  labord <- Permut[maxid, ]
  # Order according to the label switching
  MtM.ord <- model.MtM
  MtM.ord$ordcomp <- labord
  MtM.ord$pi <- model.MtM$pi[labord]
  MtM.ord$mu <- model.MtM$mu[labord, ]
  MtM.ord$sigma <- array_tree(model.MtM$sigma, 3)[labord]
  MtM.ord$df <- model.MtM$df[labord]
  MtM.ord$z_tilde <- model.MtM$z_tilde[, labord]
  return(MtM.ord)
}

result.MtM <- vector(mode = "list", length = nsample)
cl <- makeCluster(4)
registerDoParallel(cl)
result.MtM <- foreach (i = 1:nsample, .packages = c("MixtureMissing", "mclust", "purrr"), .errorhandling = "pass") %dopar% {
  internal.MtM(Xout[[i]], G, Uth, Permut)  
}
stopCluster(cl)

for (i in 1:nsample) {
  UU <- diag(G)
  UU <- UU[result.MtM[[i]]$clusters, ] 
  result.MtM[[i]]$U <- UU[, result.MtM[[i]]$ordcomp]
  result.MtM[[i]]$W <- matrix(1, n.obs, p)
  result.MtM[[i]]$W[result.MtM[[i]]$outliers, ] <- rep(0, p)
}

################################################################################
# POSTERIOR
post.th <- vector(mode = "list", length = nsample)
for (samp in 1:nsample) {
  post.th[[samp]] <- post(X[[samp]], pp, mu, Sigma)
}

################################################################################
## SINGLE-POPULATION METHODS
## cellMCD (Raymaekers and Rousseeuw, 2023) [R package: cellWise]
result.cellMCD <- vector(mode = "list", length = nsample)
cl <- makeCluster(4)
registerDoParallel(cl)
result.cellMCD <- foreach (i = 1:nsample, .packages = c("cellWise"), .errorhandling = "pass") %dopar% {
  cellMCD(Xout[[i]], alpha = 0.50, crit = 1e-6, noCits = 500, checkPars = list(coreOnly = TRUE, numDiscrete = 0))       
}
stopCluster(cl)

for (i in 1:nsample) {
  abs.error.cellMCD <- abs((result.cellMCD[[i]]$Ximp - X[[i]]))
  result.cellMCD[[i]]$error <- c(mean(abs.error.cellMCD), mean(abs.error.cellMCD / abs(X[[i]])), sqrt(mean(abs.error.cellMCD^2)))
}

## DI (Raymaekers and Rousseeuw, 2021) [R package: cellWise]
result.DI <- vector(mode = "list", length = nsample)
cl <- makeCluster(4)
registerDoParallel(cl)
result.DI <- foreach (i = 1:nsample, .packages = c("cellWise"), .errorhandling = "pass") %dopar% {
  DI(Xout[[i]], crit = 1e-6, maxits = 500, checkPars = list(coreOnly = TRUE, numDiscrete = 0))
}
stopCluster(cl)

for (i in 1:nsample) {
  result.DI[[i]]$W <- matrix(1, n.obs, p)
  result.DI[[i]]$W[result.DI[[i]]$indcells] <- 0
  abs.error.DI <- abs((result.DI[[i]]$Ximp - X[[i]]))
  result.DI[[i]]$error <- c(mean(abs.error.DI), mean(abs.error.DI / abs(X[[i]])), sqrt(mean(abs.error.DI^2)))
}

################################################################################
## MODELS' EVALUATION
ARI <- matrix(as.double(NA), nsample, 8)
Missclass <- matrix(as.double(NA), nsample, 8)
MSE.post <- matrix(as.double(NA), nsample, 8)
MSE.pi <- matrix(as.double(NA), nsample, 8)
MSE.mu <- matrix(as.double(NA), nsample, G*8)
KL.sigma <- matrix(as.double(NA), nsample, G*8)
Ws <- matrix(as.double(NA), nsample, 2*9)
param.stat <- matrix(as.double(NA), 8, 8)
stat.W <- matrix(as.double(NA), 9, 5)
count.sample <- matrix(0, 8, 1)
error <- matrix(as.double(NA), nsample, 3*4)
stat.error <- matrix(as.double(NA), 4, 3)
for (samp in 1:nsample) {
  # ARI
  ARI[samp, 1] <- adjustedRandIndex(label.th, result.cellGMM[[samp]]$label)
  count.sample[1, ] <- count.sample[1, ] + 1
  ARI[samp, 2] <- adjustedRandIndex(label.th, result.cellGMM.penopt[[samp]]$label)
  count.sample[2, ] <- count.sample[2, ] + 1
  ARI[samp, 3] <- adjustedRandIndex(label.th, result.tclust[[samp]]$label.compl)
  count.sample[3, ] <- count.sample[3, ] + 1
  if (!is.character(unlist(result.sclust.25[[samp]][1]))) {
    ARI[samp, 4] <- adjustedRandIndex(label.th, result.sclust.25[[samp]]$label.compl)
    count.sample[4, ] <- count.sample[4, ] + 1
  }
  if (!is.character(unlist(result.sclust[[samp]][1]))) {
    ARI[samp, 5] <- adjustedRandIndex(label.th, result.sclust[[samp]]$label.compl)
    count.sample[5, ] <- count.sample[5, ] + 1
  }
  if (length(result.MNM[[samp]]) > 4) {
    ARI[samp, 6] <- adjustedRandIndex(label.th, result.MNM[[samp]]$clusters)
    count.sample[6, ] <- count.sample[6, ] + 1
  }
  if (length(result.MCNM[[samp]]) > 4) {
    ARI[samp, 7] <- adjustedRandIndex(label.th, result.MCNM[[samp]]$clusters)
    count.sample[7, ] <- count.sample[7, ] + 1
  }
  if (length(result.MtM[[samp]]) > 4) {
    ARI[samp, 8] <- adjustedRandIndex(label.th, result.MtM[[samp]]$clusters)
    count.sample[8, ] <- count.sample[8, ] + 1
  }
  # Missclass
  Missclass[samp, 1] <- sum(t(result.cellGMM[[samp]]$U) %*% Uth - diag(diag(t(result.cellGMM[[samp]]$U) %*% Uth)))/n.obs
  Missclass[samp, 2] <- sum(t(result.cellGMM.penopt[[samp]]$U) %*% Uth - diag(diag(t(result.cellGMM.penopt[[samp]]$U) %*% Uth)))/n.obs
  Missclass[samp, 3] <- sum(t(result.tclust[[samp]]$U) %*% Uth - diag(diag(t(result.tclust[[samp]]$U) %*% Uth)))/n.obs
  if (!is.character(unlist(result.sclust.25[[samp]][1]))) {
    Missclass[samp, 4] <- sum(t(result.sclust.25[[samp]]$U) %*% Uth - diag(diag(t(result.sclust.25[[samp]]$U) %*% Uth)))/n.obs
  }
  if (!is.character(unlist(result.sclust[[samp]][1]))) {
    Missclass[samp, 5] <- sum(t(result.sclust[[samp]]$U) %*% Uth - diag(diag(t(result.sclust[[samp]]$U) %*% Uth)))/n.obs
  }
  if (length(result.MNM[[samp]]) > 4) {
    Missclass[samp, 6] <-sum(t(result.MNM[[samp]]$U) %*% Uth - diag(diag(t(result.MNM[[samp]]$U) %*% Uth)))/n.obs
  }
  if (length(result.MCNM[[samp]]) > 4) {
    Missclass[samp, 7] <- sum(t(result.MCNM[[samp]]$U) %*% Uth - diag(diag(t(result.MCNM[[samp]]$U) %*% Uth)))/n.obs
  }
  if (length(result.MtM[[samp]]) > 4) {
    Missclass[samp, 8] <- sum(t(result.MtM[[samp]]$U) %*% Uth - diag(diag(t(result.MtM[[samp]]$U) %*% Uth)))/n.obs
  }
  for (g in 1:G) {
    if (g == 1) {
      # MSE.pi
      MSE.pi[samp, g] <- (pp[g] - result.cellGMM[[samp]]$pp[g])^2
      MSE.pi[samp, g + 1] <- (pp[g] - result.cellGMM.penopt[[samp]]$pp[g])^2
      MSE.pi[samp, g + 2] <- (pp[g] - result.tclust[[samp]]$weights[g])^2
      if (!is.character(unlist(result.sclust.25[[samp]][1]))) {
        MSE.pi[samp, g + 3] <- (pp[g] - result.sclust.25[[samp]]$pi[g])^2
      }
      if (!is.character(unlist(result.sclust[[samp]][1]))) {
        MSE.pi[samp, g + 4] <- (pp[g] - result.sclust[[samp]]$pi[g])^2
      }
      if (length(result.MNM[[samp]]) > 4) {
        MSE.pi[samp, g + 5] <- (pp[g] - result.MNM[[samp]]$pi[g])^2
      }
      if (length(result.MCNM[[samp]]) > 4) {
        MSE.pi[samp, g + 6] <- (pp[g] - result.MCNM[[samp]]$pi[g])^2
      }
      if (length(result.MtM[[samp]]) > 4) {
        MSE.pi[samp, g + 7] <- (pp[g] - result.MtM[[samp]]$pi[g])^2
      }
      # MSE.post
      MSE.post[samp, g] <- sqrt(mean((post.th[[samp]][, g] - result.cellGMM[[samp]]$post[, g])^2))
      MSE.post[samp, g + 1] <- sqrt(mean((post.th[[samp]][, g] - result.cellGMM.penopt[[samp]]$post[, g])^2))
      MSE.post[samp, g + 2] <- sqrt(mean((post.th[[samp]][, g] - result.tclust[[samp]]$post[, g])^2))
      if (!is.character(unlist(result.sclust.25[[samp]][1]))) {
        MSE.post[samp, g + 3] <- sqrt(mean((post.th[[samp]][, g] - result.sclust.25[[samp]]$post[, g])^2))
      }
      if (!is.character(unlist(result.sclust[[samp]][1]))) {
        MSE.post[samp, g + 4] <- sqrt(mean((post.th[[samp]][, g] - result.sclust[[samp]]$post[, g])^2))
      }
      if (length(result.MNM[[samp]]) > 4) {
        MSE.post[samp, g + 5] <- sqrt(mean((post.th[[samp]][, g] - result.MNM[[samp]]$z_tilde[, g])^2))
      }
      if (length(result.MCNM[[samp]]) > 4) {
        MSE.post[samp, g + 6] <- sqrt(mean((post.th[[samp]][, g] - result.MCNM[[samp]]$z_tilde[, g])^2))
      }
      if (length(result.MtM[[samp]]) > 4) {
        MSE.post[samp, g + 7] <- sqrt(mean((post.th[[samp]][, g] - result.MtM[[samp]]$z_tilde[, g])^2))
      }
    }
    # MSE.mu
    MSE.mu[samp, g] <- mean((mu[g, ] - result.cellGMM[[samp]]$mu[g, ])^2)
    MSE.mu[samp, g + G] <- mean((mu[g, ] - result.cellGMM.penopt[[samp]]$mu[g, ])^2)
    MSE.mu[samp, g + 2*G] <- mean((mu[g, ] - result.tclust[[samp]]$centers[g, ])^2)
    if (!is.character(unlist(result.sclust.25[[samp]][1]))) {
      MSE.mu[samp, g + 3*G] <- mean((mu[g, ] - result.sclust.25[[samp]]$mu[g, ])^2)
    }
    if (!is.character(unlist(result.sclust[[samp]][1]))) {
      MSE.mu[samp, g + 4*G] <- mean((mu[g, ] - result.sclust[[samp]]$mu[g, ])^2)
    }
    if (length(result.MNM[[samp]]) > 4) {
      MSE.mu[samp, g + 5*G] <- mean((mu[g, ] - result.MNM[[samp]]$mu[g, ])^2)
    }
    if (length(result.MCNM[[samp]]) > 4) {
      MSE.mu[samp, g + 6*G] <- mean((mu[g, ] - result.MCNM[[samp]]$mu[g, ])^2)
    }
    if (length(result.MtM[[samp]]) > 4) {
      MSE.mu[samp, g + 7*G] <- mean((mu[g, ] - result.MtM[[samp]]$mu[g, ])^2)
    }
    # KL.sigma
    KL.sigma[samp, g] <- tr(result.cellGMM[[samp]]$sigma[[g]] %*% solve(Sigma[[g]])) - p - log(det(result.cellGMM[[samp]]$sigma[[g]] %*% solve(Sigma[[g]])))
    KL.sigma[samp, g + G] <- tr(result.cellGMM.penopt[[samp]]$sigma[[g]] %*% solve(Sigma[[g]])) - p - log(det(result.cellGMM.penopt[[samp]]$sigma[[g]] %*% solve(Sigma[[g]])))
    KL.sigma[samp, g + 2*G] <- tr(result.tclust[[samp]]$cov[[g]] %*% solve(Sigma[[g]])) - p - log(det(result.tclust[[samp]]$cov[[g]] %*% solve(Sigma[[g]])))
    if (!is.character(unlist(result.sclust.25[[samp]][1]))) {
      KL.sigma[samp, g + 3*G] <- tr(result.sclust.25[[samp]]$sigma[[g]] %*% solve(Sigma[[g]])) - p - log(det(result.sclust.25[[samp]]$sigma[[g]] %*% solve(Sigma[[g]])))
    }
    if (!is.character(unlist(result.sclust[[samp]][1]))) {
      KL.sigma[samp, g + 4*G] <- tr(result.sclust[[samp]]$sigma[[g]] %*% solve(Sigma[[g]])) - p - log(det(result.sclust[[samp]]$sigma[[g]] %*% solve(Sigma[[g]])))
    }
    if (length(result.MNM[[samp]]) > 4) {
      KL.sigma[samp, g + 5*G] <- tr(result.MNM[[samp]]$sigma[[g]] %*% solve(Sigma[[g]])) - p - log(det(result.MNM[[samp]]$sigma[[g]] %*% solve(Sigma[[g]])))
    }
    if (length(result.MCNM[[samp]]) > 4) {
      KL.sigma[samp, g + 6*G] <- tr(result.MCNM[[samp]]$sigma[[g]] %*% solve(Sigma[[g]])) - p - log(det(result.MCNM[[samp]]$sigma[[g]] %*% solve(Sigma[[g]])))
    }
    if (length(result.MtM[[samp]]) > 4) {
      KL.sigma[samp, g + 7*G] <- tr(result.MtM[[samp]]$sigma[[g]] %*% solve(Sigma[[g]])) - p - log(det(result.MtM[[samp]]$sigma[[g]] %*% solve(Sigma[[g]])))
    }
  }
  # Ws
  Ws[samp, 1] <- sum(which(result.cellGMM[[samp]]$W==0)%in%which(Wth[[samp]]==0), na.rm = TRUE)/(n.obs*p *alpha.out)*100
  Ws[samp, 2] <- sum(which(result.cellGMM[[samp]]$W==0)%in%which(Wth[[samp]]==1), na.rm = TRUE)/(n.obs*p - n.obs*p *alpha.out)*100
  Ws[samp, 3] <- sum(which(result.cellGMM.penopt[[samp]]$W==0)%in%which(Wth[[samp]]==0), na.rm = TRUE)/(n.obs*p *alpha.out)*100
  Ws[samp, 4] <- sum(which(result.cellGMM.penopt[[samp]]$W==0)%in%which(Wth[[samp]]==1), na.rm = TRUE)/(n.obs*p - n.obs*p *alpha.out)*100
  Ws[samp, 5] <- sum(which(result.tclust[[samp]]$W==0)%in%which(Wth[[samp]]==0), na.rm = TRUE)/(n.obs*p *alpha.out)*100
  Ws[samp, 6] <- sum(which(result.tclust[[samp]]$W==0)%in%which(Wth[[samp]]==1), na.rm = TRUE)/(n.obs*p - n.obs*p *alpha.out)*100
  if (!is.character(unlist(result.sclust.25[[samp]][1]))) {
    Ws[samp, 7] <- sum(which(result.sclust.25[[samp]]$V==0)%in%which(Wth[[samp]]==0), na.rm = TRUE)/(n.obs*p *alpha.out)*100
    Ws[samp, 8] <- sum(which(result.sclust.25[[samp]]$V==0)%in%which(Wth[[samp]]==1), na.rm = TRUE)/(n.obs*p - n.obs*p *alpha.out)*100
  }
  if (!is.character(unlist(result.sclust[[samp]][1]))) {
    Ws[samp, 9] <- sum(which(result.sclust[[samp]]$V==0)%in%which(Wth[[samp]]==0), na.rm = TRUE)/(n.obs*p *alpha.out)*100
    Ws[samp, 10] <- sum(which(result.sclust[[samp]]$V==0)%in%which(Wth[[samp]]==1), na.rm = TRUE)/(n.obs*p - n.obs*p *alpha.out)*100
  }
  if (length(result.MCNM[[samp]]) > 4) {
    Ws[samp, 11] <- sum(which(result.MCNM[[samp]]$W==0)%in%which(Wth[[samp]]==0), na.rm = TRUE)/(n.obs*p *alpha.out)*100
    Ws[samp, 12] <- sum(which(result.MCNM[[samp]]$W==0)%in%which(Wth[[samp]]==1), na.rm = TRUE)/(n.obs*p - n.obs*p *alpha.out)*100
  }
  if (length(result.MtM[[samp]]) > 4) {
    Ws[samp, 13] <- sum(which(result.MtM[[samp]]$W==0)%in%which(Wth[[samp]]==0), na.rm = TRUE)/(n.obs*p *alpha.out)*100
    Ws[samp, 14] <- sum(which(result.MtM[[samp]]$W==0)%in%which(Wth[[samp]]==1), na.rm = TRUE)/(n.obs*p - n.obs*p *alpha.out)*100
  }
  Ws[samp, 15] <- sum(which(result.cellMCD[[samp]]$W==0)%in%which(Wth[[samp]]==0), na.rm = TRUE)/(n.obs*p *alpha.out)*100
  Ws[samp, 16] <- sum(which(result.cellMCD[[samp]]$W==0)%in%which(Wth[[samp]]==1), na.rm = TRUE)/(n.obs*p - n.obs*p *alpha.out)*100
  Ws[samp, 17] <- sum(which(result.DI[[samp]]$W==0)%in%which(Wth[[samp]]==0), na.rm = TRUE)/(n.obs*p *alpha.out)*100
  Ws[samp, 18] <- sum(which(result.DI[[samp]]$W==0)%in%which(Wth[[samp]]==1), na.rm = TRUE)/(n.obs*p - n.obs*p *alpha.out)*100
  # Imputation
  error[samp, 1:3] <- result.cellGMM[[samp]]$error
  error[samp, 4:6] <- result.cellGMM.penopt[[samp]]$error
  error[samp, 7:9] <- result.cellMCD[[samp]]$error
  error[samp, 10:12] <- result.DI[[samp]]$error
}
if (any(KL.sigma == Inf)) {
  KL.sigma[KL.sigma == Inf] <- NA
}
KL.sigma.adj <- matrix(apply(KL.sigma, 2, function(x){mean(x, na.rm = TRUE)}), nrow = 8, G, byrow = TRUE)
param.stat <- data.frame(count.sample, matrix(apply(ARI, 2, function(x){mean(x, na.rm = TRUE)})), 
                         matrix(apply(Missclass, 2, function(x){mean(x, na.rm = TRUE)})),
                         matrix(apply(MSE.post, 2, function(x){mean(x, na.rm = TRUE)})), 
                         matrix(apply(MSE.pi, 2, function(x){mean(x, na.rm = TRUE)})), 
                         matrix(apply(MSE.mu, 2, function(x){mean(x, na.rm = TRUE)}), 8, G, byrow = TRUE),
                         KL.sigma.adj)
colnames(param.stat) <- c("#samples", "mARI", "MR", "MSE.post", "MSE.pi", "MSE.mu1", "MSE.mu2", "KL.sigma1", "KL.sigma2")
rownames(param.stat) <- c("cellGMM.pen0", "cellGMM.penb", "TCLUST", "sclust_25", "sclust_5", "MNM", "MCNM", "MtM")
stat.error <- data.frame(matrix(apply(error, 2, function(x){mean(x, na.rm = TRUE)}), 4, 3, byrow = TRUE))
stat.W <-data.frame(matrix(apply(Ws, 2, function(x){mean(x, na.rm = TRUE)}), 9, 2, byrow = TRUE),
                    rbind(stat.error[1:2, ], data.frame(matrix(as.double(NA), 5, 3)), stat.error[3:4, ]))
colnames(stat.W) <- c("%TP", "%FP", "MAE", "MARE", "RMSE")
rownames(stat.W) <- c("cellGMM.pen0", "cellGMM.penb", "TCLUST", "sclust_25", "sclust_5", "MCNM", "MtM", "cellMCD", "DI")

################################################################################
## REPRESENTATION

## TABLE [Table 1, Scenario 1, 5% outliers - Main Article]
# install.packages("xtable") # If necessary
library(xtable)
xtable(stat.W)

