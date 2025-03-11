## Internal functions for cellGMM ###########

################################################################################
tr <-
  function(x)
    sum(diag(as.matrix(x)))

################################################################################
norm <-
  function(x,
           type = c("standard", "center", "range")) {
    type <- match.arg(type,
                      choices = eval(formals(norm)$type),
                      several.ok = FALSE)
    x <- as.matrix(x)
    switch(
      type,
      "standard" = scale(x, center = apply(x, 2, mean, na.rm = TRUE), scale = apply(x, 2, sd, na.rm = TRUE)),
      "center"   = scale(x, center = apply(x, 2, mean, na.rm = TRUE), scale = FALSE),
      "range"    = apply(x, 2, function(x)
        (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
    )
    return(x)
  }

################################################################################

rand.member <-
  function(n.obs,
           G) {
    if (G > n.obs)
      stop('The number of groups is larger than the number of observations')
    if (G == n.obs) {
      U <- diag(n.obs)
    } else if (G == 1) {
      U <- c(rep(1, n.obs))
    }
    else {
      U <- matrix(0, n.obs, G)
      U[1:G,] = diag(G)
      U[(G + 1):n.obs, 1] <- 1
      for (i in (G + 1):n.obs) {
        U[i,] <- U[i, sample(G)]
      }
      U <- U[sample(n.obs),]
    }
    return(as.matrix(U))
  }

################################################################################
initfunc_new <- function(X, 
                         G,
                         maxfact,
                         zero.tol,
                         alpha_tclust, 
                         alpha_1,
                         alpha_2,
                         alpha.A1, 
                         alpha.A2,
                         nrep,
                         nstart,
                         niter) {
  p <- ncol(X)
  W <- w_initial(y = X, alpha_tclust = alpha_tclust, alpha_1 = alpha_1, alpha_2 = alpha_2, K = G, 1)$ww
  x <- X
  x[W==0] <- NA
  prpar <- prepars_initial(x = x, K = G, alpha = alpha.A1, q = floor(p/2) + 1, nrep = nrep, nstart = nstart, niter = niter)
  par <- pars_initial(x = x, K = G, alpha = alpha.A2, nstart = nstart, niter = niter, prpar$ctr_, prpar$cts_, prpar$rws)
  group <- informal_eval(X = x, K = G, par$init_mean)
  Sigma <- purrr::array_tree(par$init_cov, 3)
  if (any(unlist(lapply(lapply(Sigma, function(x){eigen(x)$values}), is.complex)))) {
    refg <- which(unlist(lapply(lapply(Sigma, function(x) {eigen(x)$values}), is.complex)) == TRUE)
    for (g in refg) {
      decSigma <- eigen(Sigma[[g]])
      Sigma[[g]] <- Re(decSigma$vectors)%*%diag(Re(decSigma$values))%*%t(Re(decSigma$vectors))
    }
  }
  par$init_cov_constr <- vector(mode = "list", length = G)
  par$init_cov_constr <- restr_diffax(Sigma, maxfact = maxfact, zero.tol = zero.tol, label = group)
  return(list(mu = par$init_mean,
              Sigma = par$init_cov_constr,
              pp = as.vector(par$init_pi_),
              W = W,
              label = group))
}

################################################################################
getini <- 
  function (n.obs,
            G) {
    if (G == 1) {
      return (n.obs)
    } else if (G > 1) {
      prob.init <- runif(G, min = 5/n.obs)   
      csize <- sample(G, n.obs, replace = TRUE, prob = prob.init/sum(prob.init)) 
      return (tabulate(csize, nbins = G))   
    }
  }

################################################################################
impute_mis <-
  function(X,
           difpat,
           pat,
           pp,
           mu,
           Sigma) {
    n.obs <- dim(X)[1]
    p <- dim(X)[2]
    G <- dim(mu)[1]
    Ximp <- array(rep(X, G), dim = c(n.obs, p, G))
    sigma.imp <- array(0, dim = c(p, p, n.obs, G))
    for (t in 1:nrow(difpat)) {
      var.mis <- difpat[t, , drop = FALSE] == 0
      for (g in 1:G) {
        if (any(var.mis)) {
          if (all(var.mis)) {
            Ximp[pat == t, , g] <- mu[g, , drop = FALSE]
            un.t <- which(pat==t)
            for (i in 1:length(un.t)) {
              sigma.imp[, , un.t[i], g] <- Sigma[[g]]
            }
          } else {
            var.obs <- !var.mis
            mu_m <- mu[g, var.mis, drop = FALSE]
            mu_o <- mu[g, var.obs, drop = FALSE]
            sigma_mo <- Sigma[[g]][var.mis, var.obs, drop = FALSE]
            sigma_om <- Sigma[[g]][var.obs, var.mis, drop = FALSE]
            sigma_mm <- Sigma[[g]][var.mis, var.mis, drop = FALSE]
            sigma_oo_inv <- solve(Sigma[[g]][var.obs, var.obs], tol = .Machine$double.xmin)
            Ximp[pat == t, var.mis, g] <- 
              sweep(t(sigma_mo %*% sigma_oo_inv %*% t(sweep(X[pat == t, var.obs, drop = FALSE], 2, mu_o, "-", check.margin = FALSE))), 2, mu_m, "+")
            sigma.mis <- sigma_mm - sigma_mo %*% sigma_oo_inv %*% sigma_om
            sigma.imp[var.mis, var.mis, pat==t, g] <- sigma.mis
          }
        }
      }
    }
    return(list(Ximp = Ximp,
                sigma.imp = sigma.imp))
  }

################################################################################
impute_mis_best <-
  function(X,
           Ximp,
           W,
           label) {
    n.obs <- dim(X)[1]
    p <- dim(X)[2]
    Ximp.final <- X
    for (i in which(rowSums(W) != p)) {
      var.mis <- W[i, , drop = FALSE] == 0
      Ximp.final[i, var.mis] <- Ximp[i, var.mis, label[i], drop = FALSE]
    }
    return(Ximp.final)
  }

################################################################################
update_W <- 
  function(X,
           G,
           W,
           w,
           pp,
           mu,
           Sigma,
           alpha,
           penalty) {
    n.obs <- dim(W)[1]
    p <- dim(W)[2]
    ordering <- 1:p # Alternatively, for instance, sample(1:p, p)
    for (j in 1:p) {
      unit.obs <- which(!is.na(X[, ordering[j]]))
      Delta <- rep(0, length(unit.obs))
      for (i in 1:length(unit.obs)) {
        if (sum(penalty == 0)) {
          WW1 <- W[unit.obs[i], ]
          WW1[j] <- 1
          var.obs <- WW1 == 1
          lf1 <- rep(0, G)
          lf2 <- rep(0, G)
          for (g in 1:G) {
            mu_o <- mu[g, var.obs, drop = FALSE]
            sigma_oo <- Sigma[[g]][var.obs, var.obs, drop = FALSE]
            lf1[g] <- dmnorm(X[unit.obs[i], var.obs, drop = FALSE], drop(mu_o), sigma_oo) + log(pp[g])
            var.obs_j <- var.obs
            var.obs_j[j] <- FALSE
            if (any(var.obs_j)) {
              mu_o <- mu[g, var.obs_j, drop = FALSE]
              sigma_oo <- Sigma[[g]][var.obs_j, var.obs_j, drop = FALSE]
              lf2[g] <- dmnorm(X[unit.obs[i], var.obs_j, drop = FALSE], drop(mu_o), sigma_oo) + log(pp[g])
            } else {
              lf2[g] <- dmnorm(X[unit.obs[i], , drop = FALSE], drop(mu[g, ]), Sigma[[g]]) + log(pp[g])
            }
          }
          Delta[i] <- log(sum(exp(lf1))) - log(sum(exp(lf2)))
        } else {
          var.obs <- which(W[unit.obs[i], ] == 1)
          var.obs <- var.obs[var.obs != ordering[j]]
          xC.hat <- matrix(as.double(NA), G, 2)
          if (length(var.obs) == 0) {
            for (g in 1:G) {
              xC.hat[g, 1] <- mu[g, ordering[j]]
              xC.hat[g, 2] <- Sigma[[g]][ordering[j], ordering[j]]
              if (xC.hat[g, 2] < .Machine$double.xmin) {
                xC.hat[g, 2] <- (.Machine$double.xmin)
              }
            }
          } else {
            for (g in 1:G) {
              xC.hat[g, 1] <- mu[g, ordering[j]] + Sigma[[g]][ordering[j], var.obs, drop = FALSE]%*%solve(Sigma[[g]][var.obs, var.obs,  drop = FALSE], tol = .Machine$double.xmin)%*%(X[unit.obs[i], var.obs] - mu[g, var.obs])
              xC.hat[g, 2] <- Sigma[[g]][ordering[j], ordering[j]] - Sigma[[g]][ordering[j], var.obs, drop = FALSE]%*%solve(Sigma[[g]][var.obs, var.obs, drop = FALSE], tol = .Machine$double.xmin)%*%Sigma[[g]][var.obs, ordering[j], drop = FALSE]
              if (xC.hat[g, 2] < .Machine$double.xmin) {
                xC.hat[g, 2] <- (.Machine$double.xmin)
              }
            }
          }
          Delta[i] <- (-0.5 * sum(w[unit.obs[i], ]*(log(xC.hat[, 2]) +  ((X[unit.obs[i], ordering[j]] - xC.hat[, 1])^2)/xC.hat[, 2]))) + penalty[i, ordering[j]]
        }
      }
      if (sum(penalty) == 0) {
        cutoff <- sort(Delta, decreasing = TRUE)[ceiling(alpha * n.obs)]
      } else {
        cutoff <- min(sort(Delta, decreasing = TRUE)[ceiling(alpha * n.obs)], 0)
      }
      clean <- unit.obs[Delta >= cutoff]
      out <- unit.obs[Delta < cutoff]
      W[clean, ordering[j]] <- 1
      W[out, ordering[j]] <- 0
    }
    return(W)
  }

################################################################################
update_post <-
  function(X,
           difpat,
           pat,
           pp,
           mu,
           Sigma) {
    G <- dim(mu)[1]
    w <- matrix(as.double(NA), nrow = dim(X)[1], ncol = G)
    for (t in 1:nrow(difpat)){
      var.obs <- difpat[t, , drop = FALSE] == 1
      if (any(var.obs)) {
        for (g in 1:G) {
          mu_o <- mu[g, var.obs, drop = FALSE]
          sigma_oo <- Sigma[[g]][var.obs, var.obs, drop = FALSE]
          w[pat == t, g] <- log(pp[g]) + dmnorm(X[pat == t, var.obs, drop = FALSE], drop(mu_o), sigma_oo)
        }
      } else {
        w[pat == t, ] <-  t(replicate(length(which(pat == t)), log(pp)))
      }
    }
    wnorm <- apply(w, 1, max)
    w <- exp(sweep(w, 1, wnorm, "-"))
    w <- w / rowSums(w)
    w[which(w < sqrt(.Machine$double.eps), arr.ind = TRUE)] <- sqrt(.Machine$double.eps)
    return(w)
  }

################################################################################
update_prior <-
  function(w) {
    pp <- colSums(w) / dim(w)[1]
    return(pp)
  }

################################################################################
update_param <-
  function(Ximp,
           w, 
           sigma.imp,
           maxfact,
           zero.tol,
           label) {
    n.obs <- dim(Ximp)[1]
    p <- dim(Ximp)[2]
    G <- dim(w)[2]
    mu.new <- matrix(0, G, p)
    Sigma.new <- list()
    for (g in 1:G){
      mu.new[g, ] <- sweep(t(w[, g]) %*% Ximp[, , g], 2, 1 / sum(w[, g]), "*")
      Xo <- array(0, dim = c(p, p, n.obs))
      for (i in 1:n.obs) {
        Xo[, , i] <- w[i, g]*(tcrossprod((Ximp[i, , g] - mu.new[g, ])) + sigma.imp[, , i, g])
      }
      Xo <- apply(Xo, 1:2, sum)
      Sigma.temp <-  Xo / sum(w[, g])
      Sigma.new[[g]] <- Sigma.temp
    }
    constr <- restr_diffax(Sigma.new, maxfact = maxfact, zero.tol = zero.tol, label = label)
    Sigma.new <- constr
    return(list(mu = mu.new,
                Sigma = Sigma.new))
  }

################################################################################
restr_diffax <- 
  function (Sigma, 
            maxfact, 
            zero.tol,
            label) {  
    G <- length(Sigma)
    p <- dim(Sigma[[1]])[1]
    csize <- as.numeric(table(label))
    u <- array (NA, dim = c(p, p, G))
    d <- array (NA, dim = c(p, G))
    for (g in 1:G) {
      ev <- eigen(Sigma[[g]])
      u[, , g] <- ev$vectors
      d[, g] <- ev$values
    }
    d[d < 0] <- 0
    d <- restr2_eigenv(d, csize, maxfact, zero.tol)
    if (!(max(d) > zero.tol)) {
      return(Sigma)
    }
    for (g in 1:G) {
      Sigma[[g]] <- u[, , g] %*% diag(d[, g], nrow = p) %*% t(u[, , g])
    }
    return(Sigma)
  }

################################################################################
restr2_eigenv <- 
  function(autovalues,
           csize,
           maxfact,
           zero.tol) {
    c <- maxfact
    d <- t(autovalues)
    p <- dim(autovalues)[1]
    G <- dim(autovalues)[2]
    n.obs <- sum(csize)
    nis <- matrix(data = csize, nrow = G, ncol = p)
    d_ <- sort(c(d, d/c))
    dim <- length(d_)
    d_1 <- d_
    d_1[dim + 1] <- d_[dim]*2
    d_2 <- c(0, d_)
    ed <- (d_1 + d_2)/2
    dim <- dim + 1;
    if ((max(d[nis>0]) <= zero.tol))
      return(matrix(0, nrow = p, ncol = G))
    if (abs(max(d[nis>0])/min(d[nis>0]))<=c) {
      d[nis == 0] <- mean(d[nis>0])
      return(t(d))
    }
    t <- s <- r <- array(0, dim = c(G, dim))
    sol <- sal <- array(0, dim = c(dim))
    for (mp_ in 1:dim){
      for (i in 1:G) {
        r[i, mp_] <- sum((d[i, ]<ed[mp_])) + sum((d[i, ]>ed[mp_]*c))
        s[i, mp_] <- sum(d[i, ]*(d[i, ]<ed[mp_]))
        t[i, mp_] <- sum(d[i, ]*(d[i, ]>ed[mp_]*c))
      }
      sol[mp_] <- sum(csize/n.obs*(s[, mp_] + t[, mp_]/c))/(sum(csize/n.obs*(r[, mp_])))
      e <- sol[mp_]*(d<sol[mp_]) + d*(d>=sol[mp_])*(d<=c*sol[mp_]) + (c*sol[mp_])*(d>c*sol[mp_])
      o <- -1/2*nis/n.obs*(log(e) + d/e)
      sal[mp_] <- sum(o)
    }
    eo <- which.max(c(sal))
    m <- sol[eo]
    t(m*(d<m) + d*(d>=m)*(d<=c*m) + (c*m)*(d>c*m))
  }

################################################################################
dmnorm <- 
  function(X, 
           mu,
           sigma){
    return(log(((2 * pi)^(- length(mu) / 2)) * (det(sigma)^(-1/2)) * exp (-0.5 * mahalanobis(X, mu, sigma))))
  }

################################################################################
loglik_cellGMM <-
  function(X,
           difpat,
           pat,
           pp,
           mu,
           Sigma) {
    n.obs <- dim(X)[1]
    p <- dim(X)[2]
    G <- length(pp)
    lf <- matrix(as.double(NA), n.obs, G)
    for (t in 1:nrow(difpat)){
      var.obs <- difpat[t, , drop = FALSE] == 1
      if (any(var.obs)) {
        for (g in 1:G) {
          mu_o <- mu[g, var.obs, drop = FALSE]
          sigma_oo <- Sigma[[g]][var.obs, var.obs, drop = FALSE]
          lf[pat == t, g] <- dmnorm(X[pat == t, var.obs, drop = FALSE], drop(mu_o), sigma_oo) + log(pp[g])
        }
      } else {
        lf[pat == t, ] <- log(pp)
      }
    }
    loglik <- sum(log(rowSums(exp(lf))))
    return(loglik)
  }

################################################################################
aitken_obj <-
  function(vloglik) {
    if (length(vloglik) < 3) {
      return(Inf)
    } else {
      l.loglik <- length(vloglik)
      loglik.c <- vloglik[l.loglik]
      loglik <- vloglik[l.loglik - 1]
      loglik.p <- vloglik[l.loglik - 2]
      ait.c <- (loglik.c - loglik) / (loglik - loglik.p)
      loglik.inf <- loglik  + ((loglik.c - loglik)/(1 - ait.c))
      conv <- loglik.inf - loglik
      if (conv < 0) {
        conv <- 1
      }
      if (is.nan(conv)) {
        conv <- 0
      }
      return(conv)
    }
  }

################################################################################
restr_eigen <- 
  function(Sigma,
           ceigen) {
    if (min(eigen(Sigma, symmetric = TRUE)$values) < ceigen) {
      Sigma <- Sigma + ceigen * diag(dim(Sigma)[1])
    }
    return(Sigma)
  }

