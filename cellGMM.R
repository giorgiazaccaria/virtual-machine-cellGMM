source('InternalFunctions_cellGMM.R')
source('InitializationFunctions_cellGMM.R')

###################################################################################################################################################
# Name: cellGMM
# Gaussian Mixture Models for cellwise outlier detection with constraints and allowing missing data
#
## Description
#
## Arguments:
## X:                  (n.obs x p) numeric data matrix or data frame allowing missing data
## G:                  Number of the mixture components (numeric)
## tuning_param_init:: A data frame of tuning parameters for the initialization: alpha_tclust, alpha_1, alpha_2, alpha.A1, alpha.A2, nrep, nstart, niter (numeric)
## alpha:              Fraction of cells per variable that must remain unflagged (default: 0.75, which represents a lower bound; numeric)
## penalty:            (n.obs x p) numeric matrix representing the L0 penalty term for the update of W and the objective function (default: NULL; numeric)
## quant:              Tuning constant to flag cells (default: 0.99; numeric)
## maxfact:            Constant value greater than 1. Larger values imply larger differences of component covariance matrices; a value of 1 specifies the strongest restriction (default: 30; numeric)
## zero.tol:           Tolerance for the eigenvalues (default: 1e-16; numeric)
## normalization:      Type of normalization for the data : NULL, "standard", "center", "range" "SVD" (default: NULL; string)
## maxiter:            Maximum number of iterations (default: 500; numeric)
## tol:                Tolerance value for convergence (default: 1e-6; numeric)
## stop:              "relative" = relative log-likelihood in two sequential iterations; "aitken" = Aitken acceleration-based stopping rule (default: "aitken"; string)
## rndstart:           Number of random starts (default: 20; numeric)
## manual_initparam:   A list of initial parameters (W, pp, mu, Sigma, post, label) (default: NULL, i.e., it is ignored)
## showprogress:       Progress of the code (default: TRUE; logical)

## Values:
## call:              Matched call.
## X:                 Input data matrix
## X.imputed:         Imputed data matrix for each component
## X.imputed.best:    Imputed data matrix considering the unit-component assignment
## penalty:           (n.obs x p) penalty matrix used in the algorithm
## W:                 (n.obs x p) matrix with zeros corresponding to missing and/or outlying values
## label:             (n.obs x 1) vector of cluster membership, based on the Maximum A Posteriori
## pp:                (1 x G) vector of prior probabilities
## mu:                (G x p) matrix of the component mean vectors
## sigma:             List of G (p x p) component covariance matrices
## post:              (n.obs x G) matrix of the posterior probabilities
## loglik:            Vector of the loglikelihood computed for each iteration of the best loop
## loglik.final:      Final value of the loglikelihood
## init.param:        A list of initial parameters (W, pp, mu, Sigma, post, label) or NULL
## loop:              Loop corresponding to the best model
## iter:              Actual number of iterations needed to reach convergence
##################################################################################################################################################
cellGMM <- function(X,
                    G,
                    tuning_param_init,
                    alpha = 0.75,
                    penalty = NULL,
                    quant = 0.99,
                    maxfact = 30,
                    zero.tol = 1e-16, 
                    normalization = NULL,
                    maxiter = 500,
                    tol = 1e-6,
                    stop = "aitken",
                    rndstart = 20,
                    manual_initparam = NULL,
                    showprogress = TRUE) {
  call <- mget(names(formals())[-c(1, 14)],sys.frame(sys.nframe()))
  # PRE-PROCESSING 
  if (!is.null(normalization)) {
    X <- norm(X, normalization)
  } else {
    X <- as.matrix(X)
  }
  n.obs <- dim(X)[1]
  p <- dim(X)[2]
  if (G < 1 || G > n.obs || G%%1 != 0) {
    stop("G is not properly fixed: G must be an integer chosen into [1, dim(X)[1]]", call. = TRUE)
  }
  if (is.null(penalty) || alpha == 1) {
    penalty <- matrix(0, n.obs, p)
  }
  # STARTING
  for (loop in 1:rndstart) {
    count <- 0
    conv <- 1
    wp <- 0
    # INITIALIZATION 
    if (is.null(manual_initparam)) {
      init <- initfunc_new(X = X, G = G, maxfact = maxfact, zero.tol = zero.tol, 
                           alpha_tclust = tuning_param_init$alpha_tclust, 
                           alpha_1 = tuning_param_init$alpha_1,
                           alpha_2 = tuning_param_init$alpha_2,
                           alpha.A1 = tuning_param_init$alpha.A1, 
                           alpha.A2 = tuning_param_init$alpha.A2,
                           nrep = tuning_param_init$nrep,
                           nstart = tuning_param_init$nstart,
                           niter = tuning_param_init$niter)
    } else {
      init <- manual_initparam
      if (sum(penalty) == 0 && alpha < 1) {
        penalty <- matrix(0, n.obs, p)
        for (i in 1:n.obs) {
          logEst <- matrix(0, 1, p)
          for (g in 1:G) {
            logEst <- logEst + init$post[i, g]*log(diag(solve(init$Sigma[[g]])))
          }
          penalty[i, ] <- 0.5 * (- logEst + qchisq(p = 0.99, df = 1))
        }
      }
    }
    model <- init[-which(names(init) == "label")]
    label <- init$label
    dif.pattern <- model$W[!duplicated(model$W), , drop = FALSE]
    pat.unit <- matrix(0, n.obs, 1)
    for (t in 1:nrow(dif.pattern)) {
      pat.unit[apply(model$W, 1, identical, dif.pattern[t, ]), ] <- t
    }
    model$post <- update_post(X, dif.pattern, pat.unit, model$pp, model$mu, model$Sigma)
    if (showprogress) {
      pb <- txtProgressBar(min = 0, max = rndstart, style = 3, 
                           width = 75, char = "=")
    }
    # ITERATIONS    
    while (conv > tol && count <= maxiter) {
      if (count == 0) {
        if (sum(unique(unlist(label))) != sum(1:G)) {
          wp <- 1
          loglik <- -.Machine$double.xmax
          break
        }
        loglik <- loglik_cellGMM(X, dif.pattern, pat.unit, model$pp, model$mu, model$Sigma) - sum(penalty * (1 - model$W))
        vloglik <- loglik
      }
      count <- count + 1
      ## E-STEP 
      # UPDATE W
      if (alpha < 1) {
        model$W <- update_W(X = X, G = G, W = model$W, w = model$post, pp = model$pp, mu= model$mu, Sigma = model$Sigma, alpha = alpha, penalty = penalty)
        dif.pattern <- model$W[!duplicated(model$W), , drop = FALSE]
        pat.unit <- matrix(0, n.obs, 1)
        for (t in 1:nrow(dif.pattern)) {
          pat.unit[apply(model$W, 1, identical, dif.pattern[t, ]), ] <- t
        }
      }
      # UPDATE Z AND X[W^{c}]   
      model$post <- update_post(X, dif.pattern, pat.unit, model$pp, model$mu, model$Sigma)
      label <- apply(model$post, 1, which.max)
      if (sum(unique(unlist(label))) != sum(1:G)) {
        wp <- 1
        loglik.c <- -.Machine$double.xmax
        break
      }
      val.impute <- impute_mis(X, dif.pattern, pat.unit, model$pp, model$mu, model$Sigma)
      model$X.imputed <- val.impute$Ximp
      # M-STEP   
      model$pp <- update_prior(model$post)
      comp.param <- update_param(Ximp = val.impute$Ximp, w = model$post, sigma.imp = val.impute$sigma.imp, maxfact = maxfact, zero.tol = zero.tol, label = label)
      model$mu <- comp.param$mu
      model$Sigma <- comp.param$Sigma
      # OBJECTIVE FUNCTION 
      loglik.c <- loglik_cellGMM(X, dif.pattern, pat.unit, model$pp, model$mu, model$Sigma) - sum(penalty * (1 - model$W))
      vloglik <- c(vloglik, loglik.c)
      if (stop == "aitken") {
        conv <- aitken_obj(vloglik)
      } else if (stop == "relative") {
        conv <- (loglik.c - loglik) / abs(loglik)
        loglik <- loglik.c
      }
    }
    if (loop == 1) {
      W.best <- model$W
      label.best <- label
      model.best <-  model[-which(names(model) == "W")]
      loglik.best <- vloglik
      loglik.final.best <- loglik.c
      init.best <- init
      loop.best <- loop
      iter.best <- count
      if (wp == 0) {
        X.imputed.best <- impute_mis_best(X = X, Ximp = model.best$X.imputed, W = W.best, label = label.best)
      } else if (wp == 1) {
        X.imputed.best <- NULL
      }
    } else {
      if (loglik.c > loglik.final.best) {
        W.best <- model$W
        label.best <- label
        model.best <- model[-which(names(model) == "W")]
        loglik.best <- vloglik
        loglik.final.best <- loglik.c
        init.best <- init
        loop.best <- loop
        iter.best <- count 
        if (wp == 0) {
          X.imputed.best <- impute_mis_best(X = X, Ximp = model.best$X.imputed, W = W.best, label = label.best)
        } else if (wp == 1) {
          X.imputed.best <- NULL
        }
      }
    }
    if (G == 1) {
      break
    }
    if (showprogress) {
      setTxtProgressBar(pb, loop)
      cat("Loop", loop, "/", rndstart)
    }
  } # END LOOP
  if (showprogress) {
    close(pb)
  }
  if (loglik.final.best == -.Machine$double.xmax) {
    print("The solution has a number of cluster < G and the objective is not computed.")
  }
  return(
    list(
      call = call,
      X = X,
      X.imputed = model.best$X.imputed,
      X.imputed.best = X.imputed.best,
      penalty = penalty,
      W = W.best,
      label = label.best,
      pp = model.best$pp,
      mu = model.best$mu,
      sigma = model.best$Sigma,
      post = model.best$post,
      loglik = loglik.best,
      loglik.final = loglik.final.best,
      init.param = init.best,
      loop = loop.best,
      iter = iter.best
    )
  )
} # END FUNCTION