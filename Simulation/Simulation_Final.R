args <- (commandArgs(trailingOnly=TRUE))
cat(args[1])
if(length(args) == 1){
  # BATCH = 1, 2, ..., the value of "queue" in your .sub file
  BATCH <- as.numeric(args[1])
  set.seed(BATCH)
} else {
  stop()
}

#rm(list = ls())
# Default for testing if not running from command line
#if(!exists("BATCH")) BATCH=1

Hyperparameter <- expand.grid(
  kernel=c("Gaussian","Laplacian","Energy","IPW","CBPS","TCFD","Matern"),
  DGP.type = c(1, 2),
  N     = c(100, 200, 400, 800, 1600),
  inference = c("SS","Boot"),
  SEED  = 1:550
)

kernel <- Hyperparameter[BATCH,"kernel"]
DGP.type <- Hyperparameter[BATCH, "DGP.type"]  
N <- Hyperparameter[BATCH, "N"]
inference <- Hyperparameter[BATCH, "inference"]
SEED <- Hyperparameter[BATCH, "SEED"]


library(optiSolve)
library(WeightIt)

expit <- function(v){exp(v)/(1+exp(v))}

K.T.CFD <- function(Xmat,bandwidth){
  ## suppose (V_1,...,V_10) ~ t(5)
  V.random <- matrix(rt(n=ncol(Xmat)*10000,
                        df=5,
                        ncp=0)/bandwidth, ## Same bandwidth as the Laplacian Kernel
                     10000,
                     ncol(Xmat))
  PX <- Xmat %*% t(V.random) 
  K <- matrix(0, nrow(Xmat), nrow(Xmat))
  for (ii in 1:nrow(Xmat)) {
    diffs <- sweep(PX, 2, PX[ii,], "-")   
    K[ii, ] <- apply(cos(diffs),1,mean)
  } 
  return(K)
}

Ker <- function(Xmat, name, c=1, beta=0.5, nu=2.5, l=1, bandwidth=NULL) {
  D <- as.matrix(dist(Xmat))
  if (is.null(bandwidth)) {
    gamma <- median(D)
  } else {
    gamma <- bandwidth
  }
  if (name == "Gaussian") {
    K <- exp(-D^2 / gamma^2)
    bandwidth=gamma
  } else if (name == "Energy") {
    K <- -D
  } else if (name == "Laplacian") {
    D1 <- as.matrix(dist(Xmat, method = "manhattan"))
    sigma <- ifelse(is.null(bandwidth), median(D1), bandwidth)
    K <- exp(-D1 / sigma)
    bandwidth=sigma
  }else if (name == "TCFD") {
    D1 <- as.matrix(dist(Xmat, method = "manhattan"))
    sigma <- ifelse(is.null(bandwidth), median(D1), bandwidth)
    K <- K.T.CFD(Xmat,sigma)
    bandwidth <- sigma
  } else if (name == "inv_multiquadric") {
    K <- (D^2 + c^2)^(-beta)
  } else if (name == "Matern") {
    # eps <- 1e-8
    # D_safe <- pmax(D, eps)
    # K <- 2^(1-nu)*(sqrt(2*nu)*D_safe/l)^nu*besselK(sqrt(2*nu)*D_safe/l, nu) / base::gamma(nu)
    # diag(K) <- 0
    r <- D
    if (is.null(bandwidth)) {
      l <- median(D)
    } else {
      l <- bandwidth
    }
    K <- (1 + sqrt(5) * r / l + 5 * (r^2) / (3 * l^2)) * exp(-sqrt(5) * r / l)
    diag(K) <- 1
    bandwidth <- l
  } else {
    stop("Unsupported kernel name.")
  }
  return(list(K=K,bandwidth=bandwidth))
}

DGP <- function(DGP.type=1, N, p=10) {
  muY <- function(a, X, u) {
    coef <- rep(c(1,-0.5),each=5)
    a*(1 + X %*% coef - 1*u) +
      X %*% (-coef/2) + 0.5*u
  }
  psZ <- function(DGP.type=1,X){
    if(DGP.type==1){
      coef_z <- rep(0.15,10)
      psz <- expit( X %*% coef_z )
    }
    else if(DGP.type==2){
      coef_z <- rep(0.25,10)
      psz <- expit( abs(X) %*% coef_z - 2.25 )
    }
    return(psz)
  }
  L <- function(X){
    gamma_0 <- 0.25
    gamma_x <- rep(c(-0.5,0.5),each=5)
    L0 <- -gamma_0 + X %*% gamma_x - 0.5 * abs(U)
    L1 <-  gamma_0 + X %*% gamma_x + 0.5 * abs(U)
    output <- list(L0=L0,L1=L1)
  }
  X <- matrix(rnorm(N * p), nrow = N)
  U <- rnorm(N)
  Z <- rbinom(N, 1, psZ(DGP.type,X))
  LL <- L(X)
  
  A0 <- as.numeric(LL$L0 >= 0)
  A1 <- as.numeric(LL$L1 >= 0)
  A <- A0*(1 - Z) + A1*Z
  
  potY0 <- muY(0, X, U) + rnorm(N) * 0.25
  potY1 <- muY(1, X, U) + rnorm(N) * 0.25
  Y <- potY0*(1 - A) + potY1*A
  potYz1 <- potY1*A1 + potY0*(1-A1)
  potYz0 <- potY1*A0 + potY0*(1-A0)
  
  Denom=mean((A1-A0))                 
  Numer=mean((potYz1-potYz0))             
  LATE_TRUE=mean((potY1-potY0)[A1>A0])
  out=list(X = X, A = A, Y = Y, Z = Z,
           denominator=Denom, numerator= Numer, LATE=LATE_TRUE)
  return(out)
}

distbalance <- function(treatment, covariate, name, bandwidth, lambda, K=NULL){ 
  N=nrow(covariate)
  Z=treatment
  X=covariate
  
  if(is.null(K)){
    K <- Ker(X, name, bandwidth = bandwidth)$K 
  }
  
  i1 <- as.numeric(Z == 1)
  i0 <- as.numeric(Z == 0)
  n1 <- sum(i1)
  n0 <- sum(i0)
  
  A_diag <- diag(i1)
  B_diag <- diag(i0)
  I_N <- diag(N)
  one_vec <- rep(1, N)
  
  Q.mat <- (1 / n1^2) * (A_diag %*% K %*% A_diag) +
    (1 / n0^2) * (B_diag %*% K %*% B_diag) -
    (1 / (n1 * n0)) * (A_diag %*% K %*% B_diag) +
    lambda^2 * I_N
  
  c.vec <- - (1 / (n1 * N)) * (A_diag %*% K %*% one_vec) -
    (1 / (n0 * N)) * (B_diag %*% K %*% one_vec)
  
  Sum.Mat <- matrix(0, 2, N)
  Sum.Mat[1, which(Z == 1)] <- 1
  Sum.Mat[2, which(Z == 0)] <- 1
  d.vec <- c(n1, n0)
  
  mycop <- cop(f = quadfun(Q.mat, a = as.numeric(c.vec), d = 0),
               lb = lbcon(val = rep(0, N)),
               lc = lincon(A = Sum.Mat, dir = rep("==", 2), val = d.vec
                           ,name=c("R1","R2")))
  
  weights <- solvecop(mycop, quiet = TRUE)$x
  return(weights)
}

sim <- function(DGP.type=1, N, p=10, name, seed=1, return_K=FALSE, bandwidth=NULL, K=NULL) {
  set.seed(seed)
  lambda <- 1 / (N^2)
  data_gen=DGP(DGP.type,N,p) 
  X=data_gen$X
  A=data_gen$A
  Y=data_gen$Y
  Z=data_gen$Z
  
  if(name%in%c("Gaussian","Laplacian","Energy","TCFD","Matern")){
    # Compute kernel
    if(is.null(K)){
      K <- Ker(X, name, bandwidth = bandwidth)$K 
    }
    if(is.null(bandwidth)){
      bw <- Ker(X, name, bandwidth = bandwidth)$bandwidth
    } else {
      bw <- bandwidth
    }
    
    weight <- distbalance(Z, X, name=name, bandwidth = bw, lambda = lambda, K=K)
    
    i1 <- as.numeric(Z == 1)
    i0 <- as.numeric(Z == 0)
    n1 <- sum(i1)
    n0 <- sum(i0)
    
    num <- sum(weight * Z * Y)/n1 - sum(weight * (1 - Z) * Y)/n0
    denom <- sum(weight * Z * A)/n1 - sum(weight * (1 - Z) * A)/n0
    tau <- num / denom
  }else if(name%in%c("IPW")){
    df=data.frame(Z,X)
    logit_model <- glm(Z ~ ., df, family = binomial())
    prop_score <- predict(logit_model, type = "response")
    num <- (1/N) * (sum(Y * Z / prop_score) - sum(Y * (1 - Z) / (1 - prop_score)))
    denom <- (1/N)*(sum(A * Z / prop_score) - sum(A * (1 - Z) / (1 - prop_score)))
    tau <- num/denom
  }else{
    df=data.frame(Z,X)
    cbps_model=weightit(Z ~., df, method = "cbps", estimand = "ATE")
    prop_score <- cbps_model$ps
    num <- (1/N) * (sum(Y * Z / prop_score) - sum(Y * (1 - Z) / (1 - prop_score)))
    denom <- (1/N)*(sum(A * Z / prop_score) - sum(A * (1 - Z) / (1 - prop_score)))
    tau <- num/denom
  }
  
  out <- data.frame(DGP.type=DGP.type, N = N, seed = seed, denominator = denom,
                    numerator = num,LATE = tau)
  # Only return K/lambda if it was a kernel method, otherwise return placeholders if needed
  if (return_K) return(list(result = out, K = K, bandwidth=bw, lambda = lambda, X = X, A = A, Y = Y, Z = Z))
  return(list(result = out, K = NULL, bandwidth=NULL, lambda = NULL, X = X, A = A, Y = Y, Z = Z))
}


# Subsampling ######################

choose_best_m_ss <- function(X, A, Y, Z, name, lambda, 
                             bandwidth=NULL, reps=500, VI_bw=2, K=NULL) {
  N <- length(Z)
  range.ss <- c(sqrt(N), 3 * sqrt(N))
  grid.ss <- unique(round(seq(range.ss[1], range.ss[2], length = 20)))
  length.grid.ss <- length(grid.ss)
  
  CI.num   <- matrix(NA, length.grid.ss, 2)
  CI.denom <- matrix(NA, length.grid.ss, 2)
  CI.late  <- matrix(NA, length.grid.ss, 2)
  
  for (ss in 1:length.grid.ss) {
    m <- grid.ss[ss]
    num_temp <- numeric(reps)
    denom_temp <- numeric(reps)
    late_temp <- numeric(reps)
    
    for (subsampling.index in 1:reps) {
      idx <- sample(1:N, m, replace = FALSE)
      
      # Subsampling is only used for Kernel methods in this logic
      weight <- tryCatch(
        distbalance(Z[idx], X[idx,,drop=FALSE], name=name, bandwidth=bandwidth, lambda=lambda, K=K[idx,idx]),
        error = function(e) rep(NA, m)
      )
      if (any(is.na(weight))) next
      
      n1 <- sum(Z[idx]==1)
      n0 <- sum(Z[idx]==0)
      if (n1==0 || n0==0) next
      
      num <- sum(weight * Z[idx] * Y[idx])/n1 - sum(weight * (1-Z[idx]) * Y[idx])/n0
      denom <- sum(weight * Z[idx] * A[idx])/n1 - sum(weight * (1-Z[idx]) * A[idx])/n0
      late <- num / denom
      
      num_temp[subsampling.index]   <- num
      denom_temp[subsampling.index] <- denom
      late_temp[subsampling.index]  <- late
    }
    
    CI.num[ss,]   <- quantile(num_temp,   c(0.025,0.975), na.rm=TRUE)
    CI.denom[ss,] <- quantile(denom_temp, c(0.025,0.975), na.rm=TRUE)
    CI.late[ss,]  <- quantile(late_temp,  c(0.025,0.975), na.rm=TRUE)
  }
  
  # Variability index for each effect type
  VI.calc <- function(CI.mat, ii) {
    sd(CI.mat[max(1, ii-VI_bw):min(length.grid.ss, ii+VI_bw), 1], na.rm=TRUE) +
      sd(CI.mat[max(1, ii-VI_bw):min(length.grid.ss, ii+VI_bw), 2], na.rm=TRUE)
  }
  
  VI.num   <- sapply(1:length.grid.ss, function(ii) VI.calc(CI.num, ii))
  VI.denom <- sapply(1:length.grid.ss, function(ii) VI.calc(CI.denom, ii))
  VI.late  <- sapply(1:length.grid.ss, function(ii) VI.calc(CI.late, ii))
  
  # Best m for each
  best_idx_num   <- which.min(VI.num)
  best_idx_denom <- which.min(VI.denom)
  best_idx_late  <- which.min(VI.late)
  
  return(list(
    VI.values   = list(num=VI.num, denom=VI.denom, late=VI.late),
    subsample.size = grid.ss,
    best_m      = list(num   = grid.ss[best_idx_num],
                       denom = grid.ss[best_idx_denom],
                       late  = grid.ss[best_idx_late]),
    num_CI      = CI.num[best_idx_num,],
    denom_CI    = CI.denom[best_idx_denom,],
    late_CI     = CI.late[best_idx_late,]
  ))
}



# Bootstrap ################################

run_boot <- function(X, A, Y, Z, name, lambda=NULL, reps=500, bandwidth=NULL, K=NULL) {
  N <- length(Z)
  num_vec <- numeric(reps)
  denom_vec <- numeric(reps)
  late_vec <- numeric(reps)
  
  is_kernel <- name %in% c("Gaussian","Laplacian","Energy","TCFD","Matern")
  
  for (i in 1:reps) {
    idx <- sample(1:N, N, replace=TRUE)
    
    # Store bootstrapped data for easy access
    Z_b <- Z[idx]
    X_b <- X[idx,,drop=FALSE]
    Y_b <- Y[idx]
    A_b <- A[idx]
    
    if (is_kernel) {
      weight <- tryCatch(
        distbalance(Z_b, X_b, name=name, bandwidth=bandwidth, lambda=lambda, K=K[idx,idx]),
        error = function(e) rep(NA, N)
      )
      if (any(is.na(weight))) next
      
      n1 <- sum(Z_b==1)
      n0 <- sum(Z_b==0)
      if (n1==0 || n0==0) next
      
      num <- sum(weight * Z_b * Y_b)/n1 - sum(weight * (1-Z_b) * Y_b)/n0
      denom <- sum(weight * Z_b * A_b)/n1 - sum(weight * (1-Z_b) * A_b)/n0
      late <- num / denom
      
    } else if (name == "IPW") {
      df_b <- data.frame(Z=Z_b, X_b)
      # Check for singularity or separation
      logit_model <- tryCatch(glm(Z ~ ., df_b, family = binomial()), 
                              error = function(e) NULL, warning = function(w) NULL)
      if(is.null(logit_model)) next
      
      prop_score <- predict(logit_model, type = "response")
      
      # Handle extreme prop scores to avoid Inf
      prop_score <- pmax(pmin(prop_score, 1 - 1e-6), 1e-6)
      
      num <- (1/N) * (sum(Y_b * Z_b / prop_score) - sum(Y_b * (1 - Z_b) / (1 - prop_score)))
      denom <- (1/N) * (sum(A_b * Z_b / prop_score) - sum(A_b * (1 - Z_b) / (1 - prop_score)))
      late <- num/denom
      
    } else if (name == "CBPS") {
      df_b <- data.frame(Z=Z_b, X_b)
      # Wrap weightit in tryCatch as optimization can fail on bootstrap samples
      cbps_model <- tryCatch(
        weightit(Z ~., df_b, method = "cbps", estimand = "ATE"),
        error = function(e) NULL
      )
      if(is.null(cbps_model)) next
      
      prop_score <- cbps_model$ps
      # Handle extreme prop scores
      prop_score <- pmax(pmin(prop_score, 1 - 1e-6), 1e-6)
      
      num <- (1/N) * (sum(Y_b * Z_b / prop_score) - sum(Y_b * (1 - Z_b) / (1 - prop_score)))
      denom <- (1/N) * (sum(A_b * Z_b / prop_score) - sum(A_b * (1 - Z_b) / (1 - prop_score)))
      late <- num/denom
    }
    
    num_vec[i] <- num
    denom_vec[i] <- denom
    late_vec[i] <- late
  }
  return(list(
    num_CI   = quantile(num_vec,   probs=c(0.025,0.975), na.rm=TRUE),
    denom_CI = quantile(denom_vec, probs=c(0.025,0.975), na.rm=TRUE),
    late_CI  = quantile(late_vec,  probs=c(0.025,0.975), na.rm=TRUE)
  ))
}


# Main computation #####################################

run_all <- function(DGP.type,kernel,inference) {
  
  # Run simulation for all methods to get point estimates
  # return_K is only relevant for kernels, but sim handles it gracefully now
  sim_res <- sim(DGP.type, N, 10, kernel, seed = SEED, return_K = TRUE)
  
  # Initialize result placeholders
  best_m_num <- NA
  best_m_denom <- NA
  best_m_late <- NA
  CI.num <- c(NA, NA)
  CI.denom <- c(NA, NA)
  CI.late <- c(NA, NA)
  
  is_kernel <- kernel %in% c("Gaussian","Laplacian","Energy","TCFD","Matern")
  
  if (inference == "Boot") {
    boot_CI <- run_boot(
      X = sim_res$X, A = sim_res$A, Y = sim_res$Y, Z = sim_res$Z,
      name = kernel, lambda = sim_res$lambda, reps=500, bandwidth=sim_res$bandwidth,
      K=sim_res$K
    )
    CI.num   <- boot_CI$num_CI
    CI.denom <- boot_CI$denom_CI
    CI.late  <- boot_CI$late_CI
    # For bootstrap, m is just N
    best_m_num <- N
    best_m_denom <- N
    best_m_late <- N
    
  } else if (inference == "SS" && is_kernel) {
    # Subsampling only for kernels
    ss_res <- choose_best_m_ss(
      X = sim_res$X, A = sim_res$A, Y = sim_res$Y, Z = sim_res$Z,
      name = kernel, lambda = sim_res$lambda, reps=500, VI_bw=2, bandwidth=sim_res$bandwidth,
      K=sim_res$K
    )
    CI.num   <- ss_res$num_CI
    CI.denom <- ss_res$denom_CI
    CI.late  <- ss_res$late_CI
    best_m_num   <- ss_res$best_m$num
    best_m_denom   <- ss_res$best_m$denom
    best_m_late   <- ss_res$best_m$late
  }
  
  # Construct output dataframe
  df <- list(
    kernel=kernel,
    DGP.type=DGP.type,
    N = sim_res$result$N,
    SEED = sim_res$result$seed,
    Inference = inference,
    m_num = best_m_num,
    num=sim_res$result$numerator,
    Lower.num=CI.num[1],
    Upper.num=CI.num[2],
    m_denom = best_m_denom,
    denom=sim_res$result$denominator,
    Lower.denom=CI.denom[1],
    Upper.denom=CI.denom[2],
    m_late = best_m_late,
    late=sim_res$result$LATE,
    Lower.late=CI.late[1],
    Upper.late=CI.late[2]
  )
  
  return(df)
}

RESULT <- run_all(DGP.type, kernel, inference)

write.csv(RESULT,
          file = sprintf("Result_optiSolve_%s_DGP%0.5d_N%0.5d_%s_SEED%0.5d.csv",
                         kernel, DGP.type, N,inference, SEED), row.names = FALSE
)

