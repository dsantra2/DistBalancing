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
#BATCH=1

expit <- function(v){exp(v)/(1+exp(v))}

K.T.CFD <- function(Xmat,bandwidth) {
  ## suppose (V_1,...,V_10) ~ t(5)
  V.random <- matrix(rt(n=ncol(Xmat)*10000,
                        df=5,
                        ncp=0)/bandwidth, ## Same bandwidth as the Laplacian Kernel
                     10000,
                     ncol(Xmat))
  
  PX <- V.random %*% t(Xmat)
  
  # Want: cos(P_i - P_j) averaged over all samples
  # cos(a - b) = cos(a)cos(b) + sin(a)sin(b)
  C <- cos(PX)
  S <- sin(PX) 
  K <- (t(C) %*% C + t(S) %*% S)/10000
  
  return(K)
}

Ker <- function(Xmat, name, c=1, beta=0.5, nu=1, l=1, bandwidth=NULL) {
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
  } else if (name == "TCFD") {
    D1 <- as.matrix(dist(Xmat, method = "manhattan"))
    sigma <- ifelse(is.null(bandwidth), median(D1), bandwidth)
    K <- K.T.CFD(Xmat,sigma)
    bandwidth <- sigma
  }else if (name == "inv_multiquadric") {
    K <- (D^2 + c^2)^(-beta)
    #} else if (name == "Matern") {
    #  K <- 2^(1-nu)*(sqrt(2*nu)*D/l)^nu*besselK(sqrt(2*nu)*D/l, nu) / base::gamma(nu)
    #  diag(K) <- 0
  } else if (name == "Matern") {
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

distbalance <- function(treatment, K, lambda){ 
  N <- length(treatment)
  Z <- treatment
  
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
               lc = lincon(A = Sum.Mat, dir = rep("==", 2), val = d.vec,
                           name=c("R1","R2")))
  
  weights <- solvecop(mycop, quiet = TRUE)$x
  return(weights)
}

library(optiSolve)
library(WeightIt)
library(DoubleML)
library(mlr3)
library(mlr3learners)
library(data.table)



## Real data ##

data = fetch_401k(return_type = "data.table", instrument = TRUE)

features_base = c("age", "inc", "educ", "fsize", "marr", "twoearn", "db", "pira", "hown")
data_dml = DoubleMLData$new(data,
                            y_col = "net_tfa",
                            d_cols = "p401",
                            x_cols = features_base,
                            z_cols = "e401")

data_scaled <- as.data.frame(lapply(data_dml$data, function(x) {
  if(is.numeric(x)) (x - min(x))/(max(x)-min(x)) else x
}))

Y_full <- data_dml$data$net_tfa
A_full <- data_dml$data$p401
Z_full <- data_dml$data$e401
X_full <- data_scaled[,features_base]

N <- 9915
range.ss <- c(sqrt(N), 3*sqrt(N))
grid.ss <- unique(round(seq(range.ss[1], range.ss[2], length=20)))

Hyperparameter <- expand.grid( kernel=c("Gaussian","Laplacian","Energy","TCFD","IPW","CBPS","Matern"),
                               resample.N=c(grid.ss,N),
                               SEED=1:750)

Hyperparameter$inference <- ifelse(Hyperparameter$resample.N==N,"Boot","SS")
Remove <- which(Hyperparameter$kernel%in%c("IPW","CBPS") &
                  Hyperparameter$inference=="SS")
Hyperparameter <- Hyperparameter[-Remove,]

kernel <- Hyperparameter[BATCH,"kernel"]
inference <- Hyperparameter[BATCH, "inference"]
Subsamplesize <- Hyperparameter[BATCH, "resample.N"]
SEED <- Hyperparameter[BATCH, "SEED"]
set.seed(SEED)

#Hyperparameter <- read.csv("Index_of_Rerun.csv")
#kernel <- Hyperparameter[BATCH,"kernel"]
#inference <- Hyperparameter[BATCH, "Inference"]
#Subsamplesize <- Hyperparameter[BATCH, "m"]
#SEED <- Hyperparameter[BATCH, "SEED"]
#set.seed(SEED)



## One resample per job ##
N <- length(Y_full)
if (inference == "Boot") {
  idx <- sample.int(N, N, replace = TRUE)             
} else {
  idx <- sample.int(N, Subsamplesize, replace = FALSE)
}

Y <- Y_full[idx]; A <- A_full[idx]; Z <- Z_full[idx]
X <- X_full[idx, , drop = FALSE]
m <- length(idx)

lambda <- 1 / (N^2)

if (kernel %in% c("Gaussian", "Laplacian", "Energy", "TCFD","Matern")) {
  KK <- Ker(X, name = kernel, bandwidth = NULL)
  K  <- KK$K
  weight  <- tryCatch(distbalance(Z, K, lambda = lambda),
                      error = function(e) rep(NA, m))
  i1 <- as.numeric(Z == 1)
  i0 <- as.numeric(Z == 0)
  n1 <- sum(i1)
  n0 <- sum(i0)
  num <- sum(weight * Z * Y)/n1 - sum(weight * (1 - Z) * Y)/n0
  denom <- sum(weight * Z * A)/n1 - sum(weight * (1 - Z) * A)/n0
  late <- num / denom
} else if (kernel == "IPW") {
  df <- data.frame(Z = Z, X)
  logit_model <- glm(Z ~ ., df, family = binomial())
  prop_score <- predict(logit_model, type = "response")
  num <- (1/m) * (sum(Y * Z / prop_score) - sum(Y * (1 - Z) / (1 - prop_score)))
  denom <- (1/m)*(sum(A * Z / prop_score) - sum(A * (1 - Z) / (1 - prop_score)))
  late <- num/denom
} else if (kernel == "CBPS") {
  df <- data.frame(Z = Z, X)
  cbps_model=weightit(Z ~., df, method = "cbps", estimand = "ATE")
  prop_score <- cbps_model$ps
  num <- (1/m) * (sum(Y * Z / prop_score) - sum(Y * (1 - Z) / (1 - prop_score)))
  denom <- (1/m)*(sum(A * Z / prop_score) - sum(A * (1 - Z) / (1 - prop_score)))
  late <- num/denom
} else {
  stop("Unknown kernel.")
}

RESULT <- data.frame(
  kernel       = kernel,
  m            = m,     # Subsamplesize
  SEED         = SEED,
  Inference    = inference,
  num          = num,
  denom        = denom,
  late         = late
)

write.csv(RESULT,
          file = sprintf("Result_401k_%s_N%05d_%s_SEED%05d.csv",
                         kernel, m, inference, SEED), row.names = FALSE)
