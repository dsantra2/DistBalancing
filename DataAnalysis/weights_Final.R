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

Hyperparameter <- expand.grid(
  kernel=c("Gaussian","Laplacian","Energy","TCFD","Matern")
)

kernel <- Hyperparameter[BATCH,"kernel"]

library(DoubleML)
library(mlr3)
library(mlr3learners)
library(data.table)
library(ggplot2)
library(optiSolve)
library(WeightIt)

# suppress messages during fitting
lgr::get_logger("mlr3")$set_threshold("warn")

# load data as a data.table
data = fetch_401k(return_type = "data.table", instrument = TRUE)
# dim(data)
# str(data)

# Initialize DoubleMLData with an instrument
# Basic model
features_base = c("age", "inc", "educ", "fsize", "marr", "twoearn", "db", "pira", "hown")
data_dml_base_iv = DoubleMLData$new(data,
                                    y_col = "net_tfa",
                                    d_cols = "p401",
                                    x_cols = features_base,
                                    z_cols = "e401")

data <- as.data.frame(lapply(data_dml_base_iv$data, function(x) {
  if(is.numeric(x)) {
    (x - min(x)) / (max(x) - min(x))
  } else {
    x  
  }
}))

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



treatment <- data$p401
Xmat <- as.matrix(data[, features_base])

ker_result <- Ker(Xmat, kernel, bandwidth = NULL)
K <- ker_result$K

lambda <- 1/(nrow(Xmat))^2   

weights <- distbalance(treatment, K, lambda)

write.csv(weights, file = sprintf("weights_%s.csv", kernel), row.names = FALSE)

