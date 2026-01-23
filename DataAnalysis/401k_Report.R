library(optiSolve)
library(WeightIt)
library(DoubleML)
library(mlr3)
library(mlr3learners)
library(data.table)
library(purrr)
library(ggplot2)
library(dplyr)
library(tidyr)

data = fetch_401k(return_type = "data.table", instrument = TRUE)

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
print("A vs Z table")
table(data$p401,data$e401)



results <- read.csv("Result_401k_combined.csv",header = TRUE)

results=results%>%
  filter(denom<=1 & denom>=-1)

freq <- table(results$SEED)
results <- results[results$SEED %in% names(freq[freq == 86]), ]

top500_seeds=results %>%
  count(SEED, sort = TRUE) %>%   
  slice_head(n = 500)

results <- results %>%
  filter(SEED %in% top500_seeds$SEED)



data = fetch_401k(return_type = "data.table", instrument = TRUE)

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

methods <- c("Gaussian", "Laplacian", "Energy", "TCFD","Matern")
files <- paste0("weights_", methods, ".csv")

Y <- data_dml_base_iv$data$net_tfa  # outcome
A <- data_dml_base_iv$data$p401     # treatment
Z <- data_dml_base_iv$data$e401     # instrument

i1 <- as.numeric(Z == 1)
i0 <- as.numeric(Z == 0)
n1 <- sum(i1)
n0 <- sum(i0)

compute_est <- function(weight) {
  num <- sum(weight * Z * Y)/n1 - sum(weight * (1 - Z) * Y)/n0
  denom <- sum(weight * Z * A)/n1 - sum(weight * (1 - Z) * A)/n0
  tau <- num / denom
  c(num = num, denom = denom, tau = tau)
}

estimates <- t(sapply(seq_along(methods), function(i) {
  w <- read.csv(files[i])$x
  compute_est(w)
}))

estimates=as.data.frame(estimates)
estimates$kernel <- methods
rownames(estimates) <- NULL
estimates <- estimates[, c("kernel", "num", "denom", "tau")]
estimates

N=9915
X=data[,features_base]
df=data.frame(Z,X)
logit_model <- glm(Z ~ ., df, family = binomial())
prop_score <- predict(logit_model, type = "response")
num_ipw <- (1/N) * (sum(Y * Z / prop_score) - sum(Y * (1 - Z) / (1 - prop_score)))
denom_ipw <- (1/N)*(sum(A * Z / prop_score) - sum(A * (1 - Z) / (1 - prop_score)))
tau_ipw <- num_ipw/denom_ipw

cbps_model=weightit(Z ~., df, method = "cbps", estimand = "ATE")
prop_score <- cbps_model$ps
num_cbps <- (1/N) * (sum(Y * Z / prop_score) - sum(Y * (1 - Z) / (1 - prop_score)))
denom_cbps <- (1/N)*(sum(A * Z / prop_score) - sum(A * (1 - Z) / (1 - prop_score)))
tau_cbps <- num_cbps/denom_cbps




choose_best_m_ss <- function(method="Gaussian", VI_bw=2) {
  df <- results %>%
    filter(kernel == method, Inference == "SS")
  estimates_temp=estimates %>%
    filter(kernel == method)
  
  grid.ss <- sort(unique(df$m))
  length.grid.ss <- length(grid.ss)
  
  CI.num   <- matrix(NA, length.grid.ss, 2)
  CI.denom <- matrix(NA, length.grid.ss, 2)
  CI.late  <- matrix(NA, length.grid.ss, 2)
  
  for (ss in 1:length.grid.ss) {
    df_temp=df %>% filter(m==grid.ss[ss])
    CI.num[ss,] <- quantile(df_temp$num, c(0.025,0.975), na.rm=TRUE)
    CI.denom[ss,] <- quantile(df_temp$denom, c(0.025,0.975), na.rm=TRUE)
    CI.late[ss,] <- quantile(df_temp$late, c(0.025,0.975), na.rm=TRUE)
    
    CI.num[ss,]   <- sqrt(grid.ss[ss]/9915)*(CI.num[ss,]-estimates_temp$num)+estimates_temp$num
    CI.denom[ss,] <- sqrt(grid.ss[ss]/9915)*(CI.denom[ss,]-estimates_temp$denom)+estimates_temp$denom
    CI.late[ss,]  <- sqrt(grid.ss[ss]/9915)*(CI.late[ss,]-estimates_temp$tau)+estimates_temp$tau
  }
  
  VI.calc <- function(CI.mat, ii) {
    sd(CI.mat[max(1, ii-VI_bw):min(length.grid.ss, ii+VI_bw), 1], na.rm=TRUE) +
      sd(CI.mat[max(1, ii-VI_bw):min(length.grid.ss, ii+VI_bw), 2], na.rm=TRUE)
  }
  
  VI.num   <- sapply(1:length.grid.ss, function(ii) VI.calc(CI.num, ii))
  VI.denom <- sapply(1:length.grid.ss, function(ii) VI.calc(CI.denom, ii))
  VI.late  <- sapply(1:length.grid.ss, function(ii) VI.calc(CI.late, ii))
  
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


run_boot <- function(method="Gaussian") {
  df_temp <- results %>% filter(Inference=="Boot",kernel==method)
  
  return(list(
    num_CI   = quantile(df_temp$num,   probs=c(0.025,0.975), na.rm=TRUE),
    denom_CI = quantile(df_temp$denom, probs=c(0.025,0.975), na.rm=TRUE),
    late_CI  = quantile(df_temp$late,  probs=c(0.025,0.975), na.rm=TRUE)
  ))
}

methods <- c("Gaussian","Laplacian","Energy","TCFD","Matern","IPW","CBPS")

boot_results <- map(methods, run_boot) %>% set_names(methods)

subsampling_results <- map(methods[1:5], choose_best_m_ss) %>% set_names(methods[1:5])

