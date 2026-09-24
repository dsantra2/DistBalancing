
suppressMessages(library(DistBalancing))
 
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript 06_ate_simulation_cluster.R <BATCH>")
BATCH <- as.numeric(args[1])
 
## ---- hyperparameter grid ----------------------------------------------------
## kernel: Isotropic = Matern, ED = energy distance; IPW/CBPS are comparators.
N_GRID      <- c(500, 1000, 2000, 4000)
KERNEL_GRID <- c("Gaussian", "Laplace", "Isotropic", "ED", "IPW", "CBPS")
KAPPA_GRID  <- c(1)
N_REPS      <- 55                     # simulation replicates (SEED = 1..N_REPS)
 
## ---- balancing-objective configurations -----------------------------------
## User-specified moments phi (a function of X; cfd_weights() also accepts a
## precomputed n x m matrix), and the ratio c0 between the moment and the
## distributional term. c0 = 0 (pure distributional balancing) is always fitted.
PHI_GRID <- list(
  linear    = function(X) X,                 # first moments
  quadratic = function(X) cbind(X, X^2)      # first and second moments
)
C0_GRID <- c(0.1, 1, 10)
 
## Data-driven lambda_n = n^2 L(w_IPW)/||w_IPW||^2 * LAMBDA_CONST * n^(-LAMBDA_DECAY)
LAMBDA_CONST <- 0.1
LAMBDA_DECAY <- 0.01
 
## ---- subsampling ------------------------------------------------------------
B_SS           <- 500   # final subsampling replicates per configuration
M_SEARCH_MAX_N <- 1000  # choose m by moonboot::estimate.m() for N <= this;
                         # fix m = round(2*sqrt(N)) above it
M_SEARCH_R     <- 50    # search replicates per candidate m
 
Hyperparameter <- expand.grid(kernel = KERNEL_GRID, kappa = KAPPA_GRID,
                               N = N_GRID, SEED = seq_len(N_REPS),
                               stringsAsFactors = FALSE)
if (BATCH < 1 || BATCH > nrow(Hyperparameter)) {
  stop("BATCH must be between 1 and ", nrow(Hyperparameter))
}
 
kernel_label <- Hyperparameter[BATCH, "kernel"]
kappa        <- Hyperparameter[BATCH, "kappa"]
N            <- Hyperparameter[BATCH, "N"]
SEED         <- Hyperparameter[BATCH, "SEED"]
 
density_name <- switch(kernel_label, Gaussian = "gaussian", Laplace = "laplacian",
                        Isotropic = "matern", ED = "energy", NA_character_)
is_kernel <- !is.na(density_name)
 
time_it <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  val <- expr
  list(value = val, seconds = proc.time()[["elapsed"]] - t0)
}
 
## ---- simulate one dataset --------------------------------------------------
set.seed(SEED)
d <- dgp_ate(kappa, N, p = 10)
 
## ---- per-configuration evaluation ------------------------------------------
na_ci <- c(NA_real_, NA_real_)
 
evaluate <- function(w, fit_seconds) {
  point_est <- ate_weighting_estimate(d$Y, d$A, w)
  pi_res <- time_it(plugin_ci(d$Y, d$A, w))
  ca_res <- time_it(center_adjusted_plugin_ci(d$Y, d$A, d$X, w))
  data.frame(
    point_est = point_est, Fit.Time = fit_seconds,
    Plugin.Lower = pi_res$value$ci[1], Plugin.Upper = pi_res$value$ci[2],
    Plugin.SE = pi_res$value$se, Plugin.Time = pi_res$seconds,
    CenterAdj.Lower = ca_res$value$ci[1], CenterAdj.Upper = ca_res$value$ci[2],
    CenterAdj.SE = ca_res$value$se, CenterAdj.Time = ca_res$seconds
  )
}
 
rows <- list()
 
if (is_kernel) {
  K <- cfd_kernel_gram(d$X, density = density_name)$K
 
  configs <- c(list(list(phi = "none", c0 = 0)),
               unlist(lapply(names(PHI_GRID), function(ph)
                 lapply(C0_GRID, function(c0) list(phi = ph, c0 = c0))),
                 recursive = FALSE))
 
  for (cfg in configs) {
    feats <- if (cfg$phi == "none") NULL else PHI_GRID[[cfg$phi]]
    fit <- time_it(cfd_weights(Z = d$A, X = d$X, K = K, lambda = "ipw",
                               moment_ratio = cfg$c0, moment_features = feats,
                               lambda_const = LAMBDA_CONST, lambda_decay = LAMBDA_DECAY))
    w <- fit$value$w
    ev <- evaluate(w, fit$seconds)
 
    ## Subsampling CI: same specification (kernel, phi, c0, data-driven lambda_n)
    ## re-applied on every subsample; `full_w` avoids re-solving the full sample.
    m_fixed <- if (N <= M_SEARCH_MAX_N) NULL else round(2 * sqrt(N))
    ss_res <- time_it(tryCatch(
      subsample_ci(d$X, d$A, d$Y, K = K, lambda = "ipw", moment_ratio = cfg$c0,
                   moment_features = feats, lambda_const = LAMBDA_CONST,
                   lambda_decay = LAMBDA_DECAY, full_w = w, B = B_SS, m = m_fixed,
                   m_R = M_SEARCH_R, estimand = "ate", seed = SEED),
      error = function(e) { message("subsample_ci failed: ", conditionMessage(e)); NULL }))
    ss <- ss_res$value
    ss_ci   <- if (is.null(ss)) na_ci else unname(ss$ci$ate)
    ss_wald <- if (is.null(ss)) na_ci else unname(ss$wald_ci$ate$ci)
 
    rows[[length(rows) + 1]] <- data.frame(
      kernel = kernel_label, phi = cfg$phi, c0 = cfg$c0, kappa = kappa, N = N,
      SEED = SEED,
      lambda_n = fit$value$lambda, moment_k = fit$value$moment_scale,
      moment_gamma = fit$value$moment_gamma, loss = fit$value$loss,
      ev,
      SS.m = if (is.null(ss)) NA_real_ else ss$best_m$ate,
      SS.Lower = ss_ci[1], SS.Upper = ss_ci[2],
      SS.SE = if (is.null(ss)) NA_real_ else unname(ss$se[["ate"]]),
      SS.Wald.Lower = ss_wald[1], SS.Wald.Upper = ss_wald[2],
      SS.Time = ss_res$seconds
    )
  }
} else {
  fit <- time_it(if (kernel_label == "IPW") ipw_estimate(d$A, d$X)$w else
    cbps_estimate(d$A, d$X)$w)
  rows[[1]] <- data.frame(
    kernel = kernel_label, phi = NA_character_, c0 = NA_real_, kappa = kappa, N = N,
    SEED = SEED, lambda_n = NA_real_, moment_k = NA_real_, moment_gamma = NA_real_,
    loss = NA_real_,
    evaluate(fit$value, fit$seconds),
    SS.m = NA_real_, SS.Lower = NA_real_, SS.Upper = NA_real_, SS.SE = NA_real_,
    SS.Wald.Lower = NA_real_, SS.Wald.Upper = NA_real_, SS.Time = NA_real_
  )
}
 
RESULT <- do.call(rbind, rows)
rownames(RESULT) <- NULL
 
write.csv(RESULT,
          file = sprintf("Result_ATE_%s_kappa%.1f_N%05d_SEED%05d.csv",
                          kernel_label, kappa, N, SEED),
          row.names = FALSE)