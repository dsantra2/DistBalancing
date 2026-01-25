# Replication Files for "Distributional Balancing for Causal Inference: A Unified Framework via Characteristic Function Distance"
This supplementary file contains replication codes for "Distributional Balancing for Causal Inference: A Unified Framework via Characteristic Function Distance" [Santra, Chen, Park, 2026](https://doi.org/10.48550/arXiv.2601.15449)[^1].
## Software and Packages
* Software: R version 4.4.3
* Packages: _optiSolve_ (version 1.0); _WeightIt_ (version 1.4.0); _DoubleML_ (version 1.0.2); _mlr3_ (version 1.0.1); _mlr3learners_ (version 0.12.0); _data.table_ (version 1.17.0); _ggplot2_ (version 3.5.2); _purrr_ (version 1.0.4); _dplyr_ (version 1.1.4); _tidyr_ (version 1.3.1)
# Simulation Folder
The Simulation folder contains replication files for the simulation studies presented in Section 6 of the main paper.
## Simulation Setup
For the simulations, we have considered:

* **Dimension of the data:** $p=10$. $X_i \sim [N(0,1)]^{10}$ for all $i \in \{1, 2, \dots, N\}$.
* **Unobserved confounder:** $U \sim N(0,1)$.
* **Instrumental variable ($Z$):**
  * **Linear:** $Z_i \sim \text{Bernoulli}(\text{expit}(X_i \cdot \mathbf{v}))$, where $\mathbf{v}$ is a vector $0.15 \cdot \mathbf{1}_{10}$ and $\text{expit}(x) = \frac{e^x}{1+e^x}$.
  * **Non-Linear:** $Z_i \sim \text{Bernoulli} \left(\text{expit} \left(\frac{1}{4} |X_i|[1:5] + \frac{1}{4} X_i[6:10] - 2.25 \right)\right)$.
* **Treatment indicators ($A(0), A(1), A$):** The coefficients ($\gamma_0, \gamma_1$, and $\pm 0.5$) are chosen such that the monotonicity assumption (i.e., no defier) is satisfied.
  * $L_0 = \gamma_0 + X_i \cdot \gamma_x - 0.5 \cdot |U_i|$
  * $L_1 = \gamma_0 + X_i \cdot \gamma_x + 0.5 \cdot |U_i|$
  * $A(0) = \mathbb{I}(L_0 \ge 0)$
  * $A(1) = \mathbb{I}(L_1 \ge 0)$
  * Where $\gamma_0 = -0.01$ and $\gamma_x$ is a vector $(-0.5, -0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5, 0.5)$ for $p=10$ dimensions.
  * $A = A(0)(1 - Z) + A(1) Z$.
* **Potential Outcomes:** $Y(A=a) = \mu_Y(a, X, U) + \epsilon_a$ for $a=0, 1$. Where:
  * $\mu_Y(a, X, U) = a(1 + X \cdot \text{coef} - U) - X \cdot (\text{coef}/2) + 0.5U$ and $\text{coef} = (1, 1, 1, 1, 1, -0.5, -0.5, -0.5, -0.5, -0.5)$.
  * $\epsilon_0, \epsilon_1 \sim N(0, 1)$.
  * And $Y = A Y(A=1) + (1-A) Y(A=0)$.
* **Estimand:** We consider the following three:
  * Denominator: $E[A^{(z=1)} - A^{(z=0)}] = E_X [ E[A|Z=1, X] - E[A|Z=0, X] ]$
  * Numerator: $E[Y^{(z=1)} - Y^{(z=0)}] = E_X [ E[Y|Z=1, X] - E[Y|Z=0, X] ]$
  * LATE: $E[Y(Z=1) - Y(Z=0) \mid A(Z=1) > A(Z=0)]$
* **Kernels:** Gaussian, Laplacian, TCFD, Mat\'ern and Energy.

The bandwidth parameters of Gaussian and Laplacian were chosen from the median heuristics [Garreau, D., Jitkrittum, W., & Kanagawa, M. (2018)](https://arxiv.org/abs/1707.07269)[^3]. In addition, for TCFD, the kernel is approximated by using i.i.d. $V_1, \dots, V_T$ drawn from a $t(5)$ distribution using the following equation with $T=10^4$:

$$
\tilde{k}(X_i, X_j) \simeq \frac{1}{T} \sum_{t=1}^T \cos \{ V_t'(X_i - X_j) \}
$$

* **Regularization parameter:** $\lambda = \frac{1}{N^2}$.
* **Sample size:** $N = 100, 200, 400, 800, 1600$.
* **Repetition:** For each $N$, number of estimands calculated = 500.
## Codes
### Simulation_Final.R
#### 1. Setup and Hyperparameters
* The script is designed to run in a batch environment (e.g., HTCondor or SLURM).

* **Batch Processing:** Reads `commandArgs` to set the `BATCH` index, which selects a row from the hyperparameter grid.
* **Hyperparameter Grid:** A full factorial design covering:
    * **Kernels:** `Gaussian`, `Laplacian`, `Energy`, `Matern`, `TCFD` (t-5), `IPW`, `CBPS`.
    * **DGP Types:** `1` (Linear IV), `2` (Non-linear IV).
    * **Sample Size ($N$):** 100, 200, 400, 800, 1600.
    * **Inference:** `SS` (Subsampling) or `Boot` (Bootstrap).
    * **Seeds:** 1 to 550.

#### 2. Kernel Definitions
Functions to compute the Gram matrix for various kernels.

* **`K.T.CFD(Xmat, bandwidth)`**: Implements the **TCFD Kernel** approximation. It approximates the kernel by projecting data with random weights drawn from a t-distribution ($df=5$) and computing the cosine distance.
* **`Ker(Xmat, name, ...)`**: The main wrapper for kernel calculations.
    * **Logic:** Computes distance matrices (Euclidean or Manhattan) and applies the kernel formula (e.g., Gaussian $\exp{(-D^2/\gamma^2)}$).
    * **Bandwidth Selection:** Uses the **median heuristic** (median of pairwise distances)[^3] if no bandwidth is provided.

#### 3. Data Generation
* **`DGP(DGP.type, N, p)`**: Generates synthetic datasets ($X, Z, A, Y$).
    * **Covariates ($X$):** Standard normal noise ($p=10$).
    * **Instrument ($Z$):** Generated via a logistic model. `DGP.type=1` uses a linear combination of $X$; `DGP.type=2` uses a non-linear transformation.
    * **Treatment ($A$):** Determined by latent variables $L_0, L_1$ to ensure the monotonicity assumption (defining Compliers, Always-takers, Never-takers).
    * **Outcome ($Y$):** A function of treatment and covariates with added noise.
    * **True LATE:** Calculates the theoretical true LATE for benchmarking.

#### 4. Distribution Balancing (Optimization)
* **`distbalance(treatment, covariate, name, bandwidth, lambda)`**: The core optimization routine.
    * **Goal:** Find weights that minimize the Characteristic Function Distance (CFD) between treatment groups while satisfying constraints (weights sum to sample sizes).
    * **Method:** Solves a Quadratic Programming (QP) problem using the `optiSolve` package.
    * **Regularization:** Adds a ridge term $\lambda = 1/N^2$ to the diagonal of the $Q$ matrix for stability.

#### 5. Point Estimation
* **`sim(DGP.type, N, ...)`**: Orchestrates a single simulation run.
    * **Kernel Methods:** Calls `distbalance` to get optimal weights, then calculates the Wald estimator using weighted means.
    * **IPW / CBPS:** Fits a logistic regression or uses `WeightIt::weightit` (CBPS) to estimate propensity scores, then calculates the standard Inverse Probability Weighted estimator.

#### 6. Inference Methods
Two methods are implemented to construct Confidence Intervals (CIs).

* **`choose_best_m_ss(...)` (Subsampling):** Used primarily for kernel methods.
    * **Adaptive Subsampling:** Instead of fixing the subsample size $m$, it searches a grid of 20 values between $\sqrt{N}$ and $3\sqrt{N}$.
    * **Volatility Index (VI):** Selects the optimal $m$ by minimizing the volatility (standard deviation)[^6] of the resulting confidence intervals.
* **`run_boot(...)` (Bootstrap):** Standard non-parametric bootstrap.
    * Resamples the dataset with replacement $N$ times.
    * Re-estimates parameters (including re-fitting propensity scores for IPW/CBPS) to generate empirical distribution CIs.

#### 7. Execution and Output
* **`run_all(...)`**: The main driver function.
    1.  Generates data.
    2.  Computes the point estimate.
    3.  Runs the selected inference method (Subsampling or Bootstrap).
    4.  Returns a consolidated list of results (Point Estimate, Lower CI, Upper CI, Best $m$).
* **Output:** The script writes a `.csv` file named dynamically based on the parameters (e.g., `Result_optiSolve_Gaussian_DGP00001_N00100_SS_SEED00001.csv`).
  
### Simulation_Report.R: generates the summary table using the combined result from the previous one.
This R script processes the raw simulation output (`Result_merge.csv`) to compute performance metrics (Bias, MSE, Coverage, etc.) for the LATE estimators. 

#### 1. True Value Definitions
Hardcoded "True" values derived from the Data Generating Process (DGP) are defined for benchmarking.
* **Linear DGP (Type 1):** `Denom_Linear`, `Numer_Linear`, `LATE_Linear`.
* **Non-Linear DGP (Type 2):** `Denom_nonLinear`, `Numer_nonLinear`, `LATE_nonLinear`.

#### 2. Data Cleaning & Consistency Checks
The script ensures data integrity before analysis:
* **Mismatch Filtering:** Checks that point estimates (`num`, `denom`, `late`) are identical for both inference methods (`SS` and `Boot`) within the same Seed/N/Kernel. If the point estimates differ by more than $10^{-6}$, those seeds are discarded.
* **Deduplication:** Rounds values to 7 decimal places and removes duplicate rows.
* **Trimming:** Retains the first 500 valid repetitions per group.

#### 3. Confidence Interval Rescaling (Subsampling)
Confidence intervals from subsampling must be rescaled to account for the difference between the subsample size ($m$) and the full sample size ($N$).
* **Transformation:**
    $$\text{CI}_ {\text{scaled}} = \text{PointEst} \pm \sqrt{\frac{m}{N}} \times \{\text{Limit}_{\text{sub}} - \text{PointEst}\}$$
* This converts the width of the interval from the subsample scale to the proper $\sqrt{N}$ scale required for inference on the full dataset.

#### 4. Performance Metrics Calculation
The script aggregates results by `kernel`, `N`, and `Inference` type to calculate the following metrics for the **Numerator**, **Denominator**, and **LATE** separately for each DGP type:

* **Bias:** Scaled bias ($100 \times \text{mean}(\hat{\theta} - \theta_{\text{true}})$).
* **SE (Standard Error):** Scaled standard deviation of the estimates ($100 \times \text{SD}(\hat{\theta})$).
* **MSE (Mean Squared Error):** $\text{Bias}^2 + \text{SE}^2$.
* **Coverage:** The proportion of simulations where the True value falls within the calculated Confidence Interval.
* **Length:** The average width of the Confidence Interval.

#### 5. Output Dataframes
Six summary dataframes are created:
* **DGP 1 (Linear):** `num1`, `denom1`, `late1`
* **DGP 2 (Non-Linear):** `num2`, `denom2`, `late2`

# DataAnalysis Folder
The DataAnalysis folder contains replication files for the real-data analysis in Section 7 of the main paper. 
## Description
The data set is taken from the R package _DoubleML_ and can be used to estimate the effect of 401(k) eligibility and participation on accumulated assets [DoubleML2022Python](http://jmlr.org/papers/v23/21-0862.html)[^2].In this data,

* **Number of observations (N)** = 9915
* **The outcome variable (Y)** = Net financial assets (`net_tfa` ∈ [-502302, 1536798]), numeric
* **Treatment (A)** = Participation to 401k (`p401` ∈ {0,1})
* **Instrumental variable (Z)** = Eligibility for 401k plan (`e401` ∈ {0,1}) [3682 are eligible]
* **Estimand:** Our estimand is **the local average treatment effect (LATE)**[^4][^5], the treatment effect among compliers, defined and identified as:

$$
\tau_{\text{LATE}}^{*} = E{\{ Y _i(A_i=1) - Y _ i ( A _ i=0 ) \mid A _ i ( Z _ i=1) > A_ i(Z _i=0) \}} = \frac{E\{(2Z_i-1)Y_i/\Pr(Z_i|X_i)\}}{E\{(2Z_i-1)A_i/\Pr(Z_i|X_i)\}}
$$

* **Covariates (X)**
    * Age (`age`) $\in [25, 64]$, numeric
    * Income (`inc`) $\in [-2652, 242124]$, numeric
    * Family size (`fsize`) $\in [1, 13]$, numeric
    * Years of education (`educ`) $\in [1, 18]$, numeric
    * `db`: 1 if individual has defined benefit pension, binary
    * `marr`: 1 if married, binary
    * `twoearn`: 1 if two-earner household, binary
    * `pira`: 1 if individual participates in IRA plan, binary
    * `hown`: 1 if home owner, binary
 Each of the numeric covariates are standardized so that they fall in the $[0,1]$ interval.
* **Regularization parameter:** $\lambda=\frac{1}{N^2}$.
* **Sub-sample size:** 20 different values from $[\sqrt{N}, 3\sqrt{N}]$.
    * $m \in$ { 100, 110, 121, 131, 142, 152, 162, 173, 183, 194, 204, 215, 225, 236, 246, 257, 267, 278, 288, 299 }.
* **Repetition:** For each $m$, inference method, and kernel/method, number of estimands calculated = 500.
* **Optimal sub-sample size:** Chosen among the 20 sub-sample sizes based on volatility scores[^6].

#### Methods

* **Distribution balancing:** Gaussian, Laplacian, Matern, Energy, and TCFD kernels are used to construct the gram matrix.
    * The bandwidth parameters of Gaussian, Laplacian, and Matern were chosen from the median heuristics[^3].
    * For TCFD, the kernel is approximated by using i.i.d. $V_1, \ldots, V_T$ drawn from a $t(5)/\sigma$ distribution using the following equation with $T=10^4$:
    
    $$
    \tilde{k}(X _i,X _j) \simeq \frac{1}{T}\sum _{t=1}^T \cos \{ V_t'(X_i-X_j) \}
    $$
    
    where $\sigma =$ median[Laplacian kernel].

* **Covariate balancing propensity score (CBPS):** The propensity scores are calculated using the `WeightIt` package, and the weights are calculated as the inverse propensity score. This method was proposed in Imai (2014)[^8]. R package descriptions are given in Greifer (2020)[^9].

* **Inverse probability weighting (IPW):** The propensity scores are calculated using the logistic regression model with no interaction terms, and the weights are calculated as the inverse propensity score.

* **Double debiased machine learning (DDML):** See Chernozhukov (2018)[^10].
## Codes
The functions used in the following R files are similar to the previous section (simulation).
* **weights_Final.R** calculates the optimal weights for different kernels.
* **401k_Final.R** generates a single row as output and produce a csv file with a name structure : "_Result_401k_%s_N%05d_%s_SEED%05d.csv_".
* **401k_Report.R** generates the summary table using the combined result from _401k_Final_ and weights from _weights_Final.R_.

## References
[^1]: Santra, Chen, Park (2026). **Distributional Balancing for Causal Inference: A Unified Framework via Characteristic Function Distance**, arXiv [link](https://doi.org/10.48550/arXiv.2601.15449).
[^2]: Bach, P., Chernozhukov, V., Kurz, M. S., & Spindler, M. (2022). **DoubleML – An object-oriented implementation of double machine learning in Python**. Journal of Machine Learning Research, 23 (53), 1–6. [link](http://jmlr.org/papers/v23/21-0862.html).
[^3]: Garreau, D., Jitkrittum, W., & Kanagawa, M. (2018). **Large sample analysis of the median heuristic**. [link](https://arxiv.org/abs/1707.07269).
[^4]: Angrist, J., & Imbens, G. (1995). **Identification and estimation of local average treatment effects**. _National Bureau of Economic Research Cambridge_, Mass., USA.
[^5]: Angrist, J. D., Imbens, G. W., & Rubin, D. B. (1996). **Identification of causal effects using instrumental variables**. _Journal of the American Statistical Association_, 91 (434), 444–455.
[^6]: Politis, D. N., Romano, J. P., & Wolf, M. (1999). **Subsampling**. _Springer New York_. [link](https://books.google.com/books?id=nGu6rqjE6JoC)
[^8]: Imai, K., & Ratkovic, M. (2014). **Covariate balancing propensity score**. _Journal of the Royal Statistical Society Series B: Statistical Methodology_, 76 (1), 243–263.
[^9]: Greifer, N. (2020). **WeightIt: Weighting for Covariate Balance in Observational Studies**.[link](https://ngreifer.github.io/WeightIt/).
[^10]: Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., & Robins, J. (2018). **Double/debiased machine learning for treatment and structural parameters**. _Oxford University Press_ Oxford, UK.
