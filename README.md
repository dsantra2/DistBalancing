# Replication Files for "Distributional Balancing for Causal Inference: A Unified Framework via Characteristic Function Distance"
This supplementary file contains replication codes for "Distributional Balancing for Causal Inference: A Unified Framework via Characteristic Function Distance" [Santra, Chen, Park, 2026](https://doi.org/10.48550/arXiv.2601.15449).
## Software and Packages
* Software: R version 4.4.3
* Packages: _optiSolve_ (version 1.0); _WeightIt_ (version 1.4.0); _DoubleML_ (version 1.0.2); _mlr3_ (version 1.0.1); _mlr3learners_ (version 0.12.0); _data.table_ (version 1.17.0); _ggplot2_ (version 3.5.2); _purrr_ (version 1.0.4); _dplyr_ (version 1.1.4); _tidyr_ (version 1.3.1)
## Simulation Folder
The Simulation folder contains replication files for the simulation studies presented in Section 6 of the main paper.
### Simulation Setup
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

The bandwidth parameters of Gaussian and Laplacian were chosen from the median heuristics [Garreau, D., Jitkrittum, W., & Kanagawa, M. (2018)](https://arxiv.org/abs/1707.07269). In addition, for TCFD, the kernel is approximated by using i.i.d. $V_1, \dots, V_T$ drawn from a $t(5)$ distribution using the following equation with $T=10^4$:

$$
\tilde{k}(X_i, X_j) \simeq \frac{1}{T} \sum_{t=1}^T \cos \{ V_t'(X_i - X_j) \}
$$

* **Regularization parameter:** $\lambda = \frac{1}{N^2}$.
* **Sample size:** $N = 100, 200, 400, 800, 1600$.
* **Repetition:** For each $N$, number of estimands calculated = 500.
### Codes
* **Simulation_Final.R** generates a single row as output and produce a csv file with a name structure : "_Result_optiSolve_%s_DGP%0.5d_N%0.5d_%s_SEED%0.5d.csv_".
* **Simulation_Report.R** generates the summary table using the combined result from the previous one.

## DataAnalysis Folder
The DataAnalysis folder contains replication files for the real-data analysis in Section 7 of the main paper. 
### Description
The data set is taken from the R package _DoubleML_ and can be used to estimate the effect of 401(k) eligibility and participation on accumulated assets [DoubleML2022Python](http://jmlr.org/papers/v23/21-0862.html).In this data,

* **Number of observations (N)** = 9915
* **The outcome variable (Y)** = Net financial assets (`net_tfa` ∈ [-502302, 1536798]), numeric
* **Treatment (A)** = Participation to 401k (`p401` ∈ {0,1})
* **Instrumental variable (Z)** = Eligibility for 401k plan (`e401` ∈ {0,1}) [3682 are eligible]
* **Estimand:** Our estimand is **the local average treatment effect (LATE)**[^1][^2], the treatment effect among compliers, defined and identified as:

$$
\tau_{\text{LATE}}^{*} = E\{Y_i(A_i=1) - Y_i(A_i=0) \mid A_i(Z_i=1)> A_i(Z_i=0)\} = \frac{E\{(2Z_i-1)Y_i/\Pr(Z_i|X_i)\}}{E\{(2Z_i-1)A_i/\Pr(Z_i|X_i)\}}
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
* **Optimal sub-sample size:** Chosen among the 20 sub-sample sizes based on volatility scores[^3].

#### Methods

* **Distribution balancing:** Gaussian, Laplacian, Matern, Energy, and TCFD kernels are used to construct the gram matrix.
    * The bandwidth parameters of Gaussian, Laplacian, and Matern were chosen from the median heuristics[^4].
    * For TCFD, the kernel is approximated by using i.i.d. $V_1, \ldots, V_T$ drawn from a $t(5)/\sigma$ distribution using the following equation with $T=10^4$:
    
    $$
    \tilde{k}(X_i,X_j) \simeq \frac{1}{T}\sum_{t=1}^T \cos \{ V_t'(X_i-X_j) \}
    $$
    
    where $\sigma =$ median[Laplacian kernel].

* **Covariate balancing propensity score (CBPS):** The propensity scores are calculated using the `WeightIt` package, and the weights are calculated as the inverse propensity score. This method was proposed in Imai (2014)[^5]. R package descriptions are given in Greifer (2020)[^6].

* **Inverse probability weighting (IPW):** The propensity scores are calculated using the logistic regression model with no interaction terms, and the weights are calculated as the inverse propensity score.

* **Double debiased machine learning (DDML):** See Chernozhukov (2018)[^7].
### Codes
* **weights_Final.R** calculates the optimal weights for different kernels.
* **401k_Final.R** generates a single row as output and produce a csv file with a name structure : "_Result_401k_%s_N%05d_%s_SEED%05d.csv_".
* **401k_Report.R** generates the summary table using the combined result from _401k_Final_ and weights from _weights_Final.R_.

## References
* Santra, Chen, Park (2026). **Distributional Balancing for Causal Inference: A Unified Framework via Characteristic Function Distance**, arXiv [link](https://doi.org/10.48550/arXiv.2601.15449).
* Bach, P., Chernozhukov, V., Kurz, M. S., & Spindler, M. (2022). **DoubleML – An object-oriented implementation of double machine learning in Python**. Journal of Machine Learning Research, 23 (53), 1–6. [link](http://jmlr.org/papers/v23/21-0862.html).
* Garreau, D., Jitkrittum, W., & Kanagawa, M. (2018). **Large sample analysis of the median heuristic**. [link](https://arxiv.org/abs/1707.07269).
[^1]: Angrist, J. D., & Imbens, G. W. (1995). *Identification and Estimation of Local Average Treatment Effects*.
[^2]: Angrist, J. D., Imbens, G. W., & Rubin, D. B. (1996). *Identification of Causal Effects Using Instrumental Variables*.
[^3]: Politis, D. N., Romano, J. P., & Wolf, M. (1999). *Subsampling*.
[^4]: Garreau, D., Jitkrittum, W., & Kanagawa, M. (2018). *Large Sample Analysis of the Median Heuristic*.
[^5]: Imai, K., & Ratkovic, M. (2014). *Covariate Balancing Propensity Score*.
[^6]: Greifer, N. (2020). *WeightIt: Weighting for Covariate Balance in Observational Studies*.
[^7]: Chernozhukov, V., et al. (2018). *Double/Debiased Machine Learning for Treatment and Structural Parameters*.
