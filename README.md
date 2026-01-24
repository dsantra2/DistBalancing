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
  * **Linear:** $Z_i \sim \text{Bernoulli}(\text{expit}(X_i \cdot \mathbf{v}))$, where $\mathbf{v}$ is a vector $0.15\mathbf{1}$ and $\text{expit}(x) = \frac{e^x}{1+e^x}$.
  * **Non-Linear:** $Z_i \sim \text{Bernoulli}\left(\text{expit}\left(\frac{1}{4}|X_i|_{[1:5]} + \frac{1}{4}X_i_{[6:10]} - 2.25\right)\right)$.
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
* **Kernels:** Gaussian, Laplacian, Energy, and TCFD, compared with IPW and CBPS.

The bandwidth parameters of Gaussian and Laplacian were chosen from the median heuristics [@garreau2018largesampleanalysismedian]. In addition, for TCFD, the kernel is approximated by using i.i.d. $V_1, \dots, V_T$ drawn from a $t(5)$ distribution using the following equation with $T=10^4$:

$$
\tilde{k}(X_i, X_j) \simeq \frac{1}{T} \sum_{t=1}^T \cos \{ V_t'(X_i - X_j) \}
$$

* **Regularization parameter:** $\lambda = \frac{1}{N^2}$.
* **Sample size:** $N = 100, 200, 400, 800, 1600$.
* **Repetition:** For each $N$, number of estimands calculated = 500.
### Codes

## DataAnalysis Folder
The DataAnalysis folder contains replication files for the real-data analysis in Section 7 of the main paper. 
### Description
The data set is taken from the R package _DoubleML_ and can be used to estimate the effect of 401(k) eligibility and participation on accumulated assets [DoubleML2022Python](http://jmlr.org/papers/v23/21-0862.html).In this data,

* **Number of observations (N)** = 9915
* **The outcome variable (Y)** = Net financial assets (`net_tfa` ∈ [-502302, 1536798]), numeric
* **Treatment (A)** = Participation to 401k (`p401` ∈ {0,1})
* **Instrumental variable (Z)** = Eligibility for 401k plan (`e401` ∈ {0,1}) [3682 are eligible]
### Codes

## References
* Santra, Chen, Park (2026). **Distributional Balancing for Causal Inference: A Unified Framework via Characteristic Function Distance**, arXiv [link](https://doi.org/10.48550/arXiv.2601.15449).
* Bach, P., Chernozhukov, V., Kurz, M. S., & Spindler, M. (2022). **DoubleML – An object-oriented implementation of double machine learning in Python**. Journal of Machine Learning Research, 23 (53), 1–6. [link](http://jmlr.org/papers/v23/21-0862.html).
* Garreau, D., Jitkrittum, W., & Kanagawa, M. (2018). **Large sample analysis of the median heuristic**. [link](https://arxiv.org/abs/1707.07269).
