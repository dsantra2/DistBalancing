# Replication Files for "Distributional Balancing for Causal Inference: A Unified Framework via Characteristic Function Distance"
This supplementary file contains replication codes for "Distributional Balancing for Causal Inference: A Unified Framework via Characteristic Function Distance" [Santra, Chen, Park, 2026](https://doi.org/10.48550/arXiv.2601.15449).
## Software and Packages
* Software: R version 4.4.3
* Packages: _optiSolve_ (version 1.0); _WeightIt_ (version 1.4.0); _DoubleML_ (version 1.0.2); _mlr3_ (version 1.0.1); _mlr3learners_ (version 0.12.0); _data.table_ (version 1.17.0); _ggplot2_ (version 3.5.2); _purrr_ (version 1.0.4); _dplyr_ (version 1.1.4); _tidyr_ (version 1.3.1)
## Simulation Folder
The Simulation folder contains replication files for the simulation studies presented in Section 6 of the main paper.
### Codes

## DataAnalysis Folder
The DataAnalysis folder contains replication files for the real-data analysis in Section 7 of the main paper. The data set is taken from the R package _DoubleML_ and can be used to estimate the effect of 401(k) eligibility and participation on accumulated assets [DoubleML2022Python](http://jmlr.org/papers/v23/21-0862.html).In this data,

* **Number of observations (N)** = 9915
* **The outcome variable (Y)** = Net financial assets (`net_tfa` ∈ [-502302, 1536798]), numeric
* **Treatment (A)** = Participation to 401k (`p401` ∈ {0,1})
* **Instrumental variable (Z)** = Eligibility for 401k plan (`e401` ∈ {0,1}) [3682 are eligible]
### Codes
