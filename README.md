# Replication Files for "Distributional Balancing for Causal Inference: A Unified Framework via Characteristic Function Distance"
This supplementary file contains replication codes for "Distributional Balancing for Causal Inference: A Unified Framework via Characteristic Function Distance" [Santra, Chen, Park, 2026](https://doi.org/10.48550/arXiv.2601.15449).
## Software and Packages
1. Software: R version 4.4.3
2. Packages: optiSolve (version 1.0); WeightIt (version 1.4.0); DoubleML (version 1.0.2); mlr3 (version 1.0.1); mlr3learners (version 0.12.0); data.table (version 1.17.0); ggplot2 (version 3.5.2); purrr (version 1.0.4); dplyr (version 1.1.4); tidyr (version 1.3.1)
## Simulation Folder
The Simulation folder contains replication files for the simulation studies presented in Section 6 of the main paper.
### Codes

## DataAnalysis Folder
The DataAnalysis folder contains replication files for the real-data analysis in Section 7 of the main paper. The data set is taken from DoubleML R package and can be used to estimate the effect of 401(k) eligibility and participation on accumulated assets [@DoubleML2022Python].In this data,

\begin{itemize}
  \item \textbf{Number of observations (N)} = 9915,
  \item \textbf{The outcome variable (Y)} = Net financial assets ($net\_tfa\in[-502302,1536798]$), numeric
  \item \textbf{Treatment (A)} = Participation to 401k (p401$\in\{0,1\}$),
  \item \textbf{Instrumental variable (Z)} = Eligibility for 401k plan (e401$\in\{0,1\}$) [3682 are eligible]
\end{itemize} 
### Codes
