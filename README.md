# JUMPrcpp

Joint replicability analysis of high-throughput experiments using maximum of $p$-values across two studies.

# Installation

```R
require("devtools")
devtools::install_github("YanLi15/JUMPrcpp")
```

# A numeric example

We illustrate the usage of JUMPrcpp for the replicability of two large-scale multiple testing hidden joint problems using simulated data.

```R
## Pre-specify the number of hypotheses, m, the prior probabilities of the hidden joint states, xi's, and the alternative settings
m = 10000; xi00 = 0.9; xi01 = 0.025; xi10 = 0.025; xi11 = 0.05
mu1 = 3; mu2 = 3; sigma1 = 1; sigma2 = 1

## Generate the hidden states and corresponding p-values in two studies 
data.obj <- data_generation(m, xi00, xi01, xi10, xi11, mu1, mu2, sigma1, sigma2)
p1 = data.obj$pvals1
p2 = data.obj$pvals2
states1 = data.obj$states1
states2 = data.obj$states2

## Replicability analysis
library(JUMPrcpp)
alpha <- 0.05
jump.obj <- JUMP(pvals1, pvals2, alpha = alpha)
rep.idx <- which(jump.obj$p.max <= jump.obj$jump.thr)
```

# Data and reproducibility

R functions supporting the simulations are contained in the “funcs” folder.

R scripts to reproduce the numeric and realistic simulations can be found in the “simulations” folder.

R scripts to reproduce the data analysis results are summarized in the “data analysis” folder.

The spatially resolved transcriptomic data used in the data analysis can be downloaded from the links provided in the *Data Availability Statement* in the manuscript.

Raw count data and some intermediate results supporting the reproducibility of results can be downloaded [here](https://drive.google.com/drive/folders/1nEMBS7Nwqn6JXyRsiBSMrrqyFeB3rSD_?usp=share_link).
