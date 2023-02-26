# JUMPrcpp

Replicability analysis of high-throughput experiments, which implements JUMP with Rcpp.

# Installation

```R
require("devtools")
devtools::install_github("hongyuan-cao/JUMPrcpp")
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
