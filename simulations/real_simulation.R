## Generate synthetic data guided by two MOB datasets
## Compute FDR and Power of all methods across 10 replicates

library(SPARK)
library(JUMP)
source('./funcs/SimuFunc.R')

n.rep    = 10
tau1     = 0.3
tau2     = 0.5
ipt      = "Pattern II"
n.gene   = 10000
sig1.str = "Strong"
sig2.str = "Moderate"
xi00 = 0.7
xi01 = 0.125
xi10 = 0.125
xi11 = 0.05
alphas = 0.05

methods <- c("MaxP", "BH", "JUMP")

data.obj <- list()
res <- list()

for (i in 1: length(methods)){
  res[[methods[i]]]$fdp <- matrix(NA, n.rep, length(alphas))
  res[[methods[i]]]$pd <- matrix(NA, n.rep, length(alphas))
}

for (i in 1: n.rep){
  cat("Replicate ", i, ":\n")
  h = sample(0:3, n.gene, replace = TRUE, prob = c(xi00, xi01, xi10, xi11))
  states1 = rep(0, n.gene)
  states1[which(h==2|h==3)] = 1
  states2 = rep(0, n.gene)
  states2[which(h==1|h==3)] = 1
  data.obj[[i]] <- SimuData(tau1 = tau1,
                            tau2 = tau2,
                            ipt = ipt,
                            sig1.str = sig1.str,
                            sig2.str = sig2.str,
                            truth1   = states1, 
                            truth2   = states2)
  
  pvals1 <- data.obj[[i]]$pvals1
  pvals2 <- data.obj[[i]]$pvals2
  states1 <- data.obj[[i]]$truth1
  states2 <- data.obj[[i]]$truth2
  
  truth = states1 * states2
  
  # BH
  pvals1.bh <- p.adjust(pvals1, method = "BH")
  pvals2.bh <- p.adjust(pvals2, method = "BH")
  
  # MaxP
  maxp <- pmax(pvals1, pvals2)
  pvals.maxp <- p.adjust(maxp, method = "BH")
  
  for(j in 1:length(alphas)){
    alpha = alphas[j]
    
    # BH
    res$BH$fdp[i,j] <- sum(pvals1.bh <= alpha & pvals2.bh <= alpha & !truth)/max(sum(pvals1.bh <= alpha & pvals2.bh <= alpha), 1)
    res$BH$pd[i,j] <- sum(pvals1.bh <= alpha & pvals2.bh <= alpha & truth) / sum(truth)
    
    # MaxP
    res$MaxP$fdp[i,j] <- sum(pvals.maxp <= alpha & !truth)/ max(sum(pvals.maxp <= alpha), 1)
    res$MaxP$pd[i,j]  <- sum(pvals.maxp <= alpha & truth) / sum(truth)
    
    # JUMP
    jump.obj <- JUMP(pvals1, pvals2, alpha)
    jump.thr <- jump.obj$jump.thr
    p.max <- jump.obj$p.max
    res$JUMP$fdp[i,j] <- sum(p.max <= jump.thr & !truth)/max(sum(p.max <= jump.thr), 1)
    res$JUMP$pd[i,j] <- sum(p.max <= jump.thr & truth) / sum(truth)
  }
}

for (k in 1:length(methods)){
  res[[k]]$fdr <- colMeans(res[[k]]$fdp)
  res[[k]]$pwr <- colMeans(res[[k]]$pd)
}