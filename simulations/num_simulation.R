source("./funcs/methods.R")
library(JUMPrcpp)
library(splines)

## Data simulation
n = 100 # number of replications
m = 10000
xi00 = 0.8
xi01 = 0.075
xi10 = 0.075
xi11 = 0.05
mu1 = 2.5
mu2 = 2.5
sigma1 = 0.8
sigma2 = 0.8
alphas <- 0.05

methods <- c("JUMP", "MaxP", "BH", "radjust", "IDR", "Sidak", "Lancaster", "MaRR")
res <- list()
for (i in 1: length(methods)){
  res[[methods[i]]]$fdp <- matrix(NA, n, length(alphas))
  res[[methods[i]]]$pd <- matrix(NA, n, length(alphas))
}

for(i in 1:n){
  h = sample(0:3, m, replace = TRUE, prob = c(xi00, xi01, xi10, xi11))
  states1 = rep(0, m)
  states1[which(h==2|h==3)] = 1
  states2 = rep(0, m)
  states2[which(h==1|h==3)] = 1

  stat1 = rnorm(m, states1*mu1, sigma1)
  stat2 = rnorm(m, states2*mu2, sigma2)

  p1 = 1 - pnorm(stat1, mean = 0, sd = sigma1)
  p2 = 1 - pnorm(stat2, mean = 0, sd = sigma2)

  truth <- states1 * states2
  
  # BH
  pvals1.bh <- p.adjust(p1, method = "BH")
  pvals2.bh <- p.adjust(p2, method = "BH")
  
  # MaxP
  maxp <- pmax(p1, p2)
  pvals.maxp <- p.adjust(maxp, method = "BH")
  
  # IDR
  library(idr)
  x <- cbind(-log(p1), -log(p2))
  idr.vals <- QL.randomstarts(x)
  
  # Sidak
  pvals.Sidak <- sidak(p1, p2)
  pSidak.bh <- p.adjust(pvals.Sidak, method = "BH")

  # Lancaster's method
  p.Lancaster <- lancaster(p1, p2)
  pLancaster.bh <- p.adjust(p.Lancaster, method = "BH")
  
  # MaRR
  max.rank = apply(cbind(rank(p1), rank(p2)),1,max)
  
  for (j in 1:length(alphas)){
    alpha = alphas[j]
    
    # JUMP
    jump.obj <- JUMP(p1, p2, alpha)
    jump.thr <- jump.obj$jump.thr
    p.max <- jump.obj$p.max
    res$JUMP$fdp[i,j] <- sum(p.max <= jump.thr & !truth)/max(sum(p.max <= jump.thr), 1)
    res$JUMP$pd[i,j]  <- sum(p.max <= jump.thr & truth) / sum(truth)
    
    # BH
    res$BH$fdp[i,j] <- sum(pvals1.bh <= alpha & pvals2.bh <= alpha & !truth)/max(sum(pvals1.bh <= alpha & pvals2.bh <= alpha), 1)
    res$BH$pd[i,j] <- sum(pvals1.bh <= alpha & pvals2.bh <= alpha & truth) / sum(truth)
    
    # MaxP
    res$MaxP$fdp[i,j] <- sum(pvals.maxp <= alpha & !truth)/ max(sum(pvals.maxp <= alpha), 1)
    res$MaxP$pd[i,j]  <- sum(pvals.maxp <= alpha & truth) / sum(truth)
    
    # Bogomolov & Heller 2018 (JASA)
    library(radjust)
    pa = p1
    pb = p2
    pa[which(pa==0)] <- 1e-15
    pb[which(pb==0)] <- 1e-15
    res.rv18 <- radjust_sym(pa, pb, input_type = "all", directional_rep_claim = FALSE,
                            variant = "adaptive", alpha = alpha)
    rv18 <- rep(1, m)
    rv18[as.numeric(res.rv18$results_table$name)] <- res.rv18$results_table$r_value
    res$radjust$fdp[i] <- sum(rv18 <= alpha & !truth)/max(sum(rv18 <= alpha), 1)
    res$radjust$pd[i] <- sum(rv18 <= alpha & truth) / sum(truth)
    
    # IDR
    res$IDR$fdp[i,j] <- sum(idr.vals <= alpha & !truth)/max(sum(idr.vals <= alpha), 1)
    res$IDR$pd[i,j] <- sum(idr.vals <= alpha & truth) / sum(truth)
    
    # Sidak's method
    res$Sidak$fdp[i,j] <- sum(pSidak.bh <= alpha & !truth)/ max(sum(pSidak.bh <= alpha), 1)
    res$Sidak$pd[i,j]  <- sum(pSidak.bh <= alpha & truth) / sum(truth)

    # Lancaster's method
    res$Lancaster$fdp[i,j] <- sum(pLancaster.bh <= alpha & !truth)/max(sum(pLancaster.bh <= alpha), 1)
    res$Lancaster$pd[i,j]  <- sum(pLancaster.bh <= alpha & truth)/sum(truth)
    
    # MaRR
    res.marr <- MaRR(max.rank, alpha = alpha)
    marr.rej <- rep(0, m)
    marr.rej[res.marr$which.sig] = 1
    res$MaRR$fdp[i,j] <- sum(marr.rej & !truth)/max(sum(marr.rej), 1)
    res$MaRR$pd[i,j] <- sum(marr.rej & truth)/sum(truth)
  }
  print(i)
}

for (i in 1: length(methods)){
  res[[methods[i]]]$fdr <- colMeans(res[[methods[i]]]$fdp, na.rm = TRUE)
  res[[methods[i]]]$power <- colMeans(res[[methods[i]]]$pd, na.rm = TRUE)
}
