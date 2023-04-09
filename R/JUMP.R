JUMP <- function(pvals1, pvals2, alpha = 0.05, lambda = seq(0.01, 0.8, 0.01)){
  m = length(pvals1)

  # Storey's method to estimate the proportions
  xi00.hat = c(); pi0.hat1 = c(); pi0.hat2 = c()
  for (i in 1:length(lambda)) {
    xi00.hat[i] <- sum(pvals1[1:m]>=lambda[i] & pvals2>=lambda[i]) / (m*(1-lambda[i])^2)
    pi0.hat1[i] <- sum(pvals1>=lambda[i]) / (m * (1-lambda[i]))
    pi0.hat2[i] <- sum(pvals2>=lambda[i]) / (m * (1-lambda[i]))
  }
  # fitting a cubic spline by Storey and Tibshirani (2003)
  fit1 <- lm(xi00.hat~., data.frame(cbind(xi00.hat, bs(lambda))))
  fit2 <- lm(pi0.hat1~., data.frame(cbind(pi0.hat1, bs(lambda))))
  fit3 <- lm(pi0.hat2~., data.frame(cbind(pi0.hat2, bs(lambda))))

  pred1 = predict(fit1)
  diff1 = abs(diff(pred1))

  pred2 = predict(fit2)
  diff2 = abs(diff(pred2))

  pred3 = predict(fit3)
  diff3 = abs(diff(pred3))

  xi00.hat = as.numeric(pred1[which.min(diff1)])
  xi01.hat = max(0, as.numeric(pred2[which.min(diff2)] - xi00.hat))
  xi10.hat = max(0, as.numeric(pred3[which.min(diff3)] - xi00.hat))

  xi.hat = c(xi00.hat, xi01.hat, xi10.hat)
  res <- jump_cutoff(pvals1, pvals2, xi.hat, 0.05)

  return(list(p.max = res$p_max, jump.thr = res$thr_jump))
}


