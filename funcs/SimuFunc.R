######################################################################
#                Real data guided simulation function                #
#            Study 1: Replicate 9 + Study 2: Replicate 12            #
######################################################################
#' @param tau1: a value in (0,1), the variance of non-spatial residual error for study 1
#' @param tau2: a value in (0,1), the variance of non-spatial residual error for study 2
#' @param ipt: the spatial expression pattern on which to generate the data. 'Pattern II' as default
#' @param sig1.str: a string, strength of SE gene signals. 'Strong' as default, also support 'Weak' and 'Moderate'
#' @param sig2.str: a string, strength of SE gene signals. 'Moderate' as default, also support 'Weak' and 'Strong'
#' @param truth1: a vector of true hidden states of study 1, corresponding to the genes
#' @param truth2: a vector of true hidden states of study 2, corresponding to the genes

SimuData <-
  function(tau1 = 0.3,
           tau2 = 0.5,
           ipt = "Pattern II",
           sig1.str = "Strong",
           sig2.str = "Moderate",
           truth1, 
           truth2) {
    require(SPARK)
    source("funcs/spark.test.R")
    n.gene = length(truth1)
    
    load("simulations/MOB_pattern(9+12).RData")
    
    ## Generate the study 1 spatial count data based on the SPARK results from the real
    ## mouse olfactory bulb data
    beta1 <- sapply(1:length(spark1@res_vc), function(x) {
      spark1@res_vc[[x]]$coefficients
    })
    
    # -9.94, -9.25, -8.84, -8.55
    nb1 <- sapply(1:4, function(x) {
      log(x * exp(median(beta1)))
    })  
    info1 <- cbind(spark1@location, total_counts = spark1@lib_size)
    newN1 <- info1$total_counts
    pattern1 <- datalist1[,ipt]
    grp1 <- as.numeric(pattern1 > mean(pattern1)) + 1
    
    if (sig1.str == "Weak")
      ifc1 <- 2
    if (sig1.str == "Moderate")
      ifc1 <- 3
    if (sig1.str == "Strong")
      ifc1 <- 4
    uu1  <- c(nb1[1], nb1[ifc1], nb1[1])
    
    lambda1 <- sapply(1:n.gene, function(x) {
      truth1[x] * exp(uu1[grp1] + rnorm(length(uu1[grp1]), 0, tau1)) + (!truth1[x]) *
        exp(uu1[3] + rnorm(length(uu1[grp1]), 0, tau1))
    })
    newCt1 <- lapply(1:(n.gene), function(x) {
      rpois(length(lambda1[, x]), newN1 * lambda1[, x])
    })
    countdata1 <- data.frame(do.call(rbind, newCt1))
    rownames(countdata1) <- paste0("gene", 1:nrow(countdata1))
    colnames(countdata1) <- rownames(info1)
    
    ## Compute primary p-values with SPARK
    spark1 <-
      CreateSPARKObject(
        counts = countdata1,
        location = info1[,1:2],
        percentage = 0.1,
        min_total_counts = 10
      )
    
    spark1@lib_size <- info1$total_counts
    spark1 <-
      spark.vc(
        spark1,
        covariates = NULL,
        lib_size = spark1@lib_size,
        num_core = 5,
        verbose = T,
        fit.maxiter = 500
      )
    spark1 <- spark.test(spark1, check_positive = T, verbose = T)
    pvalue1 <- spark1@res_mtest[, "combined_pvalue"]
    closeAllConnections()

    ## Generate the study 2 spatial count data based on the SPARK results from the real
    ## mouse olfactory bulb data
    beta2 <- sapply(1:length(spark2@res_vc), function(x) {
      spark2@res_vc[[x]]$coefficients
    })
    
    # -9.93, -9.23, -8.83, -8.54
    nb2 <- sapply(1:4, function(x) {
      log(x * exp(median(beta2)))
    })

    info2 <- cbind(spark2@location, total_counts = spark2@lib_size)
    newN2 <- info2$total_counts
    pattern2 <- datalist2[,ipt]
    grp2 <- as.numeric(pattern2 > mean(pattern2)) + 1
    
    if (sig2.str == "Weak")
      ifc2 <- 2
    if (sig2.str == "Moderate")
      ifc2 <- 3
    if (sig2.str == "Strong")
      ifc2 <- 4
    uu2  <- c(nb2[1], nb2[ifc2], nb2[1])
    
    lambda2 <- sapply(1:n.gene, function(x) {
      truth2[x] * exp(uu2[grp2] + rnorm(length(uu2[grp2]), 0, tau2)) + (!truth2[x]) *
        exp(uu2[3] + rnorm(length(uu2[grp2]), 0, tau2))
    })
    newCt2 <- lapply(1:(n.gene), function(x) {
      rpois(length(lambda2[, x]), newN2 * lambda2[, x])
    })
    countdata2 <- data.frame(do.call(rbind, newCt2))
    rownames(countdata2) <- paste0("gene", 1:nrow(countdata2))
    colnames(countdata2) <- rownames(info2)
    
    ## Compute primary p-values with SPARK
    spark2 <-
      CreateSPARKObject(
        counts = countdata2,
        location = info2[, 1:2],
        percentage = 0.1,
        min_total_counts = 10
      )
    
    spark2@lib_size <- info2$total_counts
    spark2 <-
      spark.vc(
        spark2,
        covariates = NULL,
        lib_size = spark2@lib_size,
        num_core = 5,
        verbose = T,
        fit.maxiter = 500
      )
    spark2 <- spark.test(spark2, check_positive = T, verbose = T)
    pvalue2 <- spark2@res_mtest[, "combined_pvalue"]
    closeAllConnections()
    
    return(list(
      pvals1 = pvalue1,
      pvals2 = pvalue2,
      truth1 = truth1,
      truth2 = truth2))
  }