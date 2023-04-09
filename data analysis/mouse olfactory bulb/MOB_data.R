#################################################################################
#              Mouse olfactory bulb data analysis (ST + 10X Visium)             #
#################################################################################
library(SPARK)

##------------------------------------------------------------------------
## Load the ST Rep9 MOB data and apply SPARK 
##------------------------------------------------------------------------
# read the raw counts (spatial data)
counts1 <- read.table("./data analysis/mouse olfactory bulb/ST/Rep9_MOB_count_matrix-1.tsv", check.names = F)
counts1 <- t(counts1)
location1 <- cbind.data.frame(x = as.numeric(sapply(strsplit(colnames(counts1), split = "x"), "[", 1)), 
                              y = as.numeric(sapply(strsplit(colnames(counts1), split = "x"), "[", 2)))
rownames(location1) <- colnames(counts1)

# SPARK analysis
spark1 <- CreateSPARKObject(counts = counts1, location = location1, 
                            percentage = 0.1, min_total_counts = 10)
spark1@lib_size <- apply(spark1@counts, 2, sum)
spark1 <- spark.vc(spark1, covariates = NULL, lib_size = spark1@lib_size, 
                   num_core = 10, verbose = T, fit.maxiter = 500)
spark1 <- spark.test(spark1, check_positive = T, verbose = T)

##------------------------------------------------------------------------
## Load the 10X Visium MOB data and apply SPARK 
##------------------------------------------------------------------------
counts2 <- Matrix::readMM("./data analysis/mouse olfactory bulb/10X Visium/filtered_feature_bc_matrix/matrix.mtx")
genes2 <- read.delim("./data analysis/mouse olfactory bulb/10X Visium/filtered_feature_bc_matrix/features.tsv", header=FALSE)
barcodes2 <- read.table("./data analysis/mouse olfactory bulb/10X Visium/filtered_feature_bc_matrix/barcodes.tsv", quote="\"", comment.char="")
location2 <- read.csv("./data analysis/mouse olfactory bulb/10X Visium/spatial/tissue_positions.csv", header=FALSE)

barcodes2 <- as.vector(t(barcodes2))
rownames(counts2) <- genes2$V2
rownames(location2) <- location2$V1
location2 <- location2[barcodes2,5:6]
colnames(location2) <- c("x", "y")
location2$x <- as.numeric(location2$x)
location2$y <- as.numeric(location2$y)
colnames(counts2) <- paste(location2$x, location2$y, sep = "x") # 32285 genes * 1185 spots
rownames(location2) <- colnames(counts2)

mt_idx2  <- grep("mt-",rownames(counts2)) 
counts2 <- counts2[-mt_idx2,] # 32272 genes * 1185 spots
rm(mt_idx2, genes2, barcodes2)

# SPARK analysis
source('./funcs/spark.test.R')
spark2 <- CreateSPARKObject(counts = counts2, location = location2, 
                            percentage = 0.1, min_total_counts = 10)
spark2@lib_size <- apply(spark2@counts, 2, sum)
spark2 <- spark.vc(spark2, covariates = NULL, lib_size = spark2@lib_size, 
                   num_core = 10, verbose = T, fit.maxiter = 500)
spark2 <- spark.test(spark2, check_positive = T, verbose = T)

##----------------------------------------------------------------------------
save(counts1, counts2, location1, location2, spark1, spark2, file = "data analysis/mouse olfactory bulb/MOB.RData")
