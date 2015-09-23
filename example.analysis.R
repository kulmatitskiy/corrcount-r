source("algorithm.R")

data <- read.csv("example_datasets/ex1.csv")
numSegs <- 8
result_uncorr <- RunPoissonEM(counts = data, numSegs = numSegs, withRandomEffects = FALSE,
                       withInflatedZeros = TRUE, log = TRUE,
                       timeLog = TRUE, numAttempts = 10)
#result_corr = RunPoissonEM(counts = data[1:100,], numSegs = numSegs,  withRandomEffects = TRUE,
#                           withInflatedZeros = TRUE, log = TRUE,
#                           timeLog = TRUE, numAttempts = 10, MC=500)