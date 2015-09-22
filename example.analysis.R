source("algorithm.R")

data = read.csv("example_datasets/ex1.csv")
result = RunPoissonEM(counts = data, numSegs=15,
                       pMeans = 0, pSDs = 0,
                       withPayments = FALSE, withInflatedZeros = TRUE,
                       log = TRUE, useFreqs = TRUE, timeLog = TRUE, numAttempts = 10)