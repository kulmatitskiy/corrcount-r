# FUN(value, i)
MapVector <- function(vec, FUN) {
  return(mapply(FUN, vec, 1:length(vec)))
}

# FUN(value, i, j)
MapMatrix <- function(X, FUN) {
  nrows <- length(X[,1])
  ncols <- length(X[1,])
  return(matrix(mapply(FUN, as.vector(X), rep(1:nrows, times=ncols), rep(1:ncols, each=nrows)), ncol=ncols))
}

# FUN(value, i, j, k)
MapCube <- function(X, FUN) {
  nrows <- length(X[ ,1,1])
  ncols <- length(X[1, ,1])
  nmats <- length(X[1,1, ])
  resultVec <- mapply(FUN, as.vector(X), rep(1:nrows, times = nmats * ncols),
                      rep(rep(1:ncols, each=nrows), times=nmats),
                      rep(1:nmats, each = nrows * ncols))
  return(array(resultVec, dim=c(nrows, ncols, nmats)))
}


RoundTo <- function(val, f) { 
  return(f * round(val / f)) 
}


SumSquares <- function(counts, pMeans, pSDs) {
  return((pSDs ^ 2) * (counts - 1) + (counts * pMeans ^ 2) )
}

