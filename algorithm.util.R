source("visualize.R")

# Preparing Input Object
MakeInputObject <- function(inputs) {
  x <- list(
    original.counts = as.matrix(inputs$counts),
    numSegs = inputs$numSegs,
    useFreqs = inputs$useFreqs,
    negativeBinomial1 = inputs$negativeBinomial1,
    withInflatedZeros = inputs$withInflatedZeros,
    hurdle = inputs$hurdle,
    withRandomEffects = inputs$withRandomEffects,
    withPayments = inputs$withPayments,
    logNormalPayments = inputs$logNormalPayments,
    pMeans = inputs$pMeans,
    pSDs   = inputs$pSDs )
  x$numCats <- ncol(inputs$counts)
  x$names   <- colnames(inputs$counts)
  x$Nall    <- nrow(inputs$counts)
  if (!inputs$useFreqs || inputs$withPayments) {
    x$useFreqs <- FALSE
    x$counts <- as.matrix(inputs$counts)
    x$freqs <- rep(1,x$Nall)
    x$Ncmb <- x$Nall
  } else {
    tmp <- aggregate( rep(1,x$Nall), by = as.list(data.frame(inputs$counts)), FUN = sum) 
    x$counts <- as.matrix(tmp[, 1:x$numCats])
    x$freqs <- tmp[, x$numCats + 1]
    x$Ncmb <- nrow(tmp)
  }
  x$binom <- (x$counts > 0)
  x$rfreqs <- x$freqs / x$Nall
  if (inputs$withRandomEffects) {
    x$groupingMatrix <- inputs$groupingMatrix
    x$numGroups <- length(inputs$groupingMatrix[, 1])
  }
  if (inputs$withPayments)
  {
    x$pMeans <- data.matrix(inputs$pMeans)
    x$pSDs   <- data.matrix(inputs$pSDs)
    x$pMask  <- !is.na(x$pMeans) & !is.na(x$pSDs)
    x$pSS    <- sumSquares(x$counts, x$pMeans, x$pSDs)
  }
  class(x) <- "PoissonEMInput"
  return(x)
}



# creates a list to store and update various information about the algorithm's runtime
MakeProcessStatisticsObject <- function (num.attempts, threshold, T.max, log, log.time, MC) {
  proc <- list( num.attempts = num.attempts, threshold = threshold,
               T.max = T.max, log = log, log.time = log.time, MC = MC,
               total.iterations = 0
  )
  if (log.time) {
    time <- Sys.time()
    proc$time <- time
    zero = time - time
    proc$time.stats <- data.frame(init = zero, posteriors = zero, 
                                  param.updates = zero, LL = zero)
  }
  return(proc)
}



LogMessage <- function (process.stats, fstring, ...) {
  if (process.stats$log) {
    cat(sprintf( fstring, ... ))
  }
}


LogTime = function(process.stats, key) {
  if (process.stats$log.time) {
    now <- Sys.time()
    process.stats$time.stats[,key] <- process.stats$time.stats[, key] + (now - process.stats$time)
    process.stats$time <- now
  }
  return(process.stats)
}



# Prepare Output Object
FillReturnObject = function (r) {
  lambdas  <- r$params$lambdas
  segProbs <- r$params$segProbs
  counts   <- r$input$counts
  freqs    <- r$input$freqs
  Ncmb     <- r$input$Ncmb
  Nall     <- r$input$Nall
  K        <- r$input$numCats
  numSegs  <- r$input$numSegs
  
  r.new       <- list()
  r.new$Data  <- list()
  r.new$Model <- list()
  r.new$Model$NumberParams <- r$input$numSegs * (K + 1)
  if (r$input$withInflatedZeros) {
    r.new$Model$NumberParams <- r.new$Model$NumberParams + r$input$numSegs * K
  }
  if (r$input$withRandomEffects) {
    r.new$Model$NumberParams <- r.new$Model$NumberParams + 2 * r$input$numSegs*r$input$numGroups
  }
  r.new$Model$Segments <- list()
  r.new$Model$Segments$Number <- r$input$numSegs
  r.new$Model$Segments$Counts <- list()
  r.new$Proc <- list()
  r.new$Proc$LL         <- r$params$LL
  r.new$Proc$Starts     <- r$procStats$numAttempts
  r.new$Proc$Iterations <- r$procStats$total.iterations
  r.new$Proc$Time       <- r$procStats$timeStats 
  
  r.new$Data$Counts <- r$input$original.counts
  
  r.new$Model$Segments$Counts$Lambdas <- r$params$lambdas # MEANS OR LAMBDAS ?
  r.new$Model$Segments$Counts$Means <- r$params$lambdas # MEANS OR LAMBDAS ?
  if (r$input$withInflatedZeros) {
    r.new$Model$Segments$Counts$BuyProbs <- r$params$buyProbs
    r.new$Model$Segments$Counts$Means <- r.new$Model$Segments$Counts$Means *  r$params$buyProbs
  }
  if (r$input$withRandomEffects) {
    r.new$Model$Segments$Counts$LogGroupEffects <- list()
    r.new$Model$Segments$Counts$LogGroupEffects$Means <- r$params$rEffMus
    r.new$Model$Segments$Counts$LogGroupEffects$Sds <- r$params$rEffSigmas
    r.new$Model$Segments$Counts$GroupEffects <- list()
    r.new$Model$Segments$Counts$GroupEffects$Means <- exp(r$params$rEffMus + 0.5 * r$params$rEffSigmas^2)
    r.new$Model$Segments$Counts$GroupEffects$Sds <- sqrt((exp(r$params$rEffSigmas^2) - 1) * exp(2 * r$params$rEffMus + r$params$rEffSigmas^2))
    # TODO: ADD UPDATE OF COUNTS MEANS
    mult <- exp( log(r.new$Model$Segments$Counts$GroupEffects$Means) %*% t(r$input$groupingMatrix) )
    r.new$Model$Segments$Counts$Means <- r.new$Model$Segments$Counts$Means * mult
  }
  
  r.new$Model$Segments$Probs <- r$params$segProbs
  r.new$Model$Segments$Sizes <- r$params$segProbs * Nall
  
  # Analyze data-based distribution of counts in categories
  realMeans <- colMeans(r$input$original.counts)
  realCentered <- counts - matrix(rep(realMeans,Ncmb), ncol = K, nrow = Ncmb, byrow = TRUE)
  realSDs    <- sqrt(colSums((realCentered^2) * freqs) / (Nall - 1))
  realNrmzd <- realCentered / matrix(rep(realSDs, Ncmb), ncol = K, nrow = Ncmb, byrow = TRUE)
  F         <- realNrmzd * sqrt(freqs) 
  realCorr  <- t(F) %*% F/ (Nall - 1)
  
  r.new$Data$Categories <- list()
  r.new$Data$Categories$Names  <- r$input$names
  r.new$Data$Categories$Counts <- list()
  r.new$Data$Categories$Counts$Means  <- realMeans
  r.new$Data$Categories$Counts$SDs    <- realSDs
  r.new$Data$Categories$Counts$Corr   <- realCorr
  
  # Analyze model-based distribution of counts in categories
  if (r$input$withInflatedZeros) {
    terms <- lambdas * segProbs * r$params$buyProbs
    terms[is.na(terms)] <- 0
    modelMeans <- colSums(terms)
    
    terms <- terms * (lambdas + 1)
    terms[is.na(terms)] <- 0
    modelSDs <- sqrt(colSums(terms) - modelMeans^2)
    
    H <- 0
  } else {
    modelMeans <- colSums(lambdas * segProbs)
    modelSDs   <- sqrt(colSums(lambdas*(lambdas+1)*segProbs) - modelMeans^2)
    H <- ( t(lambdas) %*% (lambdas * segProbs) - modelMeans %*% t(modelMeans) ) / (modelSDs %*% t(modelSDs))
  }
  
  r.new$Model$Categories <- list()
  r.new$Model$Categories$Names  <- r$input$names
  r.new$Model$Categories$Counts <- list()
  r.new$Model$Categories$Counts$Means <- modelMeans
  r.new$Model$Categories$Counts$SDs   <- modelSDs
  r.new$Model$Categories$Counts$Corr  <- H
  
  if(r$input$withPayments) {
    r.new$Model$NumberParams <- r.new$Model$NumberParams + 2 * K * r$inputs$numSegs
    mus    <- r$params$mus
    sigmas <- r$params$sigmas 
    r.new$Model$Segments$Payments <- list()
    if (r$input$logNormalPayments) {
      r.new$Data$LogPayments <- list()
      r.new$Data$LogPayments$Means  <- r$input$pMeans
      r.new$Data$LogPayments$SDs    <- r$input$pSDs
      r.new$Model$Segments$LogPayments <- list()
      r.new$Model$Segments$LogPayments$Means <- mus
      r.new$Model$Segments$LogPayments$SDs   <- sigmas
      r.new$Model$Segments$Payments$Means    <- exp(mus) + 0.5 * (sigmas^2)
      r.new$Model$Segments$Payments$Medians  <- exp(mus)
      r.new$Model$Segments$Payments$SDs      <- sqrt((exp(sigmas^2) - 1) * exp(2 * mus + sigmas^2))
    } else {
      r.new$Data$Payments <- list()
      r.new$Data$Payments$Means <- pMeans
      r.new$Data$Payments$SDs   <- pSDs
      r.new$Model$Segments$Payments$Means    <- mus
      r.new$Model$Segments$Payments$Medians  <- mus
      r.new$Model$Segments$Payments$SDs      <- sigmas
    }
  }
  
  r.new$Fitted <- list()
  r.new$Fitted$PosteriorSegProbs <- r$posteriors$segment
  #r.new$Fitted$HProbs <- computeHypTestProbs( r$input, r$params )
  
  m <- max(r$input$counts)
  #TODO: take out to separate visualize method
  AdjustGraphics()
  for (cat in 1:K) {
    nbs <- matrix(rep(0:m,numSegs),ncol=(m+1),byrow=TRUE)
    poiMeans <- lambdas[,cat]
    if (r$input$withRandomEffects) {
      poiMeans <- poiMeans * mult[,cat]
    }
    dpoi <- dpois(nbs,poiMeans)
    dpoi[is.na(dpoi)] <- 0
    if (r$input$withInflatedZeros) {
      predicted_tmp <- colSums(segProbs*( dpoi*r$params$buyProbs[,cat] + (1-r$params$buyProbs[,cat])*0^nbs))*Nall 
    } else {
      predicted_tmp <- colSums(segProbs * dpoi) * Nall
    }
    observed_tmp <- aggregate(r$input$original.counts[,cat],list(r$input$original.counts[,cat]),FUN=length)
    observed_tmp <- merge(data.frame(Group.1=0:m), observed_tmp,all.x=TRUE)
    observed_tmp[is.na(observed_tmp[,2]), 2] <- 0
    plot(observed_tmp[,1], observed_tmp[,2],
         xlim=c(0,15), lwd=5, col=rgb(0,0,0,0.75), type="h",
         xlab=sprintf("Number of purchases in %s", r$input$names[cat]),
         ylab="Number of users")
    points(0:m, predicted_tmp, type="o",col=rgb(0,0,1,0.5),lwd=2)
    if (cat==1) {
      predicted = cbind(0:m, predicted_tmp)
      observed  = observed_tmp
    } else {
      predicted = cbind(predicted, predicted_tmp)
      observed  = cbind(observed,  observed_tmp[,-1])
    }
  }
  colnames(predicted) <- c("Count", r$input$names)
  colnames(observed)  <- c("Count", r$input$names)
  r.new$Fitted$Categories$Expected.Count.Distributions <- predicted 
  r.new$Fitted$Categories$Observed.Count.Distributions <- observed
  
  if (r$input$withPayments) {
    for (cat in 1:K) {
      hist(r$input$pMeans[,cat], breaks=200, xlim=c(4,10), freq=FALSE, main="",
           xlab=sprintf("LogPayment Means in %s", r$input$names[cat]))
      x <- seq(from=4,to=10,length.out=200)
      y <- matrix(rep(x,each=numSegs),ncol=200)
      y <- segProbs * dnorm(y, r$params$mus[,cat], r$params$sigmas[,cat])
      y <- colSums(y)
      points(x, y, type='l', lwd=2, col=rgb(0,0,1,0.5))
    }
  }
  
  class(r.new) = "PoissonEMResultSet"
  return(r.new)
}








summary.PoissonEMResultSet = function (result, roundPaymentsToFactorOf = 50)
{
  s = data.frame(SegProbs = round(result$segProbs, 2))
  
  ret = cbind(s, round(result$lambdas))
  #if(result$withPayments)
  #    ret = cbind(ret, round(result$mus, 2))
  
  if (result$logNormalPayments) {
    means   <- exp(result$mus + 0.5*result$sigmas^2)
    medians <- exp(result$mus)
    sds     <- sqrt( (exp(result$sigmas^2)-1) * exp(2*result$mus + result$sigmas^2)  )
  } else {
    means   <- result$mus    
    medians <- means;
    sds     <- result$sigmas
  }
  means   <- round(means/roundPaymentsToFactorOf)*roundPaymentsToFactorOf
  medians <- round(medians/roundPaymentsToFactorOf)*roundPaymentsToFactorOf
  sds     <- round(sds/roundPaymentsToFactorOf)*roundPaymentsToFactorOf
  
  
  m <- paste(rep("(",length(means)), paste(means, sds,sep=","), rep(")",length(means)),sep="")
  if (result$logNormalPayments) {
    m <- paste(m, medians, sep=" | ")
  }
  m <- matrix(m, nrow=nrow(means))
  colnames(m) <- colnames(result$mus)
  ret <- cbind(ret, m)
  return(ret)
}








computeHypTestProbs <- function (x, params) {
  return(0); # TODO
  Pmat <- matrix(rep(params$segProbs,x$Nall),nrow=x$Nall,byrow=TRUE)
  
  # Cell Probabilities: vector of length=numCats containing probabilities of
  # observing the count+payment of each category for Individual n, given segment s
  cellProbs <- function (n, s) {
    probs <- ppois( x$original.counts[n,], params$lambdas[s,], lower.tail=FALSE )
    if (x$withPayments) {
      useInd <- x$pMask[n,]
      probs[useInd] <- probs[useInd] * pnorm( x$pMeans[n,useInd],
                                              params$mus[s,useInd], params$sigmas[s,useInd]/sqrt(x$counts[n,useInd]), lower.tail=FALSE)
    }
    return(probs);
  }
  
  jointProbs <- mapMatrix(Pmat, function(p,n,s) {
    return(prod( cellProbs(n,s) )) # 'dpois' is the only place where conditional prob functional form has a role
  })
  
  #jointProbs[rowSums(jointProbs) == 0,] <- 1e-300; # because 0s break log-likelihood
  
  return(rowSums(jointProbs*params$segProbs))
  #     return(apply(jointProbs, MARGIN=1, prod));
}

