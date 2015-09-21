source("common.R")
source("algorithm.util.R")


MakeCandidatesObject <- function (x) {
  candidates = list()
  #     candidates = list( LL = list(), segProbs = list(), lambdas = list() )
  #     if (x$withPayments)
  #     {
  #         candidates$mus    = list()
  #         candidates$sigmas = list()
  #     }
  return(candidates)
}

UpdateCandidates = function (candidates, params)
{
  if (length(candidates) == 0 || params$LL > max(candidates[[1]]$LL))
  {
    candidates[[1]] = params
  } 
  return( candidates )
}

getBestCandidate = function (candidates)
{
  return(candidates[[1]])
}




initMusAndSigmas = function (numSegs, pMeans, pSDs)
{
  for ( c in 1:ncol(pMeans) )
  {
    mus_range    = quantile(pMeans[,c],c(0.3,0.7),na.rm=TRUE)
    sigmas_range = quantile(pSDs[,c],  c(0.3,0.7),na.rm=TRUE)
    mus_v    = runif(numSegs,min=mus_range[1],max=mus_range[2])
    sigmas_v = runif(numSegs,min=sigmas_range[1],max=sigmas_range[2])
    if ( c == 1 )
    {
      mus    = mus_v
      sigmas = sigmas_v
    } else {
      mus    = cbind(mus,    mus_v)            
      sigmas = cbind(sigmas, sigmas_v)            
    }
  }
  return( list(mus = mus, sigmas = sigmas) )
}



InitStartingParams <- function (x, procStats) {
  params <- list()
  
  # find a "good" starting point
  repeat {
    # 
    # NOTE: IF WE ARE USING FREQS THEN K-MEANS FOR COUNT DATA IS BAD
    km = kmeans( x$counts, centers=x$numSegs )
    params$lambdas  = km$centers
    params$segProbs = km$size/x$Ncmb
    
    if ( x$withInflatedZeros ) {
      binom = (x$counts > 0)
      km = kmeans( binom, centers=x$numSegs)
      params$buyProbs = km$centers
      params$segProbs = (params$segProbs + km$size/(x$Ncmb))/2
    }
    
    if ( x$withRandomEffects )
    {
      params$rEffMus    = matrix(rep(1, x$numSegs * x$numGroups), nrow=x$numSegs)
      params$rEffSigmas = matrix(rep(0.1, x$numSegs * x$numGroups), nrow=x$numSegs)
    }
    
    if ( x$withPayments )
    {
      tmp = initMusAndSigmas(x$numSegs,x$pMeans, x$pSDs)
      params$mus = tmp$mus
      params$sigmas = tmp$sigmas
    } #else {
    #             params$mus    = params$segProbs * 0
    #             params$sigmas = params$lambdas * 0 
    #         }
    procStats = logTime(procStats, "init")
    
    
    conditionals = computeConditionals(x, params, procStats$MC )
    LL = computeMarginalLL(x, params, conditionals)
    if(is.finite(LL))
    {
      params$LL = LL
      return(list(params = params, procStats = procStats, conditionals=conditionals))
    }
  }
}


# The Algorithm
RunPoissonEM = function ( counts, numSegs,
                          withInflatedZeros = FALSE, hurdle = FALSE, infZerosAsLatent = TRUE,
                          withRandomEffects = FALSE, groupingMatrix = diag(length(counts[1,])),
                          withPayments = FALSE, logNormalPayments = TRUE,
                          pMeans = counts * 0, pSDs = counts * 0,
                          numAttempts = 20, threshold = 1e-8, Tmax = 350, MC = 100,
                          useFreqs = TRUE, log = FALSE, timeLog = FALSE ) 
{
  # In addition to all the arguments in the function definition, we need to
  # use a lot of variables derived from those arguments.
  # Rather than list all these derivations/preparations at the top of 
  # our algorithm's body (i.e., right here), we move them to a separate method
  # called MakeInputObject(args) and group their results into a single list (here it will be 'x')
  x <- MakeInputObject(as.list(environment()))
  
  candidates <- MakeCandidatesObject(x);
  procStats  <- MakeProcessStatisticsObject(numAttempts, threshold, Tmax, log, timeLog, MC)
  
  # Store log-likelihood values in this vector
  LLs <- c()
  
  # Iterate through different starting points
  for (attempt in 1:numAttempts) {
    LogMessage(procStats, "Initializing for attempt %d..\n", attempt)
    
    tmp = InitStartingParams(x, procStats)
    params = tmp$params
    procStats = tmp$procStats
    conditionals = tmp$conditionals
    
    # Run EM iterations given the starting point
    for ( t in 1:Tmax )
    {
      # 1. Use current estimates to compute conditional probabilities of
      # observed data given unobserved data (either all possible values or simulated ones)
      # Note: if withRandomEffects=TRUE, conditionals actually consist of many repeated MC samples to be aggregated later on 
      if (t>1)
        conditionals = computeConditionals(x, params, MC)
      
      # 2. Compute log-likelihood and compare with previous value, check if we can terminate
      if (t>1)
        params$LL = computeMarginalLL(x, params, conditionals) # !!! switch to expected joint
      procStats = logTime(procStats,"LL")
      
      LLs = c(LLs, params$LL)
      if(t>1 && !is.finite(params$LL)) # for debugging
      {
        f = 0
      }
      if(t>1 && (params$LL - LLs[length(LLs)-1] < 0) ) # for debugging
      {
        f = 0
      }
      if ( t>1 && abs(params$LL - LLs[length(LLs)-1])/x$Nall <= threshold ) # IDEALLY ABS NOT NEEDED
      {
        logMessage(procStats, "-- Terminating with %f.\n", params$LL/x$Nall)
        break;
      }
      logMessage(procStats, "-- Iteration %d...(improving on %f)\n", t, LLs[length(LLs)]/x$Nall)
      procStats$total.iterations = procStats$total.iterations + 1
      
      # 3. Compute posterior weights and expectations
      posteriors = computePosteriors(x, params, conditionals)
      procStats = logTime(procStats, "posteriors")
      
      # 4. Update parameter estimates
      params = updateParams(x, params, posteriors)
      procStats = logTime(procStats, "paramUpdates")
    }
    
    candidates = updateCandidates( candidates, params )
  }
  
  # Choose the winner and return object with all relevant stuff in it
  return( fillReturnObject( list( input = x,
                                  conditionals = conditionals,
                                  posteriors = posteriors,
                                  procStats = procStats,
                                  params = getBestCandidate(candidates) )))
}







# -- SIDE COMPUTATIONS --

computeMarginalLL = function (x, params, conditionals) 
{
  segProbs = rep(params$segProbs, each=x$Ncmb)
  return( sum(log(rowSums(conditionals$theJoint * segProbs)) * x$freqs))
}

library(tensor)
computeConditionals = function (x, params, MC_trials=10)
{
  conditionals = list()
  
  # dimensions of a typical cube
  dim = c(x$Ncmb, x$numCats, x$numSegs)
  
  # first we prepare cube(s) with poisson probs
  if ( x$withRandomEffects )
  {
    dim  = c(x$Ncmb, x$numCats, x$numSegs, MC_trials)
    dim2 = c(x$Ncmb, x$numGroups, x$numSegs, MC_trials)
    counts = rep(x$counts, times=MC_trials*x$numSegs ) # create many cubes of counts
    mus  = rep( t(params$rEffMus), each=x$Ncmb) # one cube of mus
    mus  = rep(mus, times=MC_trials) # !!! maybe not needed due to autofill
    sgms = rep( t(params$rEffSigmas), each=x$Ncmb) # on cube of sigmas
    sgms = rep(sgms, times=MC_trials) # !!! maybe not needed due to autofill
    lambdas = rep( t(params$lambdas), each=x$Ncmb ) # one cube of lambdas
    lambdas = rep( lambdas, times=MC_trials ) # maybe not needed?
    lambdas = array( lambdas, dim=dim ) # better explicitly an array/ or not??
    
    Hs = rlnorm( prod(dim2), mus, sgms )  # generate the Hs!
    if(any(is.na(Hs)))
    {
      dfg = 34 # DEBUG PURPOSE
    }
    Hs = array(Hs,dim=dim2) # must be explicitly an array
    
    logHs = log(Hs)
    sqlogHs = logHs^2
    
    Hprods       = exp( tensor(logHs, x$groupingMatrix, 2, 1) ) 
    Hprods       = aperm(Hprods,c(1,4,2,3))
    poissonMeans = lambdas * Hprods
  }
  else
  {
    counts = rep( x$counts, times=x$numSegs ) # cube of counts
    poissonMeans = rep( t(params$lambdas), each=x$Ncmb )
  }
  probs = dpois( counts, poissonMeans )
  
  # then update the cube(s) with buyProbs if needed
  if ( x$withInflatedZeros )
  {
    buyProbs = rep( t(params$buyProbs), each=x$Ncmb) # cube of buyProbs
    probs = buyProbs * probs
    probs[is.na(probs)] = 0 # we want 0*NA = 0 instead of 0*NA = NA
    buyCube = probs # OK to have zero-valued elems
    probs = probs + (1-buyProbs)*(0^counts) 
  }
  
  if (x$withPayments)
  {
    # cubes of mus,sigmas,pMask,pMeans,pSDs
    pMeans = rep( x$pMeans, times=x$numSegs )
    pSDs   = rep( x$pSDs,   times=x$numSegs )
    pMask  = rep( x$pMask,  times=x$numSegs )
    mus    = rep( t(params$mus),    each=x$Ncmb )
    sigmas = rep( t(params$sigmas), each=x$Ncmb )
    
    probs[pMask] = probs[pMask] * dnorm( pMeans[pMask], mus[pMask], sigmas[pMask]/sqrt(counts[pMask]))
  }
  
  # wrap up
  
  conditionals$cube     = array( probs, dim=dim ) 
  if ( x$withInflatedZeros )
    conditionals$buyCube = array( buyCube , dim=dim )
  if ( x$withRandomEffects )
  {
    conditionals$logHs    = logHs
    conditionals$sqLogHs  = sqlogHs
    conditionals$Hprods   = Hprods
    conditionals$joint    = apply( conditionals$cube, MARGIN=c(1,3,4), FUN=prod )
    conditionals$theJoint = apply( conditionals$joint, FUN=mean, MARGIN=c(1,2))
  }
  else
  {
    conditionals$joint    = apply( conditionals$cube, MARGIN=c(1,3), FUN=prod )
    conditionals$theJoint = conditionals$joint 
  }
  
  return(conditionals)
}


computeHypTestProbs = function (x, params)
{
  return(0);
  Pmat = matrix(rep(params$segProbs,x$Nall),nrow=x$Nall,byrow=TRUE)
  
  # Cell Probabilities: vector of length=numCats containing probabilities of
  # observing the count+payment of each category for Individual n, given segment s
  cellProbs = function (n, s)
  {
    probs =  ppois( x$original.counts[n,], params$lambdas[s,], lower.tail=FALSE )
    if (x$withPayments)
    {
      useInd = x$pMask[n,]
      probs[useInd] <- probs[useInd] * pnorm( x$pMeans[n,useInd],
                                              params$mus[s,useInd], params$sigmas[s,useInd]/sqrt(x$counts[n,useInd]), lower.tail=FALSE)
    }
    return(probs);
  }
  
  jointProbs = mapMatrix(Pmat, function(p,n,s)
  {
    return(prod( cellProbs(n,s) )) # 'dpois' is the only place where conditional prob functional form has a role
  })
  
  #jointProbs[rowSums(jointProbs) == 0,] <- 1e-300; # because 0s break log-likelihood
  
  return(rowSums(jointProbs*params$segProbs))
  #     return(apply(jointProbs, MARGIN=1, prod));
}


# -- POSTERIORS UPDATE LOGIC --

computePosteriors = function(x, params, conditionals)
{
  posteriors = list()
  
  if ( x$withRandomEffects )
    joints = apply( conditionals$joint, FUN=mean, MARGIN=c(1,2))
  else
    joints = conditionals$joint
  
  segProbs = rep(params$segProbs, each=x$Ncmb)
  Pmat = joints * segProbs
  phis = rowSums(Pmat)
  posteriors$phis    = phis
  posteriors$segment = Pmat / phis
  posteriors$segmentSums = colSums(posteriors$segment * x$freqs)
  
  if (x$withInflatedZeros)
  {
    mask = (conditionals$buyCube != 0) # exclude zero-valued buyProbs
    if ( x$withRandomEffects )
    {
      MC_trials = dim(conditionals$cube)[4]
      bpp = array(0,dim=c(x$Ncmb,x$numCats,x$numSegs,MC_trials))
      jjj = conditionals$joint # bad naming 'bpp' and 'jjj' but no idea how to name them well :(
      jjj = aperm(array(rep(jjj, each=x$numCats), dim=c(x$numCats,x$Ncmb,x$numSegs,MC_trials)),c(2,1,3,4))
      bpp[mask] = jjj[mask] * conditionals$buyCube[mask] / conditionals$cube[mask]
      res = apply(bpp, FUN=mean, MARGIN=c(1,2,3))
      
      segProbs = rep(params$segProbs, each=x$Ncmb*x$numCats)
      res = res * segProbs / phis 
    }
    else
    {
      res = array(0,dim=c(x$Ncmb,x$numCats,x$numSegs))
      ps = posteriors$segment
      # to understand what the following line is doing, run these:
      #     x = matrix(1:12, nrow=4)
      #     array(t(matrix(rep(t(x),times=2),nrow=nrow(t(x)))), dim=c(4,2,3))
      ps = array(t(matrix(rep(t(ps),times=x$numCats),nrow=x$numSegs)), dim=c(x$Ncmb,x$numCats,x$numSegs))
      
      res[mask] = ps[mask] * conditionals$buyCube[mask] / conditionals$cube[mask]
    }
    posteriors$segmentBuy = res
  }
  
  if ( x$withRandomEffects )
  {
    jj = conditionals$joint 
    jj = aperm(array(rep(jj, each=x$numCats), dim=c(x$numGroups,x$Ncmb,x$numSegs,MC_trials)),c(2,1,3,4))
    logHs   = apply( jj * conditionals$logHs,   FUN=mean, MARGIN=c(1,2,3) )
    sqLogHs = apply( jj * conditionals$sqLogHs, FUN=mean, MARGIN=c(1,2,3) )
    
    segProbs2 = rep(params$segProbs, each=x$Ncmb*x$numGroups)
    posteriors$logHs   =   logHs * segProbs2 / phis
    posteriors$sqLogHs = sqLogHs * segProbs2 / phis
    posteriors$logHsSums = t(apply( posteriors$logHs, FUN=sum, MARGIN=c(2,3) ))
    
    if(x$withInflatedZeros)
    {
      Hprods = apply( bpp * jjj, FUN=mean, MARGIN=c(1,2,3) )
      posteriors$Hprods = Hprods * segProbs / phis
    }
    else
    {
      # TODO
    }
  }
  
  return(posteriors)
}


# -- PARAMS UPDATE LOGIC --

updateParams = function (x, params, posteriors)
{
  params$segProbs = updateSegProbs(x, posteriors)
  
  if ( x$withInflatedZeros )
    params$buyProbs = updateBuyProbs(x, posteriors)
  
  if ( x$withRandomEffects )
  {
    params$rEffMus    = updateREffMus(x, posteriors)
    params$rEffSigmas = updateREffSigmas(x, posteriors, params$rEffMus)
  }
  
  params$lambdas  = updateLambdas(x, posteriors)
  
  if ( x$withPayments )
  {
    params$mus      = updateMus(x, posteriors)
    params$sigmas   = updateSigmas(x, posteriors, params$mus)
  }
  
  return(params)
}

updateSegProbs = function (x, posteriors)
{
  return(colSums(posteriors$segment * x$rfreqs)) # TODO switch to segmentSums / x$Nall for better speed
}

updateBuyProbs = function (x, posteriors)
{
  numerator   = t(apply(posteriors$segmentBuy * x$freqs, MARGIN=c(2,3), FUN=sum))
  buyProbs = numerator / posteriors$segmentSums
  if ( any(buyProbs>1+1e-3))# for debug purposes
  {
    sfsr = 0
  }
  buyProbs[buyProbs>1]=1 # FIX BY BRUTE FORCE
  return( buyProbs )
}

updateLambdas = function (x, posteriors)
{
  if ( x$withInflatedZeros)
  {
    # cube of counts
    counts = rep( x$counts, times=x$numSegs)
    
    numerator   = apply(posteriors$segmentBuy * counts * x$freqs, MARGIN=c(2,3), FUN=sum)
    denominator = apply(posteriors$segmentBuy * x$freqs, MARGIN=c(2,3), FUN=sum)
    lambdas = t(numerator/denominator)
    lambdas[is.nan(lambdas)] = NA#1e-200
    colnames(lambdas) = x$names
  }
  else
  {
    M = posteriors$segment * x$freqs # Observations-to-segment contributions of each row of data
    SC = t(M) %*% x$counts           # Counts-to-segment total 
    # Divide counts-to-Segment by Observations-to-segment 
    lambdas = SC / posteriors$segmentSums
  }
  return(lambdas)
}

updateMus = function (x, posteriors)
{
  X = x$counts*x$pMeans # Sums of payments (gaussians)
  X[!x$pMask] = 0
  Y = x$counts
  Y[!x$pMask] = 0
  return( t(posteriors$segment) %*% X / t(posteriors$segment) %*% Y)
}

updateSigmas = function (x, posteriors, mus)
{
  X = x$pSS #sumSquares(counts,pMeans, pSDs); #(pSDs^2)*(counts-1)+(counts*pMeans^2)
  X[!x$pMask] = 0
  Y = x$counts
  Y[!x$pMask] = 0
  return(sqrt( t(posteriors$segment) %*% X / t(posteriors$segment) %*% x$counts - mus^2 ))
}

updateREffMus = function (x, posteriors) {
  # ?
  return( posteriors$logHsSums / posteriors$segmentSums )
}

updateREffSigmas = function (x, posteriors, mus)
{
  res = t(apply( posteriors$sqLogHs, FUN=sum, MARGIN=c(2,3) )) - 2 * mus * posteriors$logHsSums + mus^2
  res = res / posteriors$segmentSums + mus^2
  if (any(res<0))
  {
    dfgdfg =345;#DEBU purpose
  }
  return( sqrt(res) )
}