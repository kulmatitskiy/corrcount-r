library(car) # why?

AdjustGraphics <- function() {
  par(mfrow=c(1,1), mex=0.6,
      oma=c(0,0,0,0),
      mar=c(5,5,1,1),
      cex=0.95, cex.axis=0.8, cex.lab=1, cex.main=1,
      mgp=c(3,1,0),
      tcl=-0.75,
      las=0)
}


BubblePlot <- function( names, x, y, 
                        x.bubble, y.bubble, x.bubble.sd, y.bubble.sd, 
                        props, collapse=TRUE, log=FALSE) {
  data <- cbind(x,y)
  
  if (collapse) {
    tmp <- aggregate(rep(1, nrow(data)), by=as.list(data.frame(data)), FUN=sum)
    combos <- as.matrix(tmp[, 1:2])
    freqs  <- tmp[,3]
    Ncmb   <- nrow(tmp)
    plot( combos[,2] ~ combos[,1], 
          cex=50 * freqs / nrow(data), 
          pch=20, col="#cccccc",
          xlab=names[1],ylab=names[2],
          xlim=c(0, max(x.bubble) + 1), 
          ylim=c(0, max(y.bubble) + 1))
  } else {
    if (log) {
      plot( y ~ x, pch=20, col="#cccccc", log="xy",
            xlab=names[1], ylab=names[2], 
            xlim=quantile(x,c(0.1,0.9),na.rm=TRUE), 
            ylim=quantile(y,c(0.1,0.9),na.rm=TRUE))
    } else {
      plot( y ~ x, pch=20, col="#cccccc",
            xlab=names[1], ylab=names[2],
            xlim=c(0, max(x,na.rm=TRUE) + 1), ylim=c(0, max(y, na.rm=TRUE) + 1))
    }
  }
  
  for (s in 1:length(x.bubble)) {
    if (log) {
      break
    #             ellipse( c(x.bubble[s],y.bubble[s]),
    #                      matrix(c(x.bubble[s]*(x.bubble.sd[s]-1),0,0,y.bubble[s]*(y.bubble.sd[s]-1)),nrow=2),radius=1,log="xy",col="#bbbbff",center.pch=FALSE,lwd=1)
    } else {
      ellipse(c(x.bubble[s], y.bubble[s]),
              matrix(c(x.bubble.sd[s], 0, 0, y.bubble.sd[s]), nrow=2),
              radius=1, col="#bbbbff", center.pch=FALSE, lwd=1)
    }
  }
  
  points(y.bubble ~ x.bubble, cex=50 * props)
  points(y.bubble ~ x.bubble, pch=3)
}


LLPlot <- function(results) {
  graphics1()
  segNums <- c()
  LL <- c()
  AIC <- c()
  BIC <- c()
  for (i in 1:length(results)) {
    segNums = c(segNums, results[[i]]$Model$Segments$Number)
    LLval = results[[i]]$Proc$LL
    Nparams = results[[i]]$Model$NumberParams
    LL = c(LL, LLval)
    AIC = c(AIC, LLval - Nparams)
    BIC = c(BIC, LLval - 0.5*Nparams*log(length(results[[i]]$Data$Counts[,1])))
  }
  plot(LL ~ segNums, type='b', col="blue", xlab="Number of Segments")
  #     points(AIC~segNums, type='b', col="green")
  #     points(BIC~segNums, type='b', col="red")
}


SegmentProbabilityPlot = function (results) {
  segNumRange <- c()
  for (i in 1:length(results)) {
    segNumRange = c(segNumRange, results[[i]]$Model$Segments$Number)
  }
  m <- max(segNumRange)
  s <- c()
  p <- c()
  for (i in 1:length(results)) {
    s <- c(s, rep(i, times = m))
    p <- c(p, sort(results[[i]]$Model$Segments$Probs, decreasing = TRUE), rep(0, times = m - i))
  }
  barplot(table(s,s) * 0 + p, horiz = TRUE)
  
  #     for (segNum in 1:3)
  #     {
  #         lambdas  = results[[segNum]]$lambdas
  #         segProbs = sort(results[[segNum]]$segProbs, decreasing=TRUE)
  #         left = 0
  #         y = (round(results[[segNum]]$lambdas))[order(-results[[segNum]]$segProbs),]
  #         for(i in 1:segNum)
  #         {
  #             p = segProbs[i];
  #             if ( p < 0.1)
  #                 break;
  #             left = left+p/2
  #             text(left,segNum,paste(sprintf("%d",y[i,]), collapse=","),col="#ffffff",cex=0.5,pos=4,add=TRUE)
  #             left = left+p/2
  #         }
  #     }
}




# Plot implementation for PoissonEMResultSet
plot.PoissonEMResultSet = function (result) {
  AdjustGraphics()
  lambdas <- round(result$Model$Segments$Counts$Means);
  segProbs <- result$Model$Segments$Probs;
  
  counts <- result$Data$Counts
  names <- colnames(counts)
  if ("LogPayments.Means" %in% names(result$Data)) {
    pMeans <- result$Data$LogPayments.Means
    mus    <- result$Model.Segments.LogPayments.Means
    sigmas <- result$Model.Segments.LogPayments.SDs
  }
  numCats <- ncol(lambdas);
  catNames <- colnames(lambdas)
  
  for (i in 1:(numCats - 1)) {
    for (j in (i + 1):numCats) {
      BubblePlot(names[c(i,j)], counts[, i], counts[, j],
                 lambdas[, i], lambdas[, j], 
                 sqrt(lambdas[, i]), sqrt(lambdas[, j]), segProbs)
      if ("LogPayments.Means" %in% names(result$Data)) {
        #bubble.plot(roundTo(pMeans[,i],log(100)), roundTo(pMeans[,j], log(100)), mus[,i], mus[,j], segProbs, FALSE ) 
        BubblePlot(names[c(i, j)], exp(pMeans[, i]), exp(pMeans[, j]),
                    exp(mus[, i]), exp(mus[, j]),
                    exp(sigmas[, i]), exp(sigmas[, j]), segProbs, FALSE, TRUE) 
      }
    }
  }
  
  
  #     colPalette = c( rgb(rep(19,20),0:19,rep(0,20),maxColorValue=19),
  #                     rgb(18:0,rep(19,19),rep(0,19),maxColorValue=19) )
  #     getCol = function (x)
  #     {
  #         B = 19^(-1/3)
  #         return(colPalette[round(40/(1+B^x))])
  #     }
  #     
  #     x = lambdas[,1]
  #     xPmtMean = result$fitStats$Payments.MarginalByColumn$modelMean[1]
  #     xPmtSD   = result$fitStats$Payments.MarginalByColumn$realSD[1] # change to ModelSD
  #     xPmt = (result$Model.Payment.Means[,1]-xPmtMean)/xPmtSD
  #     y = lambdas[,2]
  #     yPmt = result$mus[,2]/result$fitStats$Payments.MarginalByColumn$modelMean[2]
  #     yPmtMean = result$fitStats$Payments.MarginalByColumn$modelMean[2]
  #     yPmtSD   = result$fitStats$Payments.MarginalByColumn$realSD[2] # change to ModelSD
  #     yPmt = (result$Model.Payment.Means[,2]-yPmtMean)/yPmtSD
  #     
  #     for(s in 1:length(segProbs))
  #     {
  #         draw.sector( center = c(x[s], y[s]),
  #                      start.degree = -135, end.degree = 45,
  #                      rou1 = 5*result$segProbs[s],
  #                      col = getCol(xPmt[s]), border = "black")
  #         draw.sector( center = c(x[s], y[s]),
  #                      start.degree = 45, end.degree = 225,
  #                      rou1 = 5*result$segProbs[s],
  #                      col = getCol(yPmt[s]), border = "black")
  #     }
}

