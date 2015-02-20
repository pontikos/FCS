library(cluster)
library(feature)

# quantile normalisation
quantile.normalize <- function(x0, x1, quantiles=c(.25,.5,.75)) {
    q.x0 <- quantile(x0, probs=quantiles)
    q.x1 <- quantile(x1, probs=quantiles)
    m <- lm(q.x0 ~ q.x1)
    return( cbind(1,x1)%*%coefficients(m) )
}

# peak normalization
peak.normalize <- function(x0, x1, k=2) {
    p.x0 <- peaks(x0, k=k)
    p.x1 <- peaks(x1, k=k)
    m <- lm(p.x0 ~ p.x1)
    return( cbind(1,x1)%*%coefficients(m) )
}


peaks <- function(x, k=2) {
  if (length(x)<10000){
    return(sort(pam(x, k)$medoids))
  } else {
    return(sort(kmeans(x, k)$center))
  }
}


inflection.points <- function(x) {

}

# mode normalization
mode.normalize <- function(x0, x1) {
    q.x0 <- c(mean(x0), Mode(x0), median(x0))
    q.x1 <- c(mean(x1), Mode(x1), median(x1))
    m <- lm(q.x0 ~ q.x1)
    return( cbind(1,x1)%*%coefficients(m) )
}



# flat gradient normalization
# try to align the mean of the null gradient regions of x0
# with those of x1
features.normalize <- function(x0, x1, scaleData=FALSE) {
    f <- function(f) {
        #plot(f, addSignifGradRegion=TRUE) 
        d <- data.frame(x=f$x, curvData=f$curvData, gradData=f$gradData)
        d <- d[order(d$x),]
        # significant curvature curvData, 2nd derivative, not very sensitive
        # significant gradient gradData, 1st derivative, much more sensitive
        r <- rle(as.numeric(d$gradData))
        cuts <- c(1, cumsum(r$lengths))
        runs <- cut(d$x, breaks=d$x[cuts], include.lowest=TRUE) 
        #signif regions
        #signif.regions <- unique(runs)[!as.logical(r$values)]
        unlist(tapply(d$x, runs, mean)[!as.logical(r$values)])
    }
    print(q.x0 <- f(featureSignif(x0, scaleData=scaleData)))
    print(q.x1 <- f(featureSignif(x1, scaleData=scaleData)))
    m <- lm(q.x0 ~ q.x1)
    return( cbind(1,x1)%*%coefficients(m) )
}

# flowPeaks method
gaussNorm.normalize2 <- function(x0, x1, max.lms=3) {
    print(q.x0 <- sort(extract.landmarks(x0, max.lms=max.lms)$lms))
    print(q.x1 <- sort(extract.landmarks(x1, max.lms=max.lms)$lms))
    if (length(q.x0)!=length(q.x1)) return(x1)
    m <- lm(q.x0 ~ q.x1)
    return( cbind(1,x1)%*%coefficients(m) )
}


# flowPeaks method
# Calculate base peaks
gaussNorm.normalize <- function(X) {
    lms <- lapply(X, extract.landmarks, max.lms=10)
    #number of base peaks
    n <- Mode(sapply(lms, function(x) length(x$lms)))
    #base peaks
    base.lms <- sort(colMedians( na.omit(do.call('rbind', lapply(lms, function(x) { if (length(x$lms)==n) return(sort(x$lms)) else return(NA)}))) ))
    f <- function(x1, lms) {
        #all combination of base landmarks and landmarks found in sample
        g <- expand.grid(base.lms, lms$lms)
        #match landmarks
        lms <- sort(lms$lms[expand.grid(1:length(base.lms), 1:length(lms$lms))[order(abs(g[,1]-g[,2]))[1:length(base.lms)],2]])
        print(q.x0 <- base.lms)
        print(q.x1 <- lms)
        if (length(q.x0)!=length(q.x1)) return(x1)
        m <- lm(q.x0 ~ q.x1)
        return( cbind(1,x1)%*%coefficients(m) )
    }
    return(mapply(f, X, lms, SIMPLIFY=FALSE))
}


colMedians <- function(x) {
    apply(x, 2, median)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


#could use splinefun instead
returny <- function(A, X) sapply(x, function(x) A$y[which.min(abs(A$x-x))])


sliding.window.peaks <- function(dens, peak.density.threshold=0.05, peak.distance.threshold=0.05, span=40) {
  y <- dens$y
  x <- dens$x
  ## sliding a window of size span and returns locations where the middle point in the window is maximum.  
  ## returns the indexes of the peaks
  ind <- c()
  for( i in 1:(length(y)-span)) {
    mid <- i+span%/%2
    if ( y[mid]==max(y[i:(i+span)]) & y[mid]!=y[i] & y[mid]!=y[i+span] ) ind <- c(ind, mid)
  }
  return( data.frame(ind=ind,x=x[ind],y=y[ind]) )
}


sliding.window.peaks <- function(y, peak.density.threshold=0.05, peak.distance.threshold=0.05, span=40) {
  ## sliding a window of size span and returns locations where the middle point in the window is maximum.  
  ## returns the indexes of the peaks
  ind <- c()
  for( i in 1:(length(y)-span)) {
    mid <- i+span%/%2
    if ( y[mid]==max(y[i:(i+span)]) & y[mid]!=y[i] & y[mid]!=y[i+span] ) ind <- c(ind, mid)
  }
  return(ind)
}


sliding.window.peaks <- function(d, span=40) {
  y <- d$y
  x <- d$x
  ## sliding a window of size span and returns locations where the middle point in the window is maximum.  
  ## returns the indexes of the peaks
  ind <- c()
  for( i in 1:(length(y)-span)) {
    mid <- i+span%/%2
    if ( y[mid]==max(y[i:(i+span)]) & y[mid]!=y[i] & y[mid]!=y[i+span] ) ind <- c(ind, mid)
  }
  return(x[ind])
}


# returns the top K sliding window peaks
top.sliding.window.peaks <- function(d, K, span=40) {
  y <- d$y
  x <- d$x
  ## sliding a window of size span and returns locations where the middle point in the window is maximum.  
  ## returns the indexes of the peaks
  ind <- c()
  for( i in 1:(length(y)-span)) {
    mid <- i+span%/%2
    if ( y[mid]==max(y[i:(i+span)]) & y[mid]!=y[i] & y[mid]!=y[i+span] ) ind <- c(ind, mid)
  }
  peaks <- cbind(x=x[ind],y=y[ind])
  top.peaks <- peaks[order(peaks[,'y'],decreasing=TRUE)[1:K],]
  top.peaks <- top.peaks[order(top.peaks[,'x']),]
  return(top.peaks)
}




  ## scores the peaks
  ## the score of a peak is a function of its sharpness and height 
  #print(peaks$score <- score.peaks(dens, peaks$ind, peak.density.threshold, span))
  #peaks <- as.data.frame(peaks)
  ### If two peaks are very close to each other then remove one with lowest score.
  #for (i in which(diff(x[peaks$ind]) < diff(range(x))*peak.distance.threshold)) {
    ##remove peak with lowest score
    #peaks <- peaks[-(i+which.min(peaks$score[(i-1):i])),]
  #}
  #return(peaks)
#}
#

## Assigns a score to each peak the score of a peak is a function of its sharpness and density value.
## The peaks with density value less than peak.density.threshold*maximum.peak.density are discarded.
## Of the peaks with distance less than peak.distance.threshold*range.data only one is considered.
score.peaks <- function(dens, peak.indexes, peak.density.threshold, span) {    
  x <- dens$x
  y <- dens$y
  score <- numeric(length(peak.indexes))
  if(length(peak.indexes) == 0) return(score)
  peak.max.y <- max(y[peak.indexes])
  #
  for (peak.ind in peak.indexes) {
    #get a window of span 
    ind <- (max(peak.ind-span%/%2, 1)) : (min(peak.ind+span%/%2, length(x)))
    if(length(ind)==0) ind <- 1
    #if less than a pct of the max then not a peak
    if(y[peak.ind] <  peak.density.threshold*peak.max.y){    
      score <- c(score, 0)
    } else{
      ## computing the sharpness
      # over the window of size 64
      w <- y[peak.ind]-y[ind]
      # if it's smaller than neighbouring
      # peaks than make even more negative by multiplying by 3?!
      w[which(w<0)] <- 3*w[which(w<0)]
      score <- c(score, max(sum(w)*y[peak.ind],0))
    }
  }
  return(score)  
}



## the transform: manipulates the data in such a way that the landmark at matched.lms[i] is moved to matched.lms[i+1] for each i.
register.channel <- function(x, matched.lms, matched.lms2){
  s <- m <- shift <- vector()
  lms <- vector()
  for(i in seq(1,length(matched.lms))) {
    shift <- append(shift, matched.lms2[i]-matched.lms[i])
    lms <- append(lms, matched.lms[i]) # 
    s <- append(s, sd(x)) # sd of landmark
  }
  r.data <- register.function(x, s, lms, shift)
  return(r.data)
}

#shift the data
register.function <- function(x, s, m, shift) {    
  if(length(m)==1) return(x+shift)
  if(length(m)==2){
    sh <- (shift[1]-shift[2])
    x <- x+gau(x, s[1], m[1])*(sh/2)
    x <- x-gau(x, s[2], m[2])*(sh/2)
    return(x+shift[1]-sh/2)
  }
  max.shift <- which.max(abs(shift))
  if(shift[max.shift]>0)
    sh=(shift[max.shift]-(shift[max.shift]-min(shift[-max.shift]))/2)
  else
    sh=(shift[max.shift]-(shift[max.shift]-max(shift[-max.shift]))/2)
  x <- x+sh
  shift <- shift-sh
  for(i in 1:length(m))
    x <- x+gau(x, s[i], m[i])*shift[i]    
  return (x)
}


## gaussian function used in shifting the data.
gau <- function(d, s, m) return(exp(-(d-m)**2/(2*s**2)))


featureSignif <- function (x, bw, gridsize, scaleData = FALSE, addSignifGrad = TRUE, addSignifCurv = TRUE, signifLevel = 0.05) {
    tau <- 5
    n <- length(x)
    names.x <- deparse(substitute(x))
    if (scaleData) x <- (x - min(x))/(max(x) - min(x))
    x <- as.matrix(x)
    if (missing(gridsize)) gridsize <- 401
    if (missing(bw)) {
        bw.range <- dfltBWrange(as.matrix(x), tau)
        bw <- matrix(unlist(bw.range), nrow = 2, byrow = FALSE)
        dfltCounts.out <- dfltCounts(x, gridsize, apply(bw, 2, max))
        h.low <- bw[1, ]
        h.upp <- bw[2, ]
        hmix.prop <- 1/4
        h.init <- h.low^(hmix.prop) * h.upp^(1 - hmix.prop)
        h <- h.init
    }
    else {
        dfltCounts.out <- dfltCounts(x, gridsize, bw)
        h <- bw
    }
    gcounts <- dfltCounts.out$counts
    range.x <- dfltCounts.out$range.x
    dest <- drvkde(gcounts, drv = rep(0, d), bandwidth = h, binned = TRUE, range.x = range.x, se = FALSE, gridsize = gridsize)
    dest$est[dest$est < 0] <- 0
    SignifFeatureRegion.mat <- SignifFeatureRegion(n, d, gcounts, gridsize, dest, h, signifLevel, range.x, grad = addSignifGrad, curv = addSignifCurv)
    ESS <- n * dest$est * prod(h) * (sqrt(2 * pi)^d)
    SigESS <- ESS >= 5
    SignifGradRegion.mat <- SignifFeatureRegion.mat$grad & SigESS
    SignifGradData.mat <- SignifFeatureData(x, d, dest, SignifGradRegion.mat)
    SignifGradDataPoints <- x[SignifGradData.mat, ]
    SignifCurvRegion.mat <- SignifFeatureRegion.mat$curv & SigESS
    SignifCurvData.mat <- SignifFeatureData(x, d, dest, SignifCurvRegion.mat)
    SignifCurvDataPoints <- x[SignifCurvData.mat, ]
    feat <- c(list(x = x, names = names.x, bw = h, fhat = dest),
              SignifFeatureRegion.mat,
              list( gradData = SignifGradData.mat, 
                    gradDataPoints = SignifGradDataPoints,
                    curvData = SignifCurvData.mat, 
                    curvDataPoints = SignifCurvDataPoints))
    class(feat) <- "fs"
    return(feat)
}


########## R function: SignifFeatureRegion ##########

# For determining the region of significant
# gradient for a particular bandwidth and
# significance level.

# Last changed: 18 JAN 2006

SignifFeatureRegion <- function(n, gcounts, gridsize, dest, bandwidth, signifLevel, range.x, grad=TRUE, curv=TRUE, neg.curv.only=TRUE) {
  h <- bandwidth
  ESS <- n*dest$est*prod(h)*(sqrt(2*pi)^d)
  SigESS <- ESS >= 5
  Sig.scalar <- array(NA, dim=gridsize)
  Sig2.scalar <- array(NA, dim=gridsize)
  dest$est[dest$est<0] <- 0  
  ## constant for variance of gradient estimate
  Sig.scalar <- 1/2*(2*sqrt(pi))^(-d)*n^(-1)*prod(h)^(-1)*dest$est
  ##  constants for variance of curvature estimate  
  Sig2.scalar <- (8*sqrt(pi)*n*prod(h))^(-1)*dest$est
  if (grad)
  {
    obj1 <- drvkde(gcounts, drv=1, bandwidth=h, binned=TRUE, range.x=range.x, se=FALSE)
    fhat1 <- obj1$est
    Sig.inv12 <- 1/sqrt(Sig.scalar * h^(-2))
    WaldGrad <- (Sig.inv12 * fhat1)^2
  }
  if (curv)
  {
    obj2 <- drvkde(gcounts,drv=2,bandwidth=h,binned=TRUE,range.x=range.x, se=FALSE)
    fhat2 <- obj2$est
    Sig2.inv12 <- 1/sqrt(Sig2.scalar * 3*h^(-4))
    lambda1 <- Sig2.inv12 * fhat2
    WaldCurv <- lambda1^2
    local.mode <- (lambda1 < 0)
  }
  ## multiple hypothesis testing - based on Hochberg's method
  ## - modified Bonferroni method using ordered p-values
  ## test statistic for gradient
  if (grad)
  {
    pval.Grad <- 1 - pchisq(WaldGrad, 1)
    pval.Grad.ord <- pval.Grad[order(pval.Grad)]
    num.test <- sum(!is.na(pval.Grad.ord))

    if (num.test>=1)
      num.test.seq <- c(1:num.test, rep(NA, prod(gridsize) - num.test))
    else
      num.test.seq <- rep(NA, prod(gridsize))
    
    reject.nonzero <- ((pval.Grad.ord <= signifLevel/(num.test + 1 - num.test.seq)) & (pval.Grad.ord > 0))  
    reject.nonzero.ind <- which(reject.nonzero)

    ## p-value == 0 => reject null hypotheses automatically
    SignifGrad <- array(FALSE, dim=gridsize)
    SignifGrad[which(pval.Grad==0, arr.ind=TRUE)] <- TRUE
    
    ## p-value > 0 then reject null hypotheses indicated in reject.nonzero.ind 
    for (i in reject.nonzero.ind)
      SignifGrad[which(pval.Grad==pval.Grad.ord[i], arr.ind=TRUE)] <- TRUE 
  }
  ## test statistic for curvature
  if (curv)
  {
    pval.Curv <- 1 - pchisq(WaldCurv, 1)
    pval.Curv.ord <- pval.Curv[order(pval.Curv)]
    num.test <- sum(!is.na(pval.Curv.ord))

    if (num.test>=1)
      num.test.seq <- c(1:num.test, rep(NA, prod(gridsize) - num.test))
    else
      num.test.seq <- rep(NA, prod(gridsize))
    reject.nonzero <- ((pval.Curv.ord <= signifLevel/(num.test + 1 - num.test.seq)) &(pval.Curv.ord > 0))  
    reject.nonzero.ind <- which(reject.nonzero)

    SignifCurv <- array(FALSE, dim=gridsize)

    ## p-value == 0 => reject null hypotheses automatically
    SignifCurv[which(pval.Curv==0, arr.ind=TRUE)] <- TRUE

    ## p-value > 0 then reject null hypotheses indicated in reject.nonzero.ind
    for (i in reject.nonzero.ind)
      SignifCurv[which(pval.Curv==pval.Curv.ord[i], arr.ind=TRUE)] <- TRUE 

    if (neg.curv.only) SignifCurv <- SignifCurv & local.mode
  }
  if (grad & !curv) return(list(grad=SignifGrad))
  else if (!grad & curv) return(list(curv=SignifCurv))
  else if (grad & curv) return(list(grad=SignifGrad, curv=SignifCurv))
}


########## End of SignifFeatureRegion ##########


load('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/transforms.RData')

peaks <- function(gate=NULL,chan='CD4', from=.5, to=4, n=512, K=3) {
    fcs.files <- list.files(path='~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/', pattern='.*.RData', full.names=TRUE)
    #fcs.files <- list.files(path='~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/join/All/', pattern='.*.RData', full.names=TRUE)[c(14,45)]
    #fcs.files <- "~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/join/All//CB01513T_2012-11-29.RData"
    length(fcs.files)
    DENS <- NORM.DENS <- matrix(0, ncol=n, nrow=length(fcs.files))
    peaks <- matrix(0, ncol=K, nrow=length(fcs.files))
    #plot(NULL, xlim=c(0.5,3), ylim=c(0,0.01))
    for (i in  1:length(fcs.files)) {
        print(f <- fcs.files[[i]])
        print(load(f))
        if (!is.null(gate)) {
            print(load(file.path('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/magnetic-manual-gates2/CLR',basename(f))))
            fcs.data <- fcs.data[which(as.logical(CLR[,gate])),]
        }
        fcs.data <- applyTransforms(fcs.data, transforms)
        x <- fcs.data[,chan]-fcs.data[,'PSTAT5.1']
        x <- x[percentile.filter(x)]
        x <- sort(x)
        d <- normalised.density(x,from=from,to=to,n=n)
        #unnormalised density
        DENS[i,] <- d$y
        dens <- splinefun(x=d$x,y=d$y)
        #pick the 3 highest peaks
        #browser()
        sw.peaks.x <- sliding.window.peaks(d)
        sw.peaks.x <- sort(sw.peaks.x[order(dens(sw.peaks.x),decreasing=TRUE)[1:K]])
        peaks[i,] <- sw.peaks.x
    }
    individual <- do.call('rbind', (strsplit(gsub('.RData','',basename(fcs.files)),'_')))[,1]
    date <- do.call('rbind', (strsplit(gsub('.RData','',basename(fcs.files)),'_')))[,2]
    return(list(individual=individual, date=date, x=d$x, dens=DENS, peaks=peaks))
} 

plot.peaks <- function(d, main='',xlab='CD4') {
    plot(NULL, xlim=range(d$x), ylim=range(d$dens), main=main, xlab=xlab, ylab='')
    #densities
    for (i in 1:nrow(d$dens)) {
        lines(d$x, d$dens[i,], lwd=.5, col=alpha('black',.5))
        dens <- splinefun(x=d$x,y=d$dens[i,])
        x.peaks <- d$peaks[i,]
        points(cbind(x.peaks,dens(x.peaks)),col=1:length(x.peaks),pch=20, cex=2)
   }
}

normalise.peaks <- function(d, chan, gate, from=.5, to=4, n=512, K=3) {
    DENS <- NORM.DENS <- matrix(0, ncol=ncol(d$dens), nrow=nrow(d$dens))
    peaks <- matrix(0, ncol=K, nrow=nrow(d$dens))
    for (i in  1:nrow(d$dens)) {
        individual <- d$individual[i]
        date <- d$date[i]
        print( f <- sprintf('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/%s.RData', paste(individual, date, sep='_')))
        load(f)
        if (!is.null(gate)) {
            print(load(file.path('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/magnetic-manual-gates2/CLR',basename(f))))
            fcs.data <- fcs.data[which(as.logical(CLR[,gate])),]
        }
        print(p <- as.numeric(d$peaks[i,1:2]))
        m <- lm(colMedians(d$peaks) ~ p)
        normalise <- function(x) cbind(1,x)%*%coefficients(m) 
        x <- normalise(fcs.data[,chan]-fcs.data[,'PSTAT5.1'])
        dens <- normalised.density(x,from=from,to=to,n=n)
        #unnormalised density
        DENS[i,] <- dens$y
        #pick the 3 highest peaks
        sw.peaks.x <- sliding.window.peaks(dens)
        dens <- splinefun(x=dens$x,y=dens$y)
        sw.peaks.x <- sort(sw.peaks.x[order(dens(sw.peaks.x),decreasing=TRUE)[1:K]])
        peaks[i,] <- sw.peaks.x
    }
    return(list(individual=d$individual, date=d$date, x=d$x, dens=DENS, peaks=peaks))
} 


ungated.pstat5.4 <- peaks(chan='PSTAT5.4',from=-1,to=3, K=2)
plot.peaks(ungated.pstat5.4,xlab='pSTAT5 at 1000U in lymphocytes')

#
lymphocytes.pstat5.4 <- peaks(gate='Lymphocytes',chan='PSTAT5.4',from=-1,to=3, K=2)
norm.lymphocytes.pstat5.4 <- normalise.peaks(d=lymphocytes.pstat5.4,chan='PSTAT5.4', gate='Lymphocytes', from=-1, to=3, K=2)

cd4positive.pstat5.3 <- peaks(gate='CD4',chan='PSTAT5.3',from=-1,to=3, K=2) 
norm.cd4positive.pstat5.3 <- normalise.peaks(d=cd4positive.pstat5.3,chan='PSTAT5.3', gate='CD4', from=-1, to=2, K=2)

pdf('~nikolas/GoogleDrive/PhD/Thesis/IL2/figures/pstat5-peak-normalisation.pdf',width=10,height=10)
par(mfrow=c(2,2))
figure.labels <- iter(paste(letters,')',sep=''))
plot.peaks(lymphocytes.pstat5.4,xlab='pSTAT5 at 1000U in lymphocytes')
title(nextElem(figure.labels), adj=0)
plot.peaks(norm.lymphocytes.pstat5.4,xlab='pSTAT5 at 1000U in lymphocytes')
title(nextElem(figure.labels), adj=0)
plot.peaks(cd4positive.pstat5.3,xlab='pSTAT5 at 10U in CD4+ lymphocytes')
title(nextElem(figure.labels), adj=0)
plot.peaks( norm.cd4positive.pstat5.3, xlab='pSTAT5 at 10U in CD4+ lymphocytes')
title(nextElem(figure.labels), adj=0)
dev.off()



cd4 <- peaks(chan='CD4',to=4.5)
plot.peaks(cd4,xlab='CD4')

cd25 <- peaks(chan='CD25',from=-1, to=3,K=2)
plot.peaks(cd25,xlab='CD25')

# the stain is really bad, the second peak sometimes doesn't exist!
cd45ra <- peaks(chan='CD45RA',from=-1, to=3,K=2)
plot.peaks(cd45ra,xlab='CD45RA')

#plot(normalised.density(x), lwd=.25)

peaks <- data.frame()
for (f in list.files(pattern='*.RData')) {
    load(f)
    fcs.data <- applyTransforms(fcs.data, transforms)
    x <- fcs.data[,'PSTAT5']
    d <- normalised.density(x)
    top.peaks <- top.sliding.window.peaks(d, 2)
    n <- unlist(strsplit(gsub('.RData','',f), '_'))
    peaks <- rbind( peaks, data.frame(individual=n[[1]], dose=n[[2]], date=n[[3]], t(top.peaks[,'x'])) )
    m <- lm(0:1 ~ top.peaks[,'x'])
    x.norm <- cbind(1,x)%*%coefficients(m)
    d.norm <- normalised.density(x.norm)
    d.norm.f <- splinefun(d.norm)
    #lines(d.norm, lwd=.25)
    #points(cbind(0:1,d.norm.f(0:1)), col=1:2, pch=20)
}


for (f in list.files(pattern='CB01510Q_1000U_.*.RData')) {
    print(f)
    load(f)
    x <- fcs.data[,'PSTAT5']
    w <- seq(.1,.9,.1)
    score <- sapply(w, function(w) {
        p <- top.sliding.window.peaks(density(logicleTransform(w=w)(x)),2)
        #neither of the peaks should be in the negatives
        if ( (p[2,'x'] < 0) || (p[1,'x'] < 0) ) return(0)
        return( ( p[2,'x'] - p[1,'x'] ) / abs(p[2,'y'] - p[1,'y']) )
    })
    names(score) <- w
    w.max <- as.numeric(names(which.max(score)))
    plot(density(logicleTransform(w=w.max)(x)), main=f)
    p <- top.sliding.window.peaks(density(logicleTransform(w=w.max)(x)),2)
    points(p, col=1:2, pch=20)
}


plot(density(logicleTransform(w=.1)(x)),2)
for (w in seq(.1,.9,.1)) lines(density(logicleTransform(w=w)(x)))

#individual <- 'KM00744H'
#date <- '2012-07-02'

#~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All-ungated-normalised


#for (f in list.files(path='~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/', pattern='CB01510Q.*.RData', full.names=TRUE)) {
for (f in list.files(path='~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/', pattern='KM00744H.*.RData', full.names=TRUE)) {
    print(f)
    load(f)
    x <- fcs.data[,'PSTAT5.1']
    w <- seq(.1,.9,.1)
    score <- sapply(w, function(w) {
        p <- top.sliding.window.peaks(density(logicleTransform(w=w)(x)),2)
        #neither of the peaks should be in the negatives
        if ( (p[2,'x'] < 0) || (p[1,'x'] < 0) ) return(0)
        return( ( p[2,'x'] - p[1,'x'] ) / abs(p[2,'y'] - p[1,'y']) )
    })
    names(score) <- w
    print(w.max <- as.numeric(names(which.max(score))))
    lgcl <- logicleTransform(w=w.max)
    chan <- paste('PSTAT5',1:4,sep='.')
    #
    fcs.data[,chan] <- apply(fcs.data[,chan],2, function(x) {
        p <- top.sliding.window.peaks(density(lgcl(x)),2) 
        invlgcl <- inverseLogicleTransform(trans=lgcl)
        return ( invlgcl(cbind(1,lgcl(x))%*%coefficients(lm(1:2 ~ p[,'x']))) )
    })
    plot(normalised.density(logicleTransform(w=w.max)(fcs.data[,'PSTAT5.1'])), col='white')
    sapply( 1:4, function(i) lines(normalised.density(logicleTransform(w=.6)(fcs.data[,paste('PSTAT5',i,sep='.')]-fcs.data[,'PSTAT5.1'])),col=blues4[[i]],lwd=i) ) 
}




w.best <- 0.6
from <- -1
to <- 2
w.best <- 1
lgcl <- logicleTransform(w=w.best)
invlgcl <- inverseLogicleTransform(trans=lgcl)
files <- list.files(path='~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/RData/pstat5-join/', pattern='*.RData', full.names=TRUE)
norm.dens <- norm.dens2 <- dens2 <- dens <- matrix(0, nrow=length(files), ncol=512)
norm.peaks <- norm.peaks2 <- peaks2 <- peaks <- list()
for (i in 1:length(files)) {
    f <- files[[i]]
    print(f)
    load(f)
    fcs.data <- apply(fcs.data, 2, lgcl)
    fcs.data <- baseline.relative.pstat5(fcs.data,TRUE)
    print(load(file.path('~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/magnetic-manual-gates2/CLR',basename(f))))
    fcs.data2 <- fcs.data
    fcs.data <- fcs.data[which(as.logical(CLR[,'Single cells'])),]
    chan <- paste('PSTAT5',1:4,sep='.')
    #
    #X <- lgcl(fcs.data[,'PSTAT5.4'])-lgcl(fcs.data[,'PSTAT5.1'])
    X <- fcs.data[,'PSTAT5.4']
    d <- normalised.density((X),n=512,from=from,to=to)
    p <- top.sliding.window.peaks(d,2)
    peaks[[i]] <- p
    dens[i,] <- d$y
    X.norm <- (cbind(1,(X))%*%coefficients(lm(c(0,1) ~ p[,'x'])))
    d.norm <- normalised.density((X.norm),n=512,from=from,to=to)
    p <- top.sliding.window.peaks(d.norm,2)
    norm.peaks[[i]] <- p
    norm.dens[i,] <- d.norm$y
    x.norm <- d.norm$x
    #
    fcs.data <- fcs.data2
    fcs.data <- fcs.data[which(as.logical(CLR[,'CD4'])),]
    X <- fcs.data[,'PSTAT5.3']
    d <- normalised.density((X),n=512,from=from,to=to)
    x <- d$x
    p <- top.sliding.window.peaks(d,2)
    peaks2[[i]] <- p
    dens2[i,] <- d$y
    X.norm <- (cbind(1,(X))%*%coefficients(lm(c(0,1) ~ p[,'x'])))
    d.norm <- normalised.density((X.norm),n=512,from=from,to=to)
    p <- top.sliding.window.peaks(d.norm,2)
    norm.peaks2[[i]] <- p
    norm.dens2[i,] <- d.norm$y
    x.norm <- d.norm$x
}


par(mfrow=c(2,2))
#a
plot(NULL, xlim=range(x), ylim=range(dens))
sapply(1:length(files), function(i) lines(x,dens[i,],lwd=.25))
lapply(peaks, function(p) points(p, col=1:2))
#b
plot(NULL, xlim=range(x), ylim=range(norm.dens))
sapply(1:length(files), function(i) lines(x.norm,norm.dens[i,],lwd=.25))
lapply(norm.peaks, function(p) points(p, col=1:2))
#c
plot(NULL, xlim=range(x), ylim=range(dens2))
sapply(1:length(files), function(i) lines(x,dens2[i,],lwd=.25))
lapply(peaks2, function(p) points(p, col=1:2))
#d
plot(NULL, xlim=range(x), ylim=range(norm.dens2))
sapply(1:length(files), function(i) lines(x.norm,norm.dens2[i,],lwd=.25))
lapply(norm.peaks2, function(p) points(p, col=1:2))







