# quantile normalisation
quantile.normalize <- function(x0, x1, quantiles=c(.25,.5,.75)) {
    q.x0 <- quantile(x0, probs=quantiles)
    q.x1 <- quantile(x1, probs=quantiles)
    m <- lm(q.x0 ~ q.x1)
    return( cbind(1,x1)%*%coefficients(m) )
}

# peak normalization
peak.normalize <- function(x0, x1, k=2) {
    q.x0 <- pam(x0[sample(x0,1000)], k=k)$medoids
    q.x1 <- pam(x1[sample(x1,1000)], k=k)$medoids
    m <- lm(q.x0 ~ q.x1)
    return( cbind(1,x1)%*%coefficients(m) )
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
    library(feature)
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
    base.lms <- sort(colMedians(na.omit(do.call('rbind', sapply(lms, function(x) { if (length(x$lms)==n) return(sort(x$lms)) else return(NA) })))))
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



Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


extract.landmarks <- function(x, peak.density.threshold=0.05, peak.distance.threshold=0.05, max.lms=3) {
  ## defining output variables.
  lms.list <- list(original=list(), score=list(), lms=list(), dens=list())    
  ## finding the landmarks
  lms <- landmarker(x)
  lms.list$original <- lms
  ## returns the max.lms top score landmarks.
  filtered <- filter.lms(lms, x, max.lms,  peak.density.threshold, peak.distance.threshold)
  lms.list$lms <- filtered$lms
  lms.list$score <- filtered$score
  A <- density(x)
  lms.list$dens <- sapply(filtered$lms, function(lms) returny(A,lms))
  return(lms.list)
}


returny <- function(A, X){
  y=vector()
  i=1;
  for( x in X){
    y[i]=A$y[which(abs(A$x-x)==min(abs(A$x-x)))[1]]
    i=i+1
  }
  return (y)
}



## returns the peaks (local maxima's) in the kernel density estimate of data.
landmarker <- function(x, span=10) {
  A=density(x)
  d=A$y
  pks=c()
  ## sliding a window of size span and returns locations where the middle point in the window is maximum.  
  for( i in 1:(length(d)-span)){
    if (!is.na(d[i+span%/%2]) & (d[i+span%/%2]==max(d[c(i:(i+span))], na.rm=T)))
      if (!is.na(d[i]) & !is.na(d[i+span]) & d[i+span%/%2]!=d[i] & d[i+span%/%2]!=d[i+span])
        pks=append(pks, i+span%/%2)
    
  }
  return(A$x[pks])
}


## returns the max.lms top score landmarks.
## the score of a landmarks is a function of its sharpness and density value
filter.lms <- function(lms, x, max.lms, peak.density.threshold, peak.distance.threshold) {
  filtered <- list()
  if(length(lms) == 0){
    filtered$lms <- vector()
    filtered$score <- vector()    
    return(filtered)
  }
  filtered$score <- score.lms(lms, x, max.lms, peak.density.threshold, peak.distance.threshold)
  lms.score <- data.frame(score=filtered$score, ind=c(1:length(lms)))
  lms.score <- lms.score[do.call(order, c(lms.score["score"], decreasing=T)), ]
  ind <- which(lms.score$score>0)
  if(length(ind)==0){
    filtered$lms <- vector()
    filtered$score <- vector()        
    return(filtered)
  }
  lms.score.ind <- lms.score$ind[ind]
  if(length(lms.score.ind) < max.lms)
    filtered$lms <- sort(lms[lms.score.ind], decreasing=F)
  else
    filtered$lms <- sort(lms[lms.score.ind[c(1:max.lms)]], decreasing=F)
  return(filtered) 
}


## assigns a score to each landmark. the score of a landmarks is a function of its sharpness and density value.
## the peaks with density value less than peak.density.thr*maximum.peak.density are discarded.
## of the peaks with distance less than peak.distance.thr*range.data only one is considered.
score.lms <- function(lms, x, max.lms, peak.density.threshold, peak.distance.threshold) {    
  bw <- 64
  score <- vector()
  height.cutoff <- peak.density.threshold
  if(length(lms) == 0) return(score)
  A <- density(x)
  bw <- min(64, length(A$x)/10)
  lms.max.height <- max(returny(A, lms), na.rm=T)
  MIN.LMS.DIST <- (max(A$x, na.rm=T)-min(A$x, na.rm=T))*peak.distance.threshold
  last.lms <- -1
  last.lms.i <- -1
  last.lms.score <- 0
  for(i in 1:length(lms)){
    lms.ind=which(na.omit(abs(A$x-lms[i]))==min(na.omit(abs(A$x-lms[i])) ) )
    ind=(max(lms.ind-bw%/%2, 1)):(min(lms.ind+bw%/%2, length(A$x)))
    if(length(ind)==0) ind=1
    if(A$y[lms.ind] <   height.cutoff*lms.max.height){    
      w=0
    }
    else{
      ## computing the sharpness
      w=A$y[lms.ind]-A$y[ind]
      w[which(w<0)]=3*w[which(w<0)]
    }
    ## computing final score
    score[i]=sum(w, na.rm=T)*A$y[lms.ind]
    if(score[i]<0)
      score[i]=0
    if(last.lms<0){
      last.lms=lms[i]
      last.lms.i=i
    } else {
      ##If two lms's are very close only choose one with the better score.
      if(lms[i]-last.lms < MIN.LMS.DIST) {
        if(score[i]>score[last.lms.i]){
          last.lms=lms[i]
          score[last.lms.i]=0
          last.lms.i=i          
        } else {
          score[i]=0
        }
      } else {
        last.lms=lms[i]
        last.lms.i=i
      }
    }
  }
  return(score)  
}



## manipulates the data in such a way that the landmark at matched.lms[i] is moved to matched.lms[i+1] for each i.
register.channel <- function(x, matched.lms){
  s <- m <- shift <- vector()
  lms <- vector()
  for(i in seq(1,length(matched.lms), by=2)) {
    shift <- append(shift, matched.lms[i+1]-matched.lms[i])
    lms <- append(lms, matched.lms[i]) # 
    s <- append(s, sd(x)) # sd of landmark
  }
  r.data <- register.function(x, s, lms, shift)
  return(r.data)
}

#shift the data
register.function <- function(x, s, m, shift) {    
  sum=0
  if(length(m)==1){
    return(x+shift)
  }
  if(length(m)==2){
    sh <- (shift[1]-shift[2])
    x <- x+gau(x, s[1], m[1])*(sh/2)
    x <- x-gau(x, s[2], m[2])*(sh/2)
    return(x+shift[1]-sh/2)
  }
  max.shift=which(abs(shift)==max(abs(shift)))[1]
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
gau <- function(d, s, m) return(2.7182^(-((d-m)^2)/(2*s*s)))


