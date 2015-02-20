#!/usr/bin/env Rscript
source('~nikolas/bin/FCS/fcs.R',chdir=T)
source('~nikolas/bin/FCS/normalise-functions.R',chdir=T)
suppressPackageStartupMessages(library("optparse"))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(flowClust, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(tools, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(feature, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(KernSmooth, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(ellipse, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
#suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(RANN, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))


#c('yellow','red','purple','green','pink')

## flowClust code
ellipsePoints <- function(a,b, alpha = 0, loc = c(0,0), n = 501)
{
    ## Purpose: ellipse points,radially equispaced, given geometric par.s
    ## -------------------------------------------------------------------------
    ## Arguments: a, b : length of half axes in (x,y) direction
    ##            alpha: angle (in degrees) for rotation
    ##            loc  : center of ellipse
    ##            n    : number of points
    ## -------------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 19 Mar 2002, 16:26
    ## modified by Kenneth to get rid of the precision problem met when there's a large difference in the length of the two axes
    small <- 0
    if (a/b > 1000) {
        ratio <- a/b
        b <- a
        if (round(alpha)==0) small <- 2 else small <- 1
    }
    B <- min(a,b)
    A <- max(a,b)
    ## B <= A
    d2 <- (A-B)*(A+B)                   #= A^2 - B^2
    phi <- 2*pi*seq(0,1, len = n)
    sp <- sin(phi)
    cp <- cos(phi)
    r <- a*b / sqrt(B^2 + d2 * sp^2)
    xy <- r * cbind(cp, sp)
    ## xy are the ellipse points for alpha = 0 and loc = (0,0)
    al <- alpha * pi/180
    ca <- cos(al)
    sa <- sin(al)
    xy.new <- xy %*% rbind(c(ca, sa), c(-sa, ca))
    if (small==2) xy.new[,2]=xy.new[,2]/ratio
    if (small==1) xy.new[,1]=xy.new[,1]/ratio
    xy.new + cbind(rep(loc[1],n), rep(loc[2],n))
}


## flowClust code
plotflowclust <- function(x, subset, py=3, ellipse=T, show.outliers=T, show.rm=F, include=1:(x@K), main=NULL, grayscale=F, col=(if (grayscale) gray(1/4) else 2:(length(include)+1)), pch=".", cex=0.6, col.outliers=gray(3/4), pch.outliers=".", cex.outliers=cex, col.rm=1, pch.rm=1, cex.rm=0.6, ecol=1, elty=1, level=NULL, u.cutoff=NULL, z.cutoff=NULL, npoints=501, prior=NULL, ...) {
    if (!is.numeric(subset)) subset <- match(subset, x@varNames)
    ecol <- matrix(ecol, length(include))
    elty <- matrix(elty, length(include)) 
    if (all(x@nu!=Inf)) {
        if (x@ruleOutliers[1]==0) {     # 0 means quantile
            if(all(is.na(x@prior))){
                cc <- py * qf(x@ruleOutliers[2], py, x@nu)
            }else{
                cc <- py * qf(x@ruleOutliers[2], py, x@nu)
            }
        }  else {     # 1 means u.cutoff
            if(all(is.na(x@prior))){
                cc <- ((x@nu+py)/x@ruleOutliers[2] - x@nu)    
            }else{
                cc <- ((x@nu+py)/x@ruleOutliers[2] -x@nu)    
            }
        }
    }  else cc <- qchisq(x@ruleOutliers[2], py)
    j <- 0
    lambda <- x@lambda
    cc <- rep(cc, length.out=x@K)
    for (i in include) {
       eigenPair <- eigen(x@sigma[i,subset,subset])
       l1 <- sqrt(eigenPair$values[1]) * sqrt(cc)
       l2 <- sqrt(eigenPair$values[2]) * sqrt(cc)
       angle <- atan(eigenPair$vectors[2,1] / eigenPair$vectors[1,1]) * 180/pi 
       points(rbox(ellipsePoints(a=l1[i], b=l2[i], alpha=angle, loc=x@mu[i,subset], n=npoints), lambda), type="l", lty=1, col=i)
       #prior
       print(prior$lambda)
       if (!is.null(prior)) {
		#points(rbox(t(as.matrix(prior$Mu0[i,subset])),prior$lambda),pch=20,col=i)
		#lines(ellipse(rbox(prior$Lambda0[i,subset,subset]/(prior$nu0[i]-ncol(prior$Mu0)-1),prior$lambda),centre=rbox(prior$Mu0[i,subset],prior$lambda)),col=i,lwd=1,lty=2)
       eigenPair <- eigen(prior$Lambda0[i,subset,subset]/(prior$nu0[i]-ncol(prior$Mu0)-1))
       l1 <- sqrt(eigenPair$values[1])
       l2 <- sqrt(eigenPair$values[2])
       angle <- atan(eigenPair$vectors[2,1] / eigenPair$vectors[1,1]) * 180/pi 
       points(rbox(ellipsePoints( a=l1[i], b=l2[i], alpha=angle, loc=prior$Mu0[i,subset], n=npoints), lambda=prior$lambda), type='l',col=i,lty=2)
       }
    }  
}


###
###
plot.lymph <- function(d, res=NULL, lymph.filter=NULL, plot.file=NULL, prior=NULL, radius=.9) {
        cat('>>',plot.file,'\n')
        if (!is.null(plot.file)) png(plot.file)
        par(mfrow=c(2,2))
        for (chan in list(c('SSCA','FSCA'), c('CD4','FSCA'), c('SSCA','CD4'))) {
            smoothScatter(d[,chan])
            lymph <- d[lymph.filter,chan]
            #points(lymph, pch=20, col='pink')
            chull.lymph <- c(chull(lymph),chull(lymph))
            lines(lymph[chull.lymph,], col='pink', lwd=2, lty=2)
            #also approximate with an ellipse
            lines(ellipse::ellipse(cov(lymph),centre=colMeans(lymph)),col='pink',lwd=2,lty=1)
            #lymph.gate <- c(chull(lymph), chull(lymph)[1])
            #lymph.gate <- c(lymph.gate, lymph.gate[1])
            #chuld <- lapply(lymph,"[",chull(lymph))
            #polygon(chuld,lty=2,border="black")
            #polygon(spline.poly(as.matrix(as.data.frame(chuld)),100),border="black",lwd=2)
            if (!is.null(res)) plotflowclust(res, subset=chan, ecol=1:length(res@K), prior=prior)
        }
        #
        image(1:3,1:3,matrix(data=0,nrow=3,ncol=3), axes=FALSE, col='white', ylab='', xlab='')
        text(2,3, 'count', font=2)
        text(3,3, '%', font=2)
        text(1,1, 'lymph', font=2)
        text(2,1, round(length(lymph.filter)) )
        text(3,1, round(100*length(lymph.filter)/dim(d)[[1]]) )
        text(1,2, 'total', font=2)
        text(2,2, dim(d)[[1]])
        text(3,2, 100)
        if (!is.null(plot.file)) dev.off()
}




### The CD4 lymphocyte cluster is the one with the CD4 mean intensity which is the closest to 2.5
### (this is not always true depending on the experiment)
### As long as we can identify the CD4 lymph cluster in the pooled data then we are good.
identify.lymph.cluster <- function(d, clusters, cd4.mfi=2.5) {
    #the cluster of lymphocytes is the one with the CD4 MFI which is the closest to 2.5
    cd4.means <- tapply(d[,'CD4'], clusters, mean)
    print(i <- sort(cd4.means))
    #cluster colours are chosen by increasing CD4 MFI
    clusters <- factor(clusters, levels=names(i))
    K <- length(i)
    # lymph cluster
    #make sure we keep higher density events which we can classify with certainty as belonging to lymph cluster
    #this assumes that less than 75% of the sample are lymphocytes which seems like a reasonable assumption
    #lymph.filter <- which( (posteriors > post.threshold) & (total.dens>median(total.dens)) )
    return(lymph.cluster <- which.min(abs(cd4.means-cd4.mfi)))
}

###
flowClust2Prior <- function (x, kappa, Nt = NULL) {
    if (is.null(Nt)) Nt <- nrow(x@z)
    p <- ncol(x@mu)
    K <- x@K
    nu0 <- Ng <- x@w * Nt
    if (all((nu0*kappa-p-1)>0)) {
        Lambda0 <- x@sigma
        for (i in 1:K) 
            Lambda0[i, , ] <- Lambda0[i, , ] * (kappa*nu0[i]-p-1)
    }
    else {
        stop("Can't proceed. Prior nu0 is negative for cluster(s) ", paste(which((nu0-p-1)>0), collapse=","), "\n(p-1) = ", p-1, ": Try increasing kappa")
    }
    Omega0 <- array(0, c(K, p, p))
    for (i in 1:K) {
        Omega0[i, , ] <- diag(1, p)
        if (p == 1) {
            dS <- x@sigma[i, , ]
            dO <- Omega0[i, , ]
        }
        else {
            dS <- det(x@sigma[i, , ])
            dO <- det(Omega0[i, , ])
        }
        k <- (dO/dS)^(1/p)
        Omega0[i, , ] <- Omega0[i, , ] * k
        Omega0[i, , ] <- solve(Omega0[i, , ]*Ng[i]*kappa)
    }
    nu0 <- nu0 * kappa
    Mu0 <- x@mu
    lambda <- x@lambda
    w0 <- x@w * Nt
    prior <- list(Mu0 = Mu0, Lambda0 = Lambda0, Omega0 = Omega0, w0 = w0, nu0 = nu0, nu = x@nu, lambda = x@lambda, K = K)
    class(prior) <- "flowClustPrior"
    attr(prior, "lambda") <- x@lambda
    prior
}


###
removeComponent <- function(res, i) {
    res@K <- res@K-1
    res@w <- res@w[-i]
    res@mu <- res@mu[-i,]
    res@sigma <- res@sigma[-i,,]
    res
}


###  make prior from pooled sample file
### If nu is Inf then a Gaussian is used instead of a t distribution.
### Do not estimate lambda as this step is too computationally expensive and lambda is usually close to 1 anyway
### which just results in substracting 1 via the boxcox transform.
### Futhermore lambda makes data look different across samples so priors don't match data as well.
make.lymph.prior <- function(fcs.data, K=4, level=.5, B=500, channels=c('FSCA','SSCA','CD4'), plot.file=NULL, cd4.mfi=2.5, nu=4, pool.file=NULL) {
    num.clusters <- 0
    d <- fcs.data[,channels]
    cat('pooled data dim', dim(d), '\n')
    print(head(d))
    print(dim(d))
    #scaled.d <- apply(d,2,scale)
    res <- flowClust::flowClust(d,K=K,B=B,level=level,trans=0,lambda=1,nu=nu)
    clusters <- Map(res)
    print(num.clusters <- length(table(clusters)))
    #}
    prior <- flowClust2Prior(res,kappa=1)
    #if (!is.null(plot.file)) plotPrior(d, prior=prior, plot.file=gsub('.fcs','',plot.file))
    if (!is.null(plot.file)) plotClustRes(d, res, plot.file=gsub('.fcs','',plot.file))
    dens <- compute.dens(d, res)
    total.dens <- rowSums(dens)
    # posterior probability of belonging to each cluster
    post <- dens/total.dens 
    post <- compute.post(d, res)
    # assignment to whichever clusters has the highest posterior probability
    clustering <- apply(post,1,which.max)
    cat('>>lymph.pop', lymph.cluster <- identify.lymph.cluster(d, clustering, cd4.mfi), '\n')
    print( length( lymph.filter <- which((clustering==lymph.cluster) & (total.dens>quantile(total.dens)[['25%']])) ) )
    if (!is.null(plot.file)) plot.lymph(d, res=res, lymph.filter=lymph.filter, plot.file=plot.file)
    print(dput(getEstimates(res)))
    #
    cat('>>initial.weights', res@w, '\n' )
    w <- table(clustering)/length(clustering)
    cat('>>updated.weights', w, '\n' )
    cat('>>subset.lymph.cluster.pct',100*res@w[lymph.cluster],'\n')
    cat('>>all.lymph.cluster.pct', 100*w[lymph.cluster], '\n') 
    cat('>>gated.lymph.cluster.weight', 100*length(lymph.filter)/length(clustering), '\n') 
    # estimates of lymph cluster
    #getEstimates(res)[[
    return( list(prior=prior,res=res,data=d,lymph.cluster=lymph.cluster,lymph.filter=lymph.filter) )
}


###
gate.singlets <- function(X,lymph.filter,k=3) {
    #filter on SSCW
    sscw <- X[lymph.filter,'SSCW']
    med <- median(sscw)
    sscw.filter <- which(med-k*median(abs(med-sscw)) < sscw & sscw < med+k*median(abs(med-sscw)))
    return(sscw.filter)
    #FSCH,SSCH,FSCW,SSCW
    #m <- covMcd(X)
    #m$mah < max(m$mah)
    #mahalanobis(X)
} 


### gate lymphocytes
### assumptions:
### 1) we can distinguish 4-5 clusters using flowClust (preferably with prior) in FSC-A x SSC-A x CD4
### 2) lymphocyte cluster is the one with the CD4 mean intensity which is the closest to 2.5
###   (this is not always true depending on the experiment and needs to be an argument to this function)
### kappa is the prior weight relative to the sample, see flowClust2Prior
gate.lymph <- function(d, fcs.name=NULL, K=4, level=.9, u.cutoff=.5, B=5000, post.threshold=.99, channels=c('FSCA','SSCA','CD4'), prior=NULL, plot.file=NULL, out.file=NULL, lymph.cluster=NULL, kappa=1) {
    d.original <- d
    #prior <- prior$prior
    if (is.null(prior)) usePrior <- 'no' else usePrior <- 'yes'
    prior$prior <- flowClust2Prior(prior$res,kappa=kappa)
    if (is(d,'flowFrame')) d<-getChannels(d, channels)
    else d <- d[,channels]
    #subset
    num.clusters <- 0
    cat('Prior\n')
    print(prior$prior)
    print(prior$res)
    #
    res <- flowClust::flowClust(d,K=K,B=B,usePrior=usePrior,prior=prior$prior,level=level,u.cutoff=u.cutoff,trans=0,lambda=1,control=list(B.lambda=0))
    result <- list(res=res, prior=prior, d=d.original, fcs.name=fcs.name)
    save(result, file=out.file, compress='xz', compression_level=9)
    clusters <- Map(res)
    cat('Number of clusters', num.clusters <- length(table(clusters)), '\n')
    # checks
    # if number of clusters different than expected BAD
    if (num.clusters != K) stop('Different number of clusters, expected: ',K, ' but found ', num.clusters, '\n' )
    #distance of lymph cluster MFI from original estimate
    #likely to be an issue if far from original estimate
    print(lymph.mfi.diff <- sqrt(sum((res@mu[lymph.cluster,]-prior$res@mu[lymph.cluster,])**2)))
    if (lymph.mfi.diff > 2) warning('Possible outlier due to MFI diff', lymph.mfi.diff, '\n')
    #most sensitive to difference in CD4
    print(lymph.cd4.mfi.diff <- abs(res@mu[lymph.cluster,which('CD4'==res@varNames)]-prior$res@mu[lymph.cluster,which('CD4'==prior$res@varNames)]))
    if (lymph.cd4.mfi.diff > 1) warning('Possible outlier due to CD4 MFI diff ', lymph.cd4.mfi.diff, '\n')
    cat('Posterior\n')
    print(summary(res))
    #compute dens in whole sample
    dens <- compute.dens(d, res)
    total.dens <- rowSums(dens)
    # posterior probability of belonging to each cluster
    post <- dens/total.dens 
    # assignment to whichever clusters has the highest posterior probability
    clustering <- apply(post,1,which.max)
    #lymph.cluster <- identify.lymph.cluster(d, clustering, cd4.mfi)
    # further filtering by total density so that we exclude low density points
    print( length( lymph.filter <- which((clustering==lymph.cluster) & (total.dens>quantile(total.dens,probs=seq(0,1,.05))[['5%']])) ) )
    #
    cat('>>initial.weights', res@w, '\n' )
    w <- table(clustering)/length(clustering)
    cat('>>updated.weights', w, '\n' )
    cat('>>subset.lymph.cluster.pct',100*res@w[lymph.cluster],'\n')
    cat('>>all.lymph.cluster.pct', 100*w[lymph.cluster], '\n') 
    cat('>>gated.lymph.cluster.weight', 100*length(lymph.filter)/length(clustering), '\n') 
    #
    if (!is.null(plot.file)) plot.lymph(d, res=res, lymph.filter=lymph.filter, plot.file=plot.file)
    #if (!is.null(plot.file)) plotPrior(d, prior, plot.file=gsub('.png','-prior.png',plot.file))
    if (!is.null(plot.file)) plotClustRes(d, prior$res, outliers=TRUE, plot.file=gsub('.png','-prior.png',plot.file))
    if (!is.null(plot.file)) plotClustRes(d, res, plot.file=gsub('.png','-posterior.png',plot.file))
    if (is.null(out.file)) return(lymph.filter)
    result$lymph.filter <- lymph.filter
    return(result)
}


### MAIN

option_list <- list( 
    make_option(c("-f","--in.file"), help = 'FCS file to parse or list of files in a csv file.'),
    make_option(c('--down.sample'), default=.1, help='A numeric value between 0 and 1 specifying the ratio by which to downsample to achieve flowClust speedup.  The default 0.1, meaning that 10% of the data will be randomly sampled. Can also be a number bigger than 1 in which case it is interpreted as the number of events to randomly sample.'),
    make_option(c('--channels'), default='FSCA,SSCA,CD4,SSCW,SSCH', help='Channels on which to do the clustering.  By default CD4 lymphocyte gating works on Forward, Side Scatter and CD4.'),
    make_option(c('-K', '--clusters'), default=5, help='Number of clusters.  The default is 5.'),
    make_option(c('--level'), default=.9, help='A numeric value between 0 and 1 specifying the threshold quantile level used to call a point an outlier.  The default is 0.9, meaning that any point outside the 90% quantile region will be called an outlier.'),
    make_option(c('-B', '--em.iterations'), default=500, help='Number of EM iterations.'),
    make_option(c('--post.threshold'), default=.99, help='The posterior cutoff.'),
    make_option(c('--plot.dir'), default=NULL, help='The directory to which to send the plots.'),
    make_option(c('--out.dir'), default=NULL, help = 'Output file which may be a subset FCS file, indexes or the flowClust result object.'),
    make_option(c('--prior'), default=NULL, help='Prior file to use.'),
    make_option(c('--pool.file'), default=NULL, help='Pool file.'),
    make_option(c('--cd4.mfi'), default=2.5, help='The expected CD4 MFI of the CD4+ lymphocyte cluster.')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

# do checks first
if (is.null(opt$in.file) || !file.exists(opt$in.file))
    stop("No FCS file specified on command line or file does not exist.")
if (!is.null(opt$channels)) {
    print(channels <- unlist(strsplit(opt$channels, ",")))
} else  {
    channels <- NULL
}
if (!is.null(opt$plot.dir)) {
    dir.create(opt$plot.dir, recursive=TRUE, showWarnings=FALSE)
    prior.plot.file <- file.path(opt$plot.dir,'prior.fcs.png')
} else {
    prior.plot.file <- NULL
    plot.file <- NULL
}
if (!is.null(opt$out.dir)) {
    dir.create(opt$out.dir, recursive=TRUE, showWarnings=FALSE)
    prior.out.file <- file.path(opt$out.dir,'prior.RData')
} else {
    prior.out.file <- NULL
    out.file <- NULL
}

## READ DATA
# takes a single fcs, csv or RData file
cat('Reading', in.file <- opt$in.file, '\n')
cat('File extension', file.ext <- tolower(file_ext(opt$in.file)), '\n')
if (file.ext=='fcs') {
    #fcs file
    print(head(fcs.data <- getChannels(read.FCS(in.file, channels=channels),channels)))
    #My version of read.FCS discards events on axis so total count is less that what is in FCS file initially.
    cat('total.count', total.count <- as.numeric(read.FCSheader(in.file,keyword='$TOT')), '\n')
} else if (file.ext=='csv') {
    #csv file
    print(head(fcs.data <- read.csv(in.file)[,channels]))
    print(total.count <- dim(fcs.data)[[1]])
} else if (file.ext=='rdata') {
    #rdata file
    if (exists('x')) rm(x)
    print(load(in.file))
    if (exists('x')) fcs.data <- x
    print(head(fcs.data <- fcs.data[,channels]))
    print(total.count <- dim(fcs.data)[[1]])
} else {
    stop('Unsupported file extension', file.ext,'\n')
}

# Two modes of operation: either compute prior or run on all files with specified prior
# If no prior file specified
# make a prior from pooled fcs file specified in in.file
# otherwise load prior from specified file and write out.file
if (is.null(opt$prior)) {
    # Compute prior from pooled data
    cat('No prior specified, will estimate prior from pooled data.\n')
    print( prior.plot.file )
    prior <- make.lymph.prior(fcs.data, K=opt$clusters, plot.file=prior.plot.file, B=opt$em.iterations, level=opt$level, cd4.mfi=opt$cd4.mfi, nu=4)
    print( prior$prior )
    print(prior.out.file)
    save(prior, file=prior.out.file)
} else {
    # Cluster using prior
    cat('Prior specified.\n')
    print(load(opt$prior))
    print(fcs.name <- in.file)
    if (!is.null(opt$plot.dir)) print(plot.file <- file.path(opt$plot.dir,paste(basename(file_path_sans_ext(fcs.name)),'.png',sep='')))
    if (!is.null(opt$out.dir))  print(out.file <- file.path(opt$out.dir, paste(basename(file_path_sans_ext(fcs.name)),'.RData',sep='')))
    # flowClust results exists already
    if (file.exists(out.file)) {
        cat(out.file, 'exists not running again!\n')
        print(load(out.file))
    } else {
        result <- gate.lymph(fcs.data,
                                   fcs.name=fcs.name,
                                   K=opt$clusters,
                                   level=opt$level,
                                   plot.file=plot.file,
                                   prior=prior,
                                   post.threshold=opt$post.threshold,
                                   out.file=out.file,
                                   lymph.cluster=prior$lymph.cluster,
                                   kappa=1)
        save(result, file=out.file, compress='xz', compression_level=9)

    }
    fcs.name <- result$fcs.name
    lymph.filter <- result$lymph.filter
    if ('SSCW' %in% channels) singlet.filter <- gate.singlets(fcs.data, lymph.filter, 4)
    result$singlet.filter <- singlet.filter
    lymph.count <- length(singlet.filter)
    lymph.pct <- 100*lymph.count/total.count
    cat('>',fcs.name,',total.count,', total.count, '\n',sep='')
    cat('>',fcs.name,',lymph.count,', lymph.count, '\n',sep='')
    cat('>',fcs.name,',lymph.pct,', lymph.pct, '\n',sep='')
    print(range(lymph.filter))
    #write lymphocytes to csv file
    cd4.lymphocytes <- fcs.data[lymph.filter,channels][singlet.filter,]
    print(head(cd4.lymphocytes))
    for (i in channels) {
        print(i)
        swp <- sliding.window.peaks(density(cd4.lymphocytes[,i]),span=40)
        print(swp <- swp[which(swp$y > .05*max(swp$y)),])
        cat(i, 'peak modality', nrow(swp),'\n')
        if (nrow(swp)>2) cat('Peak modality > 2', i, '\n')
    }
    save(result, file=out.file, compress='xz', compression_level=9)
}


