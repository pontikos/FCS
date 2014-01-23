#!/usr/bin/env Rscript
source('~nikolas/bin/FCS/fcs.R',chdir=T)
suppressPackageStartupMessages(library("optparse"))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(tools, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(flowClust, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(tools, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(feature, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(KernSmooth, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(RANN, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))


## instead of returning density in grid format
## returns 2D density at each (x,y) point
kde2D <- function(d, bw=.1) {
    # compute fast kernel density estimate
    b<-bkde2D(as.matrix(d),bw)
    # this returns you a grid
    # we will use a fast nearest neighbour method
    # to find the closest point in the grid
    grid <- expand.grid(b$x1, b$x2)
    nn <- nn2(grid,d,k=1)
    return(cbind(d, dens=as.numeric(b$fhat)[nn$nn.idx]))
}


#c('yellow','red','purple','green','pink')

# to perform reverse box-cox transformation (multivariate)
rbox <- function(data, lambda) {
    if (length(lambda)>1 || lambda!=0) data <- sign(lambda*data+1)*(sign(lambda*data+1)*(lambda*data+1))^(1/lambda) else data <- exp(data)
    data
}


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
plotflowclust <- function(x, subset, py=3, ellipse=T, show.outliers=T, show.rm=F, include=1:(x@K), main=NULL, grayscale=F, col=(if (grayscale) gray(1/4) else 2:(length(include)+1)), pch=".", cex=0.6, col.outliers=gray(3/4), pch.outliers=".", cex.outliers=cex, col.rm=1, pch.rm=1, cex.rm=0.6, ecol=1, elty=1, level=NULL, u.cutoff=NULL, z.cutoff=NULL, npoints=501, add=F,...) {
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
    if (length(x@lambda)>0){
        if (any(x@lambda!=1))
            lambda<-rep(x@lambda, length.out=x@K)
        else
        #WTF.. why set lambda to 0?
            lambda<-numeric(0)
    }else{
        lambda<-numeric(0);
    }
    cc <- rep(cc, length.out=x@K)
    for (i in include) {
        eigenPair <- eigen(x@sigma[i,subset,subset])
        l1 <- sqrt(eigenPair$values[1]) * sqrt(cc)
        l2 <- sqrt(eigenPair$values[2]) * sqrt(cc)
        angle <- atan(eigenPair$vectors[2,1] / eigenPair$vectors[1,1]) * 180/pi 
       if (length(lambda)>0&any(lambda!=1))
            points(rbox(ellipsePoints(a=l1[i], b=l2[i], alpha=angle, loc=x@mu[i,subset], n=npoints), lambda[i]), type="l", lty=elty[j <- j+1], col=ecol[j])
       else
            points(ellipsePoints(a=l1[i], b=l2[i], alpha=angle, loc=x@mu[i,subset], n=npoints), type="l", lty=elty[j <- j+1], col=ecol[j])
    }  
}


###
###
plot.lymph <- function(d, res, lymph.filter=NULL, plot.file=NULL) {
        png(plot.file)
        par(mfrow=c(2,2))
        for (chan in list(c('SSCA','FSCA'), c('FSCA','CD4'), c('SSCA','CD4'))) {
            smoothScatter(d[,chan])
            lymph <- d[lymph.filter,chan]
            points(lymph, pch=20, col='pink')
            lines(lymph[c(chull(lymph),chull(lymph)[1]),], col='black', lwd=2)
            #lymph.gate <- c(chull(lymph), chull(lymph)[1])
            #lymph.gate <- c(lymph.gate, lymph.gate[1])
            #chuld <- lapply(lymph,"[",chull(lymph))
            #polygon(chuld,lty=2,border="black")
            #polygon(spline.poly(as.matrix(as.data.frame(chuld)),100),border="black",lwd=2)
            plotflowclust(res, subset=chan)
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
        dev.off()
}


#slow as fuck
#compute.density <- function(fcs.data, channels, kernel_mult=5.0, apprx_mult=1.5, med_samples=2000) {
    #if (missing(channels)) channels <- colnames(getChannels(fcs.data))
    #SPADE.density(getChannels(fcs.data, channels), kernel_mult, apprx_mult, med_samples)
#}


density.filter <- function(x, channels, quant='25%') {
    x <- getChannels(x, channels)
    x <- apply(x, 2, scale)
    dens1 <- kde2D(x[,c('FSCA','SSCA')])[,3]
    dens2 <- kde2D(x[,c('SSCA','CD4')])[,3]
    dens3 <- kde2D(x[,c('FSCA','CD4')])[,3]
    dens1 <- dens1/sum(dens1)
    dens2 <- dens2/sum(dens2)
    dens3 <- dens3/sum(dens3)
    q1 <- quantile(dens1)[[quant]]
    q2 <- quantile(dens2)[[quant]]
    q3 <- quantile(dens3)[[quant]]
    f <- dens1 > q1 & dens2 > q2 & dens3 > q3
    return(f)
}


pool.data2 <- function(d, down.sample, channels) 
    #pool data
    #d <- getChannels( do.call('rbind', lapply(d, function(x) getChannels(x[sample(1:nrow(x), down.sample),]))), channels )
    #d <- lapply(d, function(x) getChannels(x[sample(1:nrow(x), down.sample),], channels))
    #way too slow!
    #dens <- lapply(d, compute.density, channels=channels)
    getChannels( do.call('rbind', 
            lapply(d, function(x) {
                #calculate pairwise density
                table(f <- density.filter(x, channels, '25%'))
                x <- getChannels(x[f,], channels)
                return(x)
                }) ), channels)

pool.data <- function(d, down.sample, channels) return(getChannels( do.call('rbind', lapply(d, function(x) getChannels(x[sample(1:nrow(x), down.sample),]))), channels ))

### The CD4 lymphocyte cluster is the one with the CD4 mean intensity which is the closest to 2.5
### (this is not always true depending on the experiment)
identify.lymph.cluster <- function(d, clusters) {
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
    return(lymph.cluster <- which.min(abs(cd4.means-2.5)))
}


### density at each point for each component
compute.dens <- function(d, res) sapply(1:res@K, function(i) res@w[i]*flowClust::dmvt(d, mu=res@mu[i,], sigma=res@sigma[i,,], nu=res@nu, lambda=res@lambda)$value)

###  make prior by pooling samples
### 
make.lymph.prior <- function(fcs.files, down.sample=100, K=4, level=.5, B=500, channels=c('FSC-A','SSC-A','CD4'), plot.file=NULL) {
    num.clusters <- 0
    #resample until we have the desired number of clusters
    while (num.clusters!=K) {
        d <- pool.data(fcs.files, down.sample, channels)
        cat('pooled data dim', dim(d), '\n')
        print(head(d))
        #scaled.d <- apply(d,2,scale)
        res <- flowClust::flowClust(d,K=K,B=B,level=level)
        clusters <- Map(res)
        print(num.clusters <- length(table(clusters)))
    }
    dens <- compute.dens(d, res)
    total.dens <- rowSums(dens)
    # posterior probability of belonging to each cluster
    post <- dens/total.dens 
    # assignment to whichever clusters has the highest posterior probability
    clustering <- apply(post,1,which.max)
    lymph.cluster <- identify.lymph.cluster(d, clustering)
    print( length( lymph.filter <- which((clustering==lymph.cluster) & (total.dens>quantile(total.dens)[['25%']])) ) )
    #print(dput(getEstimates(res)))
    #
    cat('>>initial.weights', res@w, '\n' )
    w <- table(clustering)/length(clustering)
    cat('>>updated.weights', w, '\n' )
    cat('>>subset.lymph.cluster.pct',100*res@w[lymph.cluster],'\n')
    cat('>>all.lymph.cluster.pct', 100*w[lymph.cluster], '\n') 
    cat('>>gated.lymph.cluster.weight', 100*length(lymph.filter)/length(clustering), '\n') 
    #
    if (!is.null(plot.file)) plot.lymph(d, res=res, lymph.filter=lymph.filter, plot.file=plot.file)
    return( list(prior=flowClust::flowClust2Prior(res,kappa=1),res=res,data=d) )
}



### gate lymphocytes
### assumptions:
### 1) we can distinguish 4 clusters using flowClust (preferably with prior) in FSC-A x SSC-A x CD4
### 2) lymphocyte cluster is the one with the CD4 mean intensity which is the closest to 2.5
###   (this is not always true depending on the experiment and needs to be an argument to this function)
gate.lymph <- function(d, fcs.name=NULL, down.sample=.1, K=4, level=.5, B=500, post.threshold=.99, channels=c('FSC-A','SSC-A','CD4'), prior=NULL, plot.file=NULL, out.file=NULL) {
    d.original <- d
    if (is(d,'flowFrame')) d<-getChannels(d, channels)
    #subset
    if (is.null(prior)) usePrior <- 'no' else usePrior <- 'yes'
    num.clusters <- 0
    #resample until we have the desired number of clusters
    while (num.clusters!=K) {
        if (down.sample > 1) d.sub <- sample(1:nrow(d),round(down.sample)) else d.sub <- sample(1:nrow(d),round(nrow(d)*down.sample))
        res <- flowClust::flowClust(d[d.sub,],K=K,B=B,usePrior=usePrior,prior=prior,level=level)
        clusters <- Map(res)
        print(num.clusters <- length(table(clusters)))
    }
    #compute dens in whole sample
    dens <- compute.dens(d, res)
    total.dens <- rowSums(dens)
    # posterior probability of belonging to each cluster
    post <- dens/total.dens 
    # assignment to whichever clusters has the highest posterior probability
    clustering <- apply(post,1,which.max)
    lymph.cluster <- identify.lymph.cluster(d, clustering)
    print( length( lymph.filter <- which((clustering==lymph.cluster) & (total.dens>quantile(total.dens)[['25%']])) ) )
    #
    cat('>>initial.weights', res@w, '\n' )
    w <- table(clustering)/length(clustering)
    cat('>>updated.weights', w, '\n' )
    cat('>>subset.lymph.cluster.pct',100*res@w[lymph.cluster],'\n')
    cat('>>all.lymph.cluster.pct', 100*w[lymph.cluster], '\n') 
    cat('>>gated.lymph.cluster.weight', 100*length(lymph.filter)/length(clustering), '\n') 
    #
    if (!is.null(plot.file)) plot.lymph(d, res=res, lymph.filter=lymph.filter, plot.file=plot.file)
    if (is.null(out.file)) return(lymph.filter)
    ext <- file_ext(out.file)
    if (ext == 'idx') write.csv(lymph.filter, file=out.file, quote=FALSE, row.names=FALSE, col.names=FALSE)
    else if (ext == 'RData')
        #lymph <- d.original[lymph.filter,]
        #save(lymph, file=out.file, compress='xz', compression_level=9)
        save(res, lymph.filter, d, fcs.name, file=out.file, compress='xz', compression_level=9)
    return(lymph.filter)
}


### MAIN

option_list <- list( 
    make_option(c("-f","--in.file"), help = 'FCS file to parse or list of files in a csv file.'),
    make_option(c('--down.sample'), default=.1, help='A numeric value between 0 and 1 specifying the ratio by which to downsample to achieve flowClust speedup.  The default 0.1, meaning that 10% of the data will be randomly sampled. Can also be a number bigger than 1 in which case it is interpreted as the number of events to randomly sample.'),
    make_option(c('--channels'), default='FSCA,SSCA,CD4', help='Channels on which to do the clustering.  By default CD4 lymphocyte gating works on Forward, Side Scatter and CD4.'),
    make_option(c('-K', '--clusters'), default=5, help='Number of clusters.  The default is 5.'),
    make_option(c('--level'), default=.5, help='A numeric value between 0 and 1 specifying the threshold quantile level used to call a point an outlier.  The default is 0.9, meaning that any point outside the 90% quantile region will be called an outlier.'),
    make_option(c('-B', '--em.iterations'), default=500, help='Number of EM iterations.'),
    make_option(c('--post.threshold'), default=.99, help='The posterior cutoff.'),
    make_option(c('--plot.dir'), default=NULL, help='The directory to which to send the plots.'),
    make_option(c('--out.dir'), default=NULL, help = 'Output file which may be a subset FCS file, indexes or the flowClust result object.'),
    make_option(c('--prior'), default=NULL, help='Prior file to use.')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

# do checks first
if (is.null(opt$in.file) || !file.exists(opt$in.file))
    stop("No FCS file specified on command line or file does not exist.")
if (!is.null(opt$channels)) {
    channels <- unlist(strsplit(opt$channels, ","))
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
if (tolower(file_ext(opt$in.file))=='fcs') {
    # expect a single fcs files
    fcs.files <- list(suppressWarnings(read.FCS(opt$in.file, channels=channels)))
    names(fcs.files) <- opt$in.file
} else if (tolower(file_ext(opt$in.file))=='csv') {
    # expect a list of fcs files in the csv file
    fcs.names <- as.character(read.csv(opt$in.file)[,1])
    fcs.files <- lapply(fcs.names, read.FCS, channels=channels)
    names(fcs.files) <- fcs.names
}

# Two modes of operation: either compute prior or run on all files with specified prior
# If no prior file specified
# make a prior from a pool of fcs files
# otherwise load prior from specified file
if (is.null(opt$prior)) {
    cat('No prior specified, will estimate prior from pooled data.\n')
    print( prior.plot.file )
    prior <- make.lymph.prior( fcs.files, K=opt$clusters, channels=channels, plot.file=prior.plot.file, B=opt$em.iterations, level=opt$level, down.sample=opt$down.sample )
    print( prior$prior )
    save(prior, file=prior.out.file)
} else {
    cat('Prior specified.\n')
    print(load(opt$prior))
    for (fcs.name in names(fcs.files)) {
        fcs.data <- fcs.files[[fcs.name]]
        print(fcs.name)
        if (!is.null(opt$plot.dir)) print(plot.file <- file.path(opt$plot.dir,paste(basename(fcs.name),'.png',sep='')))
        if (!is.null(opt$out.dir))  print(out.file <- file.path(opt$out.dir, paste(basename(fcs.name),'.RData',sep='')))
        if (file.exists(out.file)) {
            cat(out.file, 'exists not running again!\n')
            print(load(out.file))
            lymph <- getChannels(fcs.data,channels)[lymph.filter,]
            write.csv(lymph,file=file.path(opt$out.dir, paste(basename(fcs.name),'.csv',sep='')),quote=FALSE,row.names=F)
        } else {
            lymph.filter <- gate.lymph( fcs.data, fcs.name=fcs.name, channels=channels, K=opt$clusters, level=opt$level, plot.file=plot.file, prior=prior$prior, post.threshold=opt$post.threshold, out.file=out.file, down.sample=opt$down.sample )
        }
        total.count <- dim(fcs.data)[[1]]
        lymph.count <- length(lymph.filter)
        lymph.pct <- 100*lymph.count/total.count
        cat('>',fcs.name,',total.count,', total.count, '\n',sep='')
        cat('>',fcs.name,',lymph.count,', lymph.count, '\n',sep='')
        cat('>',fcs.name,',lymph.pct,', lymph.pct, '\n',sep='')
    }
}


