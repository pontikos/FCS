#!/usr/bin/env Rscript
source('~nikolas/bin/FCS/fcs.R',chdir=T)
suppressPackageStartupMessages(library("optparse"))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(flowClust, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(tools, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(feature, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(KernSmooth, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(RANN, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))


###  make prior from pooled sample file
### If nu is Inf then a Gaussian is used instead of a t distribution.
### Do not estimate lambda as this step is too computationally expensive and lambda is usually close to 1 anyway
### which just results in substracting 1 via the boxcox transform.
### Futhermore lambda makes data look different across samples so priors don't match data as well.
make.prior <- function(fcs.data, K=4, level=.5, B=1000, channels=NULL, plot.file=NULL, nu=4, pool.file=NULL) {
    num.clusters <- 0
    d <- fcs.data
    cat('pooled data dim', dim(d), '\n')
    print(head(d))
    res <- flowClust::flowClust(d,K=K,B=B,level=level,trans=0,lambda=1,nu=nu)
    clusters <- Map(res)
    print(num.clusters <- length(table(clusters)))
    prior <- flowClust2Prior(res,kappa=1)
    if (!is.null(plot.file)) plotFlowClustRes(d, res, plot.file=gsub('.fcs','',plot.file))
    dens <- compute.dens(d, res)
    total.dens <- rowSums(dens)
    # posterior probability of belonging to each cluster
    post <- dens/total.dens 
    post <- compute.post(d, res)
    # assignment to whichever clusters has the highest posterior probability
    clustering <- apply(post,1,which.max)
    print(dput(getEstimates(res)))
    #
    cat('>>initial.weights', res@w, '\n' )
    w <- table(clustering)/length(clustering)
    cat('>>updated.weights', w, '\n' )
    return( list(prior=prior,res=res,data=d) )
}


### gate 
### kappa is the prior weight relative to the sample, see flowClust2Prior
gate <- function(d, fcs.name=NULL, K=4, level=.9, u.cutoff=.5, B=10000, post.threshold=.99, channels=NULL, prior=NULL, plot.file=NULL, out.file=NULL, kappa=1) {
    d.original <- d
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
    save(res, prior, d, fcs.name, file=out.file, compress='xz', compression_level=9)
    clusters <- Map(res)
    cat('Number of clusters', num.clusters <- length(table(clusters)), '\n')
    # checks
    # if number of clusters different than expected BAD
    if (num.clusters != K) stop('Different number of clusters, expected: ',K, ' but found ', num.clusters, '\n' )
    cat('Posterior\n')
    print(summary(res))
    #compute dens in whole sample
    dens <- compute.dens(d, res)
    total.dens <- rowSums(dens)
    # posterior probability of belonging to each cluster
    post <- dens/total.dens 
    # assignment to whichever clusters has the highest posterior probability
    clustering <- apply(post,1,which.max)
    #
    cat('>>initial.weights', res@w, '\n' )
    w <- table(clustering)/length(clustering)
    cat('>>updated.weights', w, '\n' )
    if (!is.null(plot.file)) plotFlowClustRes(d, prior$res, plot.file=gsub('.png','-prior.png',plot.file))
    if (!is.null(plot.file)) plotFlowClustRes(d, res, plot.file=gsub('.png','-posterior.png',plot.file))
    save(res, d, fcs.name, file=out.file, compress='xz', compression_level=9)
}



### MAIN

option_list <- list( 
    make_option(c("-f","--in.file"), help = 'FCS file to parse or list of files in a csv file.'),
    make_option(c("--out.csv"), help = 'Out csv file.'),
    make_option(c('--down.sample'), default=.1, help='A numeric value between 0 and 1 specifying the ratio by which to downsample to achieve flowClust speedup.  The default 0.1, meaning that 10% of the data will be randomly sampled. Can also be a number bigger than 1 in which case it is interpreted as the number of events to randomly sample.'),
    make_option(c('--channels'), default='FSCA,SSCA,CD4', help='Channels on which to do the clustering.  By default CD4 lymphocyte gating works on Forward, Side Scatter and CD4.'),
    make_option(c('-K', '--clusters'), default=5, help='Number of clusters.  The default is 5.'),
    make_option(c('--level'), default=.9, help='A numeric value between 0 and 1 specifying the threshold quantile level used to call a point an outlier.  The default is 0.9, meaning that any point outside the 90% quantile region will be called an outlier.'),
    make_option(c('-B', '--em.iterations'), default=10000, help='Number of EM iterations.'),
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
    #My version of read.FCS dicards events on axis so total count is less that what is in FCS file initially.
    print(total.count <- as.numeric(read.FCSheader(in.file,keyword='$TOT')))
} else if (file.ext=='csv') {
    #csv file
    print(head(fcs.data <- read.csv(in.file)[,channels]))
    print(total.count <- dim(fcs.data)[[1]])
} else if (file.ext=='RData') {
    print(load(in.file))
    print(head(fcs.data <- d$data[,channels]))
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
    prior <- make.prior(fcs.data, K=opt$clusters, plot.file=prior.plot.file, B=opt$em.iterations, level=opt$level, nu=4)
    print( prior$prior )
    print(prior.out.file)
    save(prior, file=prior.out.file)
} else {
    # Cluster using prior
    cat('Prior specified.\n')
    print(load(opt$prior))
    print(fcs.name <- in.file)
    if (!is.null(opt$plot.dir)) print(plot.file <- file.path(opt$plot.dir,paste(basename(fcs.name),'.png',sep='')))
    if (!is.null(opt$out.dir))  print(out.file <- file.path(opt$out.dir, paste(basename(fcs.name),'.RData',sep='')))
    if (!is.null(opt$out.csv)) print(out.csv <- file.path(opt$out.csv, paste(basename(fcs.name),'.csv',sep='')))
    # flowClust results exists already
    if (file.exists(out.file)) {
        cat(out.file, 'exists not running again!\n')
        print(load(out.file))
    } else {
        gate(fcs.data, channels=channels, fcs.name=fcs.name, K=opt$clusters, level=opt$level, plot.file=plot.file, prior=prior, post.threshold=opt$post.threshold, out.file=out.file, kappa=1)
    }
    cat('>',fcs.name,',total.count,', total.count, '\n',sep='')
}


