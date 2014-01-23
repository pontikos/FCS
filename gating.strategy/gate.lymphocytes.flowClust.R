#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
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
if (is.null(opt$in.file) || !file.exists(opt$in.file)) stop("No FCS file specified on command line or file does not exist.")


suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(tools, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
source('~nikolas/bin/FCS/fcs.R')
source('~nikolas/bin/FCS/gating.strategy/gate.lymphocytes.R')

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
    #gate.lymph(d, down.sample=5000, K=4, level=.5, B=500, post.threshold=.99, channels=c('FSC-A','SSC-A','CD4'), usePrior='no', prior=NULL, plot.file=NULL)
    for (fcs.name in names(fcs.files)) {
        fcs.data <- fcs.files[[fcs.name]]
        print(fcs.name)
        if (!is.null(opt$plot.dir)) print(plot.file <- file.path(opt$plot.dir,paste(basename(fcs.name),'.png',sep='')))
        if (!is.null(opt$out.dir))  print(out.file <- file.path(opt$out.dir, paste(basename(fcs.name),'.RData',sep='')))
        if (file.exists(out.file)) {
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


