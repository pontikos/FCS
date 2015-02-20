#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("scales"))
suppressPackageStartupMessages(library("tools"))
suppressPackageStartupMessages(library("mclust"))
source('~nikolas/bin/FCS/fcs.R')
source('~nikolas/bin/FCS/normalise-functions.R')

option_list <- list( 
    make_option(c("-f","--in.file"), help = 'list of files in a csv file.'),
    make_option(c('--channels'), default='FSCA,SSCA,CD4', help='CD4 lymphocyte gating works on Forward, Side Scatter and CD4.'),
    make_option(c('--plot.file'), default=NULL, help='File to which to send the plot.'),
    make_option(c('--span'), default=NULL, help='Span of sliding window for identifying peaks.'),
    make_option(c('-G', '--groups'), default=NULL, help='Max number of clusters to consider in mclust.')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

if (!is.null(opt$channels)) {
    channels <- unlist(strsplit(opt$channels, ","))
} else  {
    channels <- NULL
}

fcs.files <- read.csv(opt$in.file,header=FALSE)[,1]

dens <- list()

for (f in fcs.files) {
    print(f)
    d <- read.flow.file(f)
    peaks.fun <- function(n) {
        x <- d[,n]
        quant <- quantile(x,probs=seq(0,1,.025))
        #x <- x[which(quant[['2.5%']] < x & x < quant[['97.5%']])]
        dens <- normalised.density(x)
        cat('channel',n,'bw',dens$bw,'\n')
        print(head(d <- data.frame(x=dens$x, y=dens$y, peak=FALSE)))
        if (!is.null(opt$span)) {
            print(span <- as.numeric(opt$span))
            #print(span <- round(1/dens$bw))
            #sliding window peaks
            p <- sliding.window.peaks(d, span=span)
            p <- p[p$score>0,]
            d[p$ind,'peak'] <- TRUE
        }
        else if (!is.null(opt$groups)) {
            groups <- as.numeric(opt$groups)
            #mixture model peaks
            m <- Mclust(x,G=1:groups)
            ind <- sapply( m$parameters$mean, function(mu) which.min(abs(dens$x-mu)) )
            d[ind,'peak'] <- TRUE
        }
        return(d)
    }
    dens[[f]] <- sapply(channels, peaks.fun)
}

#save()
print(length(dens))


#ext <- file_ext(opt$plot.file)
#do.call('ext', list(opt$plot.file))
#par(mfrow=length(channels))
pdf(opt$plot.file)
for (n in channels) {
    x <- sapply(1:length(dens), function(i) dens[[i]][,n]$x)
    y <- sapply(1:length(dens), function(i) dens[[i]][,n]$y)
    d <- dens[[1]][,n]
    plot(d$x, d$y, type='l', xlab='', ylab='', main=n, xlim=range(x), ylim=range(y), lwd=.25, col=alpha('black',.25))
    #points(d$x[d$peak],d$y[d$peak], pch=20, col='red')
    for (i in 1:length(dens)) {
        d <- dens[[i]][,n]
        lines(d$x, d$y, lwd=.25, col=alpha('black',.25))
        points(d$x[d$peak],d$y[d$peak], pch=20, col=1:sum(d$peak))
    }
}
dev.off()

