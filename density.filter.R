#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
option_list <- list( 
    make_option(c("-f","--in.file"), help = 'FCS file to parse or list of files in a csv file.'),
    make_option(c('--channels'), default='FSC-A,SSC-A,CD4', help='Channels on which to do the clustering.  By default CD4 lymphocyte gating works on Forward, Side Scatter and CD4.'),
    make_option(c('--density.threshold'), default='50%', help='The posterior cutoff.'),
    make_option(c('--plot.dir'), default=NULL, help='The directory to which to send the plots.')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

# do checks first
if (is.null(opt$in.file) || !file.exists(opt$in.file)) stop("No FCS file specified on command line or file does not exist.")


suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(mclust, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(tools, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(feature, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(KernSmooth, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(RANN, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))


source('~nikolas/bin/FCS/fcs.R')

if (!is.null(opt$channels)) {
    channels <- unlist(strsplit(opt$channels, ","))
} else  {
    channels <- NULL
}

if (file_ext(opt$in.file)=='fcs') {
    # expect a single fcs files
    fcs.files <- list(suppressWarnings(read.FCS(opt$in.file, channels=channels)))
    names(fcs.files) <- opt$in.file
} else if (file_ext(opt$in.file)=='csv') {
    # expect a list of fcs files in the csv file
    fcs.names <- as.character(read.csv(opt$in.file)[,1])
    fcs.files <- lapply(fcs.names, read.FCS, channels=channels)
    names(fcs.files) <- fcs.names
    }

print(names(fcs.files))

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



###
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



density.filter <- function(x, channels, bw=.1, quant='75%') {
    x <- getChannels(x,channels=channels)
    x <- apply(x, 2, scale)
    quantiles <- apply(combn(colnames(x),2),2,
              function(p) {
                print(p)
                dens <- kde2D(x[,p],bw=.1)[,3]
                dens <- dens/sum(dens)
                return(dens > quantile(dens,probs=seq(0,1,.05))[[quant]])
              })
    return(rowSums(quantiles)==ncol(combn(channels,2)))
}


###
###
plot.density <- function(d, dens.filter, channels, plot.file, quant) {
        d <- getChannels(d, channels=channels)
        d.sub <- d[dens.filter,]
        png(plot.file)
        par(mfrow=c(2,2))
        #
        #x<-d.sub[sample(1:nrow(d.sub),1000),]
        x<-d.sub
        #classification <- Mclust(apply(x,2,scale))$classification
        #classification <- kmeans(apply(x,2,scale),centers=)$cluster
        #i <- which.min(abs(tapply(x[,'CD4'], classification, mean)-2.5))
        apply(combn(colnames(d),2),2,
              function(p) {
                smoothScatter(d[,p])
                #points(x[,p], pch=20, col=as.numeric(classification==i))
                points(x[,p], pch=20)
              })
        #
        image(1:3,1:3,matrix(data=0,nrow=3,ncol=3), axes=FALSE, col='white', ylab='', xlab='')
        text(2,3, 'count', font=2)
        text(3,3, '%', font=2)
        text(1,1, quant, font=2)
        text(2,1, round(sum(dens.filter)) )
        text(3,1, round(100*sum(dens.filter)/dim(d)[[1]]) )
        text(1,2, 'total', font=2)
        text(2,2, dim(d)[[1]])
        text(3,2, 100)
        dev.off()
}



for (fcs.name in names(fcs.files)) {
    fcs.data <- fcs.files[[fcs.name]]
    print(fcs.name)
    table(f <- density.filter( fcs.data, channels=channels, quant=opt$density.threshold, bw=.1 ))
    print(plot.file <- file.path(opt$plot.dir,paste(basename(fcs.name),'.png',sep='')))
    plot.density(d=fcs.data, dens.filter=f, channels=channels, plot.file=plot.file, quant=opt$density.threshold)

}

