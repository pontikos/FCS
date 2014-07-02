#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
source('~nikolas/bin/FCS/fcs.R')
source('~nikolas/bin/FCS/normalise-functions.R')
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(tools, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(mvoutlier, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE)))) 

option_list <- list( 
make_option(c("--in.dir"), default='', help = 'Path to RData files.'),
make_option(c("--out.file"), default='', help = ''),
make_option(c("--plot.file"), default='', help = "")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

plot.file <- opt$plot.file
out.file <- opt$out.file

setwd(opt$in.dir)
print(load('prior.RData'))
i <- prior$lymph.cluster

files <- list.files(pattern='[^prior].*.RData$')[1:10]

fcs.names <- c()
mfi <- list()
median.fi <- list()
dens <- list()
for (f in files) {
    print(f)
    print(load(f))
    if (exists("lymph.filter")) {
    mfi[[f]] <- c(res@mu[i,], res@w[[i]])
    median.fi[[f]] <- apply(d[lymph.filter,],2,median)
    fcs.names <- c(fcs.names, fcs.name)
    dens[[f]] <- sapply(res@varNames, function(v) normalised.density(d[lymph.filter,v]))
    rm(lymph.filter)
    }
}

#if (length(mfi)!=length(median.fi)!=length(dens)) stop('lengths different', length(mfi), length(median.fi), length(dens))



###
plotOutliers <- function(x, dens, plot.file=NULL, channels=NULL, outliers=x$mean.outliers) {
    if (is.null(channels)) {
        col.names <- colnames(x)
        nc <- ncol(x)
    } else {
        col.names <- channels
        nc <- length(channels)
    }
    cat('>>',plot.file,'\n')
    if (!is.null(plot.file)) pdf(plot.file)
    par(mfrow = c(nc, nc), mai=c(.2,.2,.2,.2), oma=c(2,3,2,.5))
    for (n in 1:length(col.names)) {
            ylab <- ''
            xlab <- ''
            main <- ''
        for (n2 in 1:length(col.names)) {
            # off-diagonal scatter plot median outliers
            if (n != n2) {
                dim.inds <- c(n2,n)
                plot(x[,dim.inds],pch=20,col=ifelse(outliers,'red','black'),xlab='',ylab='',main='')
            #diagonal plot univariate densities
            } else {
                xlim<-range(sapply(dens, function(d) d[,n]$x))
                ylim<-range(sapply(dens, function(d) d[,n]$y))
                plot(dens[[1]][,n], type='l', xlab='', ylab='', col=ifelse(outliers[[1]], 'red', 'black'), main='', lwd=.1, xlim=xlim, ylim=ylim)
                for (i in 1:length(dens)) {
                    swp <- sliding.window.peaks(dens[[i]][,n],span=40)
                    swp <- swp[which(swp$y > .05*max(swp$y)),]
                    print(swp)
                    points(swp$x, swp$y, pch=20)
                    lines(dens[[i]][,n], col=ifelse(((outliers[[i]]) || (nrow(swp)>1)), 'red', 'black'), lwd=.1)
                }
            }
            if (n==1) mtext(col.names[[n2]], side=3, cex=2, line=1) 
            if (n2==1) mtext(col.names[[n]], side=2, cex=2, line=2)
       }
    }
    if (!is.null(plot.file)) dev.off()
}



print(head(mfi <- do.call('rbind', mfi)))
colnames(mfi) <- c(res@varNames, 'tau')

median.fi <- do.call('rbind', median.fi)
colnames(median.fi) <- res@varNames
print(median.fi)

print(mean.outliers <- uni.plot(mfi[,res@varNames]))
print(median.outliers <- uni.plot(median.fi[,res@varNames]))
head(mfi <- data.frame(fcs.names, mfi, mean.outliers=mean.outliers$outliers, median.outliers=median.outliers$outliers))
print(dim(mfi))
write.csv(mfi, file=out.file, quote=FALSE)

plotOutliers(mfi[,res@varNames], dens, plot.file=plot.file, channels=res@varNames, outliers=mfi$mean.outliers)

