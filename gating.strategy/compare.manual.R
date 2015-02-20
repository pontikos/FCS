#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("tools"))

option_list <- list( 
    make_option(c("--manual"), help = 'manual'),
    make_option(c("--auto"), help = 'auto'),
    make_option(c("--outliers"), default=NULL, help = 'outliers'),
    make_option(c("--out.file"), help = '')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

manual <- '~nikolas/thor/Projects/IL2RA/CellPhenotypes/manual.csv'
auto <- '~nikolas/dunwich/Projects/IL2RA/Lymphocytes5/Out/lymph.count5.csv'

manual <- opt$manual
auto <- opt$auto

xlab <- gsub( paste('.',file_ext(manual),sep=''), '', basename(manual) ) 
xlab <- 'manual'
ylab <- gsub( paste('.',file_ext(auto),sep=''), '', basename(auto) ) 
ylab <- 'auto'

manual <- read.csv(manual)
colnames(manual) <- c('fcsFile','x')
auto <- read.csv(auto)
colnames(auto) <- c('fcsFile','y')

d <- merge(manual, auto)
d$fcsFile <- as.character(d$fcsFile)

if (!is.null(opt$outliers)) {
    outliers <- read.csv(opt$outliers, stringsAsFactors=FALSE)
    outliers$fcsFile <- basename(outliers$fcs.names)
    outliers <- outliers[,c('fcsFile','median.outliers')]
    print(head(outliers))
    print(dim(d<-merge(d, outliers)))
    print(head(d))
    d$col <- ifelse(d$median.outliers,'red','black')
} else {
    d$diff <- abs(d$x-d$y)
    thresh <- quantile(d$diff, seq(0,1,.05))[['95%']]
    d$col <- ifelse(d$diff<thresh,'black','red')
}

fileformat <- file_ext(opt$out.file)

#write.csv(d,quote=FALSE,row.names=F)

lim<-c(min(c(d$x,d$y)),max(c(d$x,d$y)))
p <- ggplot(d, aes(x=x, y=y)) + geom_point(col=d$col) + geom_abline(slope=1, intercept=0) + xlab(xlab) + ylab(ylab) + coord_fixed() + ylim(lim) + xlim(lim)
do.call(fileformat, list(opt$out.file))
print(p)
invisible(dev.off())



