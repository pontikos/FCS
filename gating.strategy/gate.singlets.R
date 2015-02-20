#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(tools, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))

gate.singlets <- function(X,k=3) {
    #filter on SSCW
    sscw <- X[,'SSCW']
    med <- median(sscw)
    sscw.filter <- which(med-k*median(abs(med-sscw)) < sscw & sscw < med+k*median(abs(med-sscw)))
    X[sscw.filter,]
    #FSCH,SSCH,FSCW,SSCW
    #m <- covMcd(X)
    #m$mah < max(m$mah)
    #mahalanobis(X)
}


option_list <- list( 
    make_option(c("-f","--in.file"), help = 'FCS file to parse or list of files in a csv file.')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

X <- read.csv(opt$in.file)
X <- gate.singlets(X)
write.csv(X, file=opt$in.file, quote=FALSE, row.names=FALSE)




