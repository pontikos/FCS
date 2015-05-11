#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("mclust"))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(tools, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
source('~nikolas/bin/FCS/fcs.R')

# Identifies peaks using univariate Mclust.


option_list <- list( 
    make_option(c("-f","--in.file"),
                help = 'List of files in a csv file. Files may be csv or fcs.'),
    make_option(c('--channel'), default=NULL,
                help=''),
    make_option(c('--plot.file'), default=NULL,
                help='Plot file.'),
    make_option(c('-K','--peaks'), default=5,
                help='Max number of peaks to consider')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

channel <- opt$channel

## READ DATA
# takes a single fcs, csv or RData file
cat('Reading', in.file <- opt$in.file, '\n')
cat('File extension', file.ext <- tolower(file_ext(opt$in.file)), '\n')
if (file.ext=='fcs') {
    #fcs file
    print(head(x <- getChannels(read.FCS(in.file, channel=channel),channel)))
} else if (file.ext=='csv') {
    #csv file
    print(head(x <- read.csv(in.file)[,channel]))
    print(total.count <- dim(x)[[1]])
} else if (file.ext=='RData') {
    print(load(in.file))
    print(head(x <- d$data[,channel]))
    print(total.count <- dim(x)[[1]])
} else {
    stop('Unsupported file extension', file.ext,'\n')
}


res <- Mclust(x, G=1:opt$peaks)

plot(density(x))


cat('>>')
cat(opt$in.file,sort(res$parameters$mean),sep=',')
cat('\n')

