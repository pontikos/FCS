#!/usr/bin/env Rscript
source('~nikolas/bin/FCS/fcs.R')
suppressPackageStartupMessages(library("optparse"))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(tools, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))

option_list <- list( make_option(c('--channels'), help='Channels to print out to CSV file.') )
option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

#channels <- c('PSTAT5','CD25','CD45RA','CD4','FOXP3','FSCA','SSCA')
channels <- unlist(strsplit(opt$channels, ","))

print(load('prior.RData'))
i <- prior$lymph.cluster

files <- list.files(pattern='[^prior].*.RData$')
for (f in files) {
    print(f)
    print(load(f))
    print(lymph.csv.file <- file.path(paste(basename(fcs.name),'.csv',sep='')))
    if (!file.exists(lymph.csv.file)) {
    lymph <- getChannels(read.FCS(fcs.name, channels)[lymph.filter,], channels)
    write.csv(lymph,file=lymph.csv.file,quote=FALSE,row.names=F)
    } else {
        cat(lymph.csv.file, 'exists already, skipping.\n')
    }
}

