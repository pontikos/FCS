#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(tools, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
source('~nikolas/bin/FCS/fcs.R')

# Simply pools all files together but in the future may decide to do normalisation before pooling.


option_list <- list( 
    make_option(c("-f","--in.file"),
                help = 'List of files in a csv file. Files may be RData, csv or fcs.'),
    make_option(c('-T','--total'), default=1E5,
                help='Total number of events in resulting pooled file.'),
    make_option(c('--channels'), default=NULL,
                help='Channels on which to do the clustering. By default select all channels.'),
    make_option(c('--plot.file'), default=NULL,
                help='Plot pooled file.'),
    make_option(c('-o','--out.file'),
                help='Pooled output file.')
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


# do checks first
if (is.null(opt$in.file) || !file.exists(opt$in.file)) stop("No input file specified on command line or file does not exist.")

if (!is.null(opt$channels)) {
    channels <- unlist(strsplit(opt$channels, ","))
} else  {
    channels <- NULL
}


# expect a list of fcs files in the csv file
fcs.names <- as.character(read.csv(opt$in.file,header=FALSE)[,1])

# will have to randomly sample N from each fcs file to obtain a pooled file with a total of T events
# round down
print(N <- floor(opt$total/length(fcs.names)))

cat('Pooling', N, 'lines, from each of the following', length(fcs.names), 'files on the channels', channels, '( total=', nrows <- N*length(fcs.names), ')', '\n')

fcs.filter <- function(fcs) {
    #remove min & max scatter
    min.ssc <- min(fcs[,'SSCA'])
    max.ssc <- max(fcs[,'SSCA'])
    fcs <- fcs[min.ssc < fcs[,'SSCA'] & fcs[,'SSCA'] < max.ssc,]
    min.fsc <- min(fcs[,'FSCA'])
    max.fsc <- max(fcs[,'FSCA'])
    fcs <- fcs[min.fsc < fcs[,'FSCA'] & fcs[,'FSCA'] < max.fsc,]
    fcs <- fcs[which(fcs[,'FSCA'] > 2),]
    fcs <- fcs[which(fcs[,'SSCA'] > 2),]
    return(fcs)
}


pooled.data <- do.call('rbind', lapply(fcs.names, function(fcs.name) {
    cat('Reading', fcs.name, 'extension', file.ext <- tolower(file_ext(fcs.name)), '\n')
    if (file.ext=='fcs') {
        fcs.data <- getChannels(read.FCS(fcs.name, filter=TRUE, channels=channels, which.lines=N), channels=channels)
    } else if (file.ext=='csv') {
        fcs.data <- read.csv(fcs.name)[1:N,channels]
    } else if (file.ext=='rdata') {
      print(load(fcs.name))
      #result
      fcs.data <- data.frame(fcs.data)
      #fcs.data <- fcs.filter(fcs.data)
      fcs.data <- fcs.data[1:N,channels]
    }
    fcs.data[,'fcs.name'] <- basename(fcs.name)
    print(head(fcs.data))
    return(fcs.data)
}))

print(dim(pooled.data))
print(head(pooled.data))
print(tail(pooled.data))
fcs.data <- pooled.data
save(fcs.data, file=opt$out.file, compress='xz', compression_level=9)

#d <- list(fcs.names=fcs.names, data=pooled.data)
#save(d,file=opt$out.file)

png(opt$plot.file)
pairs.smoothScatter(pooled.data[,channels])
dev.off()


###
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


#pool.data <- function(d, down.sample, channels) return(getChannels( do.call('rbind', lapply(d, function(x) getChannels(x[sample(1:nrow(x), down.sample),]))), channels ))
