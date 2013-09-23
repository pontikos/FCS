#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

option_list <- list( 
make_option(c("-f","--fcs"), help = "fcs file")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(flowCore, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))

read.FCS(opt$fcs)[,-9] -> fcs.data
c(1,2,3,5,7,8) -> p
fcs.data[rowSums(fcs.data@exprs[,p]>1)==dim(fcs.data@exprs[,p])[[2]],]->fcs.data

exprs <- function(fcs.data) {
    fcs.data@exprs -> x
    log10(x[,3:8]) -> x[,3:8]
    return(x)
}

exprs(fcs.data)->x


### gating strategy
# pe < 2
# cd4 > 2
# cd25 < 2
# cd45ra < 2
# naive: cd45ra > 1
# memory: cd45ra < 1

fsc <- x[,1]
ssc <- x[,2]

#scatter gate
x[ssc < 200 & fsc < 2000,] -> x

x[,6] -> pe

#CD4 gate
x[,7] -> cd4


#fixed gate
x[cd4>2 & pe<2,]->x

#non tregs gate
x[,5] -> apc
#x[,3] -> cd127
x[,8] -> cd45ra

#fixed gates
x[apc < 2 & cd45ra < 2,]->x
dim(x)[[1]] -> nontregs.count

x[,8] -> cd45ra

#naive
x[cd45ra > 1,]->naive
dim(naive)[[1]] -> naive.count
naive[,5] -> apc
naive[apc > 1,]->naive.cd25pos
naive[apc < 1,]->naive.cd25neg
#memory
x[cd45ra < 1,]->memory
dim(memory)[[1]] -> memory.count
memory[,5] -> apc
memory[apc > 1,]->memory.cd25pos
dim(naive.cd25pos)[[1]] -> naive.cd25pos.count
memory[apc < 1,]->memory.cd25neg



cat('>>fcsFile,APC.mfi,memory.freqpar,memory.count,naive.count,naive.cd25pos.freqpar,naive.cd25pos.count\n')
cat('>>')
cat(strsplit(tolower(opt$fcs), '.fcs')[[1]], mean(10**memory[,5]), 100*memory.count/nontregs.count, memory.count, naive.count, 100*naive.cd25pos.count/naive.count, naive.cd25pos.count, sep=",")
cat('\n')

