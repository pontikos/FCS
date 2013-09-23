#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

# Gate lymphocytes using flowStats::lymphGate
# requires asinh transform to function properly
# file is saved with same name but with lymphocytes.fcs ending

option_list <- list( 
make_option(c("-f","--fcs"), help = "fcsfile to parse")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

if (is.null(opt$fcs)) {
stop("Not FCS file specified on command line!")
}

suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(flowCore, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
#suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(flowStats, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))

#suppressWarnings(read.FCSheader(opt$fcs, keyword=keywords)) -> fcs.header
suppressWarnings(read.FCS(opt$fcs)) -> fcs.data

#apply compensation
print(fcs.data@description[["SPILL"]])
if (!is.null(fcs.data@description[["SPILL"]])) compensate(fcs.data,fcs.data@description[["SPILL"]])->fcs.data

# find the fluorochrome name for cd4
#fcs.data@description[grep("P[[:digit:]]+S", names(fcs.data@description))]->antibodies
#fcs.data@description[grep("P[[:digit:]]+N", names(fcs.data@description))]->channels
channels <- parameters(fcs.data)@data$name
antibodies <- parameters(fcs.data)@data$desc

cd4.channel <- channels[[grep('^cd4$', antibodies, ignore.case=T)]]
print(cd4.channel)

lgcl <- logicleTransform()

cd4 <- lgcl(fcs.data@exprs[,cd4.channel])
fsc <- fcs.data@exprs[,'FSC-A']
ssc <- fcs.data@exprs[,'SSC-A']

total.count <- length(cd4)
print(alpha <- 37863/total.count)

# cd4 gate
cd4.gate <- 2 < cd4 & cd4 < 3
# fsc gate
fsc.min <- 50000 
fsc.max <- 100000
fsc.gate <- fsc.min < fsc & fsc < fsc.max
table(fsc.gate&cd4.gate)
# ssc gate
q.ssc <- quantile(ssc[fsc.gate&cd4.gate],probs=seq(0,1,.05))
#ssc.max <- q.ssc[['90%']]
ssc.max <- 80000
ssc.gate <- ssc < ssc.max
# cd4 & fsc & ssc gate
lymph.gate <- cd4.gate&fsc.gate&ssc.gate
prop.table(table(lymph.gate))

lymph.count <- sum(as.numeric(lymph.gate))
lymph.pct <- 100*lymph.count/total.count


png(paste(gsub('.fcs', '', opt$fcs), '.lymphocytes2.png', sep=''))
plot(fcs.data@exprs[,c('FSC-A','SSC-A')],col=rgb(0,0,0,alpha), pch='.', main=paste('CD4+ T cells:', round(lymph.pct),'% =',lymph.count, '/', total.count ))
#cd4 gate red
points(fcs.data@exprs[cd4.gate,c('FSC-A','SSC-A')],col=rgb(1,0,0,alpha), pch=".")
#fsc gate green
points(fcs.data@exprs[fsc.gate,c('FSC-A','SSC-A')],col=rgb(0,1,0,alpha/2), pch=".")
abline(v=c(fsc.min, fsc.max), col='green')
#ssc gate blue
points(fcs.data@exprs[ssc.gate,c('FSC-A','SSC-A')],col=rgb(0,0,1,alpha/4), pch=".")
abline(h=ssc.max, col='blue')
dev.off()


cat('>',opt$fcs,',total.count,', total.count, '\n',sep='')
cat('>',opt$fcs,',lymph.count,', lymph.count, '\n',sep='')
cat('>',opt$fcs,',lymph.pct,', lymph.pct, '\n',sep='')

out.file <- paste(gsub('.fcs', '', opt$fcs), '.nonlymphocytes2.fcs', sep='')
cat('>outfile,',out.file,'\n',sep='')

write.FCS(fcs.data[-lymph.gate,], filename=out.file)


