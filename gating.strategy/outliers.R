#!/usr/bin/env Rscript
source('~nikolas/bin/FCS/fcs.R')
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(tools, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(mvoutlier, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))

print(load('prior.RData'))
i <- prior$lymph.cluster

files <- list.files(pattern='[^prior].*.RData$')[1:10]

fcs.names <- c()
mfi <- list()
dens <- list()
for (f in files) {
    print(f)
    print(load(f))
    mfi[[f]] <- c(res@mu[i,], res@w[[i]])
    fcs.names <- c(fcs.names, fcs.name)
    dens[[f]] <- sapply(res@varNames, function(v) normalised.density(d[lymph.filter,v]))
}

mfi <- do.call('rbind', mfi)
colnames(mfi) <- c(res@varNames, 'tau')

head(mfi <- data.frame(fcs.names, mfi, outliers))
write.csv(mfi, file='lymphocyte-mfi.csv', quote=FALSE)

pdf('outliers.pdf')
print(outliers <- uni.plot(mfi[,res@varNames]))
for (v in res@varNames) {
    f.1 <- files[[1]]
    plot(dens[[f.1]][,v], main=v, type='l', xlab='', ylab='', col=ifelse(mfi[f.1,'outliers'], 'red', 'black'))
    for (f in files) lines(dens[[f]][,v], col=ifelse(mfi[f,'outliers'], 'red', 'black'))
}
dev.off()




