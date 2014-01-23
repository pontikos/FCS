#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

option_list <- list( 
make_option(c("-f","--fcs"), help = "fcsfile to parse"),
make_option(c("-o","--output"), help = "output fcsfile prefix"),
make_option(c("-K","--NumC"), default=1, type="integer", help="Number of clusters.")
)

OptionParser(option_list=option_list) -> option.parser
parse_args(option.parser) -> opt

if (is.null(opt$fcs)) {
stop("Not FCS file specified on command line!")
}

suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(flowCore, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
suppressMessages(suppressWarnings( read.FCS(opt$fcs, alter.names=TRUE, column.pattern=".A") -> fcs.data ))

dimnames(fcs.data@exprs)[[2]] -> gating.parameters
cat("Parameters:", gating.parameters, "\n")

#gate on FSC SSC with flowClust
suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(flowClust, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
flowClust(fcs.data, varNames=c("FSC.A","SSC.A"), K=opt$NumC)->res
for (n in 1:opt$NumC) {
    fcs.data@exprs[which(Map(res)==n),]
}

gate <- function(fcs.data, scatter.gate, K, algorithm, gating.parameters) {
    if (algorithm=="flowMeans") {
        suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(flowMeans, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
        #flowMeans(fcs.data, varNames=gating.parameters, nstart=nstart, NumC=NumC, MaxN=MaxN, Standardize=TRUE, Update=Update) -> res
        flowMeans(fcs.data, varNames=gating.parameters, NumC=K, Standardize=TRUE, Mahalanobis=T) -> res
        levels(factor(res@Label))->clusters
        sort(sapply(clusters, function (n) mean(fcs.data@exprs[which(res@Label==n),"APC.A"])))->apc.mfi
        sort(sapply(clusters, function (n) length(fcs.data@exprs[which(res@Label==n),"APC.A"])))->apc.length
    } else if (algorithm=="flowClust") {
        #gating.parameters<-c("APC.A","PE.A")
        flowClust(fcs.data, varNames=c("APC.A"), K=6) -> res
        #levels(factor(res@Label))->clusters
        sort(sapply(1:K, function (n) mean(fcs.data@exprs[which(Map(res)==n),"APC.A"])))->apc.mfi
        print(apc.mfi)
    } else if (algorithm=="curvHDR") {
        suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(curvHDR, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
        curvHDRfilter(asinh(fcs.data@exprs[,gating.parameters]),HDRlevel=.5)->res
        print(sapply(res$polys, function(n) sinh(n)))
    } else if (algorithm == "SamSPECTRAL"){
        #SamSPECTRAL(fcs.data@exprs, dimensions=1:dim(data.points)[2], normal.sigma, separation.factor,number.of.clusters=6, scale=rep(1,dim(data.points)[2]), talk = TRUE, precision = 6, eigenvalues.num =NA, return_only.labels=TRUE, do.sampling=TRUE, beta=4, stabilizer=1000)
    } else if (algorithm == "kmeans") {
        fcs.data[which(fcs.data@exprs[,"APC.A"]>0),]->fcs.data
        print(fcs.data)
        log(fcs.data@exprs[,"APC.A"])->log.apc
        kmeans(log.apc, centers=K, nstart=opt$nstart)->res
        levels(factor(res$cluster))->clusters
        sort(sapply(clusters, function (n) mean(fcs.data@exprs[res$cluster==n,"APC.A"])))->apc.mfi
    } else if (algorithm == "pam") {
        suppressMessages(suppressWarnings(suppressPackageStartupMessages(library(cluster, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE))))
        log(fcs.data@exprs[,gating.parameters]+100)->log.intensity
        pam(log.intensity, k=K)->res
        sort(tapply(fcs.data@exprs[,gating.parameters],res$clustering,mean))->mfi
    }
    return(list(res=res,fcs.data=fcs.data,mfi=mfi))
}


