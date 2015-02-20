#!/usr/bin/env Rscript
source('~nikolas/bin/FCS/fcs.R')
suppressPackageStartupMessages(library(sp))
suppressPackageStartupMessages(library(flowClust))

#
#
flowClust2Prior <- function (x, kappa, Nt) {
    p <- ncol(x$mu)
    K <- x$K
    nu0 <- Ng <- x$w * Nt
    if (all((nu0 * kappa - p - 1) > 0)) {
        Lambda0 <- x$sigma
        for (i in 1:K) {
            Lambda0[i, , ] <- Lambda0[i, , ] * (kappa * nu0[i] - p - 1)
        }
    } else {
        stop("Can't proceed. Prior nu0 is negative for cluster(s) ", paste(which((nu0 - p - 1) > 0), collapse = ","), "\n(p-1) = ", p - 1, ": Try increasing kappa")
    }
    Omega0 <- array(0, c(K, p, p))
    for (i in 1:K) {
        Omega0[i, , ] <- diag(1, p)
        if (p == 1) {
            dS <- x$sigma[i, , ]
            dO <- Omega0[i, , ]
        }
        else {
            dS <- det(x$sigma[i, , ])
            dO <- det(Omega0[i, , ])
        }
        k <- (dO/dS)^(1/p)
        Omega0[i, , ] <- Omega0[i, , ] * k
        Omega0[i, , ] <- solve(Omega0[i, , ] * Ng[i] * kappa)
    }
    nu0 <- nu0 * kappa
    Mu0 <- x$mu
    lambda <- x$lambda
    w0 <- x$w * Nt
    prior <- list(Mu0 = Mu0, Lambda0 = Lambda0, Omega0 = Omega0, w0 = w0, nu0 = nu0, nu = x$nu, lambda = x$lambda, K = K)
    class(prior) <- "flowClustPrior"
    attr(prior, "lambda") <- x$lambda
    #name all array dimensions for later subsetting
    channels <- colnames(prior$Mu0)
    prior$Lambda0 <- array(prior$Lambda0,dim=c(K,2,2),dimnames=list(1:K,channels,channels))
    prior$Omega0 <- array(prior$Omega0,dim=c(K,2,2),dimnames=list(1:K,channels,channels))
    prior
}


# K is the cluster labels, must be numerically ordered starting at 1
# kappa > 1
# smaller kappa means less certainty
# one value of kappa probably doesn't fit all
# Returns two-component with 
# background cluster (position 1)
# and prior for gate (position 2).
# Probably better for the background to be everything.
manual.prior <- function(d, K, kappa=100, Nt=NULL) {
    print(sort(unique(K)))
    p <- ncol(d)
    if (is.null(Nt)) Nt <- nrow(d)
    Mu <- do.call('rbind', by(d, K, colMeans))
    #Mu[1,] <- background.mean[colnames(d)]
    Sigma <- array(dim=c(length(unique(K)),p,p))
    for (i in sort(unique(K))) Sigma[i,,] <- cov(d[which(K==i),])
    #Sigma[1,,] <- background.cov[colnames(d),colnames(d)]
    w <- as.numeric(prop.table(table(K)))
    w[1] <- .0001
    theta <- list(K=as.numeric(length(unique(K))), mu=Mu, sigma=Sigma, w=w, lambda=1, nu=4)
    print(prior <- flowClust2Prior(theta, kappa=kappa, Nt=Nt))
    level=.9
    u.cutoff=.5
    B=500
    kappa=1
    print(colnames(d))
    return(prior)
}


