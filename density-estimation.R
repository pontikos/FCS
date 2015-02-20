source('~nikolas/bin/FCS/hdr.R')

library(ggplot2)
library(mvtnorm)
library(spade)
library(ks)
library(KernSmooth)
library(MASS)
library(feature)
library(KernSmooth)
library(RANN)
library(FNN)



## instead of returning density in grid format
## returns 2D density at each (x,y) point
kde2D <- function(d, bw=.1) {
    # compute fast kernel density estimate
    b<-bkde2D(as.matrix(d),bw)
    # this returns you a grid
    # we will use a fast nearest neighbour method
    # to find the closest point in the grid
    grid <- expand.grid(b$x1, b$x2)
    nn <- nn2(grid,d,k=1)
    return(cbind(d, dens=as.numeric(b$fhat)[nn$nn.idx]))
}


sigma<- function(v, r, p) {
 V<- matrix(r^2, ncol=p, nrow=p)
 diag(V)<- 1
 V*v
}

# 2 dimensional
X<- data.frame(mvtnorm::rmvnorm(1000, mean=rep(0, 2), sigma(1, .5, 2)))
real.dens <- mvtnorm::dmvnorm(X, mean=rep(0, 2), sigma(1, .5, 2))

X <- kde2D(X[,1:2],.1)
X$knn.dens <- knn.dens(X,k=10)

ggplot(X, aes(x=X1,y=X2,z=z.knn)) + geom_point() + stat_contour()

ggplot(X, aes(x=X1,y=X2,z=dens)) + geom_point() + geom_density2d()

plot(real.dens, kde2D(X[,1:2],.5)[,3])
plot(real.dens, knn.dens(X,k=100))


# kernel_mult: Multiplier of the minimum median distance within which other observations are counted towards the density
# apprx_mult: Multiplier of the minimum median distance within which observations are approximated to have the same density
# med_samples: Number of observations used to estimate the minimum median distance
spade.density <- function(X, kernel_mult=5.0, apprx_mult=0.5, med_samples=2000) SPADE.density(X, kernel_mult, apprx_mult, med_samples)

# k: Number of nearest neighbours to use in density estimation.
knn.density <- function(x,k=10) {
    fnn.dist <- FNN::knn.dist(x,k)
    Dim <- dim(x)[[2]]
    N <- dim(x)[[1]]
    k / ( apply(fnn.dist,1, function(x, Radius=max(x)) Radius**Dim * pi**(Dim/2) / gamma(Dim/2+1)) ) 
}


`%in%` <- function(x, r) in.range(x,r)

#
downsample.indexes <- function(dens, desired_samples=10) {
    # Need to find target density such there are approximately desired_samples
    # remaining after downsampling. To do so we solve for the density such that
    # the sum of samples below that density plus the expected value of
    # samples retained above that density equals approximately the desired
    # number of samples
    density_s <- sort( dens )
    cdf       <- rev(cumsum(1.0/rev(density_s)))
    # Default solution if target density smaller than any present
    boundary2 <- desired_samples/cdf[1] 
    # Boundary actually falls amongst densities present
    if (boundary2 > density_s[1]) {
      targets <- (desired_samples-1:length(density_s)) / cdf 
      boundary2 <- targets[which.min(targets-density_s > 0)]
    }
    return( which( boundary2 > dens*runif(length(dens))) ) 
    #return( which( dens %in% c(boundary1, dens*runif(length(dens))) ) )
}


# D dimensional mvt normal used for testing accuracy
# k neigbours
mse.dens.fun <- function(dens.fun, D, N, ...) {
    set.seed(1234)
    X <- mvtnorm::rmvnorm(N, mean=rep(0, D), sigma(1, .5, D))
    true.dens <- mvtnorm::dmvnorm(X, mean=rep(0, D), sigma(1, .5, D))
    est.dens <- dens.fun(X,...)
    x <- true.dens/sum(true.dens)
    y <- est.dens/sum(est.dens)
    return( mean((x-y)**2) )
}

# D dimensional mvt normal used for testing accuracy
# k neigbours
agreement.dens.fun <- function(dens.fun, D,...) {
    set.seed(1234)
    X <- mvtnorm::rmvnorm(1000, mean=rep(0, D), sigma(1, .5, D))
    real.dens <- mvtnorm::dmvnorm(X, mean=rep(0, D), sigma(1, .5, D))
    est.dens <- dens.fun(X,...)
    x <- real.dens/sum(real.dens)
    y <- est.dens/sum(est.dens)
    lim <- range(c(x,y))
    plot(x, y, xlab='real density', ylab='estimated density', xlim=lim, ylim=lim, pch=20)
    abline(b=1,a=0, lwd=2)
    abline(line(y,x))
    return(coef(line(y,x))[[2]])
}


agreement.dens.fun(spade.density, D=10, kernel_mult=1000)
benchmark.dens.fun(spade.density, D=10)

system.time( spade.r.squared <- sapply( 2:50, function(n) benchmark.dens.fun(spade.density,n)) )
system.time( spade.r.squared.10 <- sapply( 2:50, function(n) benchmark.dens.fun(spade.density,n, kernel_mult=2)) )
system.time( knn.r.squared.5 <- sapply( 2:50, function(n) benchmark.dens.fun(knn.density,n, k=5)) )
system.time( knn.r.squared.10 <- sapply( 2:50, function(n) benchmark.dens.fun(knn.density,n, k=10)) )
system.time( knn.r.squared.50 <- sapply( 2:50, function(n) benchmark.dens.fun(knn.density,n, k=20)) )
system.time( knn.r.squared.50 <- sapply( 2:50, function(n) benchmark.dens.fun(knn.density,n, k=50)) )
system.time( knn.r.squared.100 <- sapply( 2:50, function(n) benchmark.dens.fun(knn.density,n, k=100)) )
system.time( knn.r.squared.150 <- sapply( 2:50, function(n) benchmark.dens.fun(knn.density,n, k=150)) )
system.time( knn.r.squared.1000 <- sapply( 2:50, function(n) benchmark.dens.fun(knn.density,n, k=999)) )

png('~/test.png')
plot(2:50,spade.r.squared, type='l', xlab='dimensions', ylab='MSE')
lines(2:50, spade.r.squared.10, col='black', lty=2)
lines(2:50, knn.r.squared.5, col='blue')
lines(2:50, knn.r.squared.10, col='red')
lines(2:50, knn.r.squared.20, col='green')
lines(2:50, knn.r.squared.50, col='purple')
lines(2:50, knn.r.squared.100, col='pink')
lines(2:50, knn.r.squared.150, col='pink')
lines(2:50, knn.r.squared.1000, col='pink', lty=2)
legend(
dev.off()

# what is the optimal k?
k <- 10:500
system.time( knn.mse.500 <- sapply( k, function(k) mse.dens.fun(knn.density, D=2, N=600, k=k)) )
system.time( knn.mse.1000 <- sapply( k, function(k) mse.dens.fun(knn.density, D=2, N=1000, k=k)) )
system.time( knn.mse.2000 <- sapply( k, function(k) mse.dens.fun(knn.density, D=2, N=2000, k=k)) )
plot(k, knn.mse.500, ylim=range(knn.mse.500, knn.mse.1000,knn.mse.2000))
points(k, knn.mse.1000, col='red')
points(k, knn.mse.2000, col='blue')

# For a given number of dimensions
# the number of neighbours has no bearing on the accuracy of the density estimation
# REALLY?? the plot above shows the opposite!
system.time( r.squared <- sapply( 2:100, function(k) benchmark.knn.dens(k=k, 2) ) )
plot(2:100,r.squared, type='l', xlab='# of neighbours')


D <- 2
X <- mvtnorm::rmvnorm(1000, mean=rep(0, D), sigma(1, .5, D))
est.dens <- knn.dens(X, k=100)
#est.dens <- spade.dens(X, kernel_mult=2.5)
#est.dens <- spade.dens(X, kernel_mult=10)
real.dens <- mvtnorm::dmvnorm(X, mean=rep(0, D), sigma(1, .5, D))
plot(est.dens/sum(est.dens), real.dens/sum(real.dens))
cor(est.dens/sum(est.dens), real.dens/sum(real.dens))**2
abline(b=1,a=0)

plot(X, pch=20)
print(i <- downsample.indexes(real.dens, desired_samples=100))
points(X[i,], col='red', pch=20, cex=2)

plot(density(X[,1]))
lines(density(X[i,1]), col='red')
plot(density(X[,2]))
lines(density(X[i,2]), col='red')


mean(replicate(1000, length(downsample.indexes(real.dens, desired_samples=10))))

#d <- cbind(X,real.dens)
#points(d[order(d[,3]),1:2][1:20,], col='red', pch=20)

# 1 D density estimation
N <- 100 
x <- rnorm(N)
true.dens.x <- dnorm(x)

bw <- seq(.01,1,.01)
mse <- sapply(bw, function(bw) {
    d <- density(x,bw=bw)
    est.dens <- splinefun(d)
    est.dens.x <- est.dens(x) 
    mean((est.dens.x-true.dens.x)**2)
})

plot(bw, mse)
print(bw[which.min(mse)])


par(mfrow=c(2,2))
for (N in 2:5) {
set.seed(1234)
x <- rnorm(10**N)
true.dens.x <- dnorm(x) 
d <- density(x)
print(d$bw)
est.dens <- splinefun(d)
est.dens.x <- est.dens(x) 
ylim <- xlim <- range(c(est.dens.x, true.dens.x))
plot( est.dens.x, true.dens.x, xlim=xlim, ylim=ylim, pch=20, cex=.5, main=sprintf('D=1 N=10^%d',N))
abline(b=1,a=0)
}

par(mfrow=c(1,1))
N <- 2**(10:20)
mse <- sapply(N, function(N) {
    x <- rnorm(N)
    true.dens.x <- dnorm(x)
    d <- density(x)
    est.dens <- splinefun(d)
    est.dens.x <- est.dens(x) 
    mean((est.dens.x-true.dens.x)**2)
})
plot(N, mse, type='l')

plot(log2(N), sapply(N, function(N) density(rnorm(N))$bw), type='l')


N <- 1000
par(mfrow=c(2,2))
for (D in 2:5) {
set.seed(1234)
X <- mvtnorm::rmvnorm(N, mean=rep(0, D), sigma(1, .5, D))
true.dens.x <- mvtnorm::dmvnorm(X, mean=rep(0, D), sigma(1, .5, D))
d <- ks::kde(X, eval.points=X)
est.dens.x <- d$estimate
ylim <- xlim <- range(c(est.dens.x, true.dens.x))
plot(est.dens.x, true.dens.x, pch=20, cex=.5, xlim=xlim, ylim=ylim, main=sprintf('D=%d N=1000',D))
abline(b=1,a=0)
}




head(MASS::kde2D(X))

bkde2D(X)
str(bkfe(X))

## instead of returning density in grid format
## returns 2D density at each (x,y) point
kde2D <- function(d) {
    # compute fast kernel density estimate
    b<-KernSmooth::bkde2D(as.matrix(d))
    # this returns you a grid
    # we will use a fast nearest neighbour method
    # to find the closest point in the grid
    grid <- expand.grid(b$x1, b$x2)
    nn <- RANN::nn2(grid,d,k=1)
    return(cbind(d, dens=as.numeric(b$fhat)[nn$nn.idx]))
}




