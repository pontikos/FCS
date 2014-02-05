library(feature)
library(KernSmooth)
library(RANN)
library(FNN)

## instead of returning density in grid format
## returns 2D density at each (x,y) point
kde2D <- function(d, bw=.1) {
    # compute fast kernel density estimate
    b<-KernSmooth::bkde2D(as.matrix(d),bw)
    # this returns you a grid
    # we will use a fast nearest neighbour method
    # to find the closest point in the grid
    grid <- expand.grid(b$x1, b$x2)
    nn <- RANN::nn2(grid,d,k=1)
    return(cbind(d, dens=as.numeric(b$fhat)[nn$nn.idx]))
}

knn.dens <- function(x,k=10) {
    fnn.dist <- FNN::knn.dist(x,k)
    Dim <- dim(x)[[2]]
    N <- dim(x)[[1]]
    k / ( apply(fnn.dist,1, function(x, Radius=max(x)) Radius**Dim * pi**(Dim/2) / gamma(Dim/2+1)) ) 
}


