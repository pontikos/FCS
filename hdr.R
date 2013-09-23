library(feature)
library(KernSmooth)
library(RANN)

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

