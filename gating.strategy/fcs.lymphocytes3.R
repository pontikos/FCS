
source('~nikolas/bin/hdr.R')

#library(sm)
png('~/plot-density-%03d.png')
#mapply(function(d,clus,m) { plot(getChannels(d,c('SSC-A','CD4')), pch=20, cex=.5, col=clus); points(m[,2:3],pch=rownames(m), cex=2) }, fcs.eff.down[1], fcs.eff.down.clust[1], medoids[1])
#plot(getChannels(fcs.eff[[1]],c('SSC-A','CD4')), pch=20, cex=.5, col=)
mapply( function(X) {
        f <- feature::featureSignif(X)
        plot(f)
        x.dens <- kde2D(X)
        #keep the top 5% densest points
        #dens <- sm.density(X, h=c(0.1,0.1), eval.points=X, eval.grid=FALSE)
        #thresh <- quantile(dens$estimate, .95)
        thresh <- quantile(x.dens[,3], .9)
        #quantile(f$fhat$est,.95)
        x.topdens <- X[x.dens[,3] > thresh,]
        #cluster top density points
        #h <- hclust(dist(x.topdens))
        #points(x.topdens, pch=20, cex=.5, col=cutree(h, k=6))
        points(x.topdens, pch=20)
        }
       #, lapply(fcs.eff, getChannels, channels=c('FSC-A','CD4','SSC-A'))[1] )
       #, lapply(fcs.eff, getChannels, channels=c('FSC-A','CD4')))
       , lapply(fcs.eff, getChannels, channels=c('FSC-A','FSC-W')))
dev.off()

png('~/plot-density-FSC-CD4-%03d.png')
lapply( lapply(fcs.eff, getChannels, channels=c('FSC-A','CD4')),
       function(X) {
        f <- feature::featureSignif(X)
        plot(f)
        x.dens <- kde2D(X)
        #keep the top 10% densest points
        thresh <- quantile(x.dens[,3], .9)
        x.topdens <- X[x.dens[,3] > thresh,]
        points(x.topdens, pch=20)
       } )
dev.off()




png('~/donor6-lymph-%03d.png')
# crude lymph gating
lapply( lapply(fcs.eff.6, getChannels, channels=c('FSC-A','CD4')),
       function(X) {
        #print(dim(X)[1])
        cd4 <- X[,'CD4'] 
        x <- X[which(2 < cd4 & cd4 < 4),]
        #print(dim(x)[1])
        print(100*dim(x)[1]/dim(X)[1])
        smoothScatter(X)
        abline(h=c(2,4))
        X.dens <- kde2D(X)
        thresh <- quantile(X.dens[,3], .9)
        X.topdens <- X[X.dens[,3] > thresh,]
        points(X.topdens, pch=20, col='red')
       } )
dev.off()



