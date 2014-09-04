library(flowCore)

#set w=0 is important to avoid spurious peaks close to zero
#w=0 is essentially equivalent to asinh according to Moore et al 2012, Cytometry A
lgcl <- logicleTransform(w=0.1)
lgcl <- function(y,m=3) logicleTransform(w=(m-log10(max(y)/abs(min(y))))/2)(y) 

estimateTransforms <- function(x, m=4.5) {
    sapply(colnames(x), function(n) {
           if (grepl('^SSC|^FSC',n))
               return(identity)
           y <- x[,n]
           y <- y[which(percentile.filter(y,bottom=.01,top=.98))]
           print(n)
           print(t <- max(y))
           print(abs(min(y)))
           print(r <- log10(t/abs(min(y))))
           if (r < 0) w <- 2.25
           else w <- (m-r)/2
           print(w)
           return(logicleTransform(w=w,m=m,t=t))
    })
}

applyTransforms.slow <- function(x, transforms) {
    sapply(colnames(x), function(n) {
             transforms[[n]](x[,n])
            })
}

applyTransforms <- function(x, transforms) {
    n <- 1
    transforms <- transforms[colnames(x)]
    apply(x, 2, function(y) {
             y <- transforms[[n]](y)
             n <<- n+1
             return(y)
            })
}


load( '~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/pstat5-join/All/CB00366X_2012-11-07.RData' )
transforms <- estimateTransforms(fcs.data)
transforms[['FOXP3']] <- logicleTransform(w=1)
transforms[['CD45RA']] <- logicleTransform(w=.5)
transforms[['CD25']] <- logicleTransform(w=.6)
transforms[['PSTAT5.1']] <- logicleTransform(w=.5)
transforms[grep('PSTAT5',names(transforms))] <- transforms['PSTAT5.1'] 
transforms[['PSTAT5']] <- transforms[['PSTAT5.1']]
fcs.data <- applyTransforms(fcs.data, transforms) 
save(transforms, file='~nikolas/dunwich/Projects/IL2/PSTAT5-CD25-CD45RA-CD4-FOXP3/transforms.RData')




