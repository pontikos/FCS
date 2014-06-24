suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(flowCore))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(library(ellipse))

#set w=0 is important to avoid spurious peaks close to zero
lgcl <- logicleTransform(w=0.1)

###
normalised.density <- function(x) {
    d <- density(x)
    d$y <- d$y/sum(d$y)
    return(d)
}

# 
box <- function (data, lambda) {
    if (length(lambda) > 1 || lambda == 1) return(data)
    else if (length(lambda) > 1 || lambda != 0) return((sign(data) * abs(data)^lambda - 1)/lambda)
    else if (is.na(lambda) || is.null(lambda)) return(data)
    else return(log(data))
}

# to perform reverse box-cox transformation (multivariate)
rbox <- function (data, lambda) 
{
    if (length(lambda) > 1 || lambda == 1) return(data)
    else if (length(lambda) > 1 || lambda != 0) return(sign(lambda * data + 1) * (sign(lambda * data + 1) * (lambda * data + 1))^(1/lambda))
    else if (is.na(lambda) || is.null(lambda)) return(data)
    else return(exp(data))
}


prior.flowClust <- function(d, prior.res, B=500, level=.9, u.cutoff=.5, post.threshold=.99, kappa=1) {
    prior <- flowClust2Prior(prior.res,kappa=kappa) 
    return(flowClust::flowClust(d,K=prior.res@K,B=B,usePrior='yes',prior=prior,level=level,u.cutoff=u.cutoff,trans=0,lambda=1,control=list(B.lambda=0)))
}

fcs3.filter <- function(fcs) {
        #remove negative scatter
        fcs <- fcs[fcs@exprs[,'SSC-A']>1,]
        fcs <- fcs[fcs@exprs[,'FSC-A']>1,]
        #truncate side scatter at 500 
        fcs <- fcs[fcs@exprs[,'SSC-A']>500,]
        #truncate forward scatter at 100 
        fcs <- fcs[fcs@exprs[,'FSC-A']>100,]
        #truncate side scatter at 15000 
        #fcs <- fcs[fcs@exprs[,'SSC-A']>15000,]
        #fcs <- fcs[fcs@exprs[,'SSC-A']<100000,]
        #fcs <- fcs[fcs@exprs[,'FSC-A']>40000,]
        return(fcs)
}

fcs.filter <- function(fcs) {
    #remove min & max scatter
    min.ssc <- min(fcs@exprs[,'SSC-A'])
    max.ssc <- max(fcs@exprs[,'SSC-A'])
    fcs <- fcs[min.ssc < fcs@exprs[,'SSC-A'] & fcs@exprs[,'SSC-A'] < max.ssc,]
    min.fsc <- min(fcs@exprs[,'FSC-A'])
    max.fsc <- max(fcs@exprs[,'FSC-A'])
    fcs <- fcs[min.fsc < fcs@exprs[,'FSC-A'] & fcs@exprs[,'FSC-A'] < max.fsc,]
    return(fcs)
}


# Wrapper for the flowCore read.FCS function.
# Optionally applies compensation, transforms
# (this needs to be a parameter)
# 
# Returns
read.FCS <- function(file, EXPRS=TRUE, FILTER=FALSE, channels=c(), COMP=TRUE, TRANS=TRUE, verbose=FALSE, ...) {
	if (verbose) fcs <- flowCore::read.FCS(file, ...)
	else fcs <- suppressWarnings(flowCore::read.FCS(file, ...))
    params <- fcs@parameters
    pd <- pData(params)
    # Replace any null descs with names (for FSC-A, FSC-W, SSC-A)
    bad_col <- grep("^[a-zA-Z0-9]+",pd$desc,invert=TRUE)
	if (length(bad_col) > 0) {
		keyval <- keyword(fcs)
		for (i in bad_col) {
			pd$desc[i] <- pd$name[i]
			keyval[[paste("$P",i,"S",sep="")]] <- pd$name[i]
		}
		pData(params) <- pd
		fcs <- flowFrame(exprs(fcs),params,description=description(fcs))
		keyword(fcs) <- keyval
	}
#if FMO is present Marcin tells me this means he forgot a stain
    if ( 'FMO' %in% colnames(getChannels(fcs)) ) {
        warning('FMO channel! Returning NULL.\n')
        return(NULL)
    }
	# Compensate data if SPILL or SPILLOVER present, stripping compensation matrix 
	# out of the flowFrame, i.e we should only have to do this once
	apply.comp <- function(in_fcs, keyword) {
		comp_fcs <- compensate(in_fcs, description(in_fcs)[[keyword]])
		#flowFrame(exprs(comp_fcs), parameters(comp_fcs), description(comp_fcs)[grep("SPILL",names(description(comp_fcs)),invert=TRUE)])
		flowFrame(exprs(comp_fcs), comp_fcs@parameters, description(comp_fcs)[grep("SPILL",names(description(comp_fcs)),invert=TRUE)])
	} 
	if (COMP && !is.null(description(fcs)$SPILL)) {
		fcs <- apply.comp(fcs, "SPILL")
	} else if (COMP && !is.null(description(fcs)$SPILLOVER)) {
		fcs <- apply.comp(fcs, "SPILLOVER")
	}
    #some filtering happens here
    dim(fcs)
    if (fcs@description$FCSversion=='3') {
        if (FILTER) fcs <- fcs3.filter(fcs)
        for (channel in channels) {
            #i <- grep( channel, parameters(fcs)@data$desc, ignore.case=T )
            if (!TRANS) trans <- identity
            else if (channel == 'FSC-A') trans <- function(x) 5*x/262144
            else if (channel == 'SSC-A') trans <- log10
            else trans <- lgcl
            #fcs@exprs[,i] <- trans(fcs@exprs[,i])
            fcs <- setChannels(fcs, channel, trans(getChannels(fcs, channel)))
        }
    } else if (fcs@description$FCSversion=='2') {
        for (channel in channels) {
            #if (channel == 'FSC-A') trans <- function(x) 5*x/262144
            if (!TRANS) trans <- identity
            else if (channel == 'FSC-A') trans <- identity
            else if (channel == 'SSC-A') trans <-identity 
            else trans <- log10
            #fcs@exprs[,i] <- trans(fcs@exprs[,i])
            fcs <- setChannels(fcs, channel, trans(getChannels(fcs, channel)))
        }
    }
    if (FILTER) fcs <- fcs.filter(fcs)
    dim(fcs)
    if (EXPRS) return(getChannels(fcs, channels=channels))
    else return(fcs)
}



# Wrapper for the flowCore read.FCS function.
# Assumes that FSC-A and SSC-A always exits in file.
# Applies compensation
# Returns matrix
read.FCS.matrix <- function(file, channels=NULL, comp=TRUE, verbose=FALSE, ...) {
    if (!file.exists(file)) stop(paste(file, 'does not exist!\n'))
    if (verbose) print(dim(fcs <- flowCore::read.FCS(file, ...)))
    else fcs <- suppressWarnings(flowCore::read.FCS(file, ...))
    # Compensate data if SPILL or SPILLOVER present, stripping compensation matrix 
    # out of the flowFrame, i.e we should only have to do this once
    apply.comp <- function(in_fcs, keyword) {
        comp_fcs <- compensate(in_fcs, description(in_fcs)[[keyword]])
        #flowFrame(exprs(comp_fcs), parameters(comp_fcs), description(comp_fcs)[grep("SPILL",names(description(comp_fcs)),invert=TRUE)])
        flowFrame(exprs(comp_fcs), comp_fcs@parameters, description(comp_fcs)[grep("SPILL",names(description(comp_fcs)),invert=TRUE)])
    } 
    if (comp && !is.null(description(fcs)$SPILL)) {
        fcs <- apply.comp(fcs, "SPILL")
    } else if (comp && !is.null(description(fcs)$SPILLOVER)) {
        fcs <- apply.comp(fcs, "SPILLOVER")
    }

    params <- fcs@parameters
    pd <- pData(params)
    # Replace any null descs with names (for FSC-A, FSC-W, SSC-A)
    bad_col <- grep("^[a-zA-Z0-9]+",pd$desc,invert=TRUE)
    if (length(bad_col) > 0) {
        keyval <- keyword(fcs)
        for (i in bad_col) {
            pd$desc[i] <- pd$name[i]
            keyval[[paste("$P",i,"S",sep="")]] <- pd$name[i]
        }
        pData(params) <- pd
        fcs <- flowFrame(exprs(fcs),params,description=description(fcs))
        keyword(fcs) <- keyval
    }
    #if FMO is present Marcin tells me this means he forgot a stain
    if ( 'FMO' %in% colnames(getChannels(fcs)) ) {
        warning('FMO channel! Returning NULL.\n')
        return(NULL)
    }
    # only keep the data matrix
    X <- fcs@exprs
    colnames(X) <- gsub('[- ]', '.', toupper(fcs@parameters@data$desc))
    #remove min & max scatter
    min.ssc <- min(X[,'SSC.A'])
    max.ssc <- max(X[,'SSC.A'])
    X <- X[min.ssc < X[,'SSC.A'] & X[,'SSC.A'] < max.ssc,]
    min.fsc <- min(X[,'FSC.A'])
    max.fsc <- max(X[,'FSC.A'])
    X <- X[min.fsc < X[,'FSC.A'] & X[,'FSC.A'] < max.fsc,]
    # do the transforms based on the fcs file format
    if (fcs@description$FCSversion=='3') {
        #remove negative scatter
        X <- X[X[,'SSC.A']>1,]
        X <- X[X[,'FSC.A']>1,]
        #truncate side scatter at 500 
        X <- X[X[,'SSC.A']>500,]
        #truncate forward scatter at 100 
        fcs <- fcs[fcs@exprs[,'FSC.A']>100,]
        #truncate side scatter at 15000 
        for (channel in colnames(X)) {
            if (channel == 'FSC.A') trans <- function(x) 5*x/262144
            else if (channel == 'SSC.A') trans <- log10
            else if (channel == 'TIME') trans <- identity
            else trans <- lgcl
            X[,channel] <- trans(X[,channel])
        }
    } else if (fcs@description$FCSversion=='2') {
        for (channel in colnames(X)) {
            if (channel == 'FSC.A') trans <- identity
            else if (channel == 'SSC.A') trans <-identity 
            else if (channel == 'TIME') trans <- identity
            else trans <- log10
            X[,channel] <- trans(X[,channel])
        }
    }
    for (channel in colnames(X)) {
        x <- X[,channel]
        X <- X[min(x) < x & x < max(x),]
    }
    if (is.null(channels)) return(X)
    else return(X[,channels])
}


build.flowFrame <- function(x) {
	if (!is.matrix(x)) {
		stop("Input must be matrix")
  }
  # Build metadata for FCS file
  pd <- c()  # 'params' phenoData
  dl <- list()  # 'description' list
  dl[["$DATATYPE"]] <- "F"
  for (c in 1:ncol(x)) {
		c_name <- colnames(x)[c]
		c_min <- min(x[,c])
		c_max <- max(x[,c])
		c_rng <- c_max - c_min + 1
		pl <- matrix(c(c_name, c_name, c_rng, c_min, c_max),nrow=1)
		colnames(pl) <- c("name", "desc", "range", "minRange", "maxRange")
		rownames(pl) <- paste("$P",c,sep="") 
		pd <- rbind(pd, pl)
		dl[[paste("$P",c,"B",sep="")]] <- "32";	    # Number of bits
		dl[[paste("$P",c,"R",sep="")]] <- toString(c_rng); # Range
		dl[[paste("$P",c,"E",sep="")]] <- "0,0";	    # Exponent
		dl[[paste("$P",c,"N",sep="")]] <- c_name;	    # Name
		dl[[paste("$P",c,"S",sep="")]] <- c_name;	    # Desc	
	}	
  flowFrame(x, as(data.frame(pd), "AnnotatedDataFrame"), description=dl)
}



getIndividual <- function(x) x@description$INDIVIDUAL

getDose <- function(x) x@description$DOSE

getDate <- function(x) format(as.Date(dmy(x@description$`$DATE`, quiet=T)), '%Y-%m-%d')

getSSC <- function(x) x@exprs[,'SSC-A']
getFSC <- function(x) x@exprs[,'FSC-A']


matchChannelNames <- function(channels, fcs.channels) {
   matches <- c()
    channels <- tolower(channels)
    fcs.channels <- tolower(fcs.channels)
    for (chan in channels) {
        if (length(m <- which(chan == fcs.channels))) {
            matches <- c(matches, m)
        } else if (length(m <- grep(chan, fcs.channels,value=TRUE))) {
            m <- m[which.min(sapply(m, length))]
            matches <- c(matches, which(m==fcs.channels))
        } else {
            warning(paste('could not find', chan, 'in', paste(fcs.channels, collapse=',')))
        }
    }
    return(matches)
}


getChannels <- function(fcs.data, channels) {
    if (is(fcs.data,'matrix') || is(fcs.data,'data.frame')) {
        if (missing(channels)) channels <- colnames(fcs.data)
        channels <- gsub('-', '',paste(channels, ifelse(!grepl('-', channels), '$', ''), sep=''))
        fcs.channels <- gsub('-','',colnames(fcs.data))
        x <- as.matrix(fcs.data[,sapply(channels, function(x) grep(x, fcs.channels, ignore.case=T))])
    } else if (is(fcs.data, 'flowFrame')) {
        fcs.channels <- gsub('-','',fcs.data@parameters@data$desc)
        if (missing(channels)) {
            channels <- fcs.channels
            x <- as.matrix( fcs.data@exprs[,sapply(channels, function(x) match(x, fcs.channels))] )
        } else {
            matches <- matchChannelNames(channels, fcs.channels)
            x <- as.matrix(fcs.data@exprs[,matches])
        }
    } else {
        stop('Unrecognised type:', class(fcs.data))
    }
    colnames(x) <- gsub('\\$','',channels)
    x
}

exprs <- function(fcs.data) {
    if (class(fcs.data)=='matrix')
        return(fcs.data)
    else
        return(fcs.data@exprs)
}


setChannels <- function(fcs.data, channels, x) {
    #colnames(x) <- gsub('\\$','',channels)
    #channels <- gsub('-','',paste(channels, ifelse(!grepl('-', channels), '$', ''), sep=''))
    channels <- gsub('-','',channels)
    if (is(fcs.data,'matrix') || is(fcs.data,'data.frame')) {
        fcs.channels <- gsub('-','',colnames(fcs.data))
        fcs.data[,sapply(channels, function(x) grep( x, fcs.channels, ignore.case=T ))] <- x
    } else if (is(fcs.data, 'flowFrame')) {
        fcs.channels <- gsub('-','',fcs.data@parameters@data$desc)
        #fcs.data@exprs[,sapply(channels, function(x) grep( x, fcs.channels, ignore.case=T ))] <- x
        fcs.data@exprs[,matchChannelNames(channels,fcs.channels)] <- x
    }
    fcs.data
}


getChannelsf <- function(fcs.data, channels, f) {
    return( apply(getChannels(fcs.data, channels), 2, f) )
}

getMFI <- function(fcs.data, channels) mean(getChannels(fcs.data, channels))

update.FCS.description <- function(file.name, keywords, read=FALSE) {
    fcs.file <- flowCore::read.FCS(file.name)
    for (k in names(keywords)) {
        fcs.file@description[[k]] <- keywords[[k]]
    }
    write.FCS(fcs.file, filename=file.name)
    if (read) {
        read.FCS(fcs.file)
    }
}



downsample.indexes <- function(fcs.data, density, desired_samples=20000, verbose=TRUE) {
    if (is(fcs.data,'flowFrame')) fcs.data <- fcs.data@exprs
    if (missing(density)) density <- fcs.data$dens
    fcs.data <- cbind(fcs.data, density)
    # exclude_pctile
    boundary1 <- quantile(density,0.01,names=FALSE)
    out_data <- subset(fcs.data, density > boundary1) # Exclusion    
    if (desired_samples >= nrow(out_data)) return( which(density > boundary1) )
    # Need to find target density such there are approximately desired_samples
    # remaining after downsampling. To do so we solve for the density such that
    # the sum of samples below that density plus the expected value of
    # samples retained above that density equals approximately the desired
    # number of samples
    density_s <- sort( out_data[,'density'] )
    cdf       <- rev(cumsum(1.0/rev(density_s)))
    # Default solution if target density smaller than any present
    boundary2 <- desired_samples/cdf[1] 
    # Boundary actually falls amongst densities present
    if (boundary2 > density_s[1]) {
      targets <- (desired_samples-1:length(density_s)) / cdf 
      boundary2 <- targets[which.min(targets-density_s > 0)]
    }
    if (verbose) print(c(boundary1, boundary2))
    return( which( (boundary1 < density) & (boundary2/density > runif(length(density))) ) )
}



# Main scatterplot
marginal.2D.density <- function(fcs.data, channels, dir) {
    dir.create(dir, recursive=T, showWarnings=F)
    comb <- t(combn(channels, 2))
    for (i in 1:nrow(comb)) {
        channels <- comb[i,]
        print(channels)
        d <- data.frame(getChannels(fcs.data, channels))
        d$col <- ifelse('col' %in% colnames(fcs.data@exprs), as.factor(fcs.data@exprs[,'col']), 'black') 
        colnames(d) <- c('x', 'y') 
        p1 <- ggplot(d, aes(x, y)) + 
          xlab(channels[1]) +
          ylab(channels[2]) +
          scale_x_continuous(expand = c(0, 0)) +
          scale_y_continuous(expand = c(0, 0)) +
          expand_limits(y = c(min(d$x) - .1*diff(range(d$x)), max(d$x) + .1*diff(range(d$x)))) +
          expand_limits(x = c(min(d$y) - .1*diff(range(d$y)), max(d$y) + .1*diff(range(d$y)))) +
          theme(plot.margin= unit(c(.2, .2, .5, .5), "lines"), plot.background=element_rect(fill='white')) +
          #geom_point(data=d, colour='black') +
          geom_density2d(data=d, colour='black')
        # Horizontal marginal density plot - to appear at the top of the chart
        p2 <- ggplot(d, aes(x = x)) + 
          geom_density() +
          scale_x_continuous(expand = c(0, 0)) +
          expand_limits(x = c(min(d$y) - .1*diff(range(d$y)), max(d$y) + .1*diff(range(d$y)))) +
          theme(axis.text = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_blank(),
                plot.margin= unit(c(1, .2, -.5, 0.5), "lines"),
                plot.background=element_rect(fill='white')) 
        # Vertical marginal density plot - to appear at the right of the chart
        p3 <- ggplot(d, aes(x = y)) + 
          geom_density() +
          scale_x_continuous(expand = c(0, 0)) +
          expand_limits(x = c(min(d$x) - .1*diff(range(d$x)), max(d$x) + .1*diff(range(d$x)))) +
          coord_flip() +
          theme(axis.text = element_blank(),
                axis.title = element_blank(),
                axis.ticks = element_blank(),
                plot.margin= unit(c(.2, 1, 0.5, -.5), "lines"),
                plot.background=element_rect(fill='white')) 
        # Get the gtables
        gt1 <- ggplot_gtable(ggplot_build(p1))
        gt2 <- ggplot_gtable(ggplot_build(p2))
        gt3 <- ggplot_gtable(ggplot_build(p3)) 
        # Get maximum widths and heights for x-axis and y-axis title and text
        maxWidth = unit.pmax(gt1$widths[2:3], gt2$widths[2:3])
        maxHeight = unit.pmax(gt1$heights[4:5], gt3$heights[4:5]) 
        # Set the maximums in the gtables for gt1, gt2 and gt3
        gt1$widths[2:3] <- as.list(maxWidth)
        gt2$widths[2:3] <- as.list(maxWidth) 
        gt1$heights[4:5] <- as.list(maxHeight)
        gt3$heights[4:5] <- as.list(maxHeight) 
        # Combine the scatterplot with the two marginal boxplots
        # Create a new gtable
        gt <- gtable(widths = unit(c(7, 2), "null"), height = unit(c(2, 7), "null")) 
        # Instert gt1, gt2 and gt3 into the new gtable
        gt <- gtable_add_grob(gt, gt1, 2, 1)
        gt <- gtable_add_grob(gt, gt2, 1, 1)
        gt <- gtable_add_grob(gt, gt3, 2, 2) 
        # And render the plot
        pdf(file.path(dir, sprintf('%s.pdf',paste(channels, collapse='-'))))
        print(grid.newpage())
        print(grid.draw(gt))
        print(grid.rect(x = 0.5, y = 0.5, height = .995, width = .995, default.units = "npc", gp = gpar(col="black", fill = NA, lwd = 1)) ) 
        dev.off()
    }
}


### plot clusters for all channels in panel
plot.clusters <- function(x, clustering, panel, pdf.filename, width=15, height=15) {
    x <- data.frame(getChannels(x, getPanelChannels(panel)))
    x$cluster <- as.character(clustering)
    x <- na.omit(x)
    pdf(pdf.filename, width=width, height=height)
    print(ggpairs(x, columns=getPanelChannels(panel), colour='cluster', lower=list(continuous='density'), axisLabels='none', upper=list(continuous='blank')))
    dev.off()
}


pairs.smoothScatter <- function(x,hdr=FALSE) print( pairs(x, lower.panel=NULL, upper.panel = function(x,y, ...) {
                                                smoothScatter(x,y, ..., nrpoints = 0, add = TRUE)
                                                if(hdr) {
                                                    X <- cbind(x,y)
                                                    fs <- featureSignif(X)
                                                    points(X[fs$curvData,],pch=20,col='red')
                                                } }) )


#pairs.featureSignif <- function(x) print( pairs(x, panel = function(x,y,...) { par(new=TRUE); plot(featureSignif(cbind(x,y))); }) )
                                               

###
plotClusters <- function(x, plot.file=NULL, channels=NULL, classification=NULL, posteriors=NULL, posterior.cutoff=.95, uncertainty=NULL, uncertainty.cutoff=.95) {
    if (is.null(channels)) {
        col.names <- colnames(x)
        nc <- ncol(x)
    } else {
        col.names <- channels
        nc <- length(channels)
    }
    cat('>>',plot.file,'\n')
    if (!is.null(plot.file)) png(plot.file)
    par(mfrow = c(nc, nc), mai=c(.2,.2,.2,.2), oma=c(2,3,2,.5))
    for (n in 1:length(col.names)) {
            ylab <- ''
            xlab <- ''
            main <- ''
        for (n2 in 1:length(col.names)) {
            #off-diagonal terms containing bivariate densities
            if (n != n2) {
               dim.inds <- c(n2,n)
               smoothScatter(x[,dim.inds],main=main,xlab=xlab,ylab=ylab)
               if (!is.null(classification))
               for (k in sort(unique(classification))) {
                   X1 <- x[which(classification==k),dim.inds]
                   p <- X1[chull(X1),]
                   p <- rbind(p, p)
                   lines(p, col=k, lwd=.5)
               } 
               if (!is.null(posteriors))
               for (k in 1:ncol(posteriors)) {
                   X1 <- x[posteriors[,k]>posterior.cutoff,dim.inds]
                   ch <- c(chull(X1),chull(X1))
                   p <- X1[ch,]
                   #p <- rbind(p, p)
                   lines(p, col=k, lwd=.5)
                   #also approximate with an ellipse
                   lines(ellipse::ellipse(cov(X1),centre=colMeans(X1)),col=k,lwd=2,lty=1)
               } 
            #diagonal contains univariate densities
            } else {
               d <- normalised.density(x[,n])
               plot(d,main=main,xlab=xlab,ylab=ylab)
               if (!is.null(classification))
               for (k in sort(unique(classification))) {
                   X1 <- x[which(classification==k),n]
                   if (length(X1)>2) lines(normalised.density(X1), col=k, lwd=2)
               } 
               if (!is.null(posteriors))
               for (k in 1:ncol(posteriors)) {
                   X1 <- x[posteriors[,k]>posterior.cutoff,n]
                   if (length(X1)>2) lines(normalised.density(X1), col=k, lwd=2)
               } 
            }
            if (n==1) mtext(col.names[[n2]], side=3, cex=2, line=1) 
            if (n2==1) mtext(col.names[[n]], side=2, cex=2, line=2)
       }
    }
    if (!is.null(plot.file)) dev.off()
}



###
plotFlowClustRes <- function(x, res, K=1:res@K, outliers=FALSE, plot.file=NULL, channels=NULL, col=1:res@K) {
    x <- box(x, res@lambda)
    if (is.null(channels)) {
        col.names <- colnames(x)
        nc <- ncol(x)
    } else {
        col.names <- channels
        nc <- length(channels)
    }
    if (!outliers) x <- x[which(!res@flagOutliers),]
    print(table(clustering <- MAP(x, res)))
    cat('>>',plot.file,'\n')
    if (!is.null(plot.file)) png(plot.file)
    #mar <- rep(3,4)
    #opar <- par(mfrow = c(nc, nc), mar = rep.int(1/2, 4), oma = oma)
    par(mfrow = c(nc, nc), mai=c(.2,.2,.2,.2), oma=c(2,3,2,.5))
    for (n in 1:length(col.names)) {
            ylab <- ''
            xlab <- ''
            main <- ''
        for (n2 in 1:length(col.names)) {
            if (n != n2) {
                dim.inds <- c(n2,n)
                smoothScatter(x[,dim.inds],main=main,xlab=xlab,ylab=ylab)
                for (i in K) {
                    points(t(as.matrix(res@mu[i,dim.inds])),pch=20,col=col[i])
                    mu <- res@mu[i,dim.inds]
                    sigma <- res@sigma[i,dim.inds,dim.inds]
                    lines(ellipse::ellipse(sigma,centre=mu),col=col[i],lwd=2,lty=1)
                }
                #for(i in 1:prior$K){
                    #points(t(as.matrix(prior$Mu0[i,dim.inds])),pch=20,col=i)
                    #lines(ellipse(prior$Lambda0[i,dim.inds,dim.inds]/(prior$nu0[i]-nc-1),centre=prior$Mu0[i,dim.inds]),col=i,lwd=2,lty=2)
                #}
            #diagonal contains univariate densities
            } else {
                d <- normalised.density(x[,n])
                plot(d,main=main,xlab=xlab,ylab=ylab)
                for (i in K) {
                    cl <- which(clustering==i)
                    if (length(cl)>2) lines(normalised.density(x[cl,n]), col=col[i], lwd=2)
                }
            }
            if (n==1) mtext(col.names[[n2]], side=3, cex=2, line=1) 
            if (n2==1) mtext(col.names[[n]], side=2, cex=2, line=2)
            #if (n2==length(col.names)) mtext(col.names[[n]], side=4, cex=2, line=2)
       }
    }
    if (!is.null(plot.file)) dev.off()
}

###
plotMClustRes <- function(x, res, K=1:res$G, plot.file=NULL) {
    col.names <- colnames(x)
    nc <- ncol(x)
    cat('>>',plot.file,'\n')
    if (!is.null(plot.file)) png(plot.file)
    print(table(clustering <- res$classification))
    oma <- rep(1,4)
    #mar <- rep(3,4)
    #opar <- par(mfrow = c(nc, nc), mar = rep.int(1/2, 4), oma = oma)
    par(mfrow = c(nc, nc))
    for (n in 1:length(col.names)) {
        for (n2 in 1:length(col.names)) {
            if (n != n2) {
                dim.inds <- c(n2,n)
                smoothScatter(x[,dim.inds])
                print(dim.inds)
                for (i in K) {
                    points(t(as.matrix(res$parameters$mean[dim.inds,i])),pch=20,col=i)
                    lines(ellipse(res$parameters$variance$sigma[dim.inds,dim.inds,i],centre=t(res$parameters$mean[dim.inds,i])),col=i,lwd=2,lty=1)
                }
            #diagonal contains univariate densities
            } else {
                d <- normalised.density(x[,n])
                plot(d,main=col.names[[n]])
                for (i in K) lines(normalised.density(x[which(clustering==i),n]), col=i, lwd=2)
            }
       }
    }
    if (!is.null(plot.file)) dev.off()
}


### density at each point for each flowClust component
compute.dens <- function(d, res) sapply(1:res@K, function(i) res@w[i]*flowClust::dmvt(d, mu=res@mu[i,], sigma=res@sigma[i,,], nu=res@nu, lambda=res@lambda)$value)

### outliers are points which follow below a given quantile of density
outliers <- function(d, res, dens.threshold=.05) {
    dens <- compute.dens(d, res)
    quants <- quantile(dens, probs=seq(0,1,dens.threshold))
    return(dens<quants[[2]])
}

### outliers are points which lie beyond a certain distance from the 

### posterior of each point for each flowClust component
compute.post <- function(d, res) {
    dens <- compute.dens(d, res)
    total.dens <- rowSums(dens)
    dens/total.dens 
}

### return the component assignment
MAP <- function(d, res, post.threshold=NULL) {
    d <- box(d, res@lambda)
    post <- compute.post(d,res)
    if (is.null(post.threshold))
        apply(post,1,which.max)
    else
        apply(post,1,function(x)ifelse(length(which(x>post.threshold))>0,which.max(x),NA))
}


###
plotDensity <- function(x, k=NULL) {
print( pairs(x, panel = function(x,y,...) {
             X <- cbind(x,y)
             #kde2D(X
       smoothScatter(x,y, ..., nrpoints = 0, add = TRUE)
       X <- cbind(x,y)
       if (!is.null(k)) {
           for (i in 1:max(k)){
               X1 <- X[k==i,]
               p <- X1[chull(X1),]
               p <- c(p, p[[1]])
               lines(p, col=i)
            }
       }
       #points(X, col=k, pch='.')
      }))
}


#dens <- apply(t(combn(colnames(x),2)), 1, function(y) kde2D(x[,y]))

dens <- function(x,bw) {
    o<-apply(x,2,order)
    x<-apply(x,2,sort)
    d <- apply(x,2, function(x) {
        count <- c()
        i <- 1
        while (i <= length(x)) {
            if (1 < i && x[i]-x[i-1] < bw**10) {
                count <- c(count,count[length(count)])
            } else {
                #forwards
                j <- k <- i
                while( j <= length(x) && x[j]-x[i] < bw ) j<-j+1
                #backwards
                while( 1 <= k && x[i]-x[k] < bw ) k<-k-1
                count <- c(count,j-k)
            }
            (i<-i+1)
        }
        return(count)
        })
    d <- t(t(sapply(1:ncol(x), function(i) d[o[,i],i])))
    colnames(d) <- colnames(x)
    return(sweep(d, 2, colSums(d), '/'))
    #return(d)
}


dens <- function(x) {
    d <- apply(x,2, function(x) {
          d <- density(x)
          f <- splinefun(d$x,d$y)
          f(x)
        })
    m <- apply(sweep(d, 2, colSums(d), '/'), 1, function(i) colnames(x)[which.max(i)])
    #apply(sweep(d, 2, colSums(d), '/'), 1, function(i) mean(i)/var(i))
    w <- prop.table(table(m))
    w <- 1/w
    w <- w/sum(w)
    d <- t(apply(d, 1, function(x) x*w))
    sweep(d, 2, colSums(d), '/')
}



dens <- function(x, bw=.1, n=512) {
    d <- apply(x, 2, function(x) density(x,bw,n)$y)
    return(sweep(d, 2, colSums(d), '/'))
}






