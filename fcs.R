suppressPackageStartupMessages(library(sp))
suppressPackageStartupMessages(library(cluster))
suppressPackageStartupMessages(library(KernSmooth))
suppressPackageStartupMessages(library(iterators))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(flowCore))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(GGally))
suppressPackageStartupMessages(library(ellipse))
suppressPackageStartupMessages(library(mvtnorm))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(mixtools))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(flowClust))


applyTransforms <- function(x, transforms) {
    n <- 1
    transforms <- transforms[colnames(x)]
    apply(x, 2, function(y) {
             y <- transforms[[n]](y)
             n <<- n+1
             return(y)
            })
}



trans <- function(x, m=4.5) {
    cbind(x[,grep('^SSC|^FSC',colnames(x))], apply(x[,grep('^SSC|^FSC',colnames(x),invert=TRUE)], 2,
                                function(y) {
                                    r <- log10(max(y)/abs(min(y)))
                                    w <- (m-r)/2
                                    logicleTransform(w=)(y)
                                }))
}

colMedians <- function(x) apply(x,2,median)
colMin <- function(x) apply(x,2,min)
colMax <- function(x) apply(x,2,max)
colQuantile <- function(x,prob) apply(x,2,quantile,prob)

###
plot.normalised.density <- function(d, include=c(0.1,.9), ...) {
    plot(d, xlim=quantile(d$x,probs=include), ...)
}

###
normalised.density <- function(x, ...) {
    d <- density(x, ...)
    d$y <- d$y/sum(d$y)
    class(d) <- 'normalised.density'
    return(d)
}

# flowClust::box
box <- function (data, lambda) {
    if (length(lambda) > 1 || lambda == 1) return(data)
    else if (length(lambda) > 1 || lambda != 0) return((sign(data) * abs(data)^lambda - 1)/lambda)
    else if (is.na(lambda) || is.null(lambda)) return(data)
    else return(log(data))
}
#reassignInPackage('box', pkgName='flowClust', rbox)

# to perform reverse box-cox transformation (multivariate)
# flowClust::rbox
rbox <- function (data, lambda) 
{
    if (length(lambda) > 1 || lambda == 1) return(data)
    else if (length(lambda) > 1 || lambda != 0) return(sign(lambda * data + 1) * (sign(lambda * data + 1) * (lambda * data + 1))^(1/lambda))
    else if (is.na(lambda) || is.null(lambda)) return(data)
    else return(exp(data))
} 
#reassignInPackage('rbox', pkgName='flowClust', rbox)


prior.flowClust <- function(d, prior.res, B=500, level=.9, u.cutoff=.5, post.threshold=.99, kappa=1) {
    prior <- flowClust2Prior(prior.res,kappa=kappa) 
    return(flowClust::flowClust(d,K=prior.res@K,B=B,usePrior='yes',prior=prior,level=level,u.cutoff=u.cutoff,trans=0,lambda=1,control=list(B.lambda=0)))
}

# findInterval returns -1 if before; 1 if within range; 2 if greater
in.range <- function(x, interval, all.inside=TRUE, ...) findInterval(x, interval, all.inside, ...)==1

#`%in%` <- function(x, r) in.range(x,r)

#
percentile.filter <- function(fcs.data, bottom=.01, top=.99) {
    f <- function(x) in.range(x, quantile(x, probs=c(bottom,top)))
    if(is.null(ncol(fcs.data)))
        f(fcs.data)
    else 
        apply( apply( fcs.data, 2, f), 1, all )
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
read.FCS <- function(file, EXPRS=TRUE, FILTER=FALSE, channels=c(), COMP=TRUE, verbose=FALSE, TRANS=logicleTransform(w=0.1), ...) {
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
            if (is.null(TRANS)) trans <- identity
            else if (channel == 'FSC-A') trans <- function(x) 5*x/262144
            else if (channel == 'SSC-A') trans <- log10
            else trans <- TRANS
            #fcs@exprs[,i] <- trans(fcs@exprs[,i])
            fcs <- setChannels(fcs, channel, trans(getChannels(fcs, channel)))
        }
    } else if (fcs@description$FCSversion=='2') {
        for (channel in channels) {
            #if (channel == 'FSC-A') trans <- function(x) 5*x/262144
            if (is.null(TRANS)) trans <- identity
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

#returns a matrix
read.flow.file <- function(filename,channels) {
  cat('File extension', file.ext <- tolower(file_ext(filename)), '\n')
  if (file.ext=='fcs') {
      #fcs file
      print(head(fcs.data <- getChannels(read.FCS(filename, channels=channels),channels)))
      #My version of read.FCS discards events on axis if filtered, so total count is less that what is in FCS file initially.
      print(total.count <- as.numeric(read.FCSheader(filename,keyword='$TOT')))
  } else if (file.ext=='csv') {
      #csv file
      print(head(fcs.data <- read.csv(filename)[,channels]))
      print(total.count <- dim(fcs.data)[[1]])
  } else if (file.ext=='rdata') {
      print(load(filename))
      if (exists('result'))
          print(head( fcs.data <- result$d.original[result$lymph.filter,][result$singlet.filter,channels] ))
      else
          print(head( fcs.data <- fcs.data[,channels] ))
      print(total.count <- dim(fcs.data)[[1]])
  } else {
      stop('Unsupported file extension', file.ext,'\n')
  }
  return(fcs.data)
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

getDate <- function(x) {
    if (class(x)=='flowCore')
        format(as.Date(dmy(x@description$`$DATE`, quiet=T)), '%Y-%m-%d')
    else if (class(x)=='character')
        format(as.Date(dmy(flowCore::read.FCS(x)@description$`$DATE`, quiet=T)), '%Y-%m-%d')
}

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
plotClusters1D <- function(x, classification=NULL, fun=normalised.density) {
   plot(fun(x))
   if (!is.null(classification))
   for (k in sort(unique(classification))) {
       X1 <- x[which(classification==k)]
       if (length(X1)>10) lines(fun(X1), col=k, lwd=2)
       else points(cbind(X1,0), col=k, pch=20)
   } 
}

                                               

###
### Can draw the convex hull or the ellipse (based on the cluster covariance) for a cluster.
plotClusters <- function(x, plot.points=NULL, plot.points.col='black', plot.points.pch=20, plot.file=NULL, classification=NULL, posteriors=NULL, posterior.cutoff=.95, outliers=FALSE, ellipses=TRUE, clusters.col=NULL, chulls=TRUE, chull.lwd=.5, upper=smoothPlot, lower=NULL, PAR=TRUE, ...) {
    if (class(x)=='numeric') {
        par(mfrow=c(1,1))
        plotClusters1D(x, classification, ...)
        return(NULL)
    }
    col.names <- colnames(x)
    nc <- ncol(x)
    cat('>>',plot.file,'\n')
    if (!is.null(plot.file)) png(plot.file)
    if (PAR) par(mfrow = c(nc, nc), mai=c(.2,.2,.2,.2), oma=c(2,3,2,.5))
    # plots the clusters either as chulls or as ellipses
    f <- function(dim.inds) {
           if(is.null(classification)) cols <- 'black' else cols <- classification
           if (!is.null(plot.points)) points(plot.points[,col.names[dim.inds]], col=plot.points.col, pch=plot.points.pch)
           if (!is.null(classification))
           for (k in sort(unique(classification))) {
               if (is.null(clusters.col)) col <- k
               else col <- clusters.col[[k]]
               X1 <- x[which(classification==k),dim.inds]
               if (chulls) {
                   p <- X1[chull(X1),]
                   p <- rbind(p, p)
                   lines(p, col=col, lwd=chull.lwd)
               }
               #also approximate with an ellipse
               if (ellipses) lines(ellipse::ellipse(cov(X1),centre=colMeans(X1)),col=col,lwd=2,lty=1)
           } 
           if (!is.null(posteriors))
           for (k in 1:ncol(posteriors)) {
               if (is.null(clusters.col)) col <- k
               else col <- clusters.col[[k]]
               i <- which(posteriors[,k]>posterior.cutoff)
               if (length(i)>2) {
                   X1 <- x[i,dim.inds]
                   if (chulls) {
                       #make sure chull is not open
                       ch <- c(chull(X1),chull(X1))
                       p <- X1[ch,]
                       lines(p, col=col, lwd=chull.lwd)
                   }
                   #also approximate with an ellipse
                   if (ellipses) lines(ellipse::ellipse(cov(X1),centre=colMeans(X1)),col=col,lwd=2,lty=1)
               }
           } 
    }
    for (n in 1:length(col.names)) {
            ylab <- ''
            xlab <- ''
            main <- ''
        for (n2 in 1:length(col.names)) {
            dim.inds <- c(n2,n)
            #upper diagonal terms containing bivariate densities
            if ( n < n2 ) {
               if (is.null(upper)) {
                   frame()
               } else {
                   upper(x[,dim.inds],main=main,xlab=xlab,ylab=ylab, outliers=outliers, classification=classification, ellipses=ellipses, clusters.col=clusters.col, chulls=chulls, ...)
                   f(dim.inds)
               }
            #lower diagonal terms containing bivariate densities
            } else if ( n > n2 ) {
                if (is.null(lower)) {
                    frame()
                } else {
                   lower(x[,dim.inds],main=main,xlab=xlab,ylab=ylab, outliers=outliers,  clusters.col=clusters.col, ...)
                   f(dim.inds)
                }
            #diagonal contains univariate densities
            # } else if (n == n2) {
            } else {
               #diagonal(x[,n])
               q1 <- quantile(x[,n],probs=c(0.01,0.99))
               d <- normalised.density(x[,n])
               if (!outliers)
                   plot(d,main=main,xlab=xlab,ylab=ylab,xlim=c(q1[['1%']],q1[['99%']]))
               else
                   plot(d,main=main,xlab=xlab,ylab=ylab)
               if (!is.null(classification))
               for (k in sort(unique(classification))) {
                   X1 <- x[which(classification==k),n]
                   if (length(X1)>10) lines(normalised.density(X1), col=k, lwd=2)
                   else points(cbind(X1,0), col=k, pch=20)
               } 
               if (!is.null(posteriors))
               for (k in 1:ncol(posteriors)) {
                   X1 <- x[posteriors[,k]>posterior.cutoff,n]
                   if (is.null(clusters.col)) col <- k
                   else col <- clusters.col[k]
                   if (length(X1)>10) lines(normalised.density(X1), col=col, lwd=2)
                   else points(cbind(X1,0), col=col, pch=20)
               } 
            }
            if (n==1) mtext(col.names[[n2]], side=3, cex=2, line=1) 
            if (n2==1) mtext(col.names[[n]], side=2, cex=2, line=2)
       }
    }
    if (!is.null(plot.file)) dev.off()
    if (PAR) par(mfrow=c(1,1))
}

### Supports flowClust and Mclust res.
### Can apply model to new data.
### @parameter 
plotClustRes <- function(x, res, before.res=NULL, outliers=TRUE, plot.file=NULL, cex=2) {
  print(col.names <- colnames(x))
  nc <- ncol(x)
  if (class(res)=='Mclust') {
    K <- 1:res$G
    p <- predict(res, x)
    print(table(clustering <- p$classification))
  } else if (class(res)=='flowClust') {
    K <- 1:res@K
    x <- box(x, res@lambda)
    #col.names <- colnames(x)
    print(col.names <- res@varNames)
    nc <- length(col.names)
    x <- x[,col.names]
    if (!outliers) x <- x[which(!res@flagOutliers),]
    print(head(clustering <- MAP(x, res)))
    if (!is.null(before.res)) print(head(before.clustering <- MAP(x, before.res)))
  }
  cat('>>',plot.file,'\n')
  if (!is.null(plot.file)) match.fun(file_ext(plot.file))(plot.file)
  par(mfrow = c(nc, nc), mai=c(.2,.2,.2,.2), oma=c(2,3,2,.5))
  for (n in 1:length(col.names)) {
            ylab <- ''
            xlab <- ''
            main <- ''
     for (n2 in 1:length(col.names)) {
            #off-diagonal terms containing bivariate densities
            if (n < n2) {
                dim.inds <- c(n2,n)
                smoothScatter(x[,dim.inds],main=main,xlab=xlab,ylab=ylab,colramp=colorRampPalette(c('white','blue','green','yellow','orange','red')))
                if (class(res)=='Mclust') {
                    for (i in K) {
                        points(t(as.matrix(res$parameters$mean[dim.inds,i])),pch=20,col=i)
                        lines(ellipse::ellipse(res$parameters$variance$sigma[dim.inds,dim.inds,i],centre=t(res$parameters$mean[dim.inds,i])),col=i,lwd=2,lty=1)
                    }
                } else if (class(res)=='flowClust') {
                    if(!is.null(before.res))
                    for (i in K) {
                        points(t(as.matrix(before.res@mu[i,dim.inds])),pch=20,col=i)
                        mu <- before.res@mu[i,dim.inds]
                        sigma <- before.res@sigma[i,dim.inds,dim.inds]
                        lines(ellipse::ellipse(sigma,centre=mu),col=i,lwd=2,lty=2)
                    }
                    for (i in K) {
                        points(t(as.matrix(res@mu[i,dim.inds])),pch=20,col=i)
                        mu <- res@mu[i,dim.inds]
                        sigma <- res@sigma[i,dim.inds,dim.inds]
                        lines(ellipse::ellipse(sigma,centre=mu),col=i,lwd=2,lty=1)
                    }
                }
            } else if (n > n2) {
                plot(  NULL, xlim=range(x[,n2]), ylim=range(x[,n]), xaxt='n', yaxt='n', ann=FALSE, bty='n' )
            #diagonal contains univariate densities
            } else {
                d <- normalised.density(x[,n])
                plot(d,main=main,xlab=xlab,ylab=ylab)
                if (class(res)=='Mclust') {
                    for (i in K) {
                        lines(normalised.density(x[which(clustering==i),n]), col=i, lwd=2)
                    }
                } else if (class(res)=='flowClust') {
                    for (i in K) {
                        cl <- which(clustering==i)
                        if (length(cl)>2) lines(normalised.density(x[cl,n]), col=i, lwd=2)
                    }
                    if (!is.null(before.res))
                    for (i in K) {
                        cl <- which(before.clustering==i)
                        if (length(cl)>2) lines(normalised.density(x[cl,n]), col=i, lwd=2, lty=2)
                    }
                }
            }
            if (n==1) mtext(col.names[[n2]], side=3, cex=cex, line=1) 
            if (n2==1) mtext(col.names[[n]], side=2, cex=cex, line=2)
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
    d <- d[,res@varNames]
    d <- box(d, res@lambda)
    print(head(post <- compute.post(d,res)))
    if (is.null(post.threshold))
        apply(post,1,which.max)
    else
        apply(post,1,function(x)ifelse(length(which(x>post.threshold))>0,which.max(x),NA))
}



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


filter.mahalanobis <- function(x, mu, sigma, p=.999) {
    d <- mahalanobis(x,mu,sigma)
    return(d<qchisq(p, df=length(mu)))
}

filter.mahalanobis.ellipse <- function(x, e, p=.999) {
    d <- mahalanobis(x[,names(e$Mu)],e$Mu,e$Sigma)
    return(d<qchisq(p, df=length(e$Mu)))
}

mahalanobis.magnetic.ellipse <- function(x, e, p=.999) {
    #Mu is updated 3 times
    d <- mahalanobis(x[,names(e$Mu)],e$Mu,e$Sigma)
    e$Mu <- colMeans(x[which(d<qchisq(p, df=length(e$Mu))),names(e$Mu)])
    d <- mahalanobis(x[,names(e$Mu)],e$Mu,e$Sigma)
    e$Mu <- colMeans(x[which(d<qchisq(p, df=length(e$Mu))),names(e$Mu)])
    d <- mahalanobis(x[,names(e$Mu)],e$Mu,e$Sigma)
    e$Mu <- colMeans(x[which(d<qchisq(p, df=length(e$Mu))),names(e$Mu)])
    return(e)
}

kmeans.magnetic.ellipses <- function(x, ellipses) {
    #repeat 3 times
    (centers <- do.call('rbind', lapply(ellipses, function(e) e$Mu)))
    res <- kmeans(x[,colnames(centers)], centers=centers)
    for (i in 1:nrow(centers)) ellipses[[i]]$Mu <- res$centers[i,]
    return(ellipses)
}


em.magnetic.ellipses <- function(x, gates, iterations=10, background=NULL) {
    tau <- sapply(gates, function(x) x$tau)
    for (j in 1:length(gates)) gates[[j]]$tau <- tau[[j]]/sum(tau)
    if (!is.null(background)) gates[[length(gates)+1]] <- list(Mu=colMeans(x),Sigma=cov(x),tau=background)
    for (i in 1:iterations) {
        #u <- 1/mahalanobis(
        dens <- sapply(gates, function(e) e$tau*mixtools::dmvnorm(x[,names(e$Mu)], mu=e$Mu, sigma=e$Sigma))
        post <- t(apply(dens, 1, function(x) x/sum(x)))
        for (j in 1:length(gates)) gates[[j]]$Mu <- apply(x[,names(gates[[j]]$Mu)], 2, weighted.mean, w=post[,j]) 
        tau <- colSums(post)/sum(post)
        for (j in 1:length(gates)) gates[[j]]$tau <- tau[[j]]
    }
    return(gates)
}



filter.mahalanobis.ellipses <- function(x, ellipses, p=.999) {
    rowSums(sapply( ellipses, function(e) mahalanobis(x[,names(e$Mu)],e$Mu,e$Sigma)<qchisq(p, df=length(e$Mu))))==length(ellipses)
}

#fcs.data needs to be 2D
classification.to.chull <- function(fcs.data, clr) {
    clr <- which(as.logical(clr))
    d <- fcs.data[clr,]
    return(d[chull(d),])
}

classification.to.ellipse <- function(fcs.data, clr) {
    clr <- which(as.logical(clr))
    d <- fcs.data[clr,]
    return(list(Mu=colMeans(d), Sigma=cov(d), tau=nrow(d)/nrow(fcs.data)))
}

points.in.ellipse <- function(X, w, npoints=255, alpha=.05) {
    e <- classification.to.ellipse( X, w )
    mu <- e$Mu
    sigma <- e$Sigma
    es <- eigen(sigma)
    e1 <- es$vec %*% diag(sqrt(es$val))
    r1 <- sqrt(qchisq(1 - alpha, 2))
    theta <- seq(0, 2 * pi, len = npoints)
    v1 <- cbind(r1 * cos(theta), r1 * sin(theta))
    pts = t(mu - (e1 %*% t(v1)))
    point.in.polygon(X[,1],X[,2],pts[,1],pts[,2])
}


classification.to.gate <- function(fcs.data, clr) {
    e <- classification.to.ellipse(fcs.data, clr)
    return(which(filter.mahalanobis.ellipse(fcs.data, e)))
}


plot.ellipses <- function(gates, col='black') {
    for (e in gates) lines(ellipse::ellipse(e$Sigma,centre=e$Mu),col=col,lwd=2,lty=2)
    points(do.call('rbind', lapply(gates,function(x)x$Mu)), pch='X', cex=2, col=col)
}


plot.gate.chull <- function(X, classification) {
    for (k in sort(unique(classification))) {
         X1 <- X[which(classification==k),]
         p <- X1[chull(X1),]
         p <- rbind(p, p)
         lines(p, col=k)
    } 
}

#also approximate with an ellipse
plot.gate.ellipse <- function(X, classification) {
    for (k in sort(unique(classification))) {
         X1 <- X[which(classification==k),]
         lines(ellipse::ellipse(cov(X1),centre=colMeans(X1)),col=k,lwd=2,lty=1)
    } 
}


### TODO
plotManualGates <- function(fcs.data, clr)  {
    ngates <- length(gatenames <- colnames(clr))
    figure.labels <- iter(paste(letters,')',sep=''))
    par(mfrow=c(round(sqrt(ngates)),round(sqrt(ngates))))
    for (gn in gatenames) {
        print(gn)
        print(name <- unlist(strsplit(gn, '_'))[[1]])
        print(channels <- unlist(strsplit(unlist(strsplit(gn, '_'))[[2]],'\\.')))
        print(head(fcs.data[,channels]))
        xquant <- quantile(fcs.data[,channels[[1]]],probs=seq(0,1,.01))
        yquant <- quantile(fcs.data[,channels[[2]]],probs=seq(0,1,.01))
        print(xlim <- c(xquant[['1%']],xquant[['99%']]))
        print(ylim <- c(yquant[['1%']],yquant[['99%']]))
        smoothScatter( fcs.data[,channels], xlab=channels[[1]], ylab=channels[[2]], main=name, colramp=colorRampPalette(c('white','blue','green','yellow','orange','red')), nrpoints=0, xlim=xlim, ylim=ylim)
        X <- fcs.data[which(as.logical(clr[,gn])),channels]
        p <- X[chull(X),]
        p <- rbind(p, p)
        #lines(p, col='red', lwd=2)
        lines(ellipse::ellipse(cov(X),centre=colMeans(X)),col='red',lwd=2,lty=1)
        title(nextElem(figure.labels), adj=0)
    }
}

### normalised.density plot
smoothPlot1D <- function( x, outliers=FALSE, ... ) {
    xquant <- quantile(x,probs=c(0.01,0.99))
    print(xlim <- c(xquant[['1%']],xquant[['99%']]))
    if (!outliers)
        plot( normalised.density(x), cex.lab=2, cex.main=2, xlim=xlim, ... )
    else
        plot( normalised.density(x), cex.lab=2, cex.main=2, ... )
}


###
smoothPlot <- function( fcs.data, nrpoints=0, colramp=colorRampPalette(c('white','blue','green','yellow','orange','red')), plot.points=NULL, plot.points.col='black', plot.file=NULL, classification=NULL, posteriors=NULL, posterior.cutoff=.95, outliers=FALSE, ellipses=TRUE, clusters.col=NULL, chulls=TRUE, chull.lwd=.5,  ... ) {
    xquant <- quantile(fcs.data[,1],probs=c(0.01,0.99))
    yquant <- quantile(fcs.data[,2],probs=c(0.01,0.99))
    print(xlim <- c(xquant[['1%']],xquant[['99%']]))
    print(ylim <- c(yquant[['1%']],yquant[['99%']]))
    if (!outliers) {
        if (nrow(fcs.data)>5000) {
            smoothScatter( fcs.data, colramp=colramp, nrpoints=nrpoints, cex.lab=2, cex.main=2, ylim=ylim, xlim=xlim, ... )
        } else {
            #plot( fcs.data, col=densCols(fcs.data,colramp=colramp), pch=20, ... )
            if (is.null(classification))
                plot( fcs.data, col='black', pch=20, ... )
            else
                #plot( fcs.data, col=densCols(fcs.data,colramp=colramp), pch=20, ... )
                plot( fcs.data, col=classification, pch=20, ... )
        }
    } else {
        if (nrow(fcs.data)>5000) {
            smoothScatter( fcs.data, colramp=colramp, nrpoints=nrpoints, cex.lab=2, cex.main=2, ... )
        } else {
            if (is.null(classification))
                plot( fcs.data, col='black', pch=20, ... )
            else
                #plot( fcs.data, col=densCols(fcs.data,colramp=colramp), pch=20, ... )
                plot( fcs.data, col=classification, pch=20, ... )
        }
    }
    col.names <- colnames(fcs.data)
    # plots the clusters either as chulls or as ellipses
    if(is.null(classification)) cols <- 'black' else cols <- classification
    if (!is.null(plot.points)) points(plot.points[,], col=plot.points.col, pch=20)
    if (!is.null(classification))
    for (k in sort(unique(classification))) {
        if (is.null(clusters.col)) col <- k
        else col <- clusters.col[[k]]
         X1 <- fcs.data[which(classification==k),col.names]
             if (chulls) {
                p <- X1[chull(X1),]
                p <- rbind(p, p)
                lines(p, col=col, lwd=chull.lwd)
             }
         #also approximate with an ellipse
         if (ellipses) lines(ellipse::ellipse(cov(X1),centre=colMeans(X1)),col=col,lwd=2,lty=1)
    } 
    if (!is.null(posteriors))
     for (k in 1:ncol(posteriors)) {
        if (is.null(clusters.col)) col <- k
        else col <- clusters.col[[k]]
        i <- which(posteriors[,k]>posterior.cutoff)
        if (length(i)>2) {
           X1 <- fcs.data[i,col.names]
           if (chulls) {
              #make sure chull is not open
              ch <- c(chull(X1),chull(X1))
              p <- X1[ch,]
              lines(p, col=col, lwd=chull.lwd)
           }
           #also approximate with an ellipse
           if (ellipses) lines(ellipse::ellipse(cov(X1),centre=colMeans(X1)),col=col,lwd=2,lty=1)
          }
     }
}

###
contourPlot <- function( fcs.data, bw=.1, outliers=FALSE, ... ) {
    xquant <- quantile(fcs.data[,1],probs=c(0.01,0.99))
    yquant <- quantile(fcs.data[,2],probs=c(0.01,0.99))
    print(xlim <- c(xquant[['1%']],xquant[['99%']]))
    print(ylim <- c(yquant[['1%']],yquant[['99%']]))
    dens <- bkde2D(fcs.data[,1:2],bw)
    if (!outliers)
        if (nrow(fcs.data)>5000) {
            contour(x=dens$x1, y=dens$x2, z=dens$fhat, cex.lab=2, cex.main=2, ylim=ylim, xlim=xlim, ... )
        } else {
            plot( fcs.data, pch=20,  ylim=ylim, xlim=xlim, ... )
        }
    else
        if (nrow(fcs.data)>5000) {
            contour(x=dens$x1, y=dens$x2, z=dens$fhat, cex.lab=2, cex.main=2, ... )
        } else {
            plot( fcs.data, pch=20, ... )
        }
}



#baseline relative pSTAT5
baseline.relative.pstat5 <- function(fcs.data, REPLACE=FALSE) {
    n <- as.numeric(sapply(strsplit(grep('^PSTAT5.',colnames(fcs.data),value=TRUE),'\\.'),`[[`,2))
    diff.pstat5 <- sapply(n, function(i) fcs.data[,paste('PSTAT5',i,sep='.')]/fcs.data[,'PSTAT5.1'])
    if (REPLACE) {
        colnames(diff.pstat5) <- paste('PSTAT5',n,sep='.')
        return(cbind(fcs.data[,-grep('PSTAT5',colnames(fcs.data))],diff.pstat5))
    } else {
        colnames(diff.pstat5) <- paste('diff','PSTAT5',n,sep='.')
        return(cbind(fcs.data,diff.pstat5))
    }
}


#baseline relative pSTAT5
baseline.relative.pstat5 <- function(fcs.data, REPLACE=FALSE) {
    n <- as.numeric(sapply(strsplit(grep('^PSTAT5.',colnames(fcs.data),value=TRUE),'\\.'),`[[`,2))
    diff.pstat5 <- sapply(n, function(i) fcs.data[,paste('PSTAT5',i,sep='.')]-fcs.data[,'PSTAT5.1'])
    if (REPLACE) {
        colnames(diff.pstat5) <- paste('PSTAT5',n,sep='.')
        return(cbind(fcs.data[,-grep('PSTAT5',colnames(fcs.data))],diff.pstat5))
    } else {
        colnames(diff.pstat5) <- paste('diff','PSTAT5',n,sep='.')
        return(cbind(fcs.data,diff.pstat5))
    }
}



otsu <- function(x) {
    #x <- sort(x)
    m1 <- cumsum(x[-length(x)])/(1:(length(x)-1))
    m2 <- (sum(x)-cumsum(x))[-length(x)]/((length(x)-1):1)
    return((1:(length(x)-1)/length(x))*((length(x)-1):1/length(x))*(m1 - m2)**2)
}

