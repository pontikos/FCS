library(gtable)
library(scales)
library(ggplot2)
library(spade)
library(flowCore)
library(lubridate)
library(GGally)

lgcl <- logicleTransform()

# Wrapper for the flowCore read.FCS function.
# Applies compensation
read.FCS <- function(file, channels=c(), comp=TRUE, verbose=FALSE, ...) {
	if (verbose) fcs <- flowCore::read.FCS(file, ...)
	else fcs <- suppressWarnings(flowCore::read.FCS(file, ...))
    	params <- parameters(fcs)
	pd     <- pData(params)
    # Replace any null descs with names (for FSC-A, FSC-W, SSC-A)
  bad_col <- grep("^[a-zA-Z0-9]+",pd$desc,invert=TRUE)
	if (length(bad_col) > 0) {
		keyval <- keyword(fcs)
		for (i in bad_col) {
			pd$desc[i] <- pd$name[i]
			keyval[[paste("$P",i,"S",sep="")]] <- pd$name[i]
		}
		pData(params) <- pd;
		fcs <- flowFrame(exprs(fcs),params,description=description(fcs));
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
		flowFrame(exprs(comp_fcs), parameters(comp_fcs), description(comp_fcs)[grep("SPILL",names(description(comp_fcs)),invert=TRUE)])
	}
	
	if (comp && !is.null(description(fcs)$SPILL)) {
		fcs <- apply.comp(fcs, "SPILL")
	} else if (comp && !is.null(description(fcs)$SPILLOVER)) {
		fcs <- apply.comp(fcs, "SPILLOVER")
	}
    #remove negative scatter
    fcs <- fcs[fcs@exprs[,'SSC-A']>1,]
    fcs <- fcs[fcs@exprs[,'FSC-A']>1,]
    #truncate side scatter at 500 
    fcs <- fcs[fcs@exprs[,'SSC-A']>500,]
    for (channel in channels) {
        #i <- grep( channel, parameters(fcs)@data$desc, ignore.case=T )
        if (channel == 'FSC-A') trans <- function(x) 5*x/262144
        else if (channel == 'SSC-A') trans <- log10
        else trans <- logicleTransform()
        #fcs@exprs[,i] <- trans(fcs@exprs[,i])
        fcs <- setChannels(fcs, channel, trans(getChannels(fcs, channel)))
    }
    #remove min & max scatter
    min.ssc <- min(fcs@exprs[,'SSC-A'])
    max.ssc <- max(fcs@exprs[,'SSC-A'])
    fcs <- fcs[min.ssc < fcs@exprs[,'SSC-A'] & fcs@exprs[,'SSC-A'] < max.ssc,]
    min.fsc <- min(fcs@exprs[,'FSC-A'])
    max.fsc <- max(fcs@exprs[,'FSC-A'])
    fcs <- fcs[min.fsc < fcs@exprs[,'FSC-A'] & fcs@exprs[,'FSC-A'] < max.fsc,]
    fcs
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

getChannels <- function(fcs.data, channels) {
    if (is(fcs.data,'matrix') || is(fcs.data,'data.frame')) {
        if (missing(channels)) channels <- colnames(fcs.data)
        channels <- gsub('-', '',paste(channels, ifelse(!grepl('-', channels), '$', ''), sep=''))
        fcs.channels <- gsub('-','',colnames(fcs.data))
        x <- as.matrix(fcs.data[,sapply(channels, function(x) grep(x, fcs.channels, ignore.case=T))])
    } else if (is(fcs.data, 'flowFrame')) {
        if (missing(channels)) channels <- parameters(fcs.data)@data$desc
        channels <- gsub('-', '', paste(channels, ifelse(!grepl('-', channels), '$', ''), sep=''))
        fcs.channels <- gsub('-','',parameters(fcs.data)@data$desc)
        x <- as.matrix(fcs.data@exprs[,sapply(channels, function(x) grep(x, fcs.channels, ignore.case=T))])
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
    channels <- gsub('-','',paste(channels, ifelse(!grepl('-', channels), '$', ''), sep=''))
    if (is(fcs.data,'matrix') || is(fcs.data,'data.frame')) {
        fcs.channels <- gsub('-','',colnames(fcs.data))
        fcs.data[,sapply(channels, function(x) grep( x, fcs.channels, ignore.case=T ))] <- x
    } else if (is(fcs.data, 'flowFrame')) {
        fcs.channels <- gsub('-','',parameters(fcs.data)@data$desc)
        fcs.data@exprs[,sapply(channels, function(x) grep( x, fcs.channels, ignore.case=T ))] <- x
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


compute.density <- function(fcs.data, channels, kernel_mult=5.0, apprx_mult=1.5, med_samples=2000) {
    if (missing(channels)) channels <- colnames(getChannels(fcs.data))
    SPADE.density(getChannels(fcs.data, channels), kernel_mult, apprx_mult, med_samples)
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




