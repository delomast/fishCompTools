#' This function uses the input created by \code{prepMultStrataPropClip} to estimate the
#' proportion of fish that are clipped for each strata.
#' 
#' A wrapper for \code{propClipped}
#' 
#' @param prepClipData the output of prepMultStrataPropClip
#' @param NumResults the number of samples from the posterior to take
#' @param seed A positive integer to seed the random number generator. If NA, a seed will
#'   be chosen based on the time.
#' @param priors a list of numeric vectors giving the parameters of a beta distribution to use
#'   as the priors for each strata (in order of the strata). The default is a uniform Beta(1,1)
#'   prior for all strata.
#' @param verbose TRUE to print some output
#' 
#' @export

estimStrataPropClip <- function(prepClipData, NumResults, seed = NA, priors = NA, verbose = TRUE){
	
	# get seed
	if(is.na(seed)){
		seed <- ceiling(as.numeric(format(Sys.time(), "%S")) * as.numeric(format(Sys.time(), "%j")) * as.numeric(format(Sys.time(), "%M")))
	}
	
	# check priors
	if(!is.list(priors)){
		if(verbose) cat("\nUsing the default of uniform priors for all strata.\n")
		priors <- list()
		for(i in 1:length(prepClipData)){
			priors[[i]] <- c(1,1)
		}
	} else {
		if(verbose) cat("\nUsing the supplied priors.\n")
		if(length(priors) != length(prepClipData)){
			stop("The number of priors supplied is not equal to the number of strata.")
		}
	}
	
	# get samples
	clipEstimates <- list()
	for(i in 1:length(prepClipData)){
		sData <- prepClipData[[i]][[1]]
		clipEstimates[[i]] <- propClipped(sData[1], sData[2], NumResults, seed, priors[[i]])
	}
	
	return(clipEstimates)
}
