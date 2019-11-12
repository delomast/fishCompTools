#' calculate the negative log-likelihood with only PBT and GSI (optional) information
#' 
#' this function does not allow another variable. for that, see \code{flex_negllh_var}
#' 
#' @param params vector of paramaters to optimize
#' @param nPBT number of PBT groups to estimate
#' @param nGSI number of GSI groups to estimate
#' @param ohnc vector of number of observed (PBT assigned) hatchery fish in each PBT and GSI group (gsi groups should be 0)
#' @param t vector of tag rates for all PBT and GSI groups (gsi groups should be 0)
#' @param utGSI vector of number of un-PBT assigned fish in each GSI group
#' @param ohnc_gsi matrix of counts of fish GSI assigned to various groups
#'  

flex_negllh_allGSI <- function(params, nPBT, nGSI, ohnc, t, utGSI, ohnc_gsi){
	# first, unpack params
	#piTot
	piTot <- params[1:(nPBT + nGSI)]
	piTot <- piTot / sum(piTot) #transform into proportions
	# print(piTot)
	if(sum(piTot < 0 | piTot > 1) != 0) return(Inf)
	#piGSI
	# piGSItemp <- matrix(params[(nPBT + nGSI):length(params)], nrow = (nPBT), ncol = (nGSI-1), byrow = TRUE)
	subParams <- params[(nPBT + nGSI + 1):length(params)]
	piGSItemp <- matrix(0, nrow = (nPBT), ncol = (nGSI)) #initiate with zeros
	if(nPBT > 0){
		for(i in 1:nPBT){
			piGSItemp[i,] <- subParams[1:nGSI]
			subParams <- subParams[(nGSI + 1):length(subParams)] #bump entries forward
			piGSItemp[i,] <- piGSItemp[i,] / sum(piGSItemp[i,]) #normalize
			# print(piGSItemp[i,])
			if(sum(piGSItemp[i,] < 0 | piGSItemp[i,] > 1) != 0) return(Inf) #make sure all entries are valid
		}
	}
	piGSItemp <- rbind(piGSItemp, diag(nGSI)) #add GSI groups as fixed 100%
	
	# now, calculate the log likelihood
	llh <- 0
	# first ohnc part
	if(nPBT > 0) llh <- sum(ohnc[1:nPBT] * log(piTot[1:nPBT] * t[1:nPBT]))
	# then utGSI part
	untag <- 1 - t
	for(j in 1:nGSI){
		llh <- llh + log(sum(piTot * untag * piGSItemp[,j])) * utGSI[j]
	}
	# then ohnc GSI part
	if(nPBT > 0){
		for(i in 1:nPBT){
			llh <- llh + sum(ohnc_gsi[i,] * log(piGSItemp[i,]))
		}
	}
	
	# returning negative log-likelihood for minimization
	return(-llh)
}


