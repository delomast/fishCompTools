#' calculate the negative log-likelihood with PBT, GSI (optional), and one other categorical variable
#' 
#' more here?
#' 
#' @param params vector of paramaters to optimize
#' @param nPBT number of pbt groups to estimate
#' @param nGSI number of GSI groups to estimate
#' @param ohnc vector of number of observed (PBT assigned) hatchery fish in each PBT and GSI group (gsi groups should be 0)
#' @param t vector of tag rates for all PBT and GSI groups (gsi groups should be 0)
#' @param utGSI vector of number of un-PBT assigned fish in each GSI group
#' @param ohnc_gsi matrix of counts of fish GSI assigned to various groups
#' @param utVar matrix of numbers of un-PBT assigned fish in each GSI group (row) x category of the variable (column)
#' @param ohnc_var matrix of counts of PBT-assigned fish in each category (rows are PBT groups)
#' @param nCat the numbers of categories in the variable
#'  

flex_negllh_var <- function(params, nPBT, nGSI, ohnc, t, utGSI, ohnc_gsi,
									 utVar, ohnc_var, nCat){
	# first, unpack params
	#piTot
	piTot <- params[1:(nPBT + nGSI)]
	piTot <- piTot / sum(piTot) #transform into proportions
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
			if(sum(piGSItemp[i,] < 0 | piGSItemp[i,] > 1) != 0) return(Inf) #make sure all entries are valid
		}
	}
	piGSItemp <- rbind(piGSItemp, diag(nGSI)) #add GSI groups as fixed 100%
	#piVar
	piVar <- matrix(0, nrow = (nPBT + nGSI), ncol = nCat) #initiate with zeros
	for(i in 1:(nPBT +  nGSI)){
		piVar[i,] <- subParams[1:nCat]
		subParams <- subParams[(nCat + 1):length(subParams)] #bump entries forward
		piVar[i,] <- piVar[i,] / sum(piVar[i,]) #normalize
		if(sum(piVar[i,] < 0 | piVar[i,] > 1) != 0) return(Inf) #make sure all entries are valid
	}

	# now, calculate the log likelihood
	llh <- 0
	# first ohnc part
	if(nPBT > 0) llh <- sum(ohnc[1:nPBT] * log(piTot[1:nPBT] * t[1:nPBT]))
	# then ohnc GSI part
	if(nPBT > 0){
		for(i in 1:nPBT){
			llh <- llh + sum(ohnc_gsi[i,] * log(piGSItemp[i,]))
		}
	}
	#then ohnc var part
	if(nPBT > 0){
		for(i in 1:nPBT){
			llh <- llh + sum(ohnc_var[i,] * log(piVar[i,]))
		}
	}
	# then ut part
	untag <- 1 - t
	for(j in 1:nCat){
		for(k in 1:nGSI){
			if(utVar[k,j] > 0) {
				tempSumProp <- sum(piTot * untag * piGSItemp[,k] * piVar[,j])
				if(tempSumProp <= 0) return(Inf)
				llh <- llh + log(tempSumProp) * utVar[k,j]
			}
		}
	}


	
	# returning negative log-likelihood for minimization
	return(-llh)
}


