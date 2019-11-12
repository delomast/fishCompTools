#' this calculates the gradient with respect to the params being optimized when only PBT
#' and GSI information are used (no other variables)
#' 
#' this is not intended for direct user-level use, but it is exported to allow
#' the user to optionally specify it as a argument to pass to \code{optim}
#' when using \code{MLEwrapper}. I recommend using it unless you have a very specific
#' reason not to.
#' 
#' uses the output of test_gradient (chain rule for partial derivatives)

#' @param params vector of paramaters to optimize
#' @param nPBT number of PBT groups to estimate
#' @param nGSI number of GSI groups to estimate
#' @param ohnc vector of number of observed (PBT assigned) hatchery fish in each PBT and GSI group (gsi groups should be 0)
#' @param t vector of tag rates for all PBT and GSI groups (gsi groups should be 0)
#' @param utGSI vector of number of un-PBT assigned fish in each GSI group
#' @param ohnc_gsi matrix of counts of fish GSI assigned to various groups
#' 
#' @export

params_grad <- function(params, nPBT, nGSI, ohnc, t, utGSI, ohnc_gsi){
	# first, unpack params
	#piTot
	piTot <- params[1:(nPBT + nGSI)]
	piTot <- piTot / sum(piTot) #transform into proportions
	#piGSI
	# piGSItemp <- matrix(params[(nPBT + nGSI):length(params)], nrow = (nPBT), ncol = (nGSI-1), byrow = TRUE)
	subParams <- params[(nPBT + nGSI + 1):length(params)]
	piGSItemp <- matrix(0, nrow = (nPBT), ncol = (nGSI)) #initiate with zeros
	if(nPBT > 0){
		for(i in 1:nPBT){
			piGSItemp[i,] <- subParams[1:nGSI]
			subParams <- subParams[(nGSI + 1):length(subParams)] #bump entries forward
			piGSItemp[i,] <- piGSItemp[i,] / sum(piGSItemp[i,]) #normalize
		}
	}
	piGSItemp <- rbind(piGSItemp, diag(nGSI)) #add GSI groups as fixed 100%
	
	#get partial derivs with respect to piTot and piGSI
	firstLevel <- test_gradient(piTot, piGSItemp, nPBT, nGSI, ohnc, t, utGSI, ohnc_gsi)
	
	gradient <- rep(NA, length(params))
	#calc gradient of piTot with respect to params
	par_piTot <- params[1:(nPBT + nGSI)]
	sumPar <- sum(par_piTot)
	for(i in 1:(nPBT+nGSI)){
		temp <- 0
		for(j in 1:(nPBT+nGSI)){
			if(i == j){
				temp <- temp + (firstLevel[j] * ((sumPar - par_piTot[j]) / (sumPar^2)))
			} else {
				temp <- temp + (firstLevel[j] * (-par_piTot[j] / (sumPar^2)))
			}
		}
		gradient[i] <- temp
	}
	
	#now gradient of piGSI
	if (nPBT > 0){
		pos <- (nPBT + nGSI + 1) #position in the gradient vector to assign values
		for(i in 1:nPBT){
			paramsTemp <- params[pos:(pos + nGSI - 1)]
			firstLevelTemp <- firstLevel[pos:(pos + nGSI - 1)]
			sumPar <- sum(paramsTemp)
			for(j in 1:nGSI){
				temp <- 0
				for(k in 1:nGSI){
					if(j == k){
						temp <- temp + (firstLevelTemp[k] * ((sumPar - paramsTemp[k]) / (sumPar^2)))
					} else {
						temp <- temp + (firstLevelTemp[k] * (-paramsTemp[k] / (sumPar^2)))
					}
				}
				gradient[pos] <- temp
				pos <- pos + 1
			}
		}
	}
	# print(gradient) #testing
	return(gradient)
}


#this calculates the gradient with respect to piTot and piGSI
# 
# @param piTot vector of piTot proportions
# @param piGSItemp matrix of piGSI proportions for each PBT group and GSI group (GSI groups fixed at 1)
# @param nPBT number of PBT groups to estimate
# @param nGSI number of GSI groups to estimate
# @param ohnc vector of number of observed (PBT assigned) hatchery fish in each PBT and GSI group (gsi groups should be 0)
# @param t vector of tag rates for all PBT and GSI groups (gsi groups should be 0)
# @param utGSI vector of number of un-PBT assigned fish in each GSI group
# @param ohnc_gsi matrix of counts of fish GSI assigned to various groups
# 

test_gradient <- function(piTot, piGSItemp, nPBT, nGSI, ohnc, t, utGSI, ohnc_gsi){
	untag <- 1 - t
	
	#now calculate the gradient
	gradient <- rep(0, length(piTot) + nPBT*nGSI)
	
	#piTot
	#first, PBT groups observed part
	if(nPBT > 0) gradient[1:nPBT] <- ohnc[1:nPBT] / piTot[1:nPBT]
	#then all the groups unobserved part
	for(i in 1:(nPBT + nGSI)){
		temp <- 0
		for(j in 1:nGSI){
			if(utGSI[j] == 0) next
			temp <- temp + 
				(
					(utGSI[j] * untag[i] * piGSItemp[i,j]) /
					sum(untag * piGSItemp[,j] * piTot)	
				)
		}
		gradient[i] <- gradient[i] + temp
	}
	
	#piGSI
	if (nPBT > 0){
		pos <- (nPBT + nGSI + 1) #position in the gradient vector to assign values
		for(i in 1:nPBT){
			for(j in 1:nGSI){
				temp <- 0
				temp <- ((utGSI[j] * piTot[i] * untag[i]) /
					sum(untag * piGSItemp[,j] * piTot))
				if(ohnc_gsi[i,j] > 0) temp <- temp + (ohnc_gsi[i,j] / piGSItemp[i,j])
				gradient[pos] <- temp
				pos <- pos + 1
			}
		}
	}
	# print(gradient)
	return(-gradient)
}
