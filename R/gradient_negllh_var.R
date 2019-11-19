#' this calculates the gradient with respect to the params being optimized when PBT,
#' GSI and another variable are used
#'
#' this is not intended for direct user-level use, but it is exported to allow
#' the user to optionally specify it as a argument to pass to \code{optim}
#' when using \code{MLEwrapper}. I recommend using it unless you have a specific
#' reason not to.
#'
#' @param params vector of paramaters to optimize
#' @param nPBT number of PBT groups to estimate
#' @param nGSI number of GSI groups to estimate
#' @param ohnc vector of number of observed (PBT assigned) hatchery fish in each PBT and GSI group (gsi groups should be 0)
#' @param t vector of tag rates for all PBT and GSI groups (gsi groups should be 0)
#' @param ohnc_gsi matrix of counts of fish GSI assigned to various groups
#' @param utVar matrix of numbers of un-PBT assigned fish in each GSI group (row) x category of the variable (column)
#' @param ohnc_var matrix of counts of PBT-assigned fish in each category (rows are PBT groups)
#' @param nCat the numbers of categories in the variable
#'
#' @export
#'
# uses the output of test_gradient_var
params_grad_var <- function(params, nPBT, nGSI, ohnc, t, ohnc_gsi,
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

	#get partial derivs with respect to piTot and piGSI
	firstLevel <- test_gradient_var(piTot, piGSItemp, piVar, nPBT, nGSI, ohnc, t, ohnc_gsi,
								utVar, ohnc_var, nCat)

	gradient <- rep(NA, length(params))

	########
	### copied and pasted
	########

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

	#now gradient of piVar
	for(i in 1:(nPBT+nGSI)){
		paramsTemp <- params[pos:(pos + nCat - 1)]
		firstLevelTemp <- firstLevel[pos:(pos + nCat - 1)]
		sumPar <- sum(paramsTemp)
		for(j in 1:nCat){
			temp <- 0
			for(k in 1:nCat){
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

	# print(gradient) #testing
	return(gradient)


 }


#' this calculates the gradient with respect to piTot, piGSI, and piVar
#'
#' @param piTot vector of piTot proportions
#' @param piGSItemp matrix of piGSI proportions for each PBT group and GSI group (GSI groups fixed at 1)
#' @param piVar matrix of piVar proportions for each PBT and GSI group
#' @param nPBT number of PBT groups to estimate
#' @param nGSI number of GSI groups to estimate
#' @param ohnc vector of number of observed (PBT assigned) hatchery fish in each PBT and GSI group (gsi groups should be 0)
#' @param t vector of tag rates for all PBT and GSI groups (gsi groups should be 0)
#' @param ohnc_gsi matrix of counts of fish GSI assigned to various groups
#' @param utVar matrix of numbers of un-PBT assigned fish in each GSI group (row) x category of the variable (column)
#' @param ohnc_var matrix of counts of PBT-assigned fish in each category (rows are PBT groups)
#' @param nCat the numbers of categories in the variable
#'
#' @keywords internal
#' @noRd
#'

test_gradient_var <- function(piTot, piGSItemp, piVar, nPBT, nGSI, ohnc, t, ohnc_gsi,
						utVar, ohnc_var, nCat){

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
			for(k in 1:nCat){
				if(utVar[j,k] == 0) next
				temp <- temp +
					(
						(utVar[j,k] * untag[i] * piGSItemp[i,j] * piVar[i,k]) /
						sum(untag * piGSItemp[,j] * piTot * piVar[,k])
					)
			}
		}
		gradient[i] <- gradient[i] + temp
	}

	#piGSI
	pos <- (nPBT + nGSI + 1) #position in the gradient vector to assign values
	if (nPBT > 0){
		for(i in 1:nPBT){
			for(j in 1:nGSI){
				temp <- 0
				for(k in 1:nCat){
					temp <- temp +
						(
							(utVar[j,k] * piTot[i] * untag[i] * piVar[i,k]) /
							sum(untag * piGSItemp[,j] * piTot * piVar[,k])
						)
				}
				#checking if greater than 0 to prevent dividing 0/0 if 0 is suggested for a piGSI value
				if(ohnc_gsi[i,j] > 0) temp <- temp + (ohnc_gsi[i,j] / piGSItemp[i,j])

				gradient[pos] <- temp
				pos <- pos + 1
			}
		}
	}

	#piVar
	for(i in 1:(nPBT+nGSI)){
		for(k in 1:nCat){
			temp <- 0
			for(j in 1:nGSI){
				temp <- temp +
					(
						(utVar[j,k] * piTot[i] * untag[i] * piGSItemp[i,j]) /
						sum(untag * piGSItemp[,j] * piTot * piVar[,k])
					)
			}
			#checking if greater than 0 to prevent dividing 0/0 if 0 is suggested for a piVar value
			if(i <= nPBT && ohnc_var[i,k] > 0) temp <- temp + (ohnc_var[i,k] / piVar[i,k])
			gradient[pos] <- temp
			pos <- pos + 1
		}
	}

	# print(gradient)
	return(-gradient)
}
