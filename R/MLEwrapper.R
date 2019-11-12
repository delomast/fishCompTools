#' Find maximum likelihood estimates for composition
#' 
#' The likelihood is maximized using \code{optim} and some typically close, but quick to calculate, starting values.
#' If estimating composition of a separate variable, all samples with missing data for this variable must be removed
#' prior to running this function.
#' 
#' @param trapData a dataframe with a rwo for each individual and columns for GSI assignment, PBT assignment, etc.
#' @param tags a dataframe with the first column containing names of PBT groups, and the second column containing
#'   tag rates
#' @param GSIcol name of column containing GSI assignments. If you have no GSI information, create a column with
#'   the same value for all samples.
#' @param PBTcol name of column containing PBT assignments
#' @param strataCol name of column indicating the strata the observation belongs to
#' @param adFinCol name of column containing adipose fin status - TRUE (or AI) being intact FALSE (or AD) being clipped, NA missing
#' @param AI TRUE to analyze ad-intact fish, FALSE to analyze ad-clipped fish
#' @param optimMethod the method to use first with \code{optim}. If this method fails, "Nelder-Mead" is attempted
#'   and a warning is issued.
#' @param variableCols name of column containing the variable to estimate composition for (optional)
#' @param ... other arguments to pass to \code{optim}
#' 
#' @export

MLEwrapper <- function(trapData, tags, GSIcol, PBTcol, strataCol, adFinCol, AI = TRUE, optimMethod = "Nelder-Mead",
							  variableCols = NULL, ...){
	
	# determine if variables used or not
	varBool <- !is.null(variableCols) #TRUE if variables are used
	
	# input checking
	if(varBool && length(variableCols) > 1){
		stop("variableCols must be either NULL or only one variable")
	}
	
	#don't need all the output from this, but it includes most things we need
	allInput <- prepStrata(trapData, tags, GSIcol, PBTcol, strataCol, variableCols = variableCols, variableColsOth = c(), adFinCol,
								AI = TRUE, GSIgroups = NA,
									 variableValues = NA, variableValuesOth = NA, verbose = FALSE, symPrior = .5)
	
	estimates <- list()
	#get estimates for each strata
	for(input in allInput){
		#pull values out of input
		ohnc <- input$ohnc
		nPBT <- input$nPBT
		t <- input$t
		ohnc_gsi <- input$ohnc_gsi
		utGSI <- c()
		for(g in input$GSI_values){
			utGSI <- c(utGSI, sum(input$gsiUT == g))
		}

		#define some values used in the function
		nGSI <- length(input$groups) - nPBT
		
		
		#pull variables values out of input
		ohnc_var <- data.frame() #define in case no PBT groups
		if(varBool){
			nCat <- length(input$values[[1]])
			if(nPBT > 0) ohnc_var <- input$pi_Vohnc[[1]][1:nPBT,]
			utVar <- matrix(0, nrow = nGSI, ncol = nCat)
			for(v in 1:nCat){
				val <- input$values[[1]][v]
				for(g in 1:nGSI){
					gVal <- input$GSI_values[g]
					utVar[g,v] <- sum(input$v_ut[,1] == val & input$gsiUT == gVal)
				}
			}
		}

		
		# determine reasonable starting values
		#for piTot
		start_piTot <- c()
		if(nPBT > 0) start_piTot <- ohnc[1:nPBT] / t[1:nPBT] #scale PBT by tag rates
		start_piTot <- c(start_piTot, utGSI) #just use observed GSI

		#for piGSI
		start_piGSI <- c()
		if(nPBT > 0){
			for(i in 1:nrow(ohnc_gsi)){
				temp <- ohnc_gsi[i,] #just use ohnc assignments
				temp[temp < 1] <- .1
				start_piGSI <- c(start_piGSI, temp)
			}
		}
		#for piVar
		start_piVar <- c()
		if(varBool){
			for(i in 1:nrow(ohnc_var)){
				temp <- ohnc_var[i,] #just use ohnc assignments
				temp[temp < 1] <- .1
				start_piVar <- c(start_piVar, temp)
			}
			# then ut fish
			for(i in 1:nGSI){
				#count un-tagged fish with this GSI assignment and each category
				temp <- c()
				for(val in input$values[[1]]){
					temp <- c(temp, sum(input$gsiUT == input$GSI_values[i] & input$v_ut[,1] == val))
				}
				start_piVar <- c(start_piVar, temp)
			}
		}


		# find mle
		if(varBool){
			tempFit <- optim(c(start_piTot, start_piGSI, start_piVar), flex_negllh_var, method = optimMethod, ...,
					 #arguments to pass to flex_negllh
					 nPBT = nPBT, nGSI = nGSI, ohnc = ohnc, t = t, utGSI = utGSI, ohnc_gsi = ohnc_gsi,
					 utVar = utVar, ohnc_var = ohnc_var, nCat = nCat)
		} else {
			
			# print(numDeriv::grad(function(u) flex_negllh_allGSI(u, nPBT = nPBT, nGSI = nGSI, ohnc = ohnc, t = t, utGSI = utGSI, ohnc_gsi = ohnc_gsi), 
			# 					c(start_piTot, start_piGSI))) #testing gradient function
			
			tempFit <- optim(c(start_piTot, start_piGSI), flex_negllh_allGSI, method = optimMethod, ...,
								 #arguments to pass to flex_negllh
								 nPBT = nPBT, nGSI = nGSI, ohnc = ohnc, t = t, utGSI = utGSI, ohnc_gsi = ohnc_gsi)
			# print(tempFit) #testing
		}
		if(tempFit$convergence != 0) cat("\nOptimizer gave convergence code of", tempFit$convergence, "in strata", input$strataName, "\n")
		
		# return(tempFit) #testing


		#unpack values
		# unpack piTot
		ptestim <- tempFit$par[1:(nPBT + nGSI)]
		ptestim <- ptestim / sum(ptestim)
		names(ptestim) <- input$groupsKey[,1]
		
		# piGSI
		subParams <- tempFit$par[(nPBT + nGSI + 1):length(tempFit$par)]
		piGSItemp <- matrix(0, nrow = (nPBT), ncol = (nGSI)) #initiate with zeros
		if(nPBT > 0){
			for(i in 1:nPBT){
				piGSItemp[i,] <- subParams[1:nGSI]
				subParams <- subParams[(nGSI + 1):length(subParams)]
				piGSItemp[i,] <- piGSItemp[i,] / sum(piGSItemp[i,])
			}
		}
		piGSItemp <- rbind(piGSItemp, diag(nGSI)) #add GSI groups as fixed 100%
		colnames(piGSItemp) <- input$GSIkey[,1]
		rownames(piGSItemp) <- input$groupsKey[,1]
		
		#piVar
		piVar <- matrix(nrow = 0, ncol = 0)
		if(varBool){
			piVar <- matrix(0, nrow = (nPBT + nGSI), ncol = nCat) #initiate with zeros
			for(i in 1:(nPBT +  nGSI)){
				piVar[i,] <- subParams[1:nCat]
				subParams <- subParams[(nCat + 1):length(subParams)]
				piVar[i,] <- piVar[i,] / sum(piVar[i,])
			}
			colnames(piVar) <- input$variKey[[1]][,1]
			rownames(piVar) <- input$groupsKey[,1]
		}
		
		estimates[[input$strataName]] <- list(piTot = ptestim, piGSI = piGSItemp, piVar = piVar, strataName = input$strataName)
	}
	
	return(estimates)
}
