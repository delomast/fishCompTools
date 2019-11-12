#' This function uses the input created by \code{prepStrata} to estimate the
#' composistion of each strata.
#'
#' A wrapper for \code{MCpbt}. This function is used to run the MCMC chain and produce estimates.
#'
#' @param prepData The output of \code{prepStrat} with priors and initial values modified as you see fit
#' @param iter The total number of iterations to run each strata, including burn-in iterations
#' @param burnIn The number of burn-in iterations
#' @param thin Thinning parameter for the MCMC chain
#' @param seed a positive integer to seed the random number generator. If not specified, chosen based on the current time.
#'
#' @importFrom Rcpp evalCpp
#' @useDynLib fishCompTools, .registration=TRUE
#'
#' @export

estimStrataMCpbt <- function(prepData, iter, burnIn, thin = 1, seed = NA){

	# get seed
	if(is.na(seed)){
		seed <- ceiling(as.numeric(format(Sys.time(), "%S")) * as.numeric(format(Sys.time(), "%j")) * as.numeric(format(Sys.time(), "%M")))
	}


	# get samples
	strataEstimates <- list()
	for(i in 1:length(prepData)){
		cStrata <- prepData[[i]]
		strataEstimates[[i]] <- MCpbt(iter = iter, burnIn = burnIn , thin = thin, seed = seed, #overall parameters for the chain
					     piTotPrior = cStrata$piTotPrior, ohnc = cStrata$ohnc, piTotInitial = cStrata$piTotInitial, #piTotal parameters
					     oUTInitial = cStrata$oUTInitial, groups = cStrata$groups,
					     nPBT = cStrata$nPBT, GSI_values = cStrata$GSI_values, gsiUT = cStrata$gsiUT, #pi_gsi parameters
					     pi_gsiInitial = cStrata$pi_gsiInitial, prior_pi_gsi = cStrata$prior_pi_gsi,
					     ohnc_gsi = cStrata$ohnc_gsi,
					     values = cStrata$values, pi_VInitial = cStrata$pi_VInitial, pi_Vohnc = cStrata$pi_Vohnc, pi_Vprior = cStrata$pi_Vprior, #pi_V parameters
					     v_ut = cStrata$v_ut,
					     initZ = cStrata$initZ, t = cStrata$t, #z parameters
					     valuesOth = cStrata$valuesOth, pi_VInitialOth = cStrata$pi_VInitialOth, pi_VohncOth = cStrata$pi_VohncOth, #pi_VOth parameters
					     pi_VpriorOth = cStrata$pi_VpriorOth,
                 v_utOth = cStrata$v_utOth
					     )

		#add dimension names to things
		#piTot
		colnames(strataEstimates[[i]]$piTot) <- cStrata$groupsKey[match(cStrata$groups, cStrata$groupsKey[,2]),1]

		#piGSI
		cur <- 1 # to iterate appropriately across columns
		end <- length(cStrata$GSI_values)
		colnames(strataEstimates[[i]]$piGSI) <- paste0("temp", 1:ncol(strataEstimates[[i]]$piGSI))
		for(j in 1:ncol(strataEstimates[[i]]$piTot)){
			# assign column names
			colnames(strataEstimates[[i]]$piGSI)[cur:end] <- paste0(cStrata$groupsKey[cStrata$groupsKey[,2] == cStrata$groups[j],1],
																 "_",
																 cStrata$GSIkey[match(cStrata$GSI_values, cStrata$GSIkey[,2]),1])
			cur <- end + 1
			end <- end + length(cStrata$GSI_values)
		}
		rm(cur)
		rm(end)

		# piV
		pos <- 1
		names(strataEstimates[[i]]$piV) <- names(cStrata$values)
		for(v in names(cStrata$values)){
			colnames(strataEstimates[[i]]$piV[[v]]) <- paste0("temp", 1:ncol(strataEstimates[[i]]$piV[[v]]))
			cur <- 1 # to iterate appropriately across columns
			end <- length(cStrata$values[[pos]])
			for(j in 1:ncol(strataEstimates[[i]]$piTot)){
				# assign column names
				colnames(strataEstimates[[i]]$piV[[v]])[cur:end] <- paste0(cStrata$groupsKey[cStrata$groupsKey[,2] == cStrata$groups[j],1],
																	 "_",
																	 cStrata$variKey[[pos]][match(cStrata$values[[pos]], cStrata$variKey[[pos]][,2]),1])
				cur <- end + 1
				end <- end + length(cStrata$values[[pos]])
			}
			rm(cur)
			rm(end)
			pos <- pos + 1
		}
		rm(pos)

		# piVOth
		pos <- 1
		names(strataEstimates[[i]]$piVOth) <- names(cStrata$valuesOth)
		for(v in names(cStrata$valuesOth)){
			colnames(strataEstimates[[i]]$piVOth[[v]]) <- paste0("temp", 1:ncol(strataEstimates[[i]]$piVOth[[v]]))
			cur <- 1 # to iterate appropriately across columns
			end <- length(cStrata$valuesOth[[pos]])
			for(j in 1:ncol(strataEstimates[[i]]$piTot)){
				# assign column names
				colnames(strataEstimates[[i]]$piVOth[[v]])[cur:end] <- paste0(cStrata$groupsKey[cStrata$groupsKey[,2] == cStrata$groups[j],1],
																	 "_",
																	 cStrata$variKeyOth[[pos]][match(cStrata$valuesOth[[pos]], cStrata$variKeyOth[[pos]][,2]),1])
				cur <- end + 1
				end <- end + length(cStrata$valuesOth[[pos]])
			}
			rm(cur)
			rm(end)
			pos <- pos + 1
		}
		rm(pos)
	}

	return(strataEstimates)

}
