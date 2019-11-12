#' Multiplies escapement estimates by composition estimates
#' 
#' @param prepData the same \code{prepData} used as input for \code{estimStrataMCpbt}
#' @param mainRes the output of \code{estimStrataMCpbt}
#' @param clipRes the output of \code{estimStrataPropClip}
#' @param popSizes a list of population size estimates for each strata, in order, 
#'   with each strata being an entry in the list
#' @param verbose FALSE to repress output written to the screen
#' @param writeSummary FALSE to repress summary files being written
#' @param prefix the prefix to use when namign the summary files
#' 
#' @export
#' 


multByEscapement <- function(prepData, mainRes, clipRes = NA, popSizes, verbose = TRUE, writeSummary = TRUE, prefix = "", alpha = .1){

	# have optional summary of posterior means, medians, and CI's written to file
	
	## check inputs
	numStrata <- length(prepData)
	if(numStrata != length(popSizes)){
		stop("prepData and popSizes do not represent the same number of strata.")
	}
	if(numStrata != length(mainRes)){
		stop("prepData and mainRes do not represent the same number of strata.")
	}
	if(is.list(clipRes) && (numStrata != length(clipRes))){
		stop("prepData and clipRes do not represent the same number of strata.")
	}
	stratSamp <- c()
	escapeSamp <- c()
	for(i in 1:numStrata){
		stratSamp <- c(stratSamp, nrow(mainRes[[i]]$piTot))
		escapeSamp <- c(escapeSamp, length(popSizes[[i]]))
	}
	stratSamp <- unique(stratSamp)
	escapeSamp <- unique(escapeSamp)
	if(length(stratSamp) > 1){
		stop("not all strata in mainRes have the same number of iterations")
	}
	if(length(escapeSamp) > 1){
		stop("not all strata in popSizes have the same length")
	}
	
	if(stratSamp != escapeSamp){
		stop("The number of iterations for each strata in mainRes does not equal the length of each strata in popSizes")
	}
	

	if(is.list(clipRes)){
		for(i in 1:numStrata){
		if(length(clipRes[[i]]) != length(popSizes[[i]])){
			err <- paste("In strata number", i, "the number of samples in clipRes is not equal", 
							 "to the number of samples in popSizes")
			stop(err)
			}
		}
	}
	
	
	## identify AI
	AI <- prepData[[1]]$AI
	
	## create clipRes if clipRes not given
	if(!is.list(clipRes)){
		if (AI){
			if(verbose) cat("\nNo input given for clipRes, assuming all fish are ad-intact\n")
			#make input for all ad-intact
			clipRes <- list()
			for(i in 1:length(prepData)){
				clipRes[[i]] <- rep(0, length(popSizes[[i]]))
			}
		} else {
			if(verbose) cat("\nNo input given for clipRes, assuming all fish are ad-clipped\n")
			#make input for all ad-clipped
			clipRes <- list()
			for(i in 1:length(prepData)){
				clipRes[[i]] <- rep(1, length(popSizes[[i]]))
			}
		}
	}
	
	if(verbose && AI) cat("\nDecomposing ad-intact fish\n")
	if(verbose && !AI) cat("\nDecomposing ad-clipped fish\n")
	
	strataEstimates <- list()
	#for each strata
	for(i in 1:length(prepData)){
		strataName <- prepData[[i]]$strataName
		#clipped and unclipped breakdown
		clipUnclip <- popSizes[[i]]*clipRes[[i]] ## calculate number clipped
		clipUnclip <- cbind(clipUnclip, popSizes[[i]] - clipUnclip) ##remainder are unclipped
		colnames(clipUnclip) <- c("Number_clipped", "Number_unclipped")
		#select clipped or unclipped total to use downstream
		if (AI) {
			total <- clipUnclip[,2]
		} else {
			total <- clipUnclip[,1] 
		}
		
		#save current strata to make access easier
		cStrataInput <- prepData[[i]]
		cStrata <- mainRes[[i]]
		
		# breakdown by piTot
		piTotNumbers <- cStrata$piTot * total
		colnames(piTotNumbers) <- cStrataInput$groupsKey[match(cStrataInput$groups, cStrataInput$groupsKey[,2]),1]
		
		# piGSI - this doesn't seem that useful, but why not calculate it all the same
		cur <- 1 # to iterate appropriately across columns
		end <- length(cStrataInput$GSI_values)
		GSInumbers <- matrix(0, nrow = nrow(cStrata$piGSI), ncol = ncol(cStrata$piGSI))
		colnames(GSInumbers) <- paste0("temp", 1:ncol(cStrata$piGSI))
		for(j in 1:ncol(piTotNumbers)){
			GSInumbers[,cur:end] <- cStrata$piGSI[,cur:end]*piTotNumbers[,j]
			# assign column names
			colnames(GSInumbers)[cur:end] <- paste0(cStrataInput$groupsKey[cStrataInput$groupsKey[,2] == cStrataInput$groups[j],1],
																 "_",
																 cStrataInput$GSIkey[match(cStrataInput$GSI_values, cStrataInput$GSIkey[,2]),1])
			cur <- end + 1
			end <- end + length(cStrataInput$GSI_values)
		}
		rm(cur)
		rm(end)
		
		# piV
		piVnumbers <- list()
		pos <- 1
		for(v in names(cStrataInput$values)){
			props <- cStrata$piV[[pos]]
			cur <- 1 # to iterate appropriately across columns
			end <- length(cStrataInput$values[[pos]])
			piVcounts <- matrix(0, nrow = nrow(props), ncol = ncol(props))
			colnames(piVcounts) <- paste0("temp", 1:ncol(props))
			for(j in 1:ncol(piTotNumbers)){
				piVcounts[,cur:end] <- props[,cur:end]*piTotNumbers[,j]
				# assign column names
				colnames(piVcounts)[cur:end] <- paste0(cStrataInput$groupsKey[cStrataInput$groupsKey[,2] == cStrataInput$groups[j],1],
																	 "_",
																	 cStrataInput$variKey[[pos]][match(cStrataInput$values[[pos]], cStrataInput$variKey[[pos]][,2]),1])
				cur <- end + 1
				end <- end + length(cStrataInput$values[[pos]])
			}
			rm(cur)
			rm(end)
			piVnumbers[[v]] <- piVcounts
			pos <- pos + 1
		}
		rm(pos)
		
		# piVOth
		piVnumbersOth <- list()
		pos <- 1
		for(v in names(cStrataInput$valuesOth)){
			props <- cStrata$piVOth[[pos]]
			cur <- 1 # to iterate appropriately across columns
			end <- length(cStrataInput$valuesOth[[pos]])
			piVcounts <- matrix(0, nrow = nrow(props), ncol = ncol(props))
			colnames(piVcounts) <- paste0("temp", 1:ncol(props))
			for(j in 1:ncol(piTotNumbers)){
				piVcounts[,cur:end] <- props[,cur:end]*piTotNumbers[,j]
				# assign column names
				colnames(piVcounts)[cur:end] <- paste0(cStrataInput$groupsKey[cStrataInput$groupsKey[,2] == cStrataInput$groups[j],1],
																	 "_",
																	 cStrataInput$variKeyOth[[pos]][match(cStrataInput$valuesOth[[pos]], cStrataInput$variKeyOth[[pos]][,2]),1])
				cur <- end + 1
				end <- end + length(cStrataInput$valuesOth[[pos]])
			}
			rm(cur)
			rm(end)
			piVnumbersOth[[v]] <- piVcounts
			pos <- pos + 1
		}
		rm(pos)
		
		strataEstimates[[strataName]] <- list(clipEstim = clipUnclip, 
														  piTotEstim = piTotNumbers,
														  GSIestim = GSInumbers,
														  variableEstim = piVnumbers,
														  variableEstimOth = piVnumbersOth,
														  strataName = strataName
														  )
	}
	
	# now sum all strata together
	
	
		# clipEstim = clipUnclip, 
	totClipUnclip <- matrix(0, nrow = nrow(strataEstimates[[1]]$clipEstim), ncol = 2)
	for(i in 1:numStrata){
		totClipUnclip[,1] <- totClipUnclip[,1] + strataEstimates[[i]]$clipEstim[,1]
		totClipUnclip[,2] <- totClipUnclip[,2] + strataEstimates[[i]]$clipEstim[,2]
	}
	colnames(totClipUnclip) <- c("Number_clipped", "Number_unclipped")
	
	## general strategy: 
	## get all column names across all strata
	## build matrix
	## sum up by mathcing column names across and within strata
	##		the within option allows for changing the key after running
	##		to easily pool groups
	
		# piTotEstim = piTotNumbers,
	## get all column names across all strata
	columns <- c()
	for(i in 1:numStrata){
		columns <- c(columns, colnames(strataEstimates[[i]]$piTotEstim))
	}
	columns <- unique(columns)
	## build matrix
	totPiTotEstim <- matrix(0, nrow = nrow(strataEstimates[[1]]$piTotEstim), ncol = length(columns))
	colnames(totPiTotEstim) <- columns
	## sum up by mathcing column names across and within strata
	for(c in columns){
		for(i in 1:numStrata){
			tempEstim <- strataEstimates[[i]]$piTotEstim
			for(j in which(colnames(tempEstim) == c)){
				totPiTotEstim[,c] <- totPiTotEstim[,c] + tempEstim[,j]
			}
		}
	}
	rm(columns)
	
		# GSIestim = GSInumbers,
	columns <- c()
	for(i in 1:numStrata){
		columns <- c(columns, colnames(strataEstimates[[i]]$GSIestim))
	}
	columns <- unique(columns)
	totGSIestim <- matrix(0, nrow = nrow(strataEstimates[[1]]$GSIestim), ncol = length(columns))
	colnames(totGSIestim) <- columns
	for(c in columns){
		for(i in 1:numStrata){
			tempEstim <- strataEstimates[[i]]$GSIestim
			for(j in which(colnames(tempEstim) == c)){
				totGSIestim[,c] <- totGSIestim[,c] + tempEstim[,j]
			}
		}
	}
	rm(columns)
	
		# variableEstim = piVnumbers,
	totpiVnumbers <- list()
	for(v in names(prepData[[1]]$values)){
		columns <- c()
		for(i in 1:numStrata){
			columns <- c(columns, colnames(strataEstimates[[i]]$variableEstim[[v]]))
		}
		columns <- unique(columns)
		totVarestim <- matrix(0, nrow = nrow(strataEstimates[[1]]$variableEstim[[v]]), ncol = length(columns))
		colnames(totVarestim) <- columns
		for(c in columns){
			for(i in 1:numStrata){
				tempEstim <- strataEstimates[[i]]$variableEstim[[v]]
				for(j in which(colnames(tempEstim) == c)){
					totVarestim[,c] <- totVarestim[,c] + tempEstim[,j]
				}
			}
		}
		rm(columns)
		totpiVnumbers[[v]] <- totVarestim
	}
		
		# variableEstimOth = piVnumbersOth,
	totpiVnumbersOth <- list()
	for(v in names(prepData[[1]]$valuesOth)){
		columns <- c()
		for(i in 1:numStrata){
			columns <- c(columns, colnames(strataEstimates[[i]]$variableEstimOth[[v]]))
		}
		columns <- unique(columns)
		totVarestimOth <- matrix(0, nrow = nrow(strataEstimates[[1]]$variableEstimOth[[v]]), ncol = length(columns))
		colnames(totVarestimOth) <- columns
		for(c in columns){
			for(i in 1:numStrata){
				tempEstim <- strataEstimates[[i]]$variableEstimOth[[v]]
				for(j in which(colnames(tempEstim) == c)){
					totVarestimOth[,c] <- totVarestimOth[,c] + tempEstim[,j]
				}
			}
		}
		rm(columns)
		totpiVnumbersOth[[v]] <- totVarestimOth
	}
		# strataName = strataName
	strataNames <- c()
	for(i in 1:numStrata){
		strataNames <- c(strataNames, strataEstimates[[i]]$strataName)
	}
	
	#write summary and CIs - optional
	if (writeSummary){
		#perform this for H,HNC,W
		H <- totClipUnclip[,1]
		HNC <- rep(0, length(H))
		W <- rep(0, length(H))
		if(AI){
			for(i in 1:numStrata){
				nPBT <- prepData[[i]]$nPBT
				if(nPBT == 1){
					HNC <- HNC + strataEstimates[[i]]$piTotEstim[,1]
				} else if (nPBT > 1) {
					HNC <- HNC + rowSums(strataEstimates[[i]]$piTotEstim[,1:nPBT])
				}
				if((nPBT + 1) == ncol(strataEstimates[[i]]$piTotEstim)){
					W <- W + strataEstimates[[i]]$piTotEstim[,(nPBT + 1)]
				} else {
					W <- W + rowSums(strataEstimates[[i]]$piTotEstim[,(nPBT+1):ncol(strataEstimates[[i]]$piTotEstim)])
				}
			}
		}
		output <- data.frame(Group = c("Clipped", "HatcheryNoClip", "Wild"),
							Mean = NA,
							Median = NA,
							Lower = NA,
							Upper = NA,
							stringsAsFactors = FALSE)
		#calculate means
		output$Mean <- c(mean(H), mean(HNC), mean(W))
		#calculate medians
		output$Median <- c(median(H), median(HNC), median(W))
		#calculate CI's
		output[1,4:5] <- quantile(H, c((alpha/2), (1 - (alpha/2))))
		output[2,4:5] <- quantile(HNC, c((alpha/2), (1 - (alpha/2))))
		output[3,4:5] <- quantile(W, c((alpha/2), (1 - (alpha/2))))
		colnames(output)[4:5] <- c(paste0("Lower_(", alpha/2, ")"), paste0("Upper_(", (1 - (alpha/2)), ")"))
		#write output
		write.table(output, paste0(prefix, "RearTypeSummary.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

		#piTot (stock comp)
		output <- data.frame(Group = colnames(totPiTotEstim),
					Mean = apply(totPiTotEstim,2,mean),
					Median = apply(totPiTotEstim,2,median),
					Lower = apply(totPiTotEstim,2,quantile, (alpha/2)),
					Upper = apply(totPiTotEstim,2,quantile, (1 - (alpha/2))),
					stringsAsFactors = FALSE)
		colnames(output)[4:5] <- c(paste0("Lower_(", alpha/2, ")"), paste0("Upper_(", (1 - (alpha/2)), ")"))
		write.table(output, paste0(prefix, "StockCompSummary.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
		#piV
		for (v in names(totpiVnumbers)){
			tempData <- totpiVnumbers[[v]]
			output <- data.frame(Group = colnames(tempData),
						Mean = apply(tempData,2,mean),
						Median = apply(tempData,2,median),
						Lower = apply(tempData,2,quantile, (alpha/2)),
						Upper = apply(tempData,2,quantile, (1 - (alpha/2))),
						stringsAsFactors = FALSE)
			colnames(output)[4:5] <- c(paste0("Lower_(", alpha/2, ")"), paste0("Upper_(", (1 - (alpha/2)), ")"))
			write.table(output, paste0(prefix, v, "CompSummary.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
		}
		#piVOth
		for (v in names(totpiVnumbersOth)){
			tempData <- totpiVnumbersOth[[v]]
			output <- data.frame(Group = colnames(tempData),
						Mean = apply(tempData,2,mean),
						Median = apply(tempData,2,median),
						Lower = apply(tempData,2,quantile, (alpha/2)),
						Upper = apply(tempData,2,quantile, (1 - (alpha/2))),
						stringsAsFactors = FALSE)
			colnames(output)[4:5] <- c(paste0("Lower_(", alpha/2, ")"), paste0("Upper_(", (1 - (alpha/2)), ")"))
			write.table(output, paste0(prefix, v, "CompSummary.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
		}
		
	} # end: if (writeSummary)
	
	
	# to return, per strata estimates, total estimates, and strata names
	return(list(totClipUnclip = totClipUnclip,
					totPiTotEstim = totPiTotEstim,
					totGSIestim = totGSIestim,
					totpiVestim = totpiVnumbers,
					totpiVestimOth = totpiVnumbersOth,
					strataNames = strataNames,
					strataEstimates = strataEstimates
					))
}
