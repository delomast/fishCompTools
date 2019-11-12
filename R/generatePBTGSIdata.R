#' This function simulates data for one strata.
#' 
#' To simulate multiple strata, one can run this repeatedly and add a strata variable
#' 
#' @param sampRate the probability of sampling an individual in the population
#' @param censusSize the total size of the populations
#' @param relSizePBTgroups a numeric vector giving the relative sizes of the adipose fin intact PBT tagged groups
#' @param tagRates a numeric vector giving the true PBT tag rates of the PBT tagged groups, in the
#'  same order as \code{relSizePBTgroups}
#' @param obsTagRates a numeric vector giving the PBT tag rates of the PBT tagged groups to use when creating the
#'  tag rate output, in the same order as \code{relSizePBTgroups}
#' @param physTagRates a numeric vector giving the true physical tagging rate (CWTs, PIT tags) of the PBT tagged groups, in the
#'  same order as \code{relSizePBTgroups}. Physical tagging is assumed to be independent of PBT tagging
#' @param true_clipped the true proportion of the population with a clipped adipose fin
#' @param true_noclip_H the true proportion of the population that is hatchery origin and has an intact adipose fin 
#' @param true_wild the true proportion of the population that is wild origin
#' @param relSizeGSIgroups a numeric vector giving the relative sizes of GSI groups within the wild population
#' @param PBT_GSI_calls a numeric matrix with rows equal to PBT groups, columns equal to GSI group and
#'   values giving the true proportions of the PBT fish that GSI assign to that GSI group
#' @param varMatList a list of numeric matrices, one for each variable. Rows are PBT/GSI groups, in order of PBT groups followed
#'   by GSI groups, columns being the categories for the variable, and values being the proportion of the
#'   group in that category
#'   
#'   
#' @export

generatePBTGSIdata <- function(sampRate = .1, censusSize = 3000, relSizePBTgroups = 1, tagRates = .8, obsTagRates = tagRates, physTagRates = 0,
				    true_clipped = .5, true_noclip_H = .2, true_wild = .3, relSizeGSIgroups = 1, PBT_GSI_calls = 1, varMatList = NA){
	
	# check tag rates all between 0 and 1
	if(sum(tagRates > 1) > 0 || sum(tagRates < 0) > 0){
		stop("tagRates cannot be greater than 1 or less than 0")
	}
	if(sum(obsTagRates > 1) > 0 || sum(obsTagRates < 0) > 0){
		stop("obsTagRates cannot be greater than 1 or less than 0")
	}
	if(sum(physTagRates > 1) > 0 || sum(physTagRates < 0) > 0){
		stop("physTagRates cannot be greater than 1 or less than 0")
	}
	#sample rates check
	if(sampRate > 1 || sampRate < 0){
		stop("sampRate cannot be greater than 1 or less than 0")
	}
	# true props
	propCheck <- mapply(function(x,y,z){
			return(!isTRUE(all.equal(sum(x,y,z), 1))) #floating point arithmetic error avoidance
		}, true_clipped, true_noclip_H, true_wild)
	if(sum(propCheck) > 0){
		stop("true_clipped, true_noclip_H, and true_wild must sum to 1")
	}
	
	# check varArray
	if(!is.na(varMatList) && !is.list(varMatList)){
		stop("varMatList must be a list or NA")
	}
	if(is.list(varMatList)){
		for(i in 1:length(varMatList)){
			row_sums <- apply(varMatList[[i]], 1, sum)
			for(s in row_sums){
				if(!isTRUE(all.equal(s, 1))){
					stop("rows of varMatList matrices must sum to 1")
				}
			}
		}
	}
	#check PBT_GSI_calls
	if (!is.matrix(PBT_GSI_calls)) PBT_GSI_calls <- as.matrix(PBT_GSI_calls)
	row_sums <- apply(PBT_GSI_calls, 1, sum)
	for(s in row_sums){
		if(!isTRUE(all.equal(s, 1))){
			stop("rows of PBT_GSI_calls must sum to 1")
		}
	}


	#simulate data
	# choose number sampled (trapped)
	nSampled <- rbinom(1, censusSize, sampRate)
	if(nSampled < 1){ #prevent strata with no trapped fish - if there were no trapped fish, the user would combine them with others
		nSampled <- 1
	}
	#choose numbers clipped, noclip_H, and wild
		## returns 1col matrix, turn into vector
	nRear <- as.vector(rmultinom(1, nSampled, c(true_clipped, true_noclip_H, true_wild)))
	#choose number of noclip_H that belong to each PBT group
		## returns 1col matrix, turn into vector
	PBTgroupNums <- as.vector(rmultinom(1, nRear[2], relSizePBTgroups))
	#determine number PBT tagged, physTagged, and untagged for each group
	#this mapply returns a matrix with columns as the PBT groups and rows
	#  as numBoth, numPBTonly, numPhysOnly, numUntag
	tagUntagNums <- mapply(function(x,y,z){
		# calculate probabilities of all four groups, assuming that being
		#   physically tagged and PBT tagged are independent
		#prob pbt and phys
		pBoth <- y*z
		#prob phys only
		pPhysOnly <- (1-y)*z
		#prob pbt only
		pPBTonly <- y*(1-z)
		#prob untagged
		pUntag <- 1 - pPBTonly - pPhysOnly - pBoth
		return(rmultinom(1, x, c(pBoth, pPBTonly, pPhysOnly, pUntag)))
		}, PBTgroupNums, tagRates, physTagRates)
	# build datasets
	pbtNames <- paste0("pbtGroup", 1:length(tagRates))

	#initiate dataset with clipped fish
	#columns needed are clip, phystag, pbt group (add gsi and variable columns later)
	strataData <- data.frame(clip = rep("AD", nRear[1]),
						physTag = rep("notag", nRear[1]),
						pbtGroup = rep("Unassigned", nRear[1]),
							stringsAsFactors = FALSE)
	# add GSI and variable columns
	if(is.list(varMatList)) {nCols <- (1+length(varMatList))
	} else {nCols <- 1}
	for(j in 1:nCols){
		strataData <- cbind(strataData, 
						rep(NA, nRear[1])
					)
	}
	# name GSI and variable columns
	if(nCols == 1){ colnames(strataData)[ncol(strataData)] <- "GSI"
	} else {
		colnames(strataData)[(ncol(strataData) - length(varMatList)):ncol(strataData)] <- c("GSI", paste0("Var", 1:length(varMatList)))
	}

	#now add various types of HNC fish
	#first untagged
	for(c in which(tagUntagNums[4,] > 0)){
		numAdd <- tagUntagNums[4,c]
		tempData <- data.frame(clip = rep("AI", numAdd),
								physTag = rep("notag", numAdd),
								pbtGroup = rep("Unassigned", numAdd),
									stringsAsFactors = FALSE
						)
		# add GSI column
		tempData <- cbind(tempData, 
				sample(1:ncol(PBT_GSI_calls), numAdd, replace = TRUE, prob = PBT_GSI_calls[c,])
			)
		colnames(tempData)[ncol(tempData)] <- "GSI"
		# add variable columns
		if(nCols > 1){
			for(j in 1:length(varMatList)){
				tempData <- cbind(tempData, 
								sample(1:ncol(varMatList[[j]]), numAdd, replace = TRUE, prob = varMatList[[j]][c,])
							)
			}
			colnames(tempData) <- colnames(strataData)
		}
		# add to data for the strata
		strataData <- rbind(strataData, tempData)
	}
	#Then PhysOnly
	for(c in which(tagUntagNums[3,] > 0)){
		numAdd <- tagUntagNums[3,c]
		tempData <- data.frame(clip = rep("AI", numAdd),
								physTag = rep("tag", numAdd),
								pbtGroup = rep("Unassigned", numAdd),
									stringsAsFactors = FALSE
						)
		# add GSI column
		tempData <- cbind(tempData, 
				sample(1:ncol(PBT_GSI_calls), numAdd, replace = TRUE, prob = PBT_GSI_calls[c,])
			)
		colnames(tempData)[ncol(tempData)] <- "GSI"
		# add variable columns
		if(nCols > 1){
			for(j in 1:length(varMatList)){
				tempData <- cbind(tempData, 
								sample(1:ncol(varMatList[[j]]), numAdd, replace = TRUE, prob = varMatList[[j]][c,])
							)
			}
			colnames(tempData) <- colnames(strataData)
		}
		# add to data for the strata
		strataData <- rbind(strataData, tempData)
	}
	#pbt only
	#only loop through groups with one or more tagged fish present
	for(c in which(tagUntagNums[2,] > 0)){
		numAdd <- tagUntagNums[2,c]
		tempData <- data.frame(clip = rep("AI", numAdd),
							physTag = rep("notag", numAdd),
							pbtGroup = rep(pbtNames[c], numAdd), #pbtNames is same order as columns
								stringsAsFactors = FALSE)
		# add GSI column
		tempData <- cbind(tempData, 
				sample(1:ncol(PBT_GSI_calls), numAdd, replace = TRUE, prob = PBT_GSI_calls[c,])
			)
		colnames(tempData)[ncol(tempData)] <- "GSI"
		# add variable columns
		if(nCols > 1){
			for(j in 1:length(varMatList)){
				tempData <- cbind(tempData, 
								sample(1:ncol(varMatList[[j]]), numAdd, replace = TRUE, prob = varMatList[[j]][c,])
							)
			}
			colnames(tempData) <- colnames(strataData)
		}
		# add to data for the strata
		strataData <- rbind(strataData, tempData)
	}
	#both
	#only loop through groups with one or more tagged fish present
	for(c in which(tagUntagNums[1,] > 0)){
		numAdd <- tagUntagNums[1,c]
		tempData <- data.frame(clip = rep("AI", numAdd),
							physTag = rep("tag", numAdd),
							pbtGroup = rep(pbtNames[c], numAdd), #pbtNames is same order as columns
								stringsAsFactors = FALSE)
		# add GSI column
		tempData <- cbind(tempData, 
				sample(1:ncol(PBT_GSI_calls), numAdd, replace = TRUE, prob = PBT_GSI_calls[c,])
			)
		colnames(tempData)[ncol(tempData)] <- "GSI"
		# add variable columns
		if(nCols > 1){
			for(j in 1:length(varMatList)){
				tempData <- cbind(tempData, 
								sample(1:ncol(varMatList[[j]]), numAdd, replace = TRUE, prob = varMatList[[j]][c,])
							)
			}
			colnames(tempData) <- colnames(strataData)
		}
		# add to data for the strata
		strataData <- rbind(strataData, tempData)
	}
	
	############################
	#now add wild fish
	############################
	
	#sample from GSI groups
	nWild <- as.vector(rmultinom(1, nRear[3], prob = relSizeGSIgroups))
	for(c in which(nWild > 0)){
		numAdd <- nWild[c]
		tempData <- data.frame(clip = rep("AI", numAdd),
							physTag = rep("notag", numAdd),
							pbtGroup = rep("Unassigned", numAdd), #pbtNames is same order as columns
								stringsAsFactors = FALSE)
		# add GSI column
		tempData <- cbind(tempData, c)
		colnames(tempData)[ncol(tempData)] <- "GSI"
		# add variable columns
		if(nCols > 1){
			gsiVarRow <- length(relSizePBTgroups) + c # gsi groups come after PBT groups in the matirx
			for(j in 1:length(varMatList)){
				tempData <- cbind(tempData, 
								sample(1:ncol(varMatList[[j]]), numAdd, replace = TRUE, prob = varMatList[[j]][gsiVarRow,])
							)
			}
			colnames(tempData) <- colnames(strataData)
		}
		# add to data for the strata
		strataData <- rbind(strataData, tempData)
	}
	
		
	colnames(strataData)[1:3] <- c("AdClip", "PhysTag", "GenParentHatchery")
	#genrerate tag rate input for convenience
	tagRateInput <- data.frame(pbtGroups = pbtNames, tagRates = obsTagRates, stringsAsFactors = FALSE) #give the function observed tag rates
	
	#output all three as a list
	return(list(trapData = strataData, tagRates = tagRateInput))

}
