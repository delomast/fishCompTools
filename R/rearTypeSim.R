#' Performs simulations to estimate precision, bias, and coverage of clipped, noclip_H, and wild estimator
#'
#'
#' @param numSims the number of simulations to perform
#' @param sampRates sampling rates (dam trap rates) for each strata. Order must match
#'   order of \code{censusSizes}
#' @param censusSizes population sizes (dam window counts) for each strata. order must
#'   match order of \code{sampRates}
#' @param relSizePBTgroups relative sizes of PBT groups for noclip_H fish ONLY, because
#'   PBT info is not used to determine the proportion clipped
#' @param tagRates the true tag rates for each PBT group in the same order as \code{relSizePBTgroups}
#' @param obsTagRates the tag rates to use in the estimation for each PBT group in the
#'   same order as \code{relSizePBTgroups}
##' @param physTagRates the true physical tagging rate for each PBT groups for noclip_H fish ONLY
#'   in the same order as \code{relSizePBTgroups}
#' @param true_clipped the true proportion of fish that are clipped, hatchery origin in each strata
#'   order must match the order of \code{censusSizes}
#' @param true_noclip_H the true proportion of fish that are unclipped, hatchery origin in each strata
#'   order must match the order of \code{censusSizes}
#' @param true_wild the true proportion of fish that are unclipped, wild origin in each strata
#'   order must match the order of \code{censusSizes}
#' @param alph the value of \code{alph} to use for CI estimation
#' @param B the number of bootstraps to perform for CI estimation
#'
#' @return a dataframe with estimates and CIs for all three categories from each iteration
#' @export
#'

rearTypeSim <- function(numSims, sampRates, censusSizes, relSizePBTgroups, tagRates, obsTagRates, physTagRates,
				    true_clipped, true_noclip_H, true_wild, alph = .1, B = 1000){

	# input checking
	if((numSims %% 1) != 0 || numSims < 1){
		stop("numsims must be an integer greater than 0")
	}
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
	if(sum(sampRates > 1) > 0 || sum(sampRates < 0) > 0){
		stop("sampRates cannot be greater than 1 or less than 0")
	}
	# true props
	propCheck <- mapply(function(x,y,z){
			return(!isTRUE(all.equal(sum(x,y,z), 1))) #floating point arithmetic error avoidance
		}, true_clipped, true_noclip_H, true_wild)
	if(sum(propCheck) > 0){
		stop("true_clipped, true_noclip_H, and true_wild must sum to 1 for all strata")
	}
	# alph
	if (alph >= 1 || alph <= 0){
		stop("alph must be between 0 and 1")
	}
	# check that B is greater than 0, or change function to work with no bootstraps
	if ((B %% 1) != 0 || B < 1){
		stop("B must be an integer greater than 0")
	}

	# create summary data storage variables
	#just creating blank and will rbind all output
	rearStoreAll <- data.frame(NULL, stringsAsFactors = FALSE)

	#count number of strata
	nStrata <- length(sampRates)

	#for each simulation
	for(ns in 1:numSims){
		#empty dataframe for this simulation
		simData <- data.frame(NULL, stringsAsFactors = FALSE)
		#simulate data
		for(i in 1:nStrata){
			# choose number sampled (trapped)
			nSampled <- rbinom(1, censusSizes[i], sampRates[i])
			if(nSampled < 1){ #prevent strata with no trapped fish - if there were no trapped fish, the user woudl combine them with others
				nSampled <- 1
			}
			#choose numbers clipped, noclip_H, and wild
				## returns 1col matrix, turn into vector
			nRear <- as.vector(rmultinom(1, nSampled, c(true_clipped[i], true_noclip_H[i], true_wild[i])))
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

			#initiate dataset with clipped and wild fish
			#columns needed are clip, phystag, pbt group
			strataData <- data.frame(clip = c(rep("AD", nRear[1]), rep("AI", nRear[3])),
								physTag = c(rep("notag", (nRear[1] + nRear[3]))),
								pbtGroup = c(rep("Unassigned", (nRear[1] + nRear[3]))),
									stringsAsFactors = FALSE)
			#now add various types of HNC fish
					#this mapply returns a matrix with columns as the PBT groups and rows
					#  as numBoth, numPBTonly, numPhysOnly, numUntag
					# tagUntagNums
			#first untagged
			numAdd <- sum(tagUntagNums[4,])
			strataData <- rbind(strataData,
							data.frame(clip = rep("AI", numAdd),
									physTag = rep("notag", numAdd),
									pbtGroup = rep("Unassigned", numAdd),
										stringsAsFactors = FALSE)
							)
			#Then PhysOnly
			numAdd <- sum(tagUntagNums[3,])
			strataData <- rbind(strataData,
							data.frame(clip = rep("AI", numAdd),
									physTag = rep("tag", numAdd),
									pbtGroup = rep("Unassigned", numAdd),
										stringsAsFactors = FALSE)
							)
			#pbt only
			#only loop through groups with one or more tagged fish present
			for(c in which(tagUntagNums[2,] > 0)){
				numAdd <- tagUntagNums[2,c]
				strataData <- rbind(strataData,
							data.frame(clip = rep("AI", numAdd),
									physTag = rep("notag", numAdd),
									pbtGroup = rep(pbtNames[c], numAdd), #pbtNames is same order as columns
										stringsAsFactors = FALSE)
							)
			}
			#both
			#only loop through groups with one or more tagged fish present
			for(c in which(tagUntagNums[1,] > 0)){
				numAdd <- tagUntagNums[1,c]
				strataData <- rbind(strataData,
							data.frame(clip = rep("AI", numAdd),
									physTag = rep("tag", numAdd),
									pbtGroup = rep(pbtNames[c], numAdd), #pbtNames is same order as columns
										stringsAsFactors = FALSE)
							)
			}
			#add simulated data for the strata to the main data frame for this simulations
			simData <- rbind(simData, cbind(i, strataData)) #columns are: strata, clip, physTag, pbtGroup
		}
		colnames(simData) <- c("WeekNumber", "AdClip", "PhysTag", "GenParentHatchery")

		windowData <- data.frame(week = 1:nStrata,
							count = censusSizes,
							collapse = 1:nStrata, stringsAsFactors = FALSE) #data is simulated by strata
		tagRateInput <- data.frame(pbtGroups = pbtNames, tagRates = obsTagRates, stringsAsFactors = FALSE) #give the function observed tag rates

		#calling SCOBI_deux_fast while repressing output to the screen
		invisible(capture.output(SCOBI_deux_fast(adultData = simData, windowData = windowData, Run = "tempSim", RTYPE = "wild", Hierarch_variables = NULL,
	                  SizeCut = NULL, alph = alph, B = B, writeBoot = FALSE, pbtRates = tagRateInput,
					screenOutput = "screen.txt")
		))

		#read output and add to summary list
		rearOutput <- read.table("tempSim_CI_Rearing.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
		rearStoreAll <- rbind(rearStoreAll, cbind(ns, rearOutput))

	}

	#remove SCOBI_deux output files - they will be written over with each simulation, so only need to explicitly delete at the end
	toRemove <- c("tempSim_CI_Rearing.txt", "tempSim_Rearing.txt", "tempSim_screen.txt")
	invisible(lapply(toRemove, file.remove))

	return(rearStoreAll)
}
