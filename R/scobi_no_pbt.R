#' generates estiamtes and CIs usings census count and trap data without any PBT
#'
#' \code{scobiNoPBT} generates estiamtes and CIs usings census count and trap data
#' without any PBT data. This function uses parametric bootstrapping to generate CIs
#'
#' Some longer explanation here... someday... probably right before publishing...
#'
#' @return something
#' @export
#'
#'
scobiNoPBT <- function(adultData = NULL, windowData = NULL, Run = "output", RTYPE = "wild", Hierarch_variables = NULL,
                  SizeCut = NULL, alph = 0.1, B = 5000, writeBoot = FALSE, adClipVariable = "AdClip",
				physTagsVariable = "PhysTag", lengthVariable = "LGDFLmm", dataGroupVariable = "WeekNumber"){

	# read in data
	if (is.character(adultData)){
		Fishdata <- read.csv(file = adultData, header = TRUE, na.strings = c("NA",""), stringsAsFactors = FALSE)
	} else {
		Fishdata <- adultData
	}
	if (is.character(windowData)){
		Windata  <- read.csv(file = windowData, header = TRUE, stringsAsFactors = FALSE)
	} else {
		Windata <- windowData
	}
	if (is.character(pbtRates)){
		pbtRate <- read.csv(file = pbtRates, header = TRUE, stringsAsFactors = FALSE)
	} else {
		pbtRate <- pbtRates
	}

	collaps <- Windata[,3]

	# combine (collapse) counts for each group according to user specifies strata
	cat("\n", dataGroupVariable, "is collapsed according to: \n")
	temp <- rbind(Windata[,1],collaps)
	rownames(temp) <- c("Week","Strata")
	print(temp)
	u_collaps <- sort(unique(collaps)) # These are the strata used for the analysis
	WinCounts <- c()
	for(i in u_collaps){
		WinCounts <- c(WinCounts, sum(Windata[Windata[,3] == i,2]))
	}
	WinData <- data.frame(u_collaps,WinCounts) # This is the window data reduced to strata defined by collaps

	if("Strata" %in% colnames(Fishdata)){
		stop("\nStrata is a reserved column name. Please rename this column in your input adultData dataset. Exiting.\n")
	}

	# Recode dataGroupVariable in fish data to collapsed strata
	nullstrat <- c()#strata without any trapped fish
	first <- TRUE
	for(strat in Windata[,1]) {  # Note that the original Windata is used here
		juststrat <- Fishdata[Fishdata[,dataGroupVariable] ==  strat,]	#trap data with only fish in the current strata
		if( nrow(juststrat) == 0 ) {#if no fish in the strata
			nullstrat <- c(nullstrat,strat)
		} else {
			juststrat$Strata <- Windata[Windata[,1] == strat,3]	#add in the strata that matches that week
			if(first) {
				FishData <- juststrat
				first <- FALSE
			} else {
				FishData <- rbind(FishData,juststrat)
			}
		}

	}
	##### input data is Fishdata, working data with Strata field is FishData, whith a capital D

	# Check to see if any of the original weeks are missing trapped fish and if any recoded week are missing fish
	if(length(nullstrat) == 0) {
		cat("\nThere are trapped fish for every week\n")
	} else if (sum(table(FishData$Strata) == 0) == 0) {
		cat("\nWeeks ", nullstrat, " have no trapped fish but it appears that \n")
		cat(" these weeks are collapsed with weeks having trapped fish.\n")
	} else{
		cat("\nWeeks ",nullstrat," have no trapped fish and it appears that \n")
		cat("some of these weeks are NOT collapsed with weeks having trapped fish.\n")
		cat("Here are the totals of trapped fish for each collapsed strata: \n" )
		print(table(FishData$Strata))
		stop("\nExiting.")
	}

	rm(Windata)  # save memory, new window data has capital D
	rm(Fishdata) #  save memory, new fish data has capital D

	# estimate number of H, HNC, and W
	# clip variable determines H vs (HNC + W)
	# phystag variable determines HNC vs W and assumes tagging and detection rates of phystags is 100%

	#discard observations with:
	##no entry for ad-clip
	num_obs <- nrow(FishData)
	FishData <- FishData[!is.na(FishData[,adClipVariable]),]
	cat("\nRemoved", num_obs - nrow(FishData), "observations for missing ad-clip information.\n")

	# Define large and small if one of the factors is fork length
	if(lengthVariable %in% Hierarch_variables) {
		FishData[!is.na(FishData[,lengthVariable]) & FishData[,lengthVariable] < SizeCut, lengthVariable] <- "Sm"
		FishData[!is.na(FishData[,lengthVariable]) & FishData[,lengthVariable] != "Sm", lengthVariable] <- "Lg"
	}

	h_hnc_w_estimates <- matrix(nrow = nrow(WinData), ncol = 4)	## this will have proportion of each type trapped each strata
	colnames(h_hnc_w_estimates) <- c("strata", "clipped_hatchery", "unclipped_hatchery", "wild")
	h_hnc_w_estimates[,1] <- WinData[,1]

	#for each strata
	for(i in 1:nrow(h_hnc_w_estimates)){
		prop_all_data <- FishData[FishData$Strata == h_hnc_w_estimates[i,1],]
		#calculate number clipped
		clip_count <- sum(prop_all_data[,adClipVariable] == "AD")

		#unclipped only data
		ai_data <- prop_all_data[prop_all_data[,adClipVariable] == "AI",]
		# calculate total number unclipped
		aiCount <- nrow(ai_data)

		#calculate number HNC and wild, if possible
		if(aiCount > 0 && sum(!is.na(ai_data[,physTagsVariable])) == 0){
			cat("Warning. No AI fish with data for", physTagsVariable, "in strata", h_hnc_w_estimates[i,1], "\n")
			cat("Assuming all fish in this strata are wild.\n")
			ai_data[,physTagsVariable] <- "notag"
		}

		# calculate number AI with physical tag - these are considered HNC
		phystag_count <- sum(!is.na(ai_data[,physTagsVariable]) & ai_data[,physTagsVariable] == "tag")
		#calculate number AI without phys tag - these are considered W
		noPhystag_count <- sum(!is.na(ai_data[,physTagsVariable]) & ai_data[,physTagsVariable] == "notag")
		# expand to account for fish with missing phys tag information
		phystag_count <- aiCount * (phystag_count / (phystag_count + noPhystag_count)) #denom can't be zero b/c if all NA, assume all wild
		noPhystag_count <- aiCount - phystag_count

		total <- sum(noPhystag_count, phystag_count, clip_count)

		h_hnc_w_estimates[i,2] <- clip_count/total
		h_hnc_w_estimates[i,3] <- phystag_count/total
		h_hnc_w_estimates[i,4] <- noPhystag_count/total
		rm(ai_data)
	}
	#use proportions combined with census (window) data to estimate counts of each category
	h_hnc_w_count_est <- h_hnc_w_estimates[,2:4]*WinData[,2]
	if(!is.matrix(h_hnc_w_count_est)){
		h_hnc_w_count_est <- t(as.matrix(h_hnc_w_count_est))
	}
	category_totals <- apply(h_hnc_w_count_est,2,sum)
	h_hnc_w_count_est <- cbind(h_hnc_w_estimates[,1], h_hnc_w_count_est)
	colnames(h_hnc_w_count_est)[1] <- "strata"

	#### status of variables here
	# h_hnc_w_estimates : proportions of each type in each strata
	# h_hnc_w_count_est : Estimates of each type in each strata after multiplying by census values
	# category_totals : sum of each type across strata, after multiplying by census values

	cat("Estimated totals for\n\tClipped:\t", category_totals[1], "\n\tHatchery No-clip:\t", category_totals[2], "\n\tWild:\t", category_totals[3], "\n")

	# now calculate hierarchical variable estimates





	# now bootstrap to obtain CIs








}
