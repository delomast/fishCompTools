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

	# check PhysTag column to make sure all values are recognized
	checkInput <- unique(FishData[!is.na(FishData[,physTagsVariable]), physTagsVariable])
	if (sum(!(checkInput %in% c("tag", "notag"))) > 0){
		errorMessage <- paste0("Unrecognized entry \"", checkInput[!(checkInput %in% c("tag", "notag"))], "\" in the physTagsVariable column.")
		stop(errorMessage)
	}

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

	#save some information for bootstrapping
	save_for_boot_h_hnc_w <- matrix(nrow = nrow(WinData), ncol = 5) #rows are strata
	colnames(save_for_boot_h_hnc_w) <- c("numADAI", "propAD", "propAIAndTagInfo", "numTagInfo", "propTag") #column names to help debug

	#for each strata
	for(i in 1:nrow(h_hnc_w_estimates)){
		prop_all_data <- FishData[FishData$Strata == h_hnc_w_estimates[i,1],]
		save_for_boot_h_hnc_w[i,1] <- nrow(prop_all_data) #save number of fish with ad-clip status
		#calculate number clipped
		clip_count <- sum(prop_all_data[,adClipVariable] == "AD")
		save_for_boot_h_hnc_w[i,2] <- clip_count / save_for_boot_h_hnc_w[i,1] #save proportion that are clipped

		#unclipped only data
		ai_data <- prop_all_data[prop_all_data[,adClipVariable] == "AI",]
		# calculate total number unclipped
		aiCount <- nrow(ai_data)

		if(aiCount == 0){ #if no unclipped fish, set all to zero
			save_for_boot_h_hnc_w[i,3] <- 0 #save prop of AI with phystag Info
			save_for_boot_h_hnc_w[i,4] <- 0 #save number of AI with phystag Info
			save_for_boot_h_hnc_w[i,5] <- 0 #save prop of AI that are HNC (tagged)

			h_hnc_w_estimates[i,2] <- 1 #clip proportion is 1 - all fish
			h_hnc_w_estimates[i,3] <- 0
			h_hnc_w_estimates[i,4] <- 0
			rm(ai_data)
		}

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

		save_for_boot_h_hnc_w[i,3] <- (phystag_count + noPhystag_count) / aiCount #save prop of AI with phystag Info
		save_for_boot_h_hnc_w[i,4] <- (phystag_count + noPhystag_count) #save number of AI with phystag Info
		save_for_boot_h_hnc_w[i,5] <- phystag_count / (phystag_count + noPhystag_count) #save prop of AI that are HNC (tagged)

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
	#H, HNC, and W
	boot_h_hnc_w <- matrix(data = 0, nrow = B, ncol = 3) #storage for bootstrap counts, this is total for the whole run, so adding each strata to is as calculated
	colnames(boot_h_hnc_w) <- c("h", "hnc", "w")
	for(i in 1:nrow(WinData)){ # for each strata
		# first AD vs AI
		clip_boot <- rbinom(B, save_for_boot_h_hnc_w[i,1], save_for_boot_h_hnc_w[i,2]) # this is number of clipped "observed" in each B
		boot_num_AI <- save_for_boot_h_hnc_w[i,1] - clip_boot # this is the number AI "observed" in each B
		# second phystag vs not
		phystag_boot <- rbinom(B, save_for_boot_h_hnc_w[i,4], save_for_boot_h_hnc_w[i,5]) # this is nmber of phystag fish "observed" in each B
		# expand phystag to apply to all AI "observed" in each B
		exp_HNC_boot <- boot_num_AI * (phystag_boot / save_for_boot_h_hnc_w[i,4]) # this is expanded number of HNC "trapped" fish
		exp_w_boot <- boot_num_AI - exp_HNC_boot # this is the expanded number of W "trapped" fish
		# calculate proportions
		prop_clip <- clip_boot / save_for_boot_h_hnc_w[i,1]
		prop_HNC <- exp_HNC_boot / save_for_boot_h_hnc_w[i,1]
		prop_W <- exp_w_boot / save_for_boot_h_hnc_w[i,1]
		# multply by window count and add to running total
		boot_h_hnc_w[,1] <- boot_h_hnc_w[,1] + (prop_clip * WinData[i,2])
		boot_h_hnc_w[,2] <- boot_h_hnc_w[,2] + (prop_HNC * WinData[i,2])
		boot_h_hnc_w[,3] <- boot_h_hnc_w[,3] + (prop_W * WinData[i,2])
	}
	#make output matrix and write to file
	CI_h_hnc_w <- matrix(nrow = 3, ncol = 8)
	colnames(CI_h_hnc_w) <- c("Group", "Estimate", paste0("Lower_(", alph/2, ")"), paste0("Upper_(", 1-(alph/2), ")"), "Percen_half_width", "Lower_simul", "Upper_simul", "Percen_half_width_simul")
	CI_h_hnc_w[,1] <- c("clipped_hatchery", "unclipped_hatchery", "wild")

	sim_ci <- simulConfInt(boot_h_hnc_w, alph) #calculate simultaneous CIs

	for(i in 1:3){
		CI_h_hnc_w[i,2] <- category_totals[i]
		CI_h_hnc_w[i,3:4] <- quantile(boot_h_hnc_w[,i],c(alph/2, 1-(alph/2)))
		CI_h_hnc_w[i,5] <- round(100 * ( (as.numeric(CI_h_hnc_w[i,4]) - as.numeric(CI_h_hnc_w[i,3]) )/ (category_totals[i]*2)),2)
		CI_h_hnc_w[i,6] <- sim_ci[i,1]
		CI_h_hnc_w[i,7] <- sim_ci[i,2]
		CI_h_hnc_w[i,8] <- round(100 * ( (sim_ci[i,2] - sim_ci[i,1] )/ (category_totals[i]*2)),2)
	}

	write.table(CI_h_hnc_w, paste0(Run, "_CI_Rearing.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)







}
