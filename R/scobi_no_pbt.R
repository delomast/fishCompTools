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
	rm(adultData) # remove in case passed as a dataframe and is large
	rm(windowData)


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


	##################################################################
	# estimate number of H, HNC, and W
	# clip variable determines H vs (HNC + W)
	# phystag variable determines HNC vs W and assumes tagging and detection rates of phystags is 100%
	##################################################################

	h_hnc_w_estimates <- matrix(nrow = nrow(WinData), ncol = 4)	## this will have proportion of each type trapped each strata
	colnames(h_hnc_w_estimates) <- c("strata", "clipped_hatchery", "unclipped_hatchery", "wild")
	h_hnc_w_estimates[,1] <- WinData[,1]
	trap_counts <- h_hnc_w_estimates # this will have counts of each type trapped

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
		trap_counts[i,2] <- clip_count #save trapped clipped count

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

			trap_counts[i,3] <- 0 #save trapped hnc count
			trap_counts[i,4] <- 0 #save trapped w count

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

		trap_counts[i,3] <- phystag_count #save trapped hnc count
		trap_counts[i,4] <- noPhystag_count #save trapped w count

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

	#write output to files
	sink(paste0(Run,"_Rearing.txt"))
	cat("Counts of fish of each type in the trap data\n")
	write.table(trap_counts, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
	cat("\n\nEstimates of each type in the run\n")
	write.table(h_hnc_w_count_est, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
	cat("Totals", category_totals[1], category_totals[2], category_totals[3], sep = "\t")
	sink()



	########################################################
	# now calculate hierarchical variable estimates
	########################################################

	if(length(Hierarch_variables) > 0){
		#select rearing type of interest
		hierarch_data <- FishData
		#select rearing type of interest(clipped, noclip_h, wild)
		if (RTYPE == "clipped"){
			hierarch_data <- hierarch_data[hierarch_data[,adClipVariable] == "AD",]
		} else if (RTYPE == "noclip_H"){
			hierarch_data <- hierarch_data[hierarch_data[,adClipVariable] == "AI" & !is.na(hierarch_data[,physTagsVariable]) & hierarch_data[,physTagsVariable] == "tag",]
		} else if (RTYPE == "wild"){
			hierarch_data <- hierarch_data[hierarch_data[,adClipVariable] == "AI" & !is.na(hierarch_data[,physTagsVariable]) & hierarch_data[,physTagsVariable] == "notag",]
		} else {
			stop("Unrecognized RTYPE. Must be one of: \"clipped\", \"noclip_h\", \"wild\"\n")
		}

		#check that there is data
		if(nrow(hierarch_data) < 1){
			stop("No fish present in the dataset of the RTYPE that you chose for the hierarchical analysis. Exiting.")
		}


		hierarch_output_list <- list()	#storage of output, flexible for different numbers of factors
		column_cats <- c()	#this is the list of factors analyzed so far plus the current
							#it is therefore the list of columns that the output for the current factor will have
		#for each factor
		for (factor in Hierarch_variables){
			#remove NAs
			removed <- sum(is.na(hierarch_data[,factor]))
			hierarch_data <- hierarch_data[!is.na(hierarch_data[,factor]),]
			cat("\nRemoved ", removed, " observations for missing data for ", factor, ".\n", sep="")

			column_cats <- c(column_cats, factor)	#list columns for this hierarchical output variable
			if(nrow(hierarch_data) < 1){
				stop("There are no fish in any strata with complete data for factors", column_cats, "   Exiting.")
			}

			if(length(column_cats) < 2){
				unique_categories <- list(unique(hierarch_data[,column_cats]))
				names(unique_categories) <- c(column_cats)
			} else {
				unique_categories <- apply(hierarch_data[,column_cats],2,unique)	# get lists of all unique categories for each variable in this analysis, combinations will be the entries for each strata in the output
				# if in rare case all column_cats have the same number of categories, a matrix is returned
				# it must be turned into a list for later functionality
				if(!is.list(unique_categories)){
					unique_categoriesMatrix <- unique_categories
					unique_categories <- list()
					for(i in 1:ncol(unique_categoriesMatrix)){
						unique_categories[[colnames(unique_categoriesMatrix)[i]]] <- unique_categoriesMatrix[,i]
					}
				}
			}

			hierarch_estimates <- as.data.frame(matrix(nrow = 0, ncol = (length(column_cats) + 2)))	## this will have columns for: strata, categories, proportion of each type trapped
			hierarch_estimates_count <- as.data.frame(matrix(nrow = 0, ncol = (length(column_cats) + 2)))	## this will have columns for: strata, categories, count of each type trapped

			replace_NA <- FALSE #this will record whether or not the function needs to check for and replace NA entries due to lack of data for one or more strata
			#for each strata
			for(s in h_hnc_w_count_est[,1]){
				strata_data <- hierarch_data[hierarch_data$Strata == s,]	#select strata
				#make sure there is data
				if(nrow(strata_data) == 0){
					cat("\nWarning. No fish with complete data for factors", column_cats, "in strata", s, "\n Using the mean proportions of the entire run for this strata.\n")
					cat("This will cause an error during bootstrapping. Please consider combining this strata with an adjacent strata.\n")
					replace_NA <- TRUE
					#write all NA as temporary measure to be fixed at the end
					#build matrix with all combinations of variables
					temp_output <- build_all_combos(unique_categories)
					temp_output <- cbind(s, temp_output, NA)
					colnames(temp_output) <- c("Strata", column_cats, "Proportion_of_strata")
					hierarch_estimates <- rbind(hierarch_estimates, temp_output, stringsAsFactors = FALSE)
					colnames(temp_output) <- c("Strata", column_cats, "Count_trapped")
					hierarch_estimates_count <- rbind(hierarch_estimates_count, temp_output, stringsAsFactors = FALSE)
					next #go to next strata
				}
				#calculate proportions of each category

				#build matrix with all combinations of variables
				temp_output <- build_all_combos(unique_categories)
				temp_output <- cbind(s, temp_output, -9) #put -9 in as a placeholder
				colnames(temp_output) <- c("Strata", column_cats, "Proportion_of_strata")

				col_num <- ncol(temp_output)
				matching_data <- strata_data[ , c("Strata", names(unique_categories))] #get data with just columns that you need, in the order you need them
				for(i in 1:nrow(temp_output)){
					#calculate number of observations that are an exact match for the category
					temp_output[i,col_num] <- sum(find_matching_rows(matching_data, temp_output[i,1:ncol(matching_data)]))
				}
				colnames(temp_output) <- c("Strata", column_cats, "Count_trapped")
				hierarch_estimates_count <- rbind(hierarch_estimates_count, temp_output, stringsAsFactors = FALSE)
				colnames(temp_output) <- c("Strata", column_cats, "Proportion_of_strata")
				temp_output[,col_num] <- as.numeric(temp_output[,col_num]) / nrow(strata_data)
				hierarch_estimates <- rbind(hierarch_estimates, temp_output, stringsAsFactors = FALSE)
			}#end of for each strata loop

			#replace NAs written to strata that did not have any data
			if(replace_NA){
				#calculate the means
				all_categories <- as.matrix(build_all_combos(unique_categories))
				#calculate means
				means <- rep(-9, nrow(all_categories)) #-9 as a placeholder
				bool_mat <- as.matrix(hierarch_estimates[,2:(ncol(all_categories) + 1)])
				for(i in 1:nrow(all_categories)){
					#calculate number of observations that are an exact match for the category
					bool_temp <- find_matching_rows(bool_mat, all_categories[i,])
					means[i] <- mean(as.numeric(hierarch_estimates[bool_temp,col_num]), na.rm = TRUE)
				}
				for(i in which(is.na(hierarch_estimates[,col_num]))){
					#find the matching entry in all_categories, order of means is the same
					bool_temp <- find_matching_rows(all_categories, unlist(hierarch_estimates[i,2:(col_num - 1)]))
					#make proportion for this strata the mean over strata for the run
					hierarch_estimates[i,col_num] <- means[bool_temp]
				}
			}

			#save output
			hierarch_output_list[[factor]] <- list(hierarch_estimates, hierarch_estimates_count)
		}# end of hierachical factor loop


		# apply window counts to generate estimates of total counts
		# apply proportions in a hierarchical fashion
		# so that numbers are estimated using as much data as possible
		####get column number to use when pulling out estimates from the H HNC W decomposition
		if (RTYPE == "clipped") {
			col_n <- 2
		} else if (RTYPE == "noclip_H") {
			col_n <- 3
		} else {
			col_n <- 4
		}

		for(i in 1:length(hierarch_output_list)){
			props <- hierarch_output_list[[i]][[1]]
			props_col <- ncol(props)
			colnames(props)[props_col] <- "Total_number_in_run"

			if(i == 1){
				#first factor, so pull numbers from the H HNC W estimates
				for (j in 1:nrow(h_hnc_w_count_est)){
					props[props[,1] == h_hnc_w_count_est[j,1],props_col] <- as.numeric(props[props[,1] == h_hnc_w_count_est[j,1],props_col]) * as.numeric(h_hnc_w_count_est[j,col_n])
				}
				hierarch_output_list[[i]][[3]] <- props
			} else {
				replace_NA_2 <- c() #keep track of row number for which there are NAs to replace with mean
				props_col_2 <- (props_col-2) # minus 2 because we are expressing everything in terms of the last estimates groups
				for (j in 1:nrow(last_estim)){
					#find matching categories
					## have to use unlist here b/c last_estim is a dataframe, and it needs to be passed as a one dimensional vector, not a one row dataframe
					bool_temp <- find_matching_rows(props[,1:props_col_2], unlist(last_estim[j,1:props_col_2]))
					if (sum(as.numeric(props[bool_temp,props_col])) != 0){
						# take proportions of strata and turn into proportions of group as analyzed by the previous factor
						props[bool_temp,props_col] <- (as.numeric(props[bool_temp,props_col]) / sum(as.numeric(props[bool_temp,props_col])))
						# multiply by estimate from analysis for the previous factor
						props[bool_temp,props_col] <- as.numeric(props[bool_temp,props_col]) * as.numeric(last_estim[j,(props_col-1)])
					} else if (as.numeric(last_estim[j,(props_col-1)]) != 0) { #in cases where there is data to estimate a non-zero value for the previous variable, but there is not data to estimate the current category conditional on the previous variables
						replace_NA_2 <- c(replace_NA_2, j)
						props[bool_temp,props_col] <- NA
						cat("\nWarning. There are no fish in")
						print(last_estim[j,1:(props_col-2)], row.names = FALSE)
						cat("with data for", factor, "\nUsing the mean composition for all other fish in this strata to estimate this group.\n")
					}
				}

				#replace NAs with means for all other categories if necessary
				if (length(replace_NA_2) > 0){
					#calculate means for all categories
					cats_means <- unique(props[,c(1,(props_col-1))]) #this will be the order of means
					cats_means <- cbind(cats_means, -9)
					colnames(cats_means)[3] <- "Total_number_in_run"
					for(j in 1:nrow(cats_means)){
						#have to unlist cat_means to prevent being passed as a one row dataframe
						bool_temp <- find_matching_rows(props[,c(1, (props_col-1))], unlist(cats_means[j,1:2]))
						cats_means[j,3] <- mean(as.numeric(props[bool_temp,props_col]), na.rm = TRUE)
					}
					for(j in unique(cats_means[,1])){
						bool_temp <- cats_means[,1] == j
						cats_means[bool_temp,3] <- as.numeric(cats_means[bool_temp,3]) / sum(as.numeric(cats_means[bool_temp,3]))
					}
					#apply means to categories with no data to estimate
					for(j in 1:length(replace_NA_2)){
						#build matrix to add on to props with the missing categories
						add_on <- cbind(last_estim[replace_NA_2[j],1:(props_col-2)], cats_means[cats_means[,1] == last_estim[replace_NA_2[j],1],2:3], row.names = NULL)
						add_on[,props_col] <- as.numeric(add_on[,props_col]) * as.numeric(last_estim[replace_NA_2[j],(props_col-1)])
						props <- rbind(props, add_on)
					}
					# remove lines with NA, they have been replaced by adding on rows
					props <- props[!is.na(props[,props_col]),]
				}
				#save as a third matrix in list for output
				hierarch_output_list[[i]][[3]] <- props
			}

			#use estimated number in the last estimate to estimate numbers in the next estimate
			last_estim <- props
		}

		# Write hierarchical estimate outputs
		for(i in 1:length(hierarch_output_list)){
			## the propor_hier output file can be helpful for troublshooting but is not very informative for the end-user
			#write.table(hierarch_output_list[[i]][[1]], paste0(Run, "_Propor_Hier_", names(hierarch_output_list)[i], ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
			write.table(hierarch_output_list[[i]][[2]], paste0(Run, "_Trap_Counts_Hier_", names(hierarch_output_list)[i], ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
			write.table(hierarch_output_list[[i]][[3]], paste0(Run, "_Estim_Totals_Hier_", names(hierarch_output_list)[i], ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
		}

		#calculate total estimates over the entire run
		for(i in 1:length(hierarch_output_list)){
			estimates_by_strata <- hierarch_output_list[[i]][[3]]
			categories_to_sum <- estimates_by_strata[,2:(ncol(estimates_by_strata) - 1)]
			categories_to_sum <- unique(categories_to_sum)
			if(!is.matrix(categories_to_sum) && !is.data.frame(categories_to_sum)){
				categories_to_sum <- as.matrix(categories_to_sum)
				colnames(categories_to_sum) <- colnames(hierarch_output_list[[i]][[3]])[2]
				#this happens when there is only variable
				#causes problems below if categories to sum is not passed as a matrix
			}
			totals_for_run <- rep(-9, nrow(categories_to_sum))
			est_by_s_mat_comp <- as.matrix(estimates_by_strata[,colnames(categories_to_sum)])
			for (j in 1:nrow(categories_to_sum)){
				bool_temp <- find_matching_rows(est_by_s_mat_comp, unlist(categories_to_sum[j,]))
				totals_for_run[j] <- sum(as.numeric(estimates_by_strata[bool_temp,"Total_number_in_run"]))
			}
			totals_for_run <- cbind(categories_to_sum, totals_for_run)
			colnames(totals_for_run)[ncol(totals_for_run)] <- "Estimated_total_for_run"
			write.table(totals_for_run, paste0(Run, "_Estim_Grand_Totals_Hier_", names(hierarch_output_list)[i], ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
			hierarch_output_list[[i]][[4]] <- totals_for_run #save for use in bootstrapping
		}
	} # end of if length(Hierarch_variable) > 0



	##################################################################
	# now bootstrap to obtain CIs
	##################################################################

	if(is.null(B) || is.na(B) || B < 1){
		return("B is less than 1. Not performing any bootstrapping. ScobiNoPBT is complete.")
	}

	#H, HNC, and W
	boot_h_hnc_w <- matrix(data = 0, nrow = B, ncol = 3) #storage for bootstrap counts, this is total for the whole run, so adding each strata to is as calculated
	colnames(boot_h_hnc_w) <- c("h", "hnc", "w")
	hier_boot_type_total <- matrix(nrow = nrow(WinData), ncol = B) # storage of type of interest for hierarchical bootstrapping

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
		# multply by window count and add type of interest to running total
		boot_h_hnc_w[,1] <- boot_h_hnc_w[,1] + (prop_clip * WinData[i,2])
		boot_h_hnc_w[,2] <- boot_h_hnc_w[,2] + (prop_HNC * WinData[i,2])
		boot_h_hnc_w[,3] <- boot_h_hnc_w[,3] + (prop_W * WinData[i,2])
		# save by strata for use in hier boot
		if (RTYPE == "clipped") {
			hier_boot_type_total[i,] <- (prop_clip * WinData[i,2])
		} else if (RTYPE == "noclip_H") {
			hier_boot_type_total[i,] <- (prop_HNC * WinData[i,2])
		} else {
			hier_boot_type_total[i,] <- (prop_W * WinData[i,2])
		}
	}

	#write bootstrap estimates
	if (writeBoot){
		write.table(boot_h_hnc_w, paste0(Run, "_Boot_Rear.txt"),
				  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
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

	# Hierarchical variables
	if(length(Hierarch_variables) > 0){
		bootHier <- list() #this will have bootstrap grand total "observations" for each set of factors
		for (i in 1:length(hierarch_output_list)){ #for each factor
			hierTrapCounts <- hierarch_output_list[[i]][[2]] # trap counts - this is used to calculate sample size
			hierEstimTotals <- hierarch_output_list[[i]][[3]] # this is Estim Totals
			countCol <- ncol(hierTrapCounts) # this is column of hierTrapCounts that has the counts
			bootStrata <- matrix(nrow = 0, ncol = B) #to save bootstrap "observed proportions" per strata
			bootMeta <- matrix(nrow = 0, ncol = (countCol - 1)) # to save strata/category order
			for(j in 1:nrow(WinData)){ #for each strata
				strataTrapCounts <- hierTrapCounts[hierTrapCounts$Strata == WinData[j,1],] #trap counts for that strata
				strataEstimTotals <- hierEstimTotals[hierEstimTotals$Strata == WinData[j,1],] #estim totals for that strata
				# calculate sample sizes from trap counts: sum(strataTrapCounts[,countCol])
				# calculate probabilities from Estim Totals: strataEstimTotals[,countCol]
					#this takes into account proportiosn as estimated using observations with missing data points
					#can then feed into one multinomial sample, preventing incompatibilites when
					#one sample is zero for a category and non-zero for a subcategory (or vice versa)
				# multinomial sample
				#then turn into proportion of the strata by dividing by sample size
				sampSize <- sum(as.numeric(strataTrapCounts[,countCol]))
				bootStrata <- rbind(bootStrata, (rmultinom(B, sampSize, strataEstimTotals[,countCol]) / sampSize))
				bootMeta <- rbind(bootMeta, strataEstimTotals[,1:(countCol - 1)])
			}
			colnames(bootMeta) <- colnames(hierEstimTotals)[1:(countCol - 1)]
			colnames(bootStrata) <- paste0("B_", 1:B)
			# add strata and category labels
			bootStrata <- cbind(bootMeta, bootStrata, stringsAsFactors = FALSE)

			#calculate per strata totals for each B
			for(j in 1:nrow(WinData)){
				strataTotal <- hier_boot_type_total[j,] #totals for that strata
				boots <- bootStrata[bootStrata$Strata == WinData[j,1], countCol:ncol(bootStrata)] # boot props for that strata
				for(k in 1:nrow(boots)){ #for each group instead of each boot iter, b/c fewer groups than boot iter expected
					boots[k,] <- boots[k,] * strataTotal
				}
				# now boots is total estimate for that group/strata/boot, so assign to corresponding entries in bootStrata
				bootStrata[bootStrata$Strata == WinData[j,1], countCol:ncol(bootStrata)] <- boots
			}
			#now calculate grand totals
			#pull categories to sum from the original data estimate to ensure the same order
			categories_to_sum <- hierarch_output_list[[i]][[4]]
			categories_to_sum <- categories_to_sum[,1:(ncol(categories_to_sum) - 1)]
			if(!is.matrix(categories_to_sum) && !is.data.frame(categories_to_sum)){
				categories_to_sum <- as.matrix(categories_to_sum)
				colnames(categories_to_sum) <- colnames(hierarch_output_list[[i]][[3]])[2]
				#this happens when there is only variable
				#causes problems below if categories to sum is not passed as a matrix
			}
			totals_for_run <- matrix(nrow = nrow(categories_to_sum), ncol = B)

			est_by_s_mat_comp <- as.matrix(bootStrata[,colnames(categories_to_sum)])
			for (j in 1:nrow(categories_to_sum)){
				bool_temp <- find_matching_rows(est_by_s_mat_comp, unlist(categories_to_sum[j,]))
				totals_for_run[j,] <- colSums(bootStrata[bool_temp,countCol:ncol(bootStrata)])
			}
			bootHier[[i]] <- totals_for_run # will have same order (of rows) as "Grand Total" hierarch_output_list[[i]][[4]]
		}

		#now calculate and write out CIs
		for(i in 1:length(bootHier)){

			#write bootstrap estimates
			if (writeBoot){
				write.table(bootHier[[i]], paste0(Run, "_Boot_Hier_", names(hierarch_output_list)[i], ".txt"),
						  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
			}

			#one less column, will stitch on categories at end
			sim_ci_hier <- simulConfInt(t(bootHier[[i]]), alph)
			CI_hier <- matrix(nrow = nrow(bootHier[[i]]), ncol = 7)
			colnames(CI_hier) <- c("Estimate", paste0("Lower_(", alph/2, ")"), paste0("Upper_(", 1-(alph/2), ")"), "Percen_half_width", "Lower_simul", "Upper_simul", "Percen_half_width_simul")
			col_index <- ncol(hierarch_output_list[[i]][[4]])

			for (j in 1:nrow(bootHier[[i]])){
				CI_hier[j,1] <- hierarch_output_list[[i]][[4]][j,col_index]
				CI_hier[j,2:3] <- quantile(bootHier[[i]][j,],c(alph/2, 1-(alph/2)))
				CI_hier[j,4] <- round(100 * ( (as.numeric(CI_hier[j,3]) - as.numeric(CI_hier[j,2]) )/ (as.numeric(CI_hier[j,1])*2)),2)
				CI_hier[j,5] <- sim_ci_hier[j,1]
				CI_hier[j,6] <- sim_ci_hier[j,2]
				CI_hier[j,7] <- round(100 * ( (sim_ci_hier[j,2] - sim_ci_hier[j,1] )/ (as.numeric(CI_hier[j,1])*2)),2)
			}
			CI_hier <- cbind(hierarch_output_list[[i]][[4]][,1:(col_index - 1)], CI_hier)
			if(col_index - 1 == 1){
				colnames(CI_hier)[1] <- colnames(hierarch_output_list[[i]][[4]])[1]
			}
			#output CIs
			write.table(CI_hier, paste0(Run, "_CI_Hier_", names(hierarch_output_list)[i], ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
		}
	}
	return("ScobiNoPBT is complete")
}
