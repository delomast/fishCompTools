#Scobi-Deux
#re-write of SCOBI to allow more functionality and more flexible usage
#Main new addition is reduction of W category based on PBT expansion and
#incorporation of spibetr process for correcting GSI estimates - planned, not complete yet
#########	inputs
# adultData is trap data
# windowData is window count data
# Run is prefix for output files
# RTYPE is group of interst, either "clipped" for clipped hatchery, "noclip_H" for unclipped hatchery, or "wild" for wild
# Hierarch_variables are categories to estiamte adn calculate CIs for, with order meaning subcategories
# SizeCut is border between large and small categories in lengthVariable is one of the Hierarch_variables
# lengthVariable is a variable separated in to two categories based on sizeCut (previously LGDFLmm)
# alph is alpha value for CIs
# B is number of bootstrap iterations
# writeBoot is TRUE or FALSE for whether to write individual estimates for each bootstrap iteration
# pbtRates is csv file with a header (colnames don't matter), containing pbtGroup,tagRate with one group per line
# adClipVariable is variable that specifies if a fish is adclipped of not with values of "AD" for clipped, "AI" for intact, and NA for unknown
# physTagsVariable is variable that specifies if a fish has a physical tag, regardless of ad-clip stats, indicating it is of hatchery origin (ie a CWT) with values of "tag" for tagged, "notag" for not tagged, and NA for unknown
# pbtGroupVariable is variable in adultData file that is the pbt groups assigned, with tag rates specified in pbtRates, unassigned
		#fish are given values of "Unassigned", NA is for fish not analyzed with PBT (for example, if they failed genotyping)
# screenOutput is the file name you want to save the onscreen output to. If NULL, output is printed on console
# spibetr is the correction of wild Hierarch_variables estimates for untagged HNC composition
#' @export

SCOBI_deux <- function(adultData = NULL, windowData = NULL, Run = "output", RTYPE = "wild", Hierarch_variables = NULL,
                  SizeCut = NULL, alph = 0.1, B = 5000, writeBoot = FALSE, pbtRates = NULL, adClipVariable = "AdClip",
				physTagsVariable = "PhysTag", pbtGroupVariable = "GenParentHatchery", lengthVariable = "LGDFLmm", dataGroupVariable = "WeekNumber",
				screenOutput = NULL, spibetr = TRUE){
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
	if (is.character(screenOutput)){
		sink_out <- TRUE
	} else {
		sink_out <- FALSE
	}

	if(sink_out){
		sink(paste0(Run, "_", screenOutput))
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


	#check that all PBT assigned groups are in the tag rate file
	if("Unassigned" %in% pbtRate[,1]){
		if(pbtRate[(pbtRate[,1] == "Unassigned"),2] != 1){
			stop("\nError: Unassigned was found in your tag rate file with a tag rate not equal to 1. Exiting\n")
		}
	} else {
		pbtRate <- rbind(pbtRate, c("Unassigned", 1))	#Add Unassigned with tag rate of 1 to make calculations below simple to code
	}

	pbt_groups <- FishData[,pbtGroupVariable]
	pbt_groups <- unique(pbt_groups[!is.na(pbt_groups)])
	if (sum(!(pbt_groups %in% pbtRate[,1])) != 0){
		cat("\nError: the below PBT groups were not found in the PBT tag rate input.\n")
		cat(pbt_groups[!(pbt_groups %in% pbtRate[,1])], sep = "\n")
		stop("Exiting.")
	}

	###estimate proportions hatchery, hatchery-no-clip, and wild
	### with PBT expansion

	if("pbt_expansion" %in% colnames(FishData)){
		stop("pbt_expansion is a reserved column name. Please rename this column in your input adultData dataset. Exiting.")
	}

	#discard observations with:
	##no entry for ad-clip
	num_obs <- nrow(FishData)
	FishData <- FishData[!is.na(FishData[,adClipVariable]),]
	cat("\nRemoved", num_obs - nrow(FishData), "observations for missing ad-clip information.\n")

	#add pbt expansion column - there will be NA in non pbt attempted fish, but that is ok
	FishData$pbt_expansion <- 1/pbtRate[match(FishData[,pbtGroupVariable], pbtRate[,1]),2]

	# Define large and small if one of the factors is fork length
	if(lengthVariable %in% Hierarch_variables) {
		FishData[!is.na(FishData[,lengthVariable]) & FishData[,lengthVariable] < SizeCut, lengthVariable] <- "Sm"
		FishData[!is.na(FishData[,lengthVariable]) & FishData[,lengthVariable] != "Sm", lengthVariable] <- "Lg"
	}

	# point estimate of proportions hatchery, hatchery-no-clip, and wild
	#calculate numbers of fish in each category, for each strata
	#expand hatchery fish and remove fish from wild as appropriate
	# calculate expansion for fish that are only tagged with PBT (no clip, CWT)
	# subtract this from the wild count

	cat("\nBegan analysis at:\n")
	print(Sys.time())
	cat("\n")

	###Function to decompose into H, HNC and W
	## organzed estimation procedures into functions to make non-parametric bootstrapping simpler

	decompose_h_hnc_w <- function(func_FishData, func_WinData){

		func_h_hnc_w_estimates <- matrix(nrow = nrow(func_WinData), ncol = 4)	## this will have proportion of each type trapped each strata
		colnames(func_h_hnc_w_estimates) <- c("strata", "clipped_hatchery", "unclipped_hatchery", "wild")
		func_h_hnc_w_estimates[,1] <- func_WinData[,1]
		func_unexp_func_h_hnc_w_estimates <- func_h_hnc_w_estimates

		#for each strata
		for(i in 1:nrow(func_h_hnc_w_estimates)){
			prop_all_data <- func_FishData[func_FishData$Strata == func_h_hnc_w_estimates[i,1],]
			#calculate number clipped
			clip_count <- sum(prop_all_data[,adClipVariable] == "AD")
			#save "unexpanded" clipped count
			func_unexp_func_h_hnc_w_estimates[i,2] <- clip_count

			#unclipped only data
			ai_data <- prop_all_data[prop_all_data[,adClipVariable] == "AI",]
			#calculate number AI with phystical tag
			phystag_count <- sum(!is.na(ai_data[,physTagsVariable]) & ai_data[,physTagsVariable] == "tag")
			#calculate number AI without phys tag or unknown phystag status - this is the number that needs to be called HNC and W based on PBT and tag rates
			noORna_phystag_count <- sum(is.na(ai_data[,physTagsVariable]) | ai_data[,physTagsVariable] == "notag")
			# select pbt-ony fish that were successfully genotyped
			ai_data <- ai_data[(is.na(ai_data[,physTagsVariable]) | ai_data[,physTagsVariable] == "notag") & !is.na(ai_data[,pbtGroupVariable]),]
			# count pbt_only_HNC and wild with PBT expansions
			pbtOnly_HNC <- sum(ai_data[ai_data[,pbtGroupVariable] != "Unassigned", "pbt_expansion"]) #this is expanded pbtonly hnc count
			expected_untagged <- pbtOnly_HNC - sum(ai_data[,pbtGroupVariable] != "Unassigned") #this is the number expanded (expected number of untagged fish)
			wild_count <- sum(ai_data[,pbtGroupVariable] == "Unassigned") # this is the wild count BEFORE subtracting pbt expanded fish
			# save unexpanded wild count
			func_unexp_func_h_hnc_w_estimates[i,4] <- wild_count
			# save unexpanded HNC count
			func_unexp_func_h_hnc_w_estimates[i,3] <- sum(ai_data[,pbtGroupVariable] != "Unassigned") + phystag_count

			if(expected_untagged > wild_count){
				cat("\nWarning: The expansion of HNC fish that appear wild in strata", func_h_hnc_w_estimates[i,1], "is estimated to be ", round(expected_untagged - wild_count, 2), "more than",
				    "the number of wild appearing fish. Capping the expansion to make the count of wild fish in this strata 0.\n")
				expected_untagged <- wild_count
			}
			wild_count <- wild_count - expected_untagged # this is the wild count after subtracting pbt expanded fish
			#turn into proportions
			total_pbtHNCwild <- pbtOnly_HNC + wild_count
			if (total_pbtHNCwild == 0){
				### if all fish are clipped - used mainly when running clipped fish separately from AI fish
				pbtOnly_HNC <- 0
				wild_count <- 0
			} else {
				pbtOnly_HNC <- pbtOnly_HNC / total_pbtHNCwild
				wild_count <- wild_count / total_pbtHNCwild
			}
			#multiply by noORna_phystag_count to get totals
			pbtOnly_HNC <- pbtOnly_HNC * noORna_phystag_count
			wild_count <- wild_count * noORna_phystag_count
			# add hnc groups together, one from physical tags, one from expanded pbt
			hnc_count <- phystag_count + pbtOnly_HNC

			total <- sum(wild_count, hnc_count, clip_count)

			func_h_hnc_w_estimates[i,2] <- clip_count/total
			func_h_hnc_w_estimates[i,3] <- hnc_count/total
			func_h_hnc_w_estimates[i,4] <- wild_count/total
			rm(ai_data)

		}
		#use proportions combined with census (window) data to estimate counts of each category
		func_h_hnc_w_count_est <- func_h_hnc_w_estimates[,2:4]*func_WinData[,2]
		if(!is.matrix(func_h_hnc_w_count_est)){
			func_h_hnc_w_count_est <- t(as.matrix(func_h_hnc_w_count_est))
		}
		func_category_totals <- apply(func_h_hnc_w_count_est,2,sum)
		func_h_hnc_w_count_est <- cbind(func_h_hnc_w_estimates[,1], func_h_hnc_w_count_est)
		colnames(func_h_hnc_w_count_est)[1] <- "strata"

		return(list(func_unexp_func_h_hnc_w_estimates, func_h_hnc_w_estimates, func_h_hnc_w_count_est, func_category_totals))
	}

	Estimate_h_hnc_w <- decompose_h_hnc_w(FishData, WinData)

	#split output into individual variables for easy referencing later on
	unexp_h_hnc_w_estimates <- Estimate_h_hnc_w[[1]]
	h_hnc_w_estimates <- Estimate_h_hnc_w[[2]]
	h_hnc_w_count_est <- Estimate_h_hnc_w[[3]]
	category_totals <- Estimate_h_hnc_w[[4]]

	rm(Estimate_h_hnc_w) #save memory

	cat("Estimated totals for\n\tClipped:\t", category_totals[1], "\n\tHatchery No-clip:\t", category_totals[2], "\n\tWild:\t", category_totals[3], "\n")

	#write output to files
	sink(paste0(Run,"_Rearing.txt"))
	cat("Counts of fish of each type in the trap data\n")
	write.table(unexp_h_hnc_w_estimates, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
	cat("\n\nPBT expanded proportions of each type\n")
	write.table(h_hnc_w_estimates, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
	cat("\n\nPBT expanded counts of each type in the run\n")
	write.table(h_hnc_w_count_est, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
	cat("Totals", category_totals[1], category_totals[2], category_totals[3], sep = "\t")
	sink()

	#run analyses of each factor, for rearing type of interest

	#recursive function that will build all combinations of a given list of variables
	build_all_combos <- function(cat_list){
		if(length(cat_list) > 1){
			table <- build_all_combos(cat_list[2:length(cat_list)])	#recrusive to allow flexible number of variables input
		} else{
			return(cat_list[[1]])	#if only one, return that one
		}
		if(!is.matrix(table)){ #if only one was returned, it is a vector, have to make a matrix
			out <- matrix(nrow = 0, ncol = 2)
		} else {
			out <- matrix(nrow = 0, ncol = (ncol(table) + 1))
		}
		for(i in cat_list[[1]]){
			out <- rbind(out, cbind(i, table))
		}
		return(out)
	}

	## this function takes a matrix/dataframe (data) of categories corresponding to each observation
	### with each row an observation and each column a different variable with multiple categories
	### the fucntion returns a boolean vector identifying the rows that match the vector of desired categories (want)
	### so it finds all rows in "data" with all cells equal to the corresponding cell in "want"
	find_matching_rows <- function(data, want){
		bool_list <- matrix(nrow = nrow(data), ncol = ncol(data))
		for(k in 1:ncol(bool_list)){
			bool_list[,k] <- (data[,k] == want[k])
		}
		bool_temp <- rowSums(bool_list)
		return(bool_temp == ncol(bool_list))
	}


	#select rearing type of interest
	hierarch_data <- FishData
	#select rearing type of interest(clipped, noclip_h, wild)
	if (RTYPE == "clipped"){
		spibetr_data <- NULL
		hierarch_data <- hierarch_data[hierarch_data[,adClipVariable] == "AD",]
	} else if (RTYPE == "noclip_H"){
		spibetr_data <- NULL
		hierarch_data <- hierarch_data[hierarch_data[,adClipVariable] == "AI" & ((!is.na(hierarch_data[,physTagsVariable]) & hierarch_data[,physTagsVariable] == "tag") | (!is.na(hierarch_data[,pbtGroupVariable]) & hierarch_data[,pbtGroupVariable] != "Unassigned")),]
	} else if (RTYPE == "wild"){
		#spibetr data is data for fish that are only known to be hatchery origin through PBT
		if (spibetr){
				spibetr_data <- hierarch_data[hierarch_data[,adClipVariable] =="AI" & (is.na(hierarch_data[,physTagsVariable]) | hierarch_data[,physTagsVariable] == "notag") & !is.na(hierarch_data[,pbtGroupVariable]) & hierarch_data[,pbtGroupVariable] != "Unassigned",]
		} else {
			spibetr_data <- NULL
		}
		hierarch_data <- hierarch_data[hierarch_data[,adClipVariable] == "AI" & (is.na(hierarch_data[,physTagsVariable]) | hierarch_data[,physTagsVariable] == "notag") & hierarch_data[,pbtGroupVariable] == "Unassigned",]

	} else {
		stop("Unrecognized RTYPE. Must be one of: \"clipped\", \"noclip_h\", \"wild\"\n")
	}

	#check that there is data
	if(nrow(hierarch_data) < 1){
		stop("No fish present in the dataset of the RTYPE that you chose for the hierarchical analysis. Exiting.")
	}

	# For troubleshooting
	# func_hierarch_data=hierarch_data
	# func_Hierarch_variables=Hierarch_variables
	# func_pbtGroupVariable=pbtGroupVariable
	# func_RTYPE=RTYPE
	# func_h_hnc_w_count_est=h_hnc_w_count_est
	#
	#


	decompose_hierarchical <- function(func_hierarch_data, func_Hierarch_variables, func_pbtGroupVariable, func_RTYPE, func_h_hnc_w_count_est, func_spibetr_data, func_spibetr){

		hierarch_output_list <- list()	#storage of output, flexible for different numbers of factors
		column_cats <- c()	#this is the list of factors analyzed so far plus the current
						#it is therefore the list of columns that the output for the current factor will have
		#for each factor
		for (factor in func_Hierarch_variables){
			#remove NAs
			removed <- sum(is.na(func_hierarch_data[,factor]))
			func_hierarch_data <- func_hierarch_data[!is.na(func_hierarch_data[,factor]),]
			cat("\nRemoved ", removed, " observations for missing data for ", factor, ".\n", sep="")

			# remove NAs from spibetr data as appropriate
			if(func_spibetr && func_RTYPE == "wild"){
					func_spibetr_data <- func_spibetr_data[!is.na(func_spibetr_data[,factor]),]
				}

			column_cats <- c(column_cats, factor)	#list columns for this hierarchical output variable
			if(nrow(func_hierarch_data) < 1){
				stop("There are no fish in any strata with complete data for factors", column_cats, "   Exiting.")
			}

			if(length(column_cats) < 2){
				unique_categories <- list(unique(func_hierarch_data[,column_cats]))
				names(unique_categories) <- c(column_cats)
			} else {
				unique_categories <- apply(func_hierarch_data[,column_cats],2,unique)	# get lists of all unique categories for each variable in this analysis, combinations will be the entries for each strata in the output
			}

			hierarch_estimates <- as.data.frame(matrix(nrow = 0, ncol = (length(column_cats) + 2)))	## this will have columns for: strata, categories, proportion of each type trapped
			hierarch_estimates_count <- as.data.frame(matrix(nrow = 0, ncol = (length(column_cats) + 2)))	## this will have columns for: strata, categories, count of each type trapped

			replace_NA <- FALSE #this will record whether or not the function needs to check for and replace NA entries due to lack of data for one or more strata
			#for each strata
			for(s in func_h_hnc_w_count_est[,1]){
				strata_data <- func_hierarch_data[func_hierarch_data$Strata == s,]	#select strata
				if(func_spibetr && func_RTYPE == "wild"){
					strata_spibetr <- func_spibetr_data[func_spibetr_data$Strata == s,]
					if(length(dim(strata_spibetr)) < 2){
						strata_spibetr <- t(as.matrix(strata_spibetr))
					}
				}
				#make sure there is data
				if(nrow(strata_data) == 0){
					cat("\nWarning. No fish with complete data for factors", column_cats, "in strata", s, "\n Using the mean proportions of the entire run for this strata.\n")
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
				if(func_RTYPE == "wild"){
					temp_output2 <- temp_output #counts of samples trapped
					matching_data <- strata_data[ , c("Strata", names(unique_categories))]	#get data with just columns that you need, in the order you need them
					if(func_spibetr && nrow(strata_spibetr) > 0){
						matching_spibetr <- strata_spibetr[,c("Strata", names(unique_categories))]
						if(length(dim(matching_spibetr)) < 2){
							matching_spibetr <- t(as.matrix(matching_spibetr))
						}
						strata_spibetr_bool <- TRUE
					}
					else {
						strata_spibetr_bool <- FALSE
					}
					for(i in 1:nrow(temp_output)){
						#calculate number of observations that are an exact match for the category
						temp_output[i,col_num] <- sum(find_matching_rows(matching_data, temp_output[i,1:ncol(matching_data)]))

						temp_output2[i,col_num] <- temp_output[i,col_num]
						if(strata_spibetr_bool){
							bool_temp <- find_matching_rows(matching_spibetr, temp_output[i,1:ncol(matching_spibetr)])
							subtract <- sum(as.numeric(strata_spibetr[bool_temp, "pbt_expansion"])) - sum(bool_temp) #get number of expanded fish (not including fish actually observed)
							temp_output[i,col_num] <- as.numeric(temp_output[i,col_num]) - subtract
							#prevent negative estimates
							if (temp_output[i,col_num] < 0){
								temp_output[i,col_num] <- 0
							}
						}

					}
					colnames(temp_output2) <- c("Strata", column_cats, "Count_trapped")
					hierarch_estimates_count <- rbind(hierarch_estimates_count, temp_output2, stringsAsFactors = FALSE)
					colnames(temp_output) <- c("Strata", column_cats, "Proportion_of_strata")
					# if sum is not 0, divide to turn into proportions, otherwise leave all as zeroes
					if(sum(as.numeric(temp_output[,col_num])) > 0){
						temp_output[,col_num] <- as.numeric(temp_output[,col_num]) / sum(as.numeric(temp_output[,col_num]))
					}
					hierarch_estimates <- rbind(hierarch_estimates, temp_output, stringsAsFactors = FALSE)

				} else if (func_RTYPE == "clipped" && !(func_pbtGroupVariable %in% column_cats)){
					#if clipped and no pbt variable, just use counts, no need to expand b/c the fish that would be accounted for by expansion are still clipped and therefore in the dataset
					#same routine as for wild, but no spibetr
					matching_data <- strata_data[ , c("Strata", names(unique_categories))]	#get data with just columns that you need, in the order you need them
					for(i in 1:nrow(temp_output)){
						#calculate number of observations that are an exact match for the category
						temp_output[i,col_num] <- sum(find_matching_rows(matching_data, temp_output[i,1:ncol(matching_data)]))
					}
					colnames(temp_output) <- c("Strata", column_cats, "Count_trapped")
					hierarch_estimates_count <- rbind(hierarch_estimates_count, temp_output, stringsAsFactors = FALSE)
					colnames(temp_output) <- c("Strata", column_cats, "Proportion_of_strata")
					temp_output[,col_num] <- as.numeric(temp_output[,col_num]) / nrow(strata_data)
					hierarch_estimates <- rbind(hierarch_estimates, temp_output, stringsAsFactors = FALSE)

				} else if (func_RTYPE == "noclip_H" && !(func_pbtGroupVariable %in% column_cats)){
					# need to use expansion, but only for pbt-only fish
					# expansion of phys tag fish not necessary b/c they will be included in the dataset
					temp_output_pbt_only <- temp_output #this will be expansion for pbt-only fish (only the unobserved expanded fish, not counting the observed fish)
					colnames(temp_output_pbt_only) <- c("Strata", column_cats, "Expanded_Count_trapped")
					matching_data <- strata_data[ , c("Strata", names(unique_categories))]	#get data with just columns that you need, in the order you need them
					for(i in 1:nrow(temp_output)){
						#calculate number of observations that are an exact match for the category
						bool_temp <- find_matching_rows(matching_data, temp_output[i,1:ncol(matching_data)])
						temp_output_pbt_only[i,col_num] <- sum(strata_data[bool_temp & (is.na(strata_data[,physTagsVariable]) | strata_data[,physTagsVariable] == "notag") ,"pbt_expansion"]) - sum(bool_temp & (is.na(strata_data[,physTagsVariable]) | strata_data[,physTagsVariable] == "notag"))
						temp_output[i,col_num] <- sum(bool_temp)	#this is sample count
					}
					#save estimates
					colnames(temp_output) <- c("Strata", column_cats, "Count_trapped")
					hierarch_estimates_count <- rbind(hierarch_estimates_count, temp_output, stringsAsFactors = FALSE)
					colnames(temp_output) <- c("Strata", column_cats, "Proportion_of_strata")
					temp_output[,col_num] <- as.numeric(temp_output[,col_num]) + as.numeric(temp_output_pbt_only[,col_num]) #add in expanded fish (fish that would appear wild)
					temp_output[,col_num] <- as.numeric(temp_output[,col_num]) / sum(as.numeric(temp_output[,col_num])) #turn into proportions
					hierarch_estimates <- rbind(hierarch_estimates, temp_output, stringsAsFactors = FALSE)

				} else {	#if group of interest is clipped or noclip_H with pbt variable, use expanded and subtract from unassigned
					temp_output_2 <- temp_output
					colnames(temp_output_2) <- c("Strata", column_cats, "Count_trapped")
					if(func_RTYPE == "noclip_H"){
						temp_output_3 <- temp_output #this will be only phystagged fish, so that there expansions can be easily calculated and subtracted from unassigned
						temp_output_4 <- temp_output #this will be only phystagged fish, so that there expansions can be easily calculated and subtracted from unassigned
						colnames(temp_output_3) <- c("Strata", column_cats, "Proportion_of_strata")
						colnames(temp_output_4) <- c("Strata", column_cats, "Count_trapped")
					}

					matching_data <- strata_data[ , c("Strata", names(unique_categories))]	#get data with just columns that you need, in the order you need them
					#only physically tagged fish, for expansion to subtract from "Unassigned"
					if (func_RTYPE == "noclip_H"){
						matching_data_phystag <- strata_data[strata_data[,physTagsVariable] == "tag" , c("Strata", names(unique_categories), "pbt_expansion")]	#get data with just columns that you need, in the order you need them
						if(length(dim(matching_data_phystag)) < 2){
							matching_data_phystag <- t(as.matrix(matching_data_phystag))
						}
					}
					for(i in 1:nrow(temp_output)){
						#calculate number of observations that are an exact match for the category
						bool_temp <- find_matching_rows(matching_data, temp_output[i,1:ncol(matching_data)])
						temp_output[i,col_num] <- sum(as.numeric(strata_data[bool_temp,"pbt_expansion"]))	#this will be proportions
						temp_output_2[i,col_num] <- sum(bool_temp)	#this will sample count
						if (func_RTYPE == "noclip_H"){
							#only physically tagged fish
							#calculate number of observations that are an exact match for the category
							bool_temp <- find_matching_rows(matching_data_phystag[,1:(ncol(matching_data_phystag) - 1)], temp_output[i,1:(ncol(matching_data_phystag) - 1)])
							temp_output_3[i,col_num] <- sum(as.numeric(matching_data_phystag[bool_temp,"pbt_expansion"]))	#this will be proportions
							temp_output_4[i,col_num] <- sum(bool_temp)	#this will sample count
						}
					}

					#### need to adjust unassigned group, if it exists
					if("Unassigned" %in% temp_output[,func_pbtGroupVariable]){
						#need to loop through all combinations of other variables
						# count number of fish expanded, then subtract that from the
						# corresponding Unassigned category, with floor of 0
						categories_to_sum <- temp_output[,1:(ncol(temp_output) - 1)]
						categories_to_sum <- categories_to_sum[,colnames(categories_to_sum) != func_pbtGroupVariable]
						categories_to_sum <- unique(categories_to_sum)
						if(!is.matrix(categories_to_sum)){
							categories_to_sum <- as.matrix(categories_to_sum)
							colnames(categories_to_sum) <- "Strata"
							#this happens when the pbt group is the only variable
							#causes problems for algorithm below (b/c categories_to_sum needs to be passed as a matrix)
						}
						## if analyzing clipped fish, then expanded from all clipped can be subtracted from all clipped unassigned
						# if analyzing unclipped hatchery fish, then only expanded from phystaged fish can be subtracted from unassigned
						###### expanded form nonphystag fish will appear in the "wild" group and must be subtracted there
						# if analyzing wild fish, then pbt should not be a variable, and so will not get here
						if (func_RTYPE == "clipped"){
							assigned_data <- temp_output[temp_output[,func_pbtGroupVariable] != "Unassigned", ]
							assigned_data2 <- temp_output_2[temp_output_2[,func_pbtGroupVariable] != "Unassigned", ]
						} else if (func_RTYPE == "noclip_H"){
							assigned_data <- temp_output_3[temp_output_3[,func_pbtGroupVariable] != "Unassigned", ]
							assigned_data2 <- temp_output_4[temp_output_4[,func_pbtGroupVariable] != "Unassigned", ]
						} else {
							stop("The PBT group is a variable in your Hierarchical variables list, but the group of interest is not \"clipped\" or \"noclip_H\".")
						}
						assigned_data_categories <- as.matrix(assigned_data[, colnames(categories_to_sum)])
						temp_output_categories <- as.matrix(temp_output[, colnames(categories_to_sum)])
						for (i in 1:nrow(categories_to_sum)){
							bool_temp <- find_matching_rows(assigned_data_categories, categories_to_sum[i,])
							num_expanded <- sum(as.numeric(assigned_data[bool_temp, col_num])) - sum(as.numeric(assigned_data2[bool_temp, col_num])) #number of untagged fish accounted for by expansion
							#now subtract from unassigned group
							bool_temp <- find_matching_rows(temp_output_categories, categories_to_sum[i,])
							temp_output[bool_temp & temp_output[,func_pbtGroupVariable] == "Unassigned", col_num] <- as.numeric(temp_output[bool_temp & temp_output[,func_pbtGroupVariable] == "Unassigned", col_num]) - num_expanded

							#prevent negative estimates
							if (temp_output[bool_temp & temp_output[,func_pbtGroupVariable] == "Unassigned", col_num] < 0){
								temp_output[bool_temp & temp_output[,func_pbtGroupVariable] == "Unassigned", col_num] <- 0
							}
						}

					}
					hierarch_estimates_count <- rbind(hierarch_estimates_count, temp_output_2, stringsAsFactors = FALSE)
					temp_output[,col_num] <- as.numeric(temp_output[,col_num]) / sum(as.numeric(temp_output[,col_num]))
					hierarch_estimates <- rbind(hierarch_estimates, temp_output, stringsAsFactors = FALSE)
				}



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
		if (func_RTYPE == "clipped") {
			col_n <- 2
		} else if (func_RTYPE == "noclip_H") {
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
				for (j in 1:nrow(func_h_hnc_w_count_est)){
					props[props[,1] == func_h_hnc_w_count_est[j,1],props_col] <- as.numeric(props[props[,1] == func_h_hnc_w_count_est[j,1],props_col]) * as.numeric(func_h_hnc_w_count_est[j,col_n])
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
		return(hierarch_output_list)
	} # end of decompose hierarchical function


	if(length(Hierarch_variables) > 0){
		Hierarchical_decomposition <- decompose_hierarchical(hierarch_data, Hierarch_variables, pbtGroupVariable, RTYPE, h_hnc_w_count_est, spibetr_data, spibetr)

		# Write hierarchical estimate outputs
		for(i in 1:length(Hierarchical_decomposition)){
			## the propor_hier output file can be helpful for troublshooting but is not very informative for the end-user
			#write.table(Hierarchical_decomposition[[i]][[1]], paste0(Run, "_Propor_Hier_", names(Hierarchical_decomposition)[i], ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
			write.table(Hierarchical_decomposition[[i]][[2]], paste0(Run, "_Trap_Counts_Hier_", names(Hierarchical_decomposition)[i], ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
			write.table(Hierarchical_decomposition[[i]][[3]], paste0(Run, "_Estim_Totals_Hier_", names(Hierarchical_decomposition)[i], ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
		}

		#calculate total estimates over the entire run
		for(i in 1:length(Hierarchical_decomposition)){
			estimates_by_strata <- Hierarchical_decomposition[[i]][[3]]
			categories_to_sum <- estimates_by_strata[,2:(ncol(estimates_by_strata) - 1)]
			categories_to_sum <- unique(categories_to_sum)
			if(!is.matrix(categories_to_sum) && !is.data.frame(categories_to_sum)){
				categories_to_sum <- as.matrix(categories_to_sum)
				colnames(categories_to_sum) <- colnames(Hierarchical_decomposition[[i]][[3]])[2]
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
			write.table(totals_for_run, paste0(Run, "_Estim_Grand_Totals_Hier_", names(Hierarchical_decomposition)[i], ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
			Hierarchical_decomposition[[i]][[4]] <- totals_for_run #save for use in bootstrapping
		}
	}

	#### bootstrap to obtain CIs for quantities/groups of interest
	#### nonparametric b/c adjusting for tag rates can't be simulated by parametric with simple multinomial distribution

	# build data storage
	boot_h_hnc_w <- matrix(nrow = B, ncol = 3)
	colnames(boot_h_hnc_w) <- c("clipped_hatchery", "unclipped_hatchery", "wild")

	if(length(Hierarch_variables) > 0){
	boot_hier <- list()
		for (i in 1:length(Hierarch_variables)){
			#storing replicates as columns b/c can keep row order identical to the rownames of the Grand Totals output
			#and then just pull rownames from the Grand Totals output to label the CI output
			boot_hier[[i]] <- matrix(nrow = nrow(Hierarchical_decomposition[[i]][[4]]), ncol = B)
		}
	}

	#bootstrap
	if(B < 1){
		cat("\nEnd analysis at:\n")
		print(Sys.time())
		cat("\n")
		if(sink_out){
			sink()
		}
		return("B is less than 1. Not performing any bootstrapping. SCOBI Deux is complete.")
	}

	for (b in 1:B){
		cat("Beginning bootstrap iteration", b, "\n")
		#select data for each strata
		data_boot <- FishData #to get R to assign appropriate column types to variables, start with original data and then overwrite
		for(s in WinData[,1]){
			strata_bool <- data_boot$Strata == s
			num_obs <- sum(strata_bool)
			strata_data <- data_boot[data_boot$Strata == s,]
			data_boot[data_boot$Strata == s,] <- strata_data[sample(1:num_obs, num_obs, replace = TRUE),]
		}

		#run algorithms
		#H, HNC, W estimation
		Estimate_h_hnc_w_boot <- decompose_h_hnc_w(data_boot, WinData)

		#save category totals - we are only interested in CIs for the estimates of the run as a whole
		boot_h_hnc_w[b,] <- Estimate_h_hnc_w_boot[[4]]

		if(length(Hierarch_variables) > 0){
			#select rearing type of interest(clipped, noclip_h, wild)
			if (RTYPE == "clipped"){
				boot_spibetr_data <- NULL
				data_boot <- data_boot[data_boot[,adClipVariable] == "AD",]
			} else if (RTYPE == "noclip_H"){
				boot_spibetr_data <- NULL
				data_boot <- data_boot[data_boot[,adClipVariable] == "AI" & ((!is.na(data_boot[,physTagsVariable]) & data_boot[,physTagsVariable] == "tag") | (!is.na(data_boot[,pbtGroupVariable]) & data_boot[,pbtGroupVariable] != "Unassigned")),]
			} else if (RTYPE == "wild"){
				#spibetr data is data for fish that are only known to be hatchery origin through PBT
				if (spibetr){
						boot_spibetr_data <- data_boot[data_boot[,adClipVariable] =="AI" & (is.na(data_boot[,physTagsVariable]) | data_boot[,physTagsVariable] == "notag") & !is.na(data_boot[,pbtGroupVariable]) & data_boot[,pbtGroupVariable] != "Unassigned",]
				} else {
					boot_spibetr_data <- NULL
				}
				data_boot <- data_boot[data_boot[,adClipVariable] == "AI" & (is.na(data_boot[,physTagsVariable]) | data_boot[,physTagsVariable] == "notag") & data_boot[,pbtGroupVariable] == "Unassigned",]
			} else {
				stop("Unrecognized RTYPE. Must be one of: \"clipped\", \"noclip_h\", \"wild\".")
			}
			if(nrow(data_boot) < 1){
				cat("\nNo fish present in the dataset of the RTYPE that you chose for the hierarchical analysis. Writing all zeros for this iteration.\n")
				boot_hier[[i]][,b] <- 0
				next
			}


			#hierarchical estimation
			Hierarchical_decomposition_boot <- decompose_hierarchical(data_boot, Hierarch_variables, pbtGroupVariable, RTYPE, Estimate_h_hnc_w_boot[[3]], boot_spibetr_data, spibetr)

			# hierarchical estimate outputs
			for(i in 1:length(Hierarchical_decomposition_boot)){
				estimates_by_strata <- Hierarchical_decomposition_boot[[i]][[3]]
				#pull categories to sum from the original data estimate, in case any group was eliminated by sampling and to maintain order
				categories_to_sum <- Hierarchical_decomposition[[i]][[4]]
				categories_to_sum <- categories_to_sum[,1:(ncol(categories_to_sum) - 1)]
				if(!is.matrix(categories_to_sum) && !is.data.frame(categories_to_sum)){
					categories_to_sum <- as.matrix(categories_to_sum)
					colnames(categories_to_sum) <- colnames(Hierarchical_decomposition_boot[[i]][[3]])[2]
					#this happens when there is only variable
					#causes problems below if categories to sum is not passed as a matrix
				}
				totals_for_run <- rep(-9, nrow(categories_to_sum))
				est_by_s_mat_comp <- as.matrix(estimates_by_strata[,colnames(categories_to_sum)])
				for (j in 1:nrow(categories_to_sum)){
					bool_temp <- find_matching_rows(est_by_s_mat_comp, unlist(categories_to_sum[j,]))
					totals_for_run[j] <- sum(as.numeric(estimates_by_strata[bool_temp,"Total_number_in_run"]))
				}
				boot_hier[[i]][,b] <- totals_for_run
			}
		}

	}

	#if write bootstrap estimates
	if (writeBoot){
		write.table(boot_h_hnc_w, paste0(Run, "_Boot_Rearing.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
		if(length(Hierarch_variables) > 0){
			for(i in 1:length(boot_hier)){
				write.table(boot_hier[[i]], paste0(Run, "_Boot_Hier_", names(Hierarchical_decomposition)[i], ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
			}
		}
	}


	#compute CIs
	cat("\nComputing CIs\n")
	#H HNC wild
	CI_h_hnc_w <- matrix(nrow = 3, ncol = 8)
	colnames(CI_h_hnc_w) <- c("Group", "Estimate", paste0("Lower_(", alph/2, ")"), paste0("Upper_(", 1-(alph/2), ")"), "Percen_half_width", "Lower_simul", "Upper_simul", "Percen_half_width_simul")
	CI_h_hnc_w[1,1] <- "clipped_hatchery"
	CI_h_hnc_w[2,1] <- "unclipped_hatchery"
	CI_h_hnc_w[3,1] <- "wild"

	#calculate simultaneous CIs
	sim_ci <- SCSrank(boot_h_hnc_w, 1-(alph/2))$conf.int

	for(i in 1:3){
		CI_h_hnc_w[i,2] <- category_totals[i]
		CI_h_hnc_w[i,3:4] <- quantile(boot_h_hnc_w[,i],c(alph/2, 1-(alph/2)))
		CI_h_hnc_w[i,5] <- round(100 * ( (as.numeric(CI_h_hnc_w[i,4]) - as.numeric(CI_h_hnc_w[i,3]) )/ (category_totals[i]*2)),2)
		CI_h_hnc_w[i,6] <- sim_ci[i,1]
		CI_h_hnc_w[i,7] <- sim_ci[i,2]
		CI_h_hnc_w[i,8] <- round(100 * ( (sim_ci[i,2] - sim_ci[i,1] )/ (category_totals[i]*2)),2)
	}
	#output CIs
	write.table(CI_h_hnc_w, paste0(Run, "_CI_Rearing.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

	#hierarchical
	if(length(Hierarch_variables) > 0){
		boot_CI <- list()
		for(i in 1:length(boot_hier)){
			#one less column, will stitch on categories at end
			sim_ci_hier <- SCSrank(t(boot_hier[[i]]), 1-(alph/2))$conf.int
			CI_hier <- matrix(nrow = nrow(boot_hier[[i]]), ncol = 7)
			colnames(CI_hier) <- c("Estimate", paste0("Lower_(", alph/2, ")"), paste0("Upper_(", 1-(alph/2), ")"), "Percen_half_width", "Lower_simul", "Upper_simul", "Percen_half_width_simul")
			col_index <- ncol(Hierarchical_decomposition[[i]][[4]])

			for (j in 1:nrow(boot_hier[[i]])){
				CI_hier[j,1] <- Hierarchical_decomposition[[i]][[4]][j,col_index]
				CI_hier[j,2:3] <- quantile(boot_hier[[i]][j,],c(alph/2, 1-(alph/2)))
				CI_hier[j,4] <- round(100 * ( (as.numeric(CI_hier[j,3]) - as.numeric(CI_hier[j,2]) )/ (as.numeric(CI_hier[j,1])*2)),2)
				CI_hier[j,5] <- sim_ci_hier[j,1]
				CI_hier[j,6] <- sim_ci_hier[j,2]
				CI_hier[j,7] <- round(100 * ( (sim_ci_hier[j,2] - sim_ci_hier[j,1] )/ (as.numeric(CI_hier[j,1])*2)),2)
			}
			CI_hier <- cbind(Hierarchical_decomposition[[i]][[4]][,1:(col_index - 1)], CI_hier)
			if(col_index - 1 == 1){
				colnames(CI_hier)[1] <- colnames(Hierarchical_decomposition[[i]][[4]])[1]
			}
			boot_CI[[i]] <- CI_hier
		}
		#output CIs
		for(i in 1:length(boot_CI)){
			write.table(boot_CI[[i]], paste0(Run, "_CI_Hier_", names(Hierarchical_decomposition)[i], ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
		}
	}

	cat("\nEnd analysis at:\n")
	print(Sys.time())
	cat("\n")

	if(sink_out){
		sink()
	}

	cat("\nScobi Deux is complete\n")


}
