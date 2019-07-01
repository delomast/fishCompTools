#marginalize_SD
#marginalize Scobi-Deux
#This function produces estimates and CIs after marginalizing out a given variable
#The idea is that it can marginalize the PBTGroup variable to estimate PBT-derived categories
#		such as PBT-age
#########	inputs
# estimatesFile is the Scobi-Deux output "Estim_Totals" or "Estim_Grand_Totals" file that you want to marginalize over
# bootHier is the "Boot_Hier" file for the corresponding variable in your "Estim_Grand_Totals" file
# 	Note that use of this option requires that you feed the program the "Estim_Grand_Totals" file, not the "Estim_Totals" file
# marginalize is the name of the variable you want to marginalize over
#Example syntax to call the function
#
# SCOBI_deux(adultData = "2018RearHNCscobi_TD_no_wk1hnc.csv", windowData = "HNC_clip_window.csv",
# 		 Run = "noclipApr10", RTYPE = "noclip_H", Hierarch_variables = c("GenParentHatchery", "swAge"),
#                   SizeCut = NULL, alph = 0.1, B = 0, writeBoot = F, pbtRates = "PBT tag rates.csv",
# 		 adClipVariable = "AdClip", physTagsVariable = "PhysTag", pbtGroupVariable = "GenParentHatchery",
#			 screenOutput = "scobi_onscreen.txt")

#' @export

marginalize_SD <- function(estimatesFile = NA, bootHier = NA, marginalize = "GenParentHatchery", alph = 0.1)
{
	genCI <- T #boolean of whether to attempt to generate CIs or not
	# input error checking
	if (is.na(estimatesFile)){
		stop("No estimatesFile given")
	}
	if (is.na(bootHier)){
		cat("\nWarning: No bootHier file given, no CI estimates will be generated\n")
		genCI <- F
	}

	# load data
	estimates <- read.table(estimatesFile, sep = "\t", header = T, stringsAsFactors = F)
	if(genCI){
		bootData <- read.table(bootHier, sep = "\t", header = T, stringsAsFactors = F)
	}

	# input error checking
	if (genCI && "Strata" %in% colnames(estimates)){
		cat("\nWarning: Strata is a column name in your estimate file. Assuming this is an \"Estim_Totals\" file. No CI estimates will be generated\n")
		genCI <- F
	}

	if (genCI && nrow(bootData) != nrow(estimates)){
		cat("\nWarning: Your bootData and estimatesFile files have different numbers of rows. Assuming this is an \"Estim_Totals\" file. No CI estimates will be generated\n")
		genCI <- F
	}

	# marginal estimates
	#get all categories except marginalized
	categories_to_sum <- estimates[,!(colnames(estimates) %in% c("Total_number_in_run", "Estimated_total_for_run", marginalize))]
	categories_to_sum <- unique(categories_to_sum)
	if(!is.matrix(categories_to_sum) && !is.data.frame(categories_to_sum)){
		categories_to_sum <- as.matrix(categories_to_sum)
		colnames(categories_to_sum) <- colnames(estimates)[!(colnames(estimates) %in% c("Total_number_in_run", "Estimated_total_for_run", marginalize))]
		#this happens when there is only one variable
		#causes problems below if categories to sum is not passed as a matrix
	}
	#standardize column name for total estimate
	colnames(estimates)[colnames(estimates) %in% c("Total_number_in_run", "Estimated_total_for_run")] <- "Total_number_in_run"
	totals_for_run <- rep(-9, nrow(categories_to_sum))
	exact_match <- ncol(categories_to_sum)
	bool_list <- matrix(nrow = nrow(estimates), ncol = exact_match)
	bool_func <- function(x){sum(x) == exact_match}
	estimates_cats <- as.matrix(estimates[,colnames(categories_to_sum)])
	for (j in 1:nrow(categories_to_sum)){
		for(k in 1:exact_match){
			bool_list[,k] <- (estimates_cats[,k] == categories_to_sum[j,k])
		}
		bool_temp <- apply(bool_list, 1, bool_func)
		totals_for_run[j] <- sum(estimates[bool_temp,"Total_number_in_run"])
	}
	totals_for_run <- cbind(categories_to_sum, totals_for_run)

	# output marginalized totals
	write.table(totals_for_run, paste0("Marginalized_", marginalize, "_", estimatesFile), col.names = T, row.names = F, quote = F, sep = "\t")

	# CIs
	if(genCI){
		#make empty matrix to hold marginalized bootstrap values
		bootDataMarginal <- matrix(nrow = nrow(categories_to_sum), ncol = ncol(bootData))
		#order of categories in bootdata is the same as the order in the Grand_Total file
		# so use the orders from the Grand_Total file to sum tbe bootstrap estimates appropriately
		bool_list <- matrix(nrow = nrow(estimates), ncol = exact_match)
		for (j in 1:nrow(categories_to_sum)){
			for(k in 1:exact_match){
				bool_list[,k] <- (estimates_cats[,k] == categories_to_sum[j,k])
			}
			bool_temp <- apply(bool_list, 1, bool_func)
			#define temporrary data in case it only has one row and converts to a vector
			data_temp <- bootData[bool_temp,]
			if(length(dim(data_temp)) < 2){
				data_temp <- t(as.matrix(data_temp))
			}
			bootDataMarginal[j,] <- colSums(data_temp)
		}

		# calculate CIs
		#one less column, will stitch on categories at end
		sim_ci_hier <- simulConfInt(t(bootDataMarginal), alph)
		CI_hier <- matrix(nrow = nrow(bootDataMarginal), ncol = 7)
		colnames(CI_hier) <- c("Estimate", paste0("Lower_(", alph/2, ")"), paste0("Upper_(", 1-(alph/2), ")"), "Percen_half_width", "Lower_simul", "Upper_simul", "Percen_half_width_simul")
		col_index <- ncol(totals_for_run)
		#add estimates
		CI_hier[,1] <- totals_for_run[,col_index]
		for (j in 1:nrow(bootDataMarginal)){
			CI_hier[j,2:3] <- quantile(bootDataMarginal[j,],c(alph/2, 1-(alph/2)))
			CI_hier[j,4] <- round(100 * ( (as.numeric(CI_hier[j,3]) - as.numeric(CI_hier[j,2]) )/ (as.numeric(CI_hier[j,1])*2)),2)
			CI_hier[j,5] <- sim_ci_hier[j,1]
			CI_hier[j,6] <- sim_ci_hier[j,2]
			CI_hier[j,7] <- round(100 * ( (sim_ci_hier[j,2] - sim_ci_hier[j,1] )/ (as.numeric(CI_hier[j,1])*2)),2)
		}
		CI_hier <- cbind(totals_for_run[,1:(col_index - 1)], CI_hier)
		if(col_index - 1 == 1){
			colnames(CI_hier)[1] <- colnames(totals_for_run)[1]
		}
		# output CIs
		write.table(CI_hier, paste0("Marginalized_", marginalize, "_CI_", gsub(".+_Hier_", "", estimatesFile)), col.names = T, row.names = F, quote = F, sep = "\t")
	}

	return("Marginalizing complete")
}
