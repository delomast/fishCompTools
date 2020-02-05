#' Produce estimates and CIs after marginalizing out a given variable from the output of \code{SCOBI_deux}
#'
#' This function requires first running \code{SCOBI_deux} with \code{writeBoot = TRUE}.
#' An example use is for marginalizing the PBTGroup variable to estimate PBT-derived categories such as PBT-age.
#'
#' description2
#'
#' @param estimatesFile is the Scobi-Deux output "Estim_Totals" or "Estim_Grand_Totals" file that you
#'  want to marginalize over.
#' @param bootHier is the "Boot_Hier" file for the corresponding variable in your "Estim_Grand_Totals" file.
#'  Note that use of this option requires that you feed the program the "Estim_Grand_Totals" file, not the
#'  "Estim_Totals" file
#' @param marginalize is the name of the variable(s) you want to marginalize over. This shoudl be a character
#'  vector and can have one or multiple values.
#' @param alph This is the alpha value to use for confidence intervals. So, an alph of .05 will give a
#'  95\% CI, with quantiles at .025 and .975.
#' @export

marginalize_SD <- function(estimatesFile = NA, bootHier = NA, marginalize = "GenParentHatchery", alph = 0.1)
{
	genCI <- TRUE #boolean of whether to attempt to generate CIs or not
	# input error checking
	if (is.na(estimatesFile)){
		stop("No estimatesFile given")
	}
	if (is.na(bootHier)){
		warning("Warning: No bootHier file given, no CI estimates will be generated")
		genCI <- FALSE
	}

	# load data
	estimates <- read.table(estimatesFile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
	if(genCI){
		bootData <- read.table(bootHier, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
	}

	# input error checking
	if (genCI && "Strata" %in% colnames(estimates)){
		warning("Warning: Strata is a column name in your estimate file. Assuming this is an \"Estim_Totals\" file. No CI estimates will be generated.")
		genCI <- FALSE
	}

	if (genCI && nrow(bootData) != nrow(estimates)){
		warning("Warning: Your bootData and estimatesFile files have different numbers of rows. Assuming this is an \"Estim_Totals\" file. No CI estimates will be generated.")
		genCI <- FALSE
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
	estimates_cats <- estimates[,colnames(categories_to_sum)]
	if(length(dim(estimates_cats)) < 2){ #if it's not a matrix or dataframe
		if(ncol(categories_to_sum) == 1){#if one category, then it should be a one column matrix
			estimates_cats <- as.matrix(estimates_cats)
		} else { #if multiple categories, then there is only one group, so should be a one row matrix
			estimates_cats <- t(as.matrix(estimates_cats))
		}
	}
	for (j in 1:nrow(categories_to_sum)){
		bool_temp <- find_matching_rows(estimates_cats, unlist(categories_to_sum[j,])) # need to unlist b/c dataframe
		totals_for_run[j] <- sum(estimates[bool_temp,"Total_number_in_run"])
	}

	totals_for_run <- cbind(categories_to_sum, totals_for_run)

	# output marginalized totals
	write.table(totals_for_run, paste0("Marginalized_", paste(marginalize, collapse = "_"), "_", estimatesFile), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

	# CIs
	if(genCI){
		#make empty matrix to hold marginalized bootstrap values
		bootDataMarginal <- matrix(nrow = nrow(categories_to_sum), ncol = ncol(bootData))
		#order of categories in bootdata is the same as the order in the Grand_Total file
		# so use the orders from the Grand_Total file to sum tbe bootstrap estimates appropriately
		for (j in 1:nrow(categories_to_sum)){
			bool_temp <- find_matching_rows(estimates_cats,unlist(categories_to_sum[j,])) # need to unlist b/c dataframe
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
		write.table(CI_hier, paste0("Marginalized_", paste(marginalize, collapse = "_"), "_CI_", gsub("Hier_", "", estimatesFile)), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
	}

	return("Marginalizing complete")
}
