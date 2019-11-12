#' This function prepares the data to run multiple strata for estiamting the proportion clipped
#' @param trapData dataframe with data for fish sampled from one strata - trap data for dam escapement
#' @param adFinCol column name of column containing adipose fin status - TRUE being intact FALSE being clipped, NA missing
#' @param strataCol column name of column indicating the strata the observation belongs to
#' @param verbose TRUE to print some messages, FALSE to not
#' 
#' @export

prepStrataPropClip <- function(trapData, strataCol, adFinCol, verbose = TRUE){
	
	#turn adFinCol into boolean if necessary
	if(!is.logical(trapData[,adFinCol])){
		nonValid <- sum(!is.na(trapData[,adFinCol]) & !(trapData[,adFinCol] %in% c("AD", "AI")))
		if(nonValid > 0){
			errMessage <- paste(nonValid, "observations that are not valid options for", adFinCol,
				"\nthe adFinCol must either be a logical variable, with TRUE for ad-intact,", 
				"or be a character variable with values of AD and AI for ad-clipped and ad-intact, respectively.", 
				"\n Missing data should have values of NA.\n")
			stop(errMessage)
		}
		trapData[,adFinCol] <- trapData[,adFinCol] == "AI"
		trapData[,adFinCol] <- as.logical(trapData[,adFinCol])
	}
	
	allStrata <- list() # list containing inputs for each strata
	trapData <- trapData[!is.na(trapData[,adFinCol]) & !is.na(trapData[,strataCol]),]
	if(verbose) cat("\nUsing", nrow(trapData), "observations with adFin status and an assigned strata\n")
	for(s in sort(unique(trapData[,strataCol]))){
		strataData <- trapData[trapData[,strataCol] == s,] #select one strata
		allStrata[[s]] <- list(c(sum(!strataData[,adFinCol]), sum(strataData[,adFinCol])), # clipped, unclipped
									  strataName = s)
	}
	
	return(allStrata)
}
