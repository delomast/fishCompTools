# functions used by several other function in the package

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
