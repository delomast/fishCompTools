#' Making run-timing graphs from output of \code{SCOBI_deux}
#'
#' This function assumes that proportions or fish of each type are constant within each strata
#' it is not meant to be a particularly robust method of looking at run timing, rather it is for exploratory use
#'
#' description2
#'
#' @param windowData This is the exact same \code{windowData} input that you used when you ran \code{SCOBI_deux}.
#' @param Run This is the exact same input for \code{Run} you used when you ran \code{SCOBI_deux}. The function looks
#'  for files with this prefix in the working directory to use to make the graph.
#' @param dailyCountFile This either a dataframe or it is the path to a csv file that has a list of counts for the days
#'  (or other timestep) that is your basic
#'  unit for graphing. This csv should have two columns, the first column is the group (matching the group - leftmost column - in your
#'  window count file) for that day and the second column is the count on that day. Each consecutive day should be a
#'  new row with the first day being the first row and the last day being the last row. The group will be repeated once for
#'  each day within that group with the count for that particular day.
#' @param makePlot If TRUE, a pdf will be saved in the working directory of a run timing plot.
#' @export

run_time_graph <- function(windowData = NULL, Run = "output", dailyCountFile = NULL, makePlot = FALSE)

{

	### function to make plots
	make_plot <- function(count_data, name){
		# convert to dataframe and assign variable types
		count_data <- as.data.frame(count_data, stringsAsFactors = F)
		count_data[,2:ncol(count_data)] <- apply(count_data[,2:ncol(count_data)], 2, as.numeric)
		# get timestep values
		days <- count_data$timestep
		# find plotting parameters for each group
		plot_data <- matrix(nrow = 0, ncol = 6)
		for(j in 3:ncol(count_data)){
			## find the 5, 25, 50, 75, 95 % of the run passing days
			x <- count_data[,j]
			total <- sum(x)
			count <- x[1]
			for(i in 2:length(x)){
				count <- c(count, count[i-1] + x[i])
			}

			lin_int <- approxfun(count)
			param_found <- 0
			param_plot <- NULL
			## look for times in .1 steps
			for(i in seq(days[1],days[length(days)],.1)){
				if(param_found == 0){
					if (lin_int(i) >= total*.05){
						param_plot <- i
						param_found <- 1
					}
				}else if(param_found == 1){
					if (lin_int(i) >= total*.25){
						param_plot <- c(param_plot, i)
						param_found <- 2
					}
				}else if(param_found == 2){
					if (lin_int(i) >= total*.5){
						param_plot <- c(param_plot, i)
						param_found <- 3
					}
				}else if(param_found == 3){
					if (lin_int(i) >= total*.75){
						param_plot <- c(param_plot, i)
						param_found <- 4
					}
				}else if(param_found == 4){
					if (lin_int(i) >= total*.95){
						param_plot <- c(param_plot, i)
						break
					}
				}
			}
			plot_data <- rbind(plot_data, c(colnames(count_data)[j], param_plot))
		}

		#created as matrix for easy building by row, change to data frame for plotting
		plot_data <- as.data.frame(plot_data, stringsAsFactors = F)
		#assign correct variable types
		plot_data[,2:6] <- apply(plot_data[2:6], 2, as.numeric)
		colnames(plot_data) <- c("Group", "p1", "p2", "p3", "p4", "p5")
		#plot
		plot_made <- ggplot2::ggplot(plot_data, ggplot2::aes(Group)) +
				ggplot2::geom_boxplot(ggplot2::aes(lower=p2, upper=p4, middle=p3, ymin =p1, ymax=p5), stat = "identity") +
				ggplot2::coord_flip() + ggplot2::ylab("timestep") + ggplot2::xlab("group") + ggplot2::scale_y_continuous(breaks = seq(min(days), max(days), by = round(max(days)/15)))
		ggplot2::ggsave(filename = paste0(name, ".pdf"), plot = plot_made, device = "pdf", units = "cm", height = (nrow(plot_data)*3.8) + 4, width = 16.5, limitsize = F)
	}

	#load daily counts
	if (is.character(dailyCountFile)){
		daily_counts <- read.csv(dailyCountFile, header = TRUE, stringsAsFactors = FALSE)
	} else {
		daily_counts <- dailyCountFile
		rm(dailyCountFile) # remove in case passed as a dataframe and is large
	}


	#strata collapsing
	if (is.character(windowData)){
		Windata  <- read.csv(file = windowData, header = TRUE, stringsAsFactors = FALSE)
	} else {
		Windata <- windowData
		rm(windowData)
	}
	for(i in 1:nrow(Windata)){
		daily_counts[daily_counts[,1] == Windata[i,1],1] <- Windata[i,3]
	}

	#normalize daily counts within strata(turn them into proporitions of the strata that were observed on that day)
	for(s in unique(daily_counts[,1])){
		daily_counts[daily_counts[,1] == s,2] <- daily_counts[daily_counts[,1] == s,2] / sum(daily_counts[daily_counts[,1] == s,2])

	}

	to_make_graphs <- c()
	#get rearing file
	files <- dir()
	rearing_file <- files[grepl(paste0("^", Run,"_Rearing.txt$"), files)]
	if(length(rearing_file) < 1){
		cat("\nError: no rearing file found, will not make a run timing distribution for rearing origin fish.\n")
	}

	#get hier_variable output files
	hier_var <- files[grepl(paste0("^", Run,"_Estim_Totals_Hier_"), files) & grepl(".txt$", files)]
	cat("\nFound", length(hier_var), "file(s) for hierarchical variables to make run timining distributions for.\n")

	# make rearing if present
	#### there shouldn't be more than one rearing file found based on the regex used above, but this makes it easy
	#### to change the regex in the future to allow more than one

	if(length(rearing_file) >= 1){
		for (file in rearing_file){
			#open file
			connection  <- file(file, open = "r")
			line <- readLines(connection, n = 1, warn = FALSE)
			# find where the data you want starts
			while (line != "PBT expanded counts of each type in the run" | length(line) < 1) {
			  line <- readLines(connection, n = 1, warn = FALSE)
			}
			if(length(line) < 1){
				cat("\nError: found rearing file, but not counts for each type in the run. Exiting.\n")
				close(connection)
				return()
			}
			line <- readLines(connection, n = 1, warn = FALSE) # read header in table
			## read strata and estimates, incorporate into a matrix
			line <- readLines(connection, n = 1, warn = FALSE) # read first strata
			rearing_table <- matrix(nrow = 0, ncol = 3)
			while (!grepl("Totals", line) | length(line) < 1) {
				rearing_table <- rbind(rearing_table, as.numeric(strsplit(line, "\t")[[1]][2:4]))
				rownames(rearing_table)[nrow(rearing_table)] <- strsplit(line, "\t")[[1]][1]
				line <- readLines(connection, n = 1, warn = FALSE)
			}
			close(connection)
			colnames(rearing_table) <- c("clipped_hatchery", "unclipped_hatchery", "wild")

			#make daily counts of each type
			daily_rearing <- matrix(nrow = nrow(daily_counts), ncol = 0)
			for (i in 1:3){
				counts_temp <- c()
				#for each strata
				for(s in rownames(rearing_table)){
					counts_temp <- c(counts_temp, daily_counts[daily_counts[,1] == s,2]*rearing_table[s,i])
				}
				daily_rearing <- cbind(daily_rearing,counts_temp)
			}
			colnames(daily_rearing) <- c("clipped_hatchery", "unclipped_hatchery", "wild")
			daily_rearing_out <- cbind(daily_counts[,1], 1:nrow(daily_rearing), daily_rearing)
			colnames(daily_rearing_out)[1:2] <- c("strata", "timestep")
			#output table of daily counts for user to graph with more flexibility
			write.table(daily_rearing_out, paste0(Run, "_daily_counts_Rearing.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
			#output graph if selected
			if (makePlot){
				make_plot(daily_rearing_out, paste0(Run, "plot_Rearing"))
			}
		}	# end of for loop for multiple rearing files
	}	#end of if length rearing file

	# make hierarch if present
	if(length(hier_var) >= 1){
		for (file in hier_var){
			# load data
			temp_data <- read.table(file, header = T, sep = "\t", colClasses = "character")
			num_cols <- ncol(temp_data)
			#make counts numeric, but keep everything else as character
			temp_data[,num_cols] <- as.numeric(temp_data[,num_cols])
			#get list of all types
			#make matrix in case there is only one column
			types <- as.matrix(unique(temp_data[,2:(num_cols - 1)]))
			#make daily counts of each type
			daily_types <- matrix(nrow = nrow(daily_counts), ncol = 0)
			for (i in 1:nrow(types)){
				type <- types[i,]
				counts_temp <- c()
				exact_match <- num_cols - 2
				#for each strata
				for(s in unique(daily_counts[,1])){
					bool_temp <- apply(as.matrix(temp_data[,2:(num_cols - 1)]), 1, function(x){sum(x == type) == exact_match})
					counts_temp <- c(counts_temp, daily_counts[daily_counts[,1] == s,2]*temp_data[temp_data[,1] == s & bool_temp, num_cols])
				}
				daily_types <- cbind(daily_types,counts_temp)
			}
			# for (i in 1:nrow(types)){
			colnames(daily_types) <- apply(types, 1, paste, collapse = "_")
			# }

			daily_types_out <- cbind(daily_counts[,1], 1:nrow(daily_types), daily_types)
			colnames(daily_types_out)[1:2] <- c("strata", "timestep")
			#output table of daily counts for user to graph with more flexibility
			var_name <- gsub("^.+Hier_", "", file)
			write.table(daily_types_out, paste0(Run, "_daily_counts_", var_name), col.names = T, row.names = F, quote = F, sep = "\t")
			#output graph
			if (makePlot){
				make_plot(daily_types_out, paste0(Run, "_plot_", gsub(".txt", "", var_name)))
			}
		}
	}
}
