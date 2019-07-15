### simultaneous ci based on Mandel and Betensky 2008 Comput Stat Data Anal. 2008 Jan 10; 52(4): 2158–2165. doi: 10.1016/j.csda.2007.07.005
### data is table of bootstrap estimates with each rwo a new bootstrap and each column a variable
### alpha is the alpha value desired for the CI
### currently, it only calculates two-sided CIs
### output is a matrix with each row a variable, first column the lower, and second column the upper CI
### note that it can give incorrect outputs for HNC and W if bootstrap iterations are too small (less than 10000)
### and particularly if no clipped fish are present
### this is because HNc adn Ware dependent, so ranks will be the same but opposite signs
### signs are given randomly in ties with different signs, so need a large amount of bootstrap iterations
### to converge on appropriate simultaneous interval (which is the same as the one-at-a-time interval under
### particular conditions of no phystags and no clipped fish)

simulConfInt <- function(data, alpha){

	# order based on each variable
	indiv_ranks <- apply(data, 2, rank, ties.method= "random")

	# calculate relative ranks
	B <- nrow(indiv_ranks) # number of bootstrap iterations, used multiple times below
	rel_ranks <- apply(indiv_ranks, 2, function(x){
		return(x - ((B + 1) / 2))
	})
	rm(indiv_ranks)
	# separate into signs and absolute ranks
	signs <- apply(rel_ranks, 2, function(x){
		s <- rep(1, B)
		s[x < 0] <- -1
		return(s)
	})
	rel_ranks <- apply(rel_ranks, 2, abs)
	# calculate max rank and sign across all variables
	max_rank <- apply(rel_ranks, 1, max)
	max_sign <- rep(0, B)
	for(i in 1:B){
		row_signs <- signs[i, rel_ranks[i,] == max_rank[i]] # get all signs of relranks equal to the max
		if (length(row_signs) == 1) { #only one, so use it
			max_sign[i] <- row_signs
		} else { #if matches, randomly sample the signs
			max_sign[i] <- sample(row_signs, 1)
		}
	}
	rm(rel_ranks)
	final_ranks <- ((B + 1) / 2) + (max_rank * max_sign)

	half_alph <- alpha / 2
	upBound <- round(B * (1 - half_alph), 0)
	lowBound <- round(B * half_alph, 0)
	if (lowBound == 0){lowBound <- 1} #prevent accessing the 0th rank, may happen when bootstrap iterations are too small
	sorted <- sort(final_ranks)
	lowBound <- sorted[lowBound]
	upBound <- sorted[upBound]

	cis <- apply(data, 2, function(x){
		sorted <- sort(x)
		return(c(sorted[lowBound], sorted[upBound]))
	})

	return(t(cis))
}
