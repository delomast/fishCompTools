#include <Rcpp.h>
#include <random>
#include <vector>
#include "misc_math.h"

using namespace std;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

//This is the main function to run the MCMC composition model
//
//
//groups is numeric vector of all the groups ie 1,2,3,...in the same order as all other group based variables
//
//nPBT is the number of PBT groups
//GSI_values are the possible values in the GSI category
//gsiUT are the values of the untagged
//
//
//values is a list of NumericVectors giving the different categories for each variable
//pi_VInitial is a list of NumericMatrices giving the composition of each variable(matrix) by group(rows) and values (column)
//pi_Vohnc is a list of NumericMatrices giving the counts of ohnc each variable(matrix) by group(rows) and values (column)
//  It includes the wild, with counts of zero
//pi_Vprior is a list of NumericMatrices giving the alphas for a Dirichlet prior of each variable(matrix) by group(rows) and values (column)
//v_ut is a matrix with a row for each untagged individual, a column for each variable, and the data being the values for each variable.
//  -9 is missing
//t is tag rates in same order as groups
//' @export
// [[Rcpp::export]]
Rcpp::List MCpbt(int iter, int burnIn, int thin, unsigned int seed, //overall parameters for the chain
                          Rcpp::NumericVector piTotPrior, Rcpp::NumericVector ohnc, Rcpp::NumericVector piTotInitial, //piTotal parameters
                          Rcpp::NumericVector oUTInitial, Rcpp::NumericVector groups,
                          int nPBT, Rcpp::NumericVector GSI_values, Rcpp::NumericVector gsiUT, //pi_gsi parameters
                          Rcpp::NumericMatrix pi_gsiInitial, Rcpp::NumericMatrix prior_pi_gsi,
                          Rcpp::NumericMatrix ohnc_gsi,
                          Rcpp::List values, Rcpp::List pi_VInitial, Rcpp::List pi_Vohnc, Rcpp::List pi_Vprior, //pi_V parameters
                          Rcpp::NumericMatrix v_ut, 
                          Rcpp::NumericVector initZ, Rcpp::NumericVector t, //z parameters
                          Rcpp::List valuesOth, Rcpp::List pi_VInitialOth, Rcpp::List pi_VohncOth, Rcpp::List pi_VpriorOth, //pi_VOth parameters
                          Rcpp::NumericMatrix v_utOth
                          )
{
	// initiate random number generator
	mt19937 rg (seed);
	mt19937 * rgPoint = &rg; //pointer to random number generator to pass to subfunctions

	if(thin < 1) Rcpp::stop("thin must be 1 or greater.");

	int NumResults = (iter - burnIn) / thin; //number of results to store
	
	if (NumResults < 1) Rcpp::stop("iter, burnIn, and thin parameters lead to no results being recorded.");

	
	//set quantities frequently used
	int nGroups = groups.size(); //number of groups in piTot, including wild/unassigned
	int nGSI = nGroups - nPBT;
	vector <int> groupsC (nGroups); //needs to be vector of ints to pass to sampleC
	vector <double> priorPlusOhnc (nGroups); //creating prior plus observed PBT assigned
	vector <double> untagRates (nGroups);
	for(int i=0; i < nGroups; i++){
		priorPlusOhnc[i] = piTotPrior[i] + ohnc[i];
		untagRates[i] = 1 - t[i];
		groupsC[i] = groups[i];
	}
	vector <vector <double>> priorPlusOhnc_pi_gsi;
	for(int i=0; i<nPBT; i++){
		vector <double> tempVec (nGSI, 0);
		for(int j=0; j<nGSI; j++){
			tempVec[j] = prior_pi_gsi(i,j) + ohnc_gsi(i,j);
		}
		priorPlusOhnc_pi_gsi.push_back(tempVec);
	}
	
	int tempCol; //stores position of vector to take a value from
	int tempInt; //temporary counting value
	
	
	//define temporary variables used in calculations
	vector <double> piTot_tempAlphas (nGroups, 0);
	vector <double> tempZprobs (nGroups, 0);
	vector <double> pi_gsi_tempAlphas (nGSI, 0);
	
	//define main variables and set initial values
	vector <double> piTot (nGroups, 0); //proportions of each group in the population
	vector <int> oUT (nGroups,0); //number of untagged fish in each group
	for(int i=0; i < nGroups; i++){
		piTot[i] = piTotInitial[i];
		oUT[i] = oUTInitial[i];
		
	}

	vector <vector <double>> pi_gsi;
	for(int i=0; i < nGroups; i++){
		vector <double> tempVec (nGSI, 0);
		for(int j=0; j < nGSI; j++){
			tempVec[j] = pi_gsiInitial(i,j);
		}
		pi_gsi.push_back(tempVec);
	}
	
	vector <int> z; //declare vector of z's
	for(int i=0, max=initZ.size(); i < max; i++){
		z.push_back(initZ[i]); //initialize with given values
	}
	int nZ = z.size(); //number of untagged fish
	
	/*now for the variables
	 we need:
	 	a structure to be able to loop through the variables - 0:nVar
	 	then, for each variable - set up as vectors for each need, with element 0 being var1, 1 being var2, etc...
		 	a list of values - valuesC
			a list of proportions for each group (a matrix of proportions within groups) - pi_V
	 		a temporary list of counts, reset for each group when calculating - could allocate this on the fly, but that's a lot of allocating

	  		a prior for each group - combination of prior and observed (ohnc) values
	 Instead of a specific structure to loop through, just keep the order the same for all
	 */ 
	int nVar = values.size(); //number of variables
	vector <vector <int>> valuesC (nVar); //possible values of variables represented as (positive) ints, (with -9 in data indicates missing)
	vector <vector <vector <double>>> pi_V (nVar); //a list each variable, with proportions of each group for that variable
	vector <vector <double>> tempCounts (nVar); //temporary counts for each variable - used in calculations - using double b/c prior added to it

	vector <vector <vector <double>>> priorPlusOhncV (nVar); //creating prior plus observed PBT assigned
	for(int i=0; i<nVar; i++){
		//add values to valuesC
		Rcpp::NumericVector tempV = values[i]; //pull NumericVector out of the list of the possible values for a variable
		valuesC[i].assign(tempV.size(), 0);
		for(int j=0, max=tempV.size(); j<max; j++){
			valuesC[i][j] = tempV[j];
		}
		//initiate temporary counts with 0's
		tempCounts[i].assign(tempV.size(), 0);
		//add initial values to pi_V and priorPlusOhncV
		Rcpp::NumericMatrix tempVmat = pi_VInitial[i]; //this is matrix of initial values for pi_V of that variable
		Rcpp::NumericMatrix tempVohnc = pi_Vohnc[i]; //this is matrix of ohnc counts for pi_V of that variable
		Rcpp::NumericMatrix tempVprior = pi_Vprior[i];//this is matrix of alphas for Dirichlet prior for pi_V of that variable
		for(int j=0; j < nGroups; j++){ //for each row - group
			vector <double> tempVec (tempV.size()); //make empty vector to accept values for one group
			vector <double> tempVecPrior (tempV.size()); //make empty vector to accept values for one group
			for(int k=0, max=tempV.size(); k < max; k++){//for each column -variable value
				tempVec[k] = tempVmat(j,k); //take value from initial values matrix and assign to vector
				tempVecPrior[k] = tempVohnc(j,k) + tempVprior(j,k); //take value from initial values matrix and assign to vector
			}
			pi_V[i].push_back(tempVec); //add initial values for variable i and group j
			priorPlusOhncV[i].push_back(tempVecPrior); //add prior + ohnc values for variable i and group j
		}
	}
	
	/* Now we are creating the same containers for the "other" variables (not used in clustering) that we did for
	 * the main pi_V variables.
	 */
	int nVarOth = valuesOth.size(); //number of variables
	vector <vector <int>> valuesCOth (nVarOth); //possible values of variables represented as (positive) ints, (with -9 in data indicates missing)
	vector <vector <vector <double>>> pi_VOth (nVarOth); //a list each variable, with proportions of each group for that variable
	vector <vector <double>> tempCountsOth (nVarOth); //temporary counts for each variable - used in calculations - using double b/c prior added to it

	vector <vector <vector <double>>> priorPlusOhncVOth (nVarOth); //creating prior plus observed PBT assigned
	for(int i=0; i<nVarOth; i++){
		//add values to valuesC
		Rcpp::NumericVector tempV = valuesOth[i]; //pull NumericVector out of the list of the possible values for a variable
		valuesCOth[i].assign(tempV.size(), 0);
		for(int j=0, max=tempV.size(); j<max; j++){
			valuesCOth[i][j] = tempV[j];
		}
		//initiate temporary counts with 0's
		tempCountsOth[i].assign(tempV.size(), 0);
		//add initial values to pi_V and priorPlusOhncV
		Rcpp::NumericMatrix tempVmat = pi_VInitialOth[i]; //this is matrix of initial values for pi_V of that variable
		Rcpp::NumericMatrix tempVohnc = pi_VohncOth[i]; //this is matrix of ohnc counts for pi_V of that variable
		Rcpp::NumericMatrix tempVprior = pi_VpriorOth[i];//this is matrix of alphas for Dirichlet prior for pi_V of that variable
		for(int j=0; j < nGroups; j++){ //for each row - group
			vector <double> tempVec (tempV.size()); //make empty vector to accept values for one group
			vector <double> tempVecPrior (tempV.size()); //make empty vector to accept values for one group
			for(int k=0, max=tempV.size(); k < max; k++){//for each column -variable value
				tempVec[k] = tempVmat(j,k); //take value from initial values matrix and assign to vector
				tempVecPrior[k] = tempVohnc(j,k) + tempVprior(j,k); //take value from initial values matrix and assign to vector
			}
			pi_VOth[i].push_back(tempVec); //add initial values for variable i and group j
			priorPlusOhncVOth[i].push_back(tempVecPrior); //add prior + ohnc values for variable i and group j
		}
	}
	
	
	
	
	//////////////////////////////////
	//allocate result storage
	int currentEntry = 0; //keep track of what entry to record next
	
	Rcpp::NumericMatrix r_PiTot(NumResults, nGroups); //pi_Tot unclipped
	Rcpp::NumericMatrix r_Z(NumResults, nZ); //z
	
	/* For r_pi_gsi structure is a matrix with rows beign interations and columns being
	 * group1 gsi1, group1 gsi2, ..., group2 gsi1, ..., groupN gsiK
	 * so number of rows is number of recorded iterations
	 * number of columns is number of groups x number of gsi categories
	 * wild groups are recorded even though currently their pi_gsi doesn't change
	 * In order to make future changes easy if we do want to update them as well as
	 * to make post-processing of the data easier when looping over the pi_gsi estimates
	 */
	
	Rcpp::NumericMatrix r_pi_gsi (NumResults, (nGSI * nGroups)); //pi_gsi
	
	/* For r_pi_V, structure is the same as for r_pi_gsi, except there
	 * is one NumericMatrix for each variable. These NumericMatrices are
	 * stored in a vector, and then converted to a list to be output back
	 * to R
	 */
	
	vector <Rcpp::NumericMatrix> r_pi_V (nVar); //pi_V
	for (int i=0; i<nVar;i++){
		r_pi_V[i] = Rcpp::NumericMatrix(NumResults, (nGroups * valuesC[i].size()));
	}
	
	//the same strategy for pi_VOth as for pi_V
	vector <Rcpp::NumericMatrix> r_pi_VOth (nVarOth); //pi_V
	for (int i=0; i<nVarOth;i++){
		r_pi_VOth[i] = Rcpp::NumericMatrix(NumResults, (nGroups * valuesCOth[i].size()));
	}
	
	// end allocate result storage
	///////////////////////////////////
	
	
	////////////////////////////////////////////////////
	//cycle through iterations
	int thin2 = 1;
	for (int r=0, maxIter = NumResults + burnIn; r < maxIter; r++){
		if(r == burnIn) thin2 = thin; //this runs the number of burnin iterations, then runs enough to get the specified number of recordings
		//cycle through thinning reps between recording values
		for (int th=0; th < thin2; th++){
			
			// sample from pi_tot
			for(int i=0; i < nGroups; i++){ //calculating the alphas of the Dirichlet
				piTot_tempAlphas[i] = priorPlusOhnc[i] + oUT[i];
			}
			piTot = randDirich(piTot_tempAlphas, rgPoint);
			
			// sample from pi_v's
			for(int v=0; v < nVar; v++){ //cycle through variables
				for(int g=0; g < nGroups; g++){ //cycle through groups
					tempCounts[v].assign(valuesC[v].size(),0); //zero the vector
					for(int i=0, max = valuesC[v].size(); i < max; i++){ //cycle through categories to count
						//now count all in that group with that value for that variable
						for(int j=0; j < nZ; j++){ //for each untagged individual
							if(z[j] == groupsC[g] && v_ut(j,v) == valuesC[v][i]) //v_ut is kept as a NumericMatrix b/c no manipulation needed. row is indiv col is variable
								tempCounts[v][i]++; //this is utVcounts equivalent in the R script
						}
					}
					//now tempCounts has the counts for that group and variable
					//so now add to ohnc counts, prior, and then sample from Dirichlet
					for(int i=0, max = valuesC[v].size(); i < max; i++){
						tempCounts[v][i] += priorPlusOhncV[v][g][i];
					}
					pi_V[v][g] = randDirich(tempCounts[v], rgPoint);
				}
			}
			
			// NOTE- when adding in accounting for GSI variability, priorPlusOhnc_pi_gsi will
			//   have to be calculated in each iteration
			// sample from pi_gsi's - for PBT groups only
			for(int i=0; i<nPBT; i++){
				for(int j=0; j < nGSI; j++){ // calculating the alphas of the Dirichlet
					//count up z's with this GSI value
					tempInt = 0;
					for(int zi=0; zi<nZ; zi++){
						if(z[zi] == groupsC[i] && gsiUT[zi] == GSI_values[j]) tempInt++;
					}
					pi_gsi_tempAlphas[j] = priorPlusOhnc_pi_gsi[i][j] + tempInt;
				}
				pi_gsi[i] = randDirich(pi_gsi_tempAlphas, rgPoint);
				
			}
			
			
			
			// sample from z's
			for(int i=0; i < nZ; i++){ //for each z
				//first assign just the piTot and tag rates to the probs
				for(int g=0; g < nGroups; g++){
					tempZprobs[g] = piTot[g] * untagRates[g];
				}
				// take into account observed GSI assignment
				tempCol = 0;
				for(int j=0; j<nGSI; j++){ //determine which column corresponds to the observed value
					if(gsiUT[i] == GSI_values[j]){ //if observed GSI for that individual matches GSI_values[j]
						tempCol = j;
						break;
					}
				}
				//now multiply by GSI proportions
				for(int g=0; g < nGroups; g++){
					tempZprobs[g] *= pi_gsi[g][tempCol];
				}
				
				//now refine probs using pi_V and the observed variables - if they are present
				for(int v=0; v < nVar; v++){
					tempCol = 0;
					if(v_ut(i,v) == -9) continue; //if value is missing, skip
					//find position of that value - they should be in order, so if speed is a big issue, can replace this just using the value as the position
					for(int k=0, max = valuesC[v].size(); k < max; k++){
						if(valuesC[v][k] == v_ut(i,v)){
							tempCol = k; //if this never happens, will use 0, which would be a problem. Only happens when input is bad.
							break;
						}
					}
					for(int g=0; g < nGroups; g++){ //for each group
						// find probability of that value for that group and multiply
						tempZprobs[g] *= pi_V[v][g][tempCol];
					}
				}
				// now sample z for that individual
				z[i] = sampleC(groupsC, tempZprobs, rgPoint);
				
			}
			
			
			// calculate oUT for ease of later steps
			oUT.assign(nGroups, 0); //clear vector
			for(int i=0; i < nZ; i++){ //for each z
				for(int g=0; g < nGroups; g++){ //for each group
					if(z[i] == groupsC[g]){
						oUT[g]++; // add 1 to that group
						break;
					}
				}
			}
			
			//sample from variables that are not used in the "clustering" of individuals to groups
			//pi_VOth
			for(int v=0; v < nVarOth; v++){ //cycle through variables
				for(int g=0; g < nGroups; g++){ //cycle through groups
					tempCountsOth[v].assign(valuesCOth[v].size(),0); //zero the vector
					for(int i=0, max = valuesCOth[v].size(); i < max; i++){ //cycle through categories to count
						//now count all in that group with that value for that variable
						for(int j=0; j < nZ; j++){ //for each untagged individual
							if(z[j] == groupsC[g] && v_utOth(j,v) == valuesCOth[v][i]) //v_utOth is kept as a NumericMatrix b/c no manipulation needed. row is indiv col is variable
								tempCountsOth[v][i]++; //this is utVcounts equivalent in the R script
						}
					}
					//now tempCounts has the counts for that group and variable
					//so now add to ohnc counts, prior, and then sample from Dirichlet
					for(int i=0, max = valuesCOth[v].size(); i < max; i++){
						tempCountsOth[v][i] += priorPlusOhncVOth[v][g][i];
					}
					pi_VOth[v][g] = randDirich(tempCountsOth[v], rgPoint);
				}
			}
			
		} //end of thinning loop
		
		//record values
		if(r >= burnIn){
			//piTot
			for(int g=0; g < nGroups; g++){
				r_PiTot(currentEntry, g) = piTot[g];
			}
			//z
			for(int i=0;i<nZ;i++){
				r_Z(currentEntry,i) = z[i];
			}
			//pi_gsi
			for(int i=0;i<nGroups;i++){
				int temp = i * nGSI;
				for(int j=0;j<nGSI;j++){
					r_pi_gsi(currentEntry, (temp + j)) = pi_gsi[i][j];
				}
			}
			//pi_V
			for (int i=0; i<nVar;i++){
				for (int j=0; j<nGroups; j++){
					int temp = j * valuesC[i].size();
					for(int k=0, max = valuesC[i].size(); k<max; k++){
						r_pi_V[i](currentEntry, (temp + k)) = pi_V[i][j][k];
					}
				}
			}
			//pi_VOth
			for (int i=0; i<nVarOth;i++){
				for (int j=0; j<nGroups; j++){
					int temp = j * valuesCOth[i].size();
					for(int k=0, max = valuesCOth[i].size(); k<max; k++){
						r_pi_VOth[i](currentEntry, (temp + k)) = pi_VOth[i][j][k];
					}
				}
			}
			currentEntry++;
		}
		
		//check user interrupt every so often
		if (r % 200 == 0) Rcpp::checkUserInterrupt();
	}
	

	//organize output to return to R
	//turn r_pi_V into a Rcpp::List of NumericMatrices instead of a vector
	Rcpp::List r_pi_V_out (nVar);
	for(int i=0; i<nVar; i++){
		r_pi_V_out[i] = r_pi_V[i];
	}
	//same strategy for pi_VOth
	Rcpp::List r_pi_VOth_out (nVarOth);
	for(int i=0; i<nVarOth; i++){
		r_pi_VOth_out[i] = r_pi_VOth[i];
	}
	
	Rcpp::List outputList = Rcpp::List::create(Rcpp::Named("piTot") = r_PiTot,
                                            Rcpp::Named("z") = r_Z,
                                            Rcpp::Named("piGSI") = r_pi_gsi,
                                            Rcpp::Named("piV") = r_pi_V_out,
                                            Rcpp::Named("piVOth") = r_pi_VOth_out);

	return outputList;
}
