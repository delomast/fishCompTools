#include <Rcpp.h>
#include <random>
#include "misc_math.h"

using namespace std;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

/*This is the function to estimate the proportion clipped.
	The model has this separate from all other parameters,
 	so estimating it separately in a separate function.
 	Probably don't really need Monte Carlo to estimate this - is simple binomial -
 	but inclusion makes calculating uncertainty in complete composition
 	easy by just multiplying within all the iterations. It is also concievable
 	that later version will treat sepearte strata as related, for example as
 	a correlated random walk, or just draws from a single distribution. in these cases,
 	having it set up to estimate through MCMC may make implementation of these
 	possibilities easier.
*/
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector propClipped(int Nclip, int Nunclip, int NumResults, unsigned int seed,
                             Rcpp::NumericVector clipPrior = Rcpp::NumericVector::create(1,1)) {
	
	// initiate random number generator
	mt19937 rg (seed);
	mt19937 * rgPoint = &rg; //pointer to random number generator to pass to subfunctions
	
	//allocate result storage
	Rcpp::NumericVector r_Propclip (NumResults);
	
	//sample posterior
	double postA = clipPrior[0] + Nclip;
	double postB = clipPrior[1] + Nunclip;
	//no need to burnin and thin, this is just sampling a random Beta
	for(int i=0; i < NumResults; i++){
		r_Propclip[i] = randBeta(postA, postB, rgPoint);
	}
	
	return r_Propclip;
}
