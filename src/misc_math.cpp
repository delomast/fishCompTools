#include <Rcpp.h>
#include <random>
#include <vector>
using namespace std;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// this function generates random numbers from a beta distribution
//not exported
double randBeta(double alpha, double beta, mt19937 * rNum) {

	//input checking
	if (alpha + beta == 0 || alpha < 0 || beta < 0){
		Rcpp::stop("invalid alpha or beta to randBeta");
	}
	
	// initiate random gammas
	gamma_distribution<double> rGamma1 (alpha,1); //alpha of beta distribution
	gamma_distribution<double> rGamma2 (beta,1); //beta of beta distribution

	double x1;
	double x2;
	x1 = rGamma1(*rNum);
	x2 = rGamma2(*rNum);

	return x1 / (x1 + x2); //random beta computed from two random gammas
}

// this function generates random numbers from a Dirichlet distribution
//not exported
vector<double> randDirich(vector<double> alphas, mt19937 * rNum) {

	int dim = alphas.size(); //get number of dimensions of requested Dirichlet
	
	//input checking
	double sumA = 0;
	for(int i=0; i < dim; i++){
		if(alphas[i] < 0){
			Rcpp::stop("invalid alphas to randDirichlet");
		}
		sumA += alphas[i];
	}
	if (sumA == 0){
		Rcpp::stop("invalid alphas to randDirichlet");
	}
	
	// allocate storage of random gammas
	vector<double> rGammaValues (dim);
	double sum = 0; //sum for normalizing
	//sample gammas
	for (int i = 0; i < dim; i++){
		gamma_distribution<double> rGamma1 (alphas[i],1);
		rGammaValues[i] = rGamma1(*rNum);
		sum += rGammaValues[i]; //calculate sum
	}
	//normalize
	for (int i = 0; i < dim; i++){
		rGammaValues[i] = rGammaValues[i] / sum;
	}

	return rGammaValues; //random Dirichlet computed from random gammas
}

// this function is replacement for "sample" in R
//items are items to select from, represented by ints
//probs are relative probabilities (normalized within the function)
//not exported
int sampleC(vector<int> items, vector<double> probs, mt19937 * rNum) {
	
	//input checking and remove items with 0 probability
	if(items.size() != probs.size()) Rcpp::stop("invalid input to sampleC, dimensions differ");
	
	vector<double> newProbs;
	vector<int> newItems;
	
	for(int i=0, max = probs.size(); i<max; i++){
		if(probs[i] < 0) Rcpp::stop("invalid probs to sampleC");
		if(probs[i] > 0){
			newProbs.push_back(probs[i]);
			newItems.push_back(items[i]);
		}
	}
	if (newProbs.size() == 0) Rcpp::stop("all 0 probs to sampleC");
	//end input check
	
	//now change to section of the interval from 0-cumulSum
	double cumulSum = 0; //cumulative sum
	for(int i=0, max=newProbs.size(); i<max; i++){
		cumulSum += newProbs[i];
		newProbs[i] = cumulSum;
	}

	//get random uniform
	uniform_real_distribution <double> rUnif(0.0,cumulSum);
	double rN = rUnif(*rNum);
	//assign result
	int result = 0;
	for(int i=0, max=newProbs.size(); i<max; i++){
		if(rN < newProbs[i]){
			result = newItems[i];
			break;
		}
	}
	
	return result;
}
