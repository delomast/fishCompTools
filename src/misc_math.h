// This is the header file for misc_math functions
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

#ifndef MISC_MATH_H
#define MISC_MATH_H

#include <Rcpp.h>
#include <random>
#include <vector>
double randBeta(double alpha, double beta, std::mt19937 * rNum);
std::vector<double> randDirich(std::vector<double> alphas, std::mt19937 * rNum);
int sampleC(std::vector<int> items, std::vector<double> probs, std::mt19937 * rNum);

#endif
