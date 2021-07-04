/*
The MIT License (MIT) 
Copyright (c) <2017> <(Ken) Shun Hang Yip> <shunyip@bu.edu>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/
#include <R_ext/Rdynload.h>
#include <R.h>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <RcppArmadillo.h>
#include <vector>
#include <limits>
#include <stdint.h>
// [[Rcpp::plugins(cpp11)]]


using namespace std;
using namespace Rcpp;

//Convert dataset into Relative Expression. Please note that even though they are called "Per million", they ceased to be per million during development.
// [[Rcpp::export]]
SEXP XPMCpp (SEXP xSEXP){
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	arma::colvec RowSums = sum(GeneExp, 1);
	int_fast32_t i=0, n=0;
	arma::mat::iterator it_end = GeneExp.end();
	//One pass linear regression with one pass variance, skewness
	for (arma::mat::iterator it = GeneExp.begin(); it != it_end; ++it) {
		//std::cout << (*it) << std::endl;
		*it = *it/RowSums.at(n);
		n++;
		if (n == int(GeneExp.n_rows)) {
			n=0;
			i++;
		}
	}
	return Rcpp::wrap(GeneExp);
}
// [[Rcpp::export]]
SEXP tXPMCpp (SEXP xSEXP){
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	arma::rowvec ColumnSums = sum(GeneExp, 0);
	int_fast32_t i=0, n=0;
	arma::mat::iterator it_end = GeneExp.end();
	//One pass linear regression with one pass variance, skewness
	for (arma::mat::iterator it = GeneExp.begin(); it != it_end; ++it) {
		//std::cout << (*it) << std::endl;
		*it = *it/ColumnSums.at(i);
		n++;
		if (n == int(GeneExp.n_rows)) {
			n=0;
			i++;
		}
	}
	return Rcpp::wrap(trans(GeneExp));
}


//Algorithm to normalize batch effect
// [[Rcpp::export]]
SEXP BatchEffectCpp(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP z2SEXP){
	arma::mat GeneExp = Rcpp::as<arma::mat>(xSEXP);
	arma::mat Output = Rcpp::as<arma::mat>(ySEXP);
	arma::vec meanvec = Rcpp::as<arma::vec>(zSEXP);
	double BE_strength = Rcpp::as<double>(z2SEXP);
	
	//Objects to store the sums of variables that will be needed to perform linear regression.

	arma::vec SumX(GeneExp.n_rows);
	SumX.fill(0);
	arma::vec SumXSq(GeneExp.n_rows);
	SumXSq.fill(0);
	arma::vec SumXY(GeneExp.n_rows);
	SumXY.fill(0);
	arma::vec SumY(GeneExp.n_rows);
	SumY.fill(0);
	arma::vec nRows(GeneExp.n_rows);
	nRows.fill(0);
	
	double store;
	arma::mat::iterator it_end = GeneExp.end();
	int_fast32_t i=0, n=0;
	//One pass linear regression on all rows
	for (arma::mat::iterator it = GeneExp.begin(); it != it_end; ++it) {
		//std::cout << (*it) << std::endl;
		if (*it != 0) {
			store = log(*it);
			SumXY.at(n) += store *  meanvec.at(i);
			SumX.at(n) += store;
			nRows.at(n)++;
			SumXSq.at(n) += pow( store,2);
			SumY.at(n) += meanvec.at(i);
		}
		n++;
		if (n == int(GeneExp.n_rows)) {
			n=0;
			i++;
		}
	}
	
	//y = Mx + C, here are the M and C results of the linear regression
	double add = (1/BE_strength) - 1;
	arma::vec Slope(GeneExp.n_rows);
	arma::vec constant(GeneExp.n_rows);
	for (int_fast32_t k = 0; k < int_fast32_t(GeneExp.n_rows); k++) {
		Slope.at(k) = ((nRows.at(k) * SumXY.at(k) - SumY.at(k) * SumX.at(k))/(nRows.at(k) * SumXSq.at(k) - pow(SumX.at(k),2) )+ add) * BE_strength;
		constant.at(k) = ((SumY.at(k) - Slope.at(k) * SumX.at(k))/nRows.at(k)) * BE_strength;
	}

	
	arma::mat::iterator it_end2 = Output.end();
	i=0;
	n=0;
	//Update Output matrix for output
	for (arma::mat::iterator it2 = Output.begin(); it2 != it_end2; ++it2) {
	  //std::cout << (*it) << std::endl;
	  if (*it2 != 0) {
	   // std::cout <<"Value of n: "<<n<<"Value of i: "<<i<<"Gene Expression :  "<< *it2 <<" Slope is: "<<Slope.at(n)<<" Constant is: "<< constant.at(n) <<std::endl;
	    *it2 = exp(Slope.at(n) * log(*it2) + constant.at(n));
	  }
	  n++;
	  if (n == int(Output.n_rows)) {
	    n=0;
	    i++;
	  }
	}
 return Rcpp::wrap(Output);
}

static const R_CallMethodDef callMethods[] = {
	{"XPMCpp", (DL_FUNC) &XPMCpp, 1},
	{"tXPMCpp", (DL_FUNC) &tXPMCpp, 1},
	{"BatchEffectCpp", (DL_FUNC) &BatchEffectCpp, 4},
	{NULL, NULL, 0}
};

extern "C" {
	void R_init_Linnorm(DllInfo *info) { 
		R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	}
}
