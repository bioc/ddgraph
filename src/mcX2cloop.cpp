#include <Rcpp.h>

using namespace Rcpp;

RcppExport SEXP mcX2CLoopC( SEXP B, SEXP numTable, SEXP rowSums, SEXP colSums ){
	Rcpp::NumericMatrix rowS(rowSums);
	Rcpp::NumericMatrix colS(colSums);
	Rcpp::NumericVector numTableVal(numTable);
	Rcpp::NumericVector BVal(B);
	int iB = (int)(BVal[0]); // number of MC replicates
	int iNumTable = (int)(numTableVal[0]); // number of conditional tables
	
	// the result vector
	Rcpp:NumericVector res(iB);
	
	// debug stuff
	bool DEBUG = false;
	bool USE_ALTERNATIVE = false;
	
	// temporary variables
	double newStat = 0; // value of chi-sqaure stat on random table
	double total = 0; // total number of observations
	double diff = 0; // used to calculate chi-square stat
	Rcpp::NumericVector vRowS(2); // row sums for the particular conditional table
	Rcpp::NumericVector vColS(2); // col sums for the particular conditional table
	Rcpp::NumericMatrix E(iNumTable,4); // table of expected values (rows are different tables)
	
	if(DEBUG){
		Rprintf("rowS: %f,%f,%f,%f\\n", rowS(0,0), rowS(1,0), rowS(0,1), rowS(1,1));
		Rprintf("colS: %f,%f,%f,%f\\n", colS(0,0), colS(1,0), colS(0,1), colS(1,1));
	}
	
	// this will have the simulated tables (rows are simulations, columns are divided
	// into serialized 2x2 contingency tables)
	Rcpp::NumericMatrix simTable(iB, iNumTable*4);
	// pre-calculate simulated tables and expected values
	for(int i=0; i < iNumTable; i++){
		// extract row/col sums for this conditional table
		vRowS[0] = rowS(i,0);
		vRowS[1] = rowS(i,1);
		vColS[0] = colS(i,0);
		vColS[1] = colS(i,1);
		
		Language call("r2dtable", iB, vRowS, vColS);
		Rcpp::List tables = call.eval();
		
		// copy all the table into simTable
		int pos = i*4;
		for(int j=0; j < iB; j++ ){
			Rcpp::NumericMatrix st = tables[j];
			simTable(j,pos) = st(0,0);
			simTable(j,pos+1) = st(1,0);
			simTable(j,pos+2) = st(0,1);
			simTable(j,pos+3) = st(1,1);
		}
		
		// calculate the expected values
		total = vRowS[0] + vRowS[1];
		E(i,0) = vRowS[0] * vColS[0] / total;
		E(i,1) = vRowS[1] * vColS[0] / total;		
		E(i,2) = vRowS[0] * vColS[1] / total;		
		E(i,3) = vRowS[1] * vColS[1] / total;
	}
	
	// start replicates
	for(int j=0; j < iB; j++ ){
		newStat = 0;
		// for each replicate do all the conditional tables
		for(int i=0; i < iNumTable; i++){
			vRowS[0] = rowS(i,0);
			vRowS[1] = rowS(i,1);
			vColS[0] = colS(i,0);
			vColS[1] = colS(i,1);
			
			// position in simTable
			int pos = i*4;
			
			if(DEBUG){
				Rprintf("\\n");
				Rprintf("simTable: %f,%f,%f,%f\\n", simTable(j,pos), simTable(j,pos+1), simTable(j,pos+2), simTable(j,pos+3));
				Rprintf("row/col Sums: %f,%f,%f,%f\\n", vRowS[0], vRowS[1], vColS[0], vColS[1]);
				Rprintf("calcula Sums: %f,%f,%f,%f\\n", simTable(j,pos)+simTable(j,pos+2), simTable(j,pos+1)+simTable(j,pos+3), 
					simTable(j,pos)+simTable(j,pos+1), simTable(j,pos+2)+simTable(j,pos+3));
				Rprintf("E: %f,%f,%f,%f\\n", E(i,0), E(i,1), E(i,2), E(i,3));
				
			}
			
			// draw independently from r2dtable
			if(USE_ALTERNATIVE){
				Language call("r2dtable", 1, vRowS, vColS);
				Rcpp::List tables = call.eval();
				Rcpp::NumericMatrix st = tables[0];
				simTable(j,pos) = st(0,0);
				simTable(j,pos+1) = st(1,0);
				simTable(j,pos+2) = st(0,1);
				simTable(j,pos+3) = st(1,1);
			}
			
			// calculate the chi-square statistics
			if(!(vRowS[0] == 0 || vRowS[1] == 0 ||
			   vColS[0] == 0 || vColS[1] == 0)){			
				
				for(int k=0; k<4; k++){
					diff = E(i,k) - simTable(j,pos+k);
					newStat += diff * diff / E(i,k);
				}				
			}
		}
		// store the results of new statistics
		res[j] = newStat;
	}
	
	
	return res;
}
